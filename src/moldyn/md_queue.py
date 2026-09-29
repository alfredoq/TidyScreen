"""
Machine-wide MD run queue.

Serializes AMBER MD executions across every TidyScreen project on this
machine — GPU contention is a single physical-machine resource, not a
per-project one, and `moldyn.py` has no CUDA device-selection logic to
otherwise keep two `run_md.sh` runs from fighting over the same GPU.

This lives as its own module (not private `MolDyn` methods) because it
must be importable both by `MolDyn` (to enqueue/inspect/cancel jobs from
an interactive session) and by the standalone, detached worker process
that actually drains the queue — spawned via `subprocess.Popen`, see
`_spawn_worker()`.

State lives in a global sqlite database, sibling to the projects registry
itself (`<site-packages>/tidyscreen/projects_db/`), so it is scoped to
this Python environment/machine, not to any one project — mirrors the
existing `projects_database.db` path convention used throughout
`tidyscreen.py` (`f"{site.getsitepackages()[0]}/tidyscreen/projects_db/projects_database.db"`).

Design mirrors two conventions already established elsewhere in this
codebase:
  - Detached background launch: `MolDock._run_docking_background()`
    (`src/moldock/moldock.py`) — `subprocess.Popen(..., start_new_session=True)`
    with stdout/stderr redirected to a log file. The worker here reuses
    those exact launch mechanics; only its logic is factored differently
    (an importable module rather than a regenerated inline script), since
    the worker loop takes no per-call parameters beyond its lane and is
    identical on every invocation.
  - PID-liveness checking: `src/actionlog/action_logger.py`'s `_pid_alive()`
    helper, duplicated here rather than imported — this codebase's
    established convention is to duplicate small private helpers between
    modules rather than share them (see `CLAUDE.md`'s note on
    `_get_tleap_bond_lines` existing separately in `MolDock` and `MolDyn`).

Lanes: every job belongs to one resource lane (the `resource` column), and
each lane has its own worker process and lock file:
  - 'gpu' — MD runs (`job_type='md'`, `run_md.sh`, pmemd.cuda).
  - 'cpu' — MM-GBSA runs (`job_type='mmgbsa'`, `run_mmgbsa.sh`,
    MMPBSA.py[.MPI]), which don't touch the GPU. Keeping them in their own
    lane means a minutes-long MM-GBSA job never waits behind a days-long MD
    run, while two (possibly MPI, many-core) MM-GBSA jobs still never
    compete with each other for the CPU.
Within a lane only one job runs at a time, machine-wide, FIFO by enqueue
time; the lanes run independently of each other. Rows predating the lane
columns are treated as 'md'/'gpu' (see `_RESOURCE_SQL`).

A crashed job (worker killed mid-run, or the machine rebooted) is marked
'crashed' and is never silently re-run: `run_md.sh`'s own `MAX_RESTARTS`
retry logic only covers a single AMBER stage failing within one *live*
script invocation — it has no notion of the wrapper process itself being
killed, so blindly re-running `run_md.sh` could clobber a partially
written stage. A crashed/failed job requires an explicit
`MolDyn.requeue_md_queue_job()` call.

An MM-GBSA job's results are parsed and stored into the project's
md_registers.db by the worker itself once the run finishes (see
`_store_mmgbsa_job_results()`), exactly as a foreground run of
`MolDyn.compute_mmgbsa_on_trajectory()` would.
"""

import json
import os
import site
import sqlite3
import subprocess
import sys
import time
from datetime import datetime

from tidyscreen.databases import DatabaseManager as dbm

GLOBAL_QUEUE_DIR = os.path.join(site.getsitepackages()[0], 'tidyscreen', 'projects_db')
GLOBAL_QUEUE_DB = os.path.join(GLOBAL_QUEUE_DIR, 'md_global_queue.db')
# 'gpu' lane lock — file name kept from before lanes existed, so a worker
# spawned by the previous version is still recognised as the gpu worker.
LOCK_FILE = os.path.join(GLOBAL_QUEUE_DIR, 'md_queue_worker.lock')
WORKER_SCRIPT_DIR = os.path.join(GLOBAL_QUEUE_DIR, 'md_queue_scripts')
WORKER_LOG_DIR = os.path.join(GLOBAL_QUEUE_DIR, 'md_queue_logs')

RESOURCES = ('gpu', 'cpu')
_LOCK_FILES = {
    'gpu': LOCK_FILE,
    'cpu': os.path.join(GLOBAL_QUEUE_DIR, 'md_queue_worker_cpu.lock'),
}
_LANE_LABELS = {'gpu': 'MD (GPU)', 'cpu': 'MM-GBSA (CPU)'}

# Legacy rows (enqueued before the lane columns were added) are GPU/MD jobs.
_RESOURCE_SQL = "COALESCE(resource, 'gpu')"

QUEUE_TABLE = 'md_queue'

# Additive-only schema — never run dbm.remove_legacy_table_columns() against
# this table (see module docstring).
COLUMNS_DICT = {
    'queue_id': 'INTEGER PRIMARY KEY AUTOINCREMENT',
    'project_name': 'TEXT NOT NULL',
    'project_path': 'TEXT NOT NULL',
    'assay_id': 'INTEGER NOT NULL',
    'assay_folder_path': 'TEXT NOT NULL',
    'md_assay_label': 'TEXT',
    'status': "TEXT NOT NULL DEFAULT 'queued'",  # queued|running|completed|failed|cancelled|crashed
    'enqueued_at': 'TEXT NOT NULL',
    'started_at': 'TEXT',
    'finished_at': 'TEXT',
    # For a queue-managed job: the worker process's own pid (it runs
    # the job's script via a blocking subprocess.run inside itself). For a
    # foreground-started job (see register_foreground_run()): the pid of
    # the script process itself, since there is no separate worker.
    # Either way, _pid_alive(worker_pid) is what "is this row still
    # actually executing" means throughout this module.
    'worker_pid': 'INTEGER',
    'run_log_file': 'TEXT',
    'error_message': 'TEXT',
    'job_type': "TEXT DEFAULT 'md'",   # md|mmgbsa
    'resource': "TEXT DEFAULT 'gpu'",  # gpu|cpu — the lane (see module docstring)
    # Script to run. NULL for 'md' jobs (always <assay_folder_path>/run_md.sh);
    # for 'mmgbsa' jobs the run_mmgbsa.sh inside the specific mmgbsa[_N] run
    # folder, since one MD assay can have several MM-GBSA runs.
    'script_path': 'TEXT',
    # JSON: for 'mmgbsa' jobs, {'mmgbsa_params': ..., 'assay_info': ...} —
    # what the worker needs to store the results once the run finishes.
    'job_params': 'TEXT',
}

_LIST_COLUMNS = [
    'queue_id', 'project_name', 'assay_id', 'md_assay_label', 'status',
    'enqueued_at', 'started_at', 'finished_at', 'error_message',
    'job_type', 'resource', 'script_path',
]
_JOB_COLUMNS = [
    'queue_id', 'project_name', 'project_path', 'assay_id',
    'assay_folder_path', 'md_assay_label', 'job_type', 'script_path', 'job_params',
]


def _now():
    return datetime.now().strftime('%Y-%m-%d %H:%M:%S')


def _pid_alive(pid):
    """Mirrors src/actionlog/action_logger.py's _pid_alive() (duplicated
    rather than imported — see module docstring)."""
    if not pid:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True  # process exists, just owned by someone else
    except OSError:
        return False
    return True


def _connect():
    os.makedirs(GLOBAL_QUEUE_DIR, exist_ok=True)
    conn = sqlite3.connect(GLOBAL_QUEUE_DB, timeout=30)
    cursor = conn.cursor()
    dbm.create_table_from_columns_dict(cursor, QUEUE_TABLE, COLUMNS_DICT, verbose=False)
    dbm.update_legacy_table_columns(cursor, QUEUE_TABLE, COLUMNS_DICT, verbose=False)
    conn.commit()
    return conn, cursor


def enqueue_job(project_name, project_path, assay_id, assay_folder_path, md_assay_label=None,
                job_type='md', resource='gpu', script_path=None, job_params=None):
    """Insert a new 'queued' row and return its queue_id."""
    conn, cursor = _connect()
    try:
        data_dict = {
            'project_name': project_name,
            'project_path': project_path,
            'assay_id': assay_id,
            'assay_folder_path': assay_folder_path,
            'md_assay_label': md_assay_label,
            'status': 'queued',
            'enqueued_at': _now(),
            'job_type': job_type,
            'resource': resource,
            'script_path': script_path,
            'job_params': json.dumps(job_params) if job_params is not None else None,
        }
        queue_id = dbm.insert_data_dinamically_into_table(cursor, QUEUE_TABLE, data_dict)
        conn.commit()
        return queue_id
    finally:
        conn.close()


def queue_position(queue_id):
    """Number of active (queued/running) jobs in the same lane at or before
    this one (1 = running or about to run next)."""
    conn, cursor = _connect()
    try:
        cursor.execute(
            f"SELECT enqueued_at, {_RESOURCE_SQL} FROM {QUEUE_TABLE} WHERE queue_id = ?", (queue_id,)
        )
        row = cursor.fetchone()
        if row is None:
            return 0
        cursor.execute(
            f"SELECT COUNT(*) FROM {QUEUE_TABLE} "
            f"WHERE status IN ('queued','running') AND enqueued_at <= ? AND {_RESOURCE_SQL} = ?",
            (row[0], row[1]),
        )
        return cursor.fetchone()[0]
    finally:
        conn.close()


def can_start_foreground(resource='gpu'):
    """
    Cheap pre-check for a foreground start: True if nothing is currently
    'running' in this lane of the shared queue (queue-managed *or* a
    previously registered foreground run). Callers should still call
    register_foreground_run() atomically right before actually starting
    the process — this function only avoids launching a process that is
    obviously going to have to be refused, it does not itself close the
    (small, human-interactive-timescale) race between this check and the
    process actually starting.
    """
    _reap_stale_state()
    conn, cursor = _connect()
    try:
        cursor.execute(
            f"SELECT COUNT(*) FROM {QUEUE_TABLE} WHERE status = 'running' AND {_RESOURCE_SQL} = ?",
            (resource,),
        )
        running = cursor.fetchone()[0]
    finally:
        conn.close()
    if running > 0:
        return False, f"Another {_LANE_LABELS[resource]} job is already running on the shared MD queue."
    return True, None


def register_foreground_run(project_name, project_path, assay_id, assay_folder_path, md_assay_label, pid,
                            job_type='md', resource='gpu', script_path=None, job_params=None):
    """
    Register a directly (foreground-)started run as a 'running' row in its
    lane so other sessions calling ensure_worker_running()/
    can_start_foreground() see it and won't start a second job in that lane
    concurrently — this is what makes 'fg' and 'queue' starts correctly
    wait on each other, not just queue starts on other queue starts.

    Must be called with the process's real pid *after* it has actually
    been started (so a live pid is always available to check liveness
    against later — see _reap_stale_state()), but the check-and-insert
    below is still atomic (BEGIN IMMEDIATE) against a concurrent claim/
    registration happening in the tiny window since can_start_foreground()
    was last called.

    Returns (queue_id, None) on success, or (None, error_message) if
    something else is already running in the lane — the caller already
    started the process in that case and must decide how to handle it (see
    MolDyn._start_md_simulation()'s fg branch).
    """
    conn, cursor = _connect()
    try:
        conn.execute("BEGIN IMMEDIATE")
        cursor.execute(
            f"SELECT COUNT(*) FROM {QUEUE_TABLE} WHERE status = 'running' AND {_RESOURCE_SQL} = ?",
            (resource,),
        )
        if cursor.fetchone()[0] > 0:
            conn.rollback()
            return None, f"Another {_LANE_LABELS[resource]} job is already running on the shared MD queue."
        data_dict = {
            'project_name': project_name,
            'project_path': project_path,
            'assay_id': assay_id,
            'assay_folder_path': assay_folder_path,
            'md_assay_label': md_assay_label,
            'status': 'running',
            'enqueued_at': _now(),
            'started_at': _now(),
            'worker_pid': pid,  # here: the script process's own pid, not a queue worker's
            'job_type': job_type,
            'resource': resource,
            'script_path': script_path,
            # Kept so the row can be requeued as a worker job later.
            'job_params': json.dumps(job_params) if job_params is not None else None,
        }
        queue_id = dbm.insert_data_dinamically_into_table(cursor, QUEUE_TABLE, data_dict)
        conn.commit()
        return queue_id, None
    finally:
        conn.close()


def finish_foreground_run(queue_id, project_path, assay_id, status, error_message=None, job_type='md'):
    """Record the outcome of a foreground run registered via
    register_foreground_run(), then release anything that queued up while
    it was running."""
    _finish_row(queue_id, status, error_message)
    if job_type == 'md':
        _mirror_status_to_project_db(project_path, assay_id, status)
    ensure_worker_running()


def _set_run_log(queue_id, run_log):
    conn, cursor = _connect()
    try:
        cursor.execute(f"UPDATE {QUEUE_TABLE} SET run_log_file = ? WHERE queue_id = ?", (run_log, queue_id))
        conn.commit()
    finally:
        conn.close()


def _finish_row(queue_id, status, error_message):
    conn, cursor = _connect()
    try:
        cursor.execute(
            f"UPDATE {QUEUE_TABLE} SET status = ?, finished_at = ?, error_message = ? WHERE queue_id = ?",
            (status, _now(), error_message, queue_id),
        )
        conn.commit()
    finally:
        conn.close()


def _worker_lock_alive(resource='gpu'):
    lock_file = _LOCK_FILES[resource]
    if not os.path.exists(lock_file):
        return False
    try:
        with open(lock_file, 'r') as f:
            content = f.read().strip()
        pid = int(content) if content else None
    except (OSError, ValueError):
        return False
    return _pid_alive(pid)


def _reap_stale_state():
    """
    Self-healing check: if a lane's worker lock is held by a dead PID
    (worker crashed, or the machine rebooted), release it and mark any
    'running' row whose worker_pid is no longer alive as 'crashed'. Called
    opportunistically by every user-facing entry point (enqueue, list,
    cancel) so the queue never gets stuck without requiring manual DB/file
    cleanup.
    """
    for resource in RESOURCES:
        lock_file = _LOCK_FILES[resource]
        if os.path.exists(lock_file) and not _worker_lock_alive(resource):
            try:
                os.remove(lock_file)
            except OSError:
                pass

    conn, cursor = _connect()
    try:
        cursor.execute(f"SELECT queue_id, worker_pid FROM {QUEUE_TABLE} WHERE status = 'running'")
        rows = cursor.fetchall()
        for queue_id, worker_pid in rows:
            if not _pid_alive(worker_pid):
                cursor.execute(
                    f"UPDATE {QUEUE_TABLE} SET status = 'crashed', finished_at = ?, error_message = ? "
                    f"WHERE queue_id = ? AND status = 'running'",
                    (_now(), 'Worker process no longer alive (crash or reboot); '
                             'requeue manually with MolDyn.requeue_md_queue_job() if appropriate.', queue_id),
                )
        conn.commit()
    finally:
        conn.close()


def ensure_worker_running(resource=None):
    """
    Make sure a worker process is actively draining the given lane
    (every lane when resource is None), spawning one if needed. Safe/cheap
    to call on every enqueue/list/cancel — it is the only place new workers
    get spawned, and it is what makes the queue self-healing after a crash
    (see _reap_stale_state()).

    Returns True if this call spawned a fresh worker, False otherwise (a
    worker was already running, nothing is queued, or another process won
    the race to spawn one).
    """
    _reap_stale_state()
    resources = RESOURCES if resource is None else (resource,)
    # List, not a generator, so every lane gets checked (no short-circuit).
    return any([_ensure_lane_worker(r) for r in resources])


def _ensure_lane_worker(resource):
    if _worker_lock_alive(resource):
        return False

    conn, cursor = _connect()
    try:
        # A 'running' row can belong either to a queue-managed job (whose
        # worker holds the lane's lock file — already checked above) or to a
        # directly foreground-started job (see register_foreground_run()),
        # which holds no worker lock at all. Checking the table directly, not
        # just the lock file, is what makes fg and queue starts correctly
        # wait on each other.
        cursor.execute(
            f"SELECT COUNT(*) FROM {QUEUE_TABLE} WHERE status = 'running' AND {_RESOURCE_SQL} = ?",
            (resource,),
        )
        if cursor.fetchone()[0] > 0:
            return False
        cursor.execute(
            f"SELECT COUNT(*) FROM {QUEUE_TABLE} WHERE status = 'queued' AND {_RESOURCE_SQL} = ?",
            (resource,),
        )
        n_queued = cursor.fetchone()[0]
    finally:
        conn.close()
    if n_queued == 0:
        return False

    lock_file = _LOCK_FILES[resource]
    os.makedirs(GLOBAL_QUEUE_DIR, exist_ok=True)
    try:
        # Atomic OS-level exclusive create — the mechanism that prevents two
        # near-simultaneous enqueues (possibly from different projects) from
        # both spawning a worker for the same lane.
        fd = os.open(lock_file, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        return False  # another process just won the race

    try:
        pid, log_file = _spawn_worker(resource)
        os.write(fd, str(pid).encode())
        os.close(fd)
    except Exception as e:
        try:
            os.close(fd)
        except OSError:
            pass
        try:
            os.remove(lock_file)
        except OSError:
            pass
        print(f"⚠️  Could not start {_LANE_LABELS[resource]} queue worker: {e}")
        return False

    print(f"   🚀 {_LANE_LABELS[resource]} queue worker started (pid={pid}); log: {log_file}")
    return True


def _spawn_worker(resource='gpu'):
    """
    Write a minimal bootstrap script that imports this module and runs the
    worker loop for one lane, then launch it fully detached — mirrors
    MolDock._run_docking_background()'s launch mechanics
    (src/moldock/moldock.py) exactly: Popen + start_new_session=True +
    output redirected to a log file. Returns (pid, log_file).
    """
    os.makedirs(WORKER_SCRIPT_DIR, exist_ok=True)
    os.makedirs(WORKER_LOG_DIR, exist_ok=True)

    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    script_file = os.path.join(WORKER_SCRIPT_DIR, f'md_queue_worker_{resource}_{timestamp}.py')
    log_file = os.path.join(WORKER_LOG_DIR, f'md_queue_worker_{resource}_{timestamp}.log')

    with open(script_file, 'w') as f:
        f.write(
            "from tidyscreen.moldyn import md_queue\n"
            f"md_queue.run_worker_loop({resource!r})\n"
        )
    os.chmod(script_file, 0o755)

    with open(log_file, 'w') as logf:
        process = subprocess.Popen(
            [sys.executable, script_file],
            stdout=logf,
            stderr=subprocess.STDOUT,
            stdin=subprocess.DEVNULL,
            start_new_session=True,
        )
    return process.pid, log_file


def _claim_next_job(resource='gpu'):
    """
    Atomically claim the oldest 'queued' row of this lane (FIFO by
    enqueued_at) for this worker process — but only if nothing else in the
    lane is currently 'running'.

    That last check matters even though ensure_worker_running() already
    refuses to *spawn* a worker while something is running: a foreground
    job (see register_foreground_run()) can start concurrently with an
    already-running worker's brief gap between finishing one job and
    claiming the next, so the claim itself must re-verify exclusivity, not
    just rely on the spawn-time check. The whole check-then-claim runs
    inside one BEGIN IMMEDIATE transaction so it is atomic against any
    other connection (another worker, or a foreground run registering)
    doing the same check concurrently.

    Returns a dict with the claimed job's fields, or None if the lane is
    empty or something else in it is currently running (the caller should
    treat both cases the same way — see run_worker_loop()).
    """
    conn, cursor = _connect()
    try:
        conn.execute("BEGIN IMMEDIATE")
        cursor.execute(
            f"SELECT COUNT(*) FROM {QUEUE_TABLE} WHERE status = 'running' AND {_RESOURCE_SQL} = ?",
            (resource,),
        )
        if cursor.fetchone()[0] > 0:
            conn.rollback()
            return None

        cursor.execute(
            f"SELECT {', '.join(_JOB_COLUMNS)} FROM {QUEUE_TABLE} "
            f"WHERE status = 'queued' AND {_RESOURCE_SQL} = ? "
            f"ORDER BY enqueued_at ASC, queue_id ASC LIMIT 1",
            (resource,),
        )
        row = cursor.fetchone()
        if row is None:
            conn.rollback()
            return None

        queue_id = row[0]
        cursor.execute(
            f"UPDATE {QUEUE_TABLE} SET status = 'running', started_at = ?, worker_pid = ? WHERE queue_id = ?",
            (_now(), os.getpid(), queue_id),
        )
        conn.commit()
        return dict(zip(_JOB_COLUMNS, row))
    finally:
        conn.close()


def _mirror_status_to_project_db(project_path, assay_id, status):
    """
    Best-effort: reflect an MD job's final status into that project's own
    md_assays table (md_registers.db), lazily adding the 'status' column if
    needed — mirrors the existing lazy ALTER TABLE convention used for
    mmgbsa_results (see MolDyn._store_mmgbsa_results). Never raises — a
    failure here must not affect the authoritative md_queue row.

    Only for job_type 'md': an MM-GBSA job shares its MD assay's assay_id,
    and must not overwrite that assay's MD status with its own.
    """
    try:
        md_registers_db = os.path.join(project_path, 'dynamics', 'md_registers', 'md_registers.db')
        if not os.path.exists(md_registers_db):
            return
        conn = sqlite3.connect(md_registers_db)
        cursor = conn.cursor()
        cursor.execute("PRAGMA table_info(md_assays)")
        existing_columns = {row[1] for row in cursor.fetchall()}
        if 'status' not in existing_columns:
            cursor.execute("ALTER TABLE md_assays ADD COLUMN status TEXT")
        cursor.execute("UPDATE md_assays SET status = ? WHERE assay_id = ?", (status, assay_id))
        conn.commit()
        conn.close()
    except Exception:
        pass


def _run_job(job):
    """Run one claimed job to completion (blocking) and record the outcome.
    Called only from inside run_worker_loop()."""
    if (job.get('job_type') or 'md') == 'mmgbsa':
        _run_mmgbsa_job(job)
        return

    queue_id = job['queue_id']
    assay_folder_path = job['assay_folder_path']
    script_path = os.path.join(assay_folder_path, 'run_md.sh')
    run_log = os.path.join(assay_folder_path, 'run_md_queue.log')

    _set_run_log(queue_id, run_log)

    error_message = None
    if not os.path.exists(script_path):
        status = 'failed'
        error_message = f'run_md.sh not found at {script_path}'
    else:
        try:
            with open(run_log, 'w') as logf:
                result = subprocess.run(
                    [script_path], cwd=assay_folder_path, stdout=logf, stderr=subprocess.STDOUT,
                )
            status = 'completed' if result.returncode == 0 else 'failed'
            if status == 'failed':
                error_message = f'run_md.sh exited with code {result.returncode}'
        except Exception as e:
            status = 'failed'
            error_message = f'Error launching run_md.sh: {e}'

    _finish_row(queue_id, status, error_message)
    _mirror_status_to_project_db(job['project_path'], job['assay_id'], status)


def _run_mmgbsa_job(job):
    """
    Run one claimed MM-GBSA job's run_mmgbsa.sh (blocking), then parse and
    store its results into the project (see _store_mmgbsa_job_results()).
    The row stays 'running' — holding the cpu lane — until the results are
    stored too.
    """
    queue_id = job['queue_id']
    script_path = job.get('script_path')
    run_folder = os.path.dirname(script_path) if script_path else job['assay_folder_path']
    run_log = os.path.join(run_folder, 'run_mmgbsa_queue.log')

    _set_run_log(queue_id, run_log)

    error_message = None
    if not script_path or not os.path.exists(script_path):
        status = 'failed'
        error_message = f'run_mmgbsa.sh not found at {script_path}'
    else:
        results_dat = os.path.join(run_folder, 'mmgbsa_results.dat')
        started = time.time()
        try:
            with open(run_log, 'w') as logf:
                result = subprocess.run(
                    [script_path], cwd=run_folder, stdout=logf, stderr=subprocess.STDOUT,
                )
            # run_mmgbsa.sh ends with an unconditional echo, so a zero exit
            # code doesn't prove MMPBSA.py succeeded — a results file written
            # by this run is the real success marker.
            if result.returncode != 0:
                status = 'failed'
                error_message = f'run_mmgbsa.sh exited with code {result.returncode}'
            elif not os.path.exists(results_dat) or os.path.getmtime(results_dat) < started - 1:
                status = 'failed'
                error_message = (f'MMPBSA.py produced no results file ({results_dat}); '
                                 f'see {os.path.join(run_folder, "mmpbsa.log")}')
            else:
                status = 'completed'
                error_message = _store_mmgbsa_job_results(job, run_folder, run_log)
        except Exception as e:
            status = 'failed'
            error_message = f'Error launching run_mmgbsa.sh: {e}'

    _finish_row(queue_id, status, error_message)


def _store_mmgbsa_job_results(job, run_folder, run_log):
    """
    Parse and store a finished MM-GBSA job's results into its project, via
    the same MolDyn._finalize_mmgbsa_run() a foreground run uses. Output
    (results summary, decomposition table) is appended to the job's run log.

    Returns None on success, or a message for the row's error_message. The
    job itself stays 'completed' either way: the calculation succeeded, and
    the results can still be stored later by re-running
    compute_mmgbsa_on_trajectory() and choosing to parse the existing run.
    """
    import contextlib

    fallback = ('re-run MolDyn.compute_mmgbsa_on_trajectory() and choose to parse '
                'the existing run to store them')
    try:
        params = json.loads(job.get('job_params') or '{}')
        with open(run_log, 'a') as logf, contextlib.redirect_stdout(logf), contextlib.redirect_stderr(logf):
            print("\n--- Storing MM-GBSA results (md_queue worker) ---")
            # Imported here, not at module level: moldyn.py imports this module.
            from tidyscreen.tidyscreen import ActivateProject
            from tidyscreen.moldyn.moldyn import MolDyn
            moldyn = MolDyn(ActivateProject(job['project_name']))
            stored = moldyn._finalize_mmgbsa_run(
                job['assay_id'], run_folder,
                params.get('mmgbsa_params', {}), params.get('assay_info', {}),
            )
        if not stored:
            return f'MM-GBSA finished but results could not be parsed/stored (see {run_log}); {fallback}.'
        return None
    except Exception as e:
        return f'MM-GBSA finished but results could not be stored ({e}); {fallback}.'


def run_worker_loop(resource='gpu'):
    """
    Entry point for the detached worker process of one lane (see
    _spawn_worker()). Repeatedly claims and runs the lane's oldest queued
    job — machine-wide, across every project — until _claim_next_job()
    returns None, then releases the lane's worker lock and exits. That
    happens either because the lane is genuinely empty, or because a
    foreground-started run (see register_foreground_run()) is currently
    occupying the lane's "one job at a time" slot; either way it is safe to
    exit here, since whichever run is still 'running' will call
    ensure_worker_running() itself when it finishes (this worker's own
    next-job claim for the queue-managed case, or finish_foreground_run()
    for the foreground case) — spawning a fresh worker then if anything is
    still queued.
    """
    lock_file = _LOCK_FILES[resource]
    while True:
        job = _claim_next_job(resource)
        if job is None:
            try:
                os.remove(lock_file)
            except OSError:
                pass
            # One more check to avoid a lost wakeup: a job may have been
            # enqueued in the instant between the empty SELECT above and
            # removing the lock. Anything that still slips through this
            # narrow window is picked up by the next ensure_worker_running()
            # call from any project — this loop does not try to be a
            # perfectly race-free protocol, only self-healing on next touch.
            job = _claim_next_job(resource)
            if job is None:
                return
        _run_job(job)


def list_jobs(project_name=None):
    """Return queue rows (as dicts) ordered by enqueue time, optionally
    filtered to one project."""
    conn, cursor = _connect()
    try:
        if project_name:
            cursor.execute(
                f"SELECT {', '.join(_LIST_COLUMNS)} FROM {QUEUE_TABLE} "
                f"WHERE project_name = ? ORDER BY enqueued_at ASC, queue_id ASC",
                (project_name,),
            )
        else:
            cursor.execute(
                f"SELECT {', '.join(_LIST_COLUMNS)} FROM {QUEUE_TABLE} "
                f"ORDER BY enqueued_at ASC, queue_id ASC"
            )
        rows = cursor.fetchall()
    finally:
        conn.close()
    return [dict(zip(_LIST_COLUMNS, row)) for row in rows]


def cancel_job(queue_id):
    """Cancel a job that has not started yet. Returns (success, error_message|None)."""
    conn, cursor = _connect()
    try:
        cursor.execute(
            f"UPDATE {QUEUE_TABLE} SET status = 'cancelled', finished_at = ? "
            f"WHERE queue_id = ? AND status = 'queued'",
            (_now(), queue_id),
        )
        conn.commit()
        if cursor.rowcount == 1:
            return True, None
        cursor.execute(f"SELECT status FROM {QUEUE_TABLE} WHERE queue_id = ?", (queue_id,))
        row = cursor.fetchone()
        if row is None:
            return False, f"No job with queue_id={queue_id}."
        return False, f"Job {queue_id} is '{row[0]}'; only 'queued' jobs can be cancelled."
    finally:
        conn.close()


def requeue_job(queue_id):
    """Re-queue a failed/crashed/cancelled job. Returns (success, error_message|None)."""
    conn, cursor = _connect()
    try:
        cursor.execute(
            f"UPDATE {QUEUE_TABLE} SET status = 'queued', started_at = NULL, finished_at = NULL, "
            f"worker_pid = NULL, error_message = NULL, enqueued_at = ? "
            f"WHERE queue_id = ? AND status IN ('failed', 'crashed', 'cancelled')",
            (_now(), queue_id),
        )
        conn.commit()
        if cursor.rowcount == 1:
            return True, None
        cursor.execute(f"SELECT status FROM {QUEUE_TABLE} WHERE queue_id = ?", (queue_id,))
        row = cursor.fetchone()
        if row is None:
            return False, f"No job with queue_id={queue_id}."
        return False, f"Job {queue_id} is '{row[0]}'; only failed/crashed/cancelled jobs can be requeued."
    finally:
        conn.close()
