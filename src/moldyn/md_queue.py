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
    the worker loop takes no per-call parameters and is identical on every
    invocation.
  - PID-liveness checking: `src/actionlog/action_logger.py`'s `_pid_alive()`
    helper, duplicated here rather than imported — this codebase's
    established convention is to duplicate small private helpers between
    modules rather than share them (see `CLAUDE.md`'s note on
    `_get_tleap_bond_lines` existing separately in `MolDock` and `MolDyn`).

Only one MD job runs at a time, machine-wide, FIFO by enqueue time. A
crashed job (worker killed mid-run, or the machine rebooted) is marked
'crashed' and is never silently re-run: `run_md.sh`'s own `MAX_RESTARTS`
retry logic only covers a single AMBER stage failing within one *live*
script invocation — it has no notion of the wrapper process itself being
killed, so blindly re-running `run_md.sh` could clobber a partially
written stage. A crashed/failed job requires an explicit
`MolDyn.requeue_md_queue_job()` call.
"""

import os
import site
import sqlite3
import subprocess
import sys
from datetime import datetime

from tidyscreen.databases import DatabaseManager as dbm

GLOBAL_QUEUE_DIR = os.path.join(site.getsitepackages()[0], 'tidyscreen', 'projects_db')
GLOBAL_QUEUE_DB = os.path.join(GLOBAL_QUEUE_DIR, 'md_global_queue.db')
LOCK_FILE = os.path.join(GLOBAL_QUEUE_DIR, 'md_queue_worker.lock')
WORKER_SCRIPT_DIR = os.path.join(GLOBAL_QUEUE_DIR, 'md_queue_scripts')
WORKER_LOG_DIR = os.path.join(GLOBAL_QUEUE_DIR, 'md_queue_logs')

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
    'worker_pid': 'INTEGER',
    'run_log_file': 'TEXT',
    'error_message': 'TEXT',
}

_LIST_COLUMNS = [
    'queue_id', 'project_name', 'assay_id', 'md_assay_label', 'status',
    'enqueued_at', 'started_at', 'finished_at', 'error_message',
]
_JOB_COLUMNS = [
    'queue_id', 'project_name', 'project_path', 'assay_id',
    'assay_folder_path', 'md_assay_label',
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


def enqueue_job(project_name, project_path, assay_id, assay_folder_path, md_assay_label=None):
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
        }
        queue_id = dbm.insert_data_dinamically_into_table(cursor, QUEUE_TABLE, data_dict)
        conn.commit()
        return queue_id
    finally:
        conn.close()


def queue_position(queue_id):
    """Number of active (queued/running) jobs at or before this one
    (1 = running or about to run next)."""
    conn, cursor = _connect()
    try:
        cursor.execute(
            f"SELECT enqueued_at FROM {QUEUE_TABLE} WHERE queue_id = ?", (queue_id,)
        )
        row = cursor.fetchone()
        if row is None:
            return 0
        cursor.execute(
            f"SELECT COUNT(*) FROM {QUEUE_TABLE} "
            f"WHERE status IN ('queued','running') AND enqueued_at <= ?",
            (row[0],),
        )
        return cursor.fetchone()[0]
    finally:
        conn.close()


def _worker_lock_alive():
    if not os.path.exists(LOCK_FILE):
        return False
    try:
        with open(LOCK_FILE, 'r') as f:
            content = f.read().strip()
        pid = int(content) if content else None
    except (OSError, ValueError):
        return False
    return _pid_alive(pid)


def _reap_stale_state():
    """
    Self-healing check: if the worker lock is held by a dead PID (worker
    crashed, or the machine rebooted), release it and mark any 'running'
    row whose worker_pid is no longer alive as 'crashed'. Called
    opportunistically by every user-facing entry point (enqueue, list,
    cancel) so the queue never gets stuck without requiring manual DB/file
    cleanup.
    """
    if os.path.exists(LOCK_FILE) and not _worker_lock_alive():
        try:
            os.remove(LOCK_FILE)
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


def ensure_worker_running():
    """
    Make sure a worker process is actively draining the queue, spawning one
    if needed. Safe/cheap to call on every enqueue/list/cancel — it is the
    only place new workers get spawned, and it is what makes the queue
    self-healing after a crash (see _reap_stale_state()).

    Returns True if this call spawned a fresh worker, False otherwise (a
    worker was already running, nothing is queued, or another process won
    the race to spawn one).
    """
    _reap_stale_state()

    if _worker_lock_alive():
        return False

    conn, cursor = _connect()
    try:
        cursor.execute(f"SELECT COUNT(*) FROM {QUEUE_TABLE} WHERE status = 'queued'")
        n_queued = cursor.fetchone()[0]
    finally:
        conn.close()
    if n_queued == 0:
        return False

    os.makedirs(GLOBAL_QUEUE_DIR, exist_ok=True)
    try:
        # Atomic OS-level exclusive create — the mechanism that prevents two
        # near-simultaneous enqueues (possibly from different projects) from
        # both spawning a worker.
        fd = os.open(LOCK_FILE, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        return False  # another process just won the race

    try:
        pid, log_file = _spawn_worker()
        os.write(fd, str(pid).encode())
        os.close(fd)
    except Exception as e:
        try:
            os.close(fd)
        except OSError:
            pass
        try:
            os.remove(LOCK_FILE)
        except OSError:
            pass
        print(f"⚠️  Could not start MD queue worker: {e}")
        return False

    print(f"   🚀 MD queue worker started (pid={pid}); log: {log_file}")
    return True


def _spawn_worker():
    """
    Write a minimal bootstrap script that imports this module and runs the
    worker loop, then launch it fully detached — mirrors
    MolDock._run_docking_background()'s launch mechanics
    (src/moldock/moldock.py) exactly: Popen + start_new_session=True +
    output redirected to a log file. Returns (pid, log_file).
    """
    os.makedirs(WORKER_SCRIPT_DIR, exist_ok=True)
    os.makedirs(WORKER_LOG_DIR, exist_ok=True)

    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    script_file = os.path.join(WORKER_SCRIPT_DIR, f'md_queue_worker_{timestamp}.py')
    log_file = os.path.join(WORKER_LOG_DIR, f'md_queue_worker_{timestamp}.log')

    with open(script_file, 'w') as f:
        f.write(
            "from tidyscreen.moldyn import md_queue\n"
            "md_queue.run_worker_loop()\n"
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


def _claim_next_job():
    """
    Atomically claim the oldest 'queued' row (FIFO by enqueued_at) for this
    worker process. Returns a dict with the claimed job's fields, or None
    if the queue is empty.
    """
    conn, cursor = _connect()
    try:
        while True:
            cursor.execute(
                f"SELECT {', '.join(_JOB_COLUMNS)} FROM {QUEUE_TABLE} "
                f"WHERE status = 'queued' ORDER BY enqueued_at ASC, queue_id ASC LIMIT 1"
            )
            row = cursor.fetchone()
            if row is None:
                return None
            queue_id = row[0]
            cursor.execute(
                f"UPDATE {QUEUE_TABLE} SET status = 'running', started_at = ?, worker_pid = ? "
                f"WHERE queue_id = ? AND status = 'queued'",
                (_now(), os.getpid(), queue_id),
            )
            conn.commit()
            if cursor.rowcount == 1:
                return dict(zip(_JOB_COLUMNS, row))
            # Someone else (e.g. a concurrent cancel) changed this row
            # between the SELECT and the UPDATE above — retry with the
            # next-oldest candidate.
    finally:
        conn.close()


def _mirror_status_to_project_db(project_path, assay_id, status):
    """
    Best-effort: reflect the job's final status into that project's own
    md_assays table (md_registers.db), lazily adding the 'status' column if
    needed — mirrors the existing lazy ALTER TABLE convention used for
    mmgbsa_results (see MolDyn._store_mmgbsa_results). Never raises — a
    failure here must not affect the authoritative md_queue row.
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
    """Run one claimed job's run_md.sh to completion (blocking) and record
    the outcome. Called only from inside run_worker_loop()."""
    queue_id = job['queue_id']
    assay_folder_path = job['assay_folder_path']
    script_path = os.path.join(assay_folder_path, 'run_md.sh')
    run_log = os.path.join(assay_folder_path, 'run_md_queue.log')

    conn, cursor = _connect()
    try:
        cursor.execute(f"UPDATE {QUEUE_TABLE} SET run_log_file = ? WHERE queue_id = ?", (run_log, queue_id))
        conn.commit()
    finally:
        conn.close()

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

    conn, cursor = _connect()
    try:
        cursor.execute(
            f"UPDATE {QUEUE_TABLE} SET status = ?, finished_at = ?, error_message = ? WHERE queue_id = ?",
            (status, _now(), error_message, queue_id),
        )
        conn.commit()
    finally:
        conn.close()

    _mirror_status_to_project_db(job['project_path'], job['assay_id'], status)


def run_worker_loop():
    """
    Entry point for the detached worker process (see _spawn_worker()).
    Repeatedly claims and runs the oldest queued job — machine-wide, across
    every project — until the queue is empty, then releases the worker
    lock and exits. A fresh worker is spawned by the next
    ensure_worker_running() call (from any project) once something new is
    enqueued.
    """
    while True:
        job = _claim_next_job()
        if job is None:
            try:
                os.remove(LOCK_FILE)
            except OSError:
                pass
            # One more check to avoid a lost wakeup: a job may have been
            # enqueued in the instant between the empty SELECT above and
            # removing the lock. Anything that still slips through this
            # narrow window is picked up by the next ensure_worker_running()
            # call from any project — this loop does not try to be a
            # perfectly race-free protocol, only self-healing on next touch.
            job = _claim_next_job()
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
