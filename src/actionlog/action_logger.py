"""
Project-scoped activity log recording every public API/GUI call made against
MolDock, MolDyn, ChemSpace and MachineLearning, including background subprocess
launches (docking, MD, fingerprinting, CSV loading, RF predictions, etc).

The log lives at ``<project_path>/project_action_log.db`` (sibling of
``project_notes.db``) and is populated transparently via the ``log_action``
decorator / ``log_all_public_methods`` class decorator applied to each domain
class. Background-job completion is never actively reported by the child
process; it is inferred at read time by ``refresh_action_log_statuses`` via
PID liveness plus a log-file tail check, mirroring the existing
``get_md_assay_status`` filesystem-heuristic convention.
"""

import functools
import os
import sqlite3
import sys
import threading
from datetime import datetime, timedelta
from typing import Any, Dict, List, Optional

from tidyscreen.databases import DatabaseManager as dbm

ACTION_LOG_DB_FILENAME = 'project_action_log.db'
ACTION_LOG_TABLE = 'action_log'
MAX_LOG_ROWS = 2000

_COLUMNS_DICT = {
    'id': 'INTEGER PRIMARY KEY AUTOINCREMENT',
    'started_at': 'TEXT NOT NULL',
    'ended_at': 'TEXT',
    'class_name': 'TEXT NOT NULL',
    'method_name': 'TEXT NOT NULL',
    'source': 'TEXT NOT NULL',
    'process_pid': 'INTEGER NOT NULL',
    'status': "TEXT NOT NULL DEFAULT 'running'",
    'background': 'INTEGER NOT NULL DEFAULT 0',
    'child_pid': 'INTEGER',
    'log_file': 'TEXT',
    'error_message': 'TEXT',
    'args_summary': 'TEXT',
    'description': 'TEXT',
}

# Reserved kwarg name a caller can pass to any logged method to attach a
# free-text description to its log row, e.g. moldock_obj.dock_table(log_description="...").
# The decorator strips it before calling the real method.
LOG_DESCRIPTION_KWARG = 'log_description'

_COMPLETION_MARKERS_SUCCESS = ('Total wall time:', 'FINAL RESULTS')
_COMPLETION_MARKERS_ERROR = ('Traceback (most recent call last):',)


def _now() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def _db_path(project_path: str) -> str:
    return os.path.join(project_path, ACTION_LOG_DB_FILENAME)


def _ensure_table(cursor: sqlite3.Cursor) -> None:
    # Additive-only: never call remove_legacy_table_columns here, it can drop
    # columns without confirmation when verbose=False.
    dbm.create_table_from_columns_dict(cursor, ACTION_LOG_TABLE, _COLUMNS_DICT, verbose=False)
    dbm.update_legacy_table_columns(cursor, ACTION_LOG_TABLE, _COLUMNS_DICT, verbose=False)


class ActionLogger:
    """Records TidyScreen API/GUI method calls for a single project."""

    def __init__(self, project_path: str):
        self.project_path = project_path
        self.db_path = _db_path(project_path)

    def start(
        self,
        class_name: str,
        method_name: str,
        source: str,
        args_summary: str = '',
        description: Optional[str] = None,
    ) -> int:
        conn = sqlite3.connect(self.db_path)
        try:
            cursor = conn.cursor()
            _ensure_table(cursor)
            data_dict = {
                'started_at': _now(),
                'class_name': class_name,
                'method_name': method_name,
                'source': source,
                'process_pid': os.getpid(),
                'status': 'running',
                'background': 0,
                'args_summary': (args_summary or '')[:500],
                'description': (description[:500] if description else None),
            }
            action_id = dbm.insert_data_dinamically_into_table(cursor, ACTION_LOG_TABLE, data_dict)

            # Opportunistic soft-cap trim; never removes a still-running row.
            cursor.execute(
                f"""
                DELETE FROM {ACTION_LOG_TABLE}
                WHERE status != 'running'
                  AND id NOT IN (
                      SELECT id FROM {ACTION_LOG_TABLE} ORDER BY started_at DESC LIMIT ?
                  )
                """,
                (MAX_LOG_ROWS,),
            )
            conn.commit()
            return action_id
        finally:
            conn.close()

    def finish(self, action_id: int, status: str, error_message: Optional[str] = None) -> None:
        conn = sqlite3.connect(self.db_path)
        try:
            cursor = conn.cursor()
            _ensure_table(cursor)
            cursor.execute(
                f"UPDATE {ACTION_LOG_TABLE} SET ended_at = ?, status = ?, error_message = ? WHERE id = ?",
                (_now(), status, (error_message[:1000] if error_message else None), action_id),
            )
            conn.commit()
        finally:
            conn.close()

    def mark_background(self, action_id: int, child_pid: Optional[int], log_file: Optional[str]) -> None:
        conn = sqlite3.connect(self.db_path)
        try:
            cursor = conn.cursor()
            _ensure_table(cursor)
            cursor.execute(
                f"UPDATE {ACTION_LOG_TABLE} SET background = 1, child_pid = ?, log_file = ? WHERE id = ?",
                (child_pid, log_file, action_id),
            )
            conn.commit()
        finally:
            conn.close()

    def clear(self, older_than_days: Optional[float] = None, statuses: Optional[List[str]] = None) -> int:
        return clear_action_log(self.project_path, older_than_days=older_than_days, statuses=statuses)


def _pid_alive(pid: Optional[int]) -> bool:
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


def _tail(file_path: Optional[str], n_bytes: int = 4096) -> str:
    if not file_path:
        return ''
    try:
        with open(file_path, 'rb') as fh:
            fh.seek(0, os.SEEK_END)
            size = fh.tell()
            fh.seek(max(0, size - n_bytes), os.SEEK_SET)
            return fh.read().decode('utf-8', errors='replace')
    except OSError:
        return ''


def refresh_action_log_statuses(project_path: str) -> None:
    """Infer completion for rows still marked 'running', via PID liveness + log tail."""
    db_path = _db_path(project_path)
    if not os.path.exists(db_path):
        return
    conn = sqlite3.connect(db_path)
    try:
        cursor = conn.cursor()
        _ensure_table(cursor)
        cursor.execute(
            f"SELECT id, process_pid, background, child_pid, log_file "
            f"FROM {ACTION_LOG_TABLE} WHERE status = 'running'"
        )
        rows = cursor.fetchall()
        for row_id, process_pid, background, child_pid, log_file in rows:
            if background:
                if child_pid and _pid_alive(child_pid):
                    continue  # subprocess still running

                if not child_pid:
                    continue  # launch hasn't reported its pid yet

                tail = _tail(log_file)
                if any(marker in tail for marker in _COMPLETION_MARKERS_ERROR):
                    new_status = 'error'
                elif any(marker in tail for marker in _COMPLETION_MARKERS_SUCCESS):
                    new_status = 'success'
                else:
                    new_status = 'completed (unconfirmed)'
                cursor.execute(
                    f"UPDATE {ACTION_LOG_TABLE} SET status = ?, ended_at = ? WHERE id = ?",
                    (new_status, _now(), row_id),
                )
            else:
                if not _pid_alive(process_pid):
                    cursor.execute(
                        f"UPDATE {ACTION_LOG_TABLE} SET status = 'interrupted', ended_at = ? WHERE id = ?",
                        (_now(), row_id),
                    )
        conn.commit()
    finally:
        conn.close()


def get_action_log(project_path: str, limit: int = 50) -> List[Dict[str, Any]]:
    refresh_action_log_statuses(project_path)
    db_path = _db_path(project_path)
    if not os.path.exists(db_path):
        return []
    conn = sqlite3.connect(db_path)
    try:
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()
        _ensure_table(cursor)
        cursor.execute(f"SELECT * FROM {ACTION_LOG_TABLE} ORDER BY started_at DESC LIMIT ?", (limit,))
        return [dict(row) for row in cursor.fetchall()]
    finally:
        conn.close()


def clear_action_log(
    project_path: str,
    older_than_days: Optional[float] = None,
    statuses: Optional[List[str]] = None,
    dry_run: bool = False,
) -> int:
    """Delete finished rows from the action log. Never deletes status='running' rows.

    With dry_run=True, returns the row count that *would* be deleted without
    actually deleting anything (used to show a confirmation prompt beforehand).
    """
    db_path = _db_path(project_path)
    if not os.path.exists(db_path):
        return 0
    conn = sqlite3.connect(db_path)
    try:
        cursor = conn.cursor()
        _ensure_table(cursor)
        clauses = ["status != 'running'"]
        params: List[Any] = []
        if older_than_days is not None:
            cutoff = (datetime.now() - timedelta(days=older_than_days)).strftime("%Y-%m-%d %H:%M:%S")
            clauses.append("started_at < ?")
            params.append(cutoff)
        if statuses:
            placeholders = ', '.join(['?'] * len(statuses))
            clauses.append(f"status IN ({placeholders})")
            params.extend(statuses)
        where_sql = ' AND '.join(clauses)

        cursor.execute(f"SELECT COUNT(*) FROM {ACTION_LOG_TABLE} WHERE {where_sql}", params)
        count = cursor.fetchone()[0]
        if dry_run:
            return count
        cursor.execute(f"DELETE FROM {ACTION_LOG_TABLE} WHERE {where_sql}", params)
        conn.commit()
        return count
    finally:
        conn.close()


# Thread-local call-depth counter. A logged public method commonly calls other
# logged public methods on the same (or another) object internally, e.g.
# ChemSpace.plot_projected_chemical_space() calling self.get_projected_data() or
# similar. Only the outermost call is the "root" user action; nested calls are
# implementation detail and must not create their own log rows / description
# prompts, or a single user action floods the log with duplicates.
_call_depth = threading.local()


def _current_depth() -> int:
    return getattr(_call_depth, 'depth', 0)


def log_action(func):
    """Decorator: records a start/finish row in the project's action log around func.

    Callers may pass a reserved ``log_description`` kwarg to any decorated method to
    attach a free-text description to its log row, e.g.
    ``moldock_obj.dock_table(log_description="testing new receptor")``. It is stripped
    before the real method is called, so it never reaches ``func``.

    If this call happens while another logged method is already running on the call
    stack (e.g. a logged method calling another logged method internally), it is
    treated as a nested implementation detail: no log row is created for it, only
    the outermost ("root") call is recorded.
    """
    @functools.wraps(func)
    def wrapper(self, *args, **kwargs):
        description = kwargs.pop(LOG_DESCRIPTION_KWARG, None)

        if _current_depth() > 0:
            _call_depth.depth = _current_depth() + 1
            try:
                return func(self, *args, **kwargs)
            finally:
                _call_depth.depth -= 1

        # From here down this is the root call: any logged method invoked by func()
        # below must see depth > 0 and skip its own logging, regardless of which of
        # the early-return branches (no project_path, logger.start() failure) is hit.
        _call_depth.depth = 1
        try:
            project_path = getattr(self, 'path', None)
            if not project_path:
                return func(self, *args, **kwargs)

            logger = ActionLogger(project_path)
            source = 'gui' if 'streamlit' in sys.modules else 'api'

            # Only prompt for a plain, interactive terminal API call that didn't already
            # supply log_description. Never prompt for GUI calls or non-interactive/
            # background processes (no stdin to read, would hang).
            if description is None and source == 'api' and sys.stdin.isatty():
                try:
                    description = input(
                        f"📝 Describe this action ({type(self).__name__}.{func.__name__}) "
                        f"[optional, Enter to skip]: "
                    ).strip() or None
                except (EOFError, KeyboardInterrupt):
                    print()
                    description = None

            try:
                action_id = logger.start(
                    type(self).__name__, func.__name__, source, repr(kwargs), description=description,
                )
            except Exception:
                # A logging failure (e.g. locked db) must never break the actual action.
                return func(self, *args, **kwargs)

            try:
                result = func(self, *args, **kwargs)
            except Exception as e:
                try:
                    logger.finish(action_id, status='error', error_message=str(e))
                except Exception:
                    pass
                raise
            else:
                try:
                    if isinstance(result, dict) and result.get('background'):
                        logger.mark_background(action_id, result.get('pid'), result.get('log_file'))
                    else:
                        logger.finish(action_id, status='success')
                except Exception:
                    pass
                return result
        finally:
            _call_depth.depth = 0
    return wrapper


def log_all_public_methods(cls):
    """Class decorator: wraps every public (non-underscore) method with log_action."""
    for name, attr in list(vars(cls).items()):
        if name.startswith('_'):
            continue
        if isinstance(attr, (staticmethod, classmethod)):
            continue
        if callable(attr):
            setattr(cls, name, log_action(attr))
    return cls
