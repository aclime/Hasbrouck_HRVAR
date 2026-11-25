"""
WRDS TAQ (millisecond) queries shared by all analyses.

The WRDS connection is opened lazily on first use, so importing this module
(e.g. for `--help`, or when re-using saved data) never prompts for a login.
"""
from datetime import date
from typing import Iterable, Optional

import pandas as pd

from replication.helpers import ensure_project_on_path

WRDS_PG_HOST = "wrds-pgdata.wharton.upenn.edu"

_db = None


def get_db():
    """Shared wrds.Connection, using the credentials in ~/.pgpass."""
    global _db
    if _db is None:
        ensure_project_on_path()
        import wrds
        from wrds_creds import read_pgpass
        username, password = read_pgpass()
        _db = wrds.Connection(wrds_host=WRDS_PG_HOST, wrds_username=username,
                              wrds_password=password)
    return _db


def close_db():
    global _db
    if _db is not None:
        _db.close()
        _db = None


def _sql_list(values: Iterable[str]) -> str:
    return ", ".join("'" + str(v).replace("'", "''") + "'" for v in values)


def cqm_quotes(day: date, symbols: Iterable[str], primary_only: bool = True,
               db=None) -> pd.DataFrame:
    """All consolidated quotes (taqm_YYYY.cqm_YYYYMMDD) for `symbols` on `day`.

    This is the input to Hasbrouck's taqExtract01.sas (getTAQq). No filters on
    quote condition or cancellations are applied, matching the SAS program.

    Timestamps come back as integer microseconds after midnight:
        time_m_us     SIP timestamp (time_m)
        part_time_us  participant (exchange) timestamp (part_time); may be null

    Rows are in the original TAQ sequence (time_m, then qu_seqnum). The NBBO
    builder relies on this order to break ties when re-sorting by part_time.

    primary_only=True keeps only the common stock (sym_suffix is null).
    """
    db = db or get_db()
    suffix_filter = "AND t.sym_suffix IS NULL" if primary_only else ""
    sql = f"""
        SELECT t.sym_root,
               t.ex,
               t.bid,
               t.bidsiz,
               t.ask,
               t.asksiz,
               t.qu_seqnum,
               (EXTRACT(EPOCH FROM t.time_m) * 1000000)::bigint    AS time_m_us,
               (EXTRACT(EPOCH FROM t.part_time) * 1000000)::bigint AS part_time_us
        FROM taqm_{day.year}.cqm_{day:%Y%m%d} t
        WHERE t.sym_root IN ({_sql_list(symbols)})
              {suffix_filter}
        ORDER BY t.sym_root, t.time_m, t.qu_seqnum
    """
    df = db.raw_sql(sql)
    for col in ("bid", "bidsiz", "ask", "asksiz", "time_m_us", "part_time_us"):
        df[col] = pd.to_numeric(df[col], errors="coerce").astype("float64")
    df["ex"] = df["ex"].astype(str)
    df["date"] = pd.Timestamp(day)
    return df


def table_exists(schema: str, table: str, db=None) -> bool:
    """True if schema.table exists on WRDS (cheap check used to skip missing days)."""
    db = db or get_db()
    sql = f"""
        SELECT 1 FROM information_schema.tables
        WHERE table_schema = '{schema}' AND table_name = '{table}'
        LIMIT 1
    """
    return not db.raw_sql(sql).empty


def cqm_exists(day: date, db: Optional[object] = None) -> bool:
    return table_exists(f"taqm_{day.year}", f"cqm_{day:%Y%m%d}", db)
