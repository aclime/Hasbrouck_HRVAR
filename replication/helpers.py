"""
General helpers shared by all analyses: blacklists, trading-day iteration,
chunked WRDS pulls, and a simple listing (log) writer.
"""
import sys
import time
from datetime import date, datetime, timedelta
from pathlib import Path
from typing import Callable, Iterator, Optional

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
BLACKLIST_DIR = PROJECT_ROOT / "blacklist_files"

# TAQ table family -> blacklist file of dates whose table is missing on WRDS
BLACKLIST_FILES = {
    "cqm": "missing_days_cqm.txt",
    "ctm": "missing_days_ctm.txt",
    "nbbom": "missing_days_nbbom.txt",
}


def load_blacklist(table: str) -> set:
    """Dates (as 'YYYY-MM-DD' strings) with no `table` (cqm, ctm or nbbom) on WRDS."""
    path = BLACKLIST_DIR / BLACKLIST_FILES[table]
    if not path.exists():
        return set()
    with open(path, encoding="utf-8") as f:
        return {line.strip() for line in f if line.strip()}


def iter_days(start: date, end: date, blacklist: Optional[set] = None) -> Iterator[date]:
    """Calendar days from start to end (inclusive), skipping blacklisted ones."""
    blacklist = blacklist or set()
    d = start
    while d <= end:
        if d.strftime("%Y-%m-%d") not in blacklist:
            yield d
        d += timedelta(days=1)


def bulk_data_pull(function: Callable, start_date: date, end_date: date, days_per_chunk: int,
                   *args, **kwargs) -> pd.DataFrame:
    """Run `function(chunk_start, chunk_end, *args, **kwargs)` over consecutive
    calendar-day chunks and concatenate the results.

    Chunks are built from a list of whole days (not start + i*step), so they never
    overlap. Empty chunks (all weekends/holidays) are skipped."""
    all_days = []
    d = start_date
    while d <= end_date:
        all_days.append(d)
        d += timedelta(days=1)

    dfs = []
    t0 = time.perf_counter()
    for i in range(0, len(all_days), days_per_chunk):
        chunk = all_days[i:i + days_per_chunk]
        sub_df = function(chunk[0], chunk[-1], *args, **kwargs)
        if sub_df is None or sub_df.empty:
            continue
        dfs.append(sub_df)
    if not dfs:
        raise RuntimeError(f"no rows returned between {start_date:%Y-%m-%d} and {end_date:%Y-%m-%d}")
    print(f"Bulk df generated in {time.perf_counter() - t0:.2f} seconds", flush=True)
    return pd.concat(dfs, ignore_index=True)


def chunk_dataframe(df: pd.DataFrame, rows_per_chunk: int) -> Iterator[pd.DataFrame]:
    """Yield consecutive row slices of df with at most rows_per_chunk rows each."""
    for i in range(0, len(df), rows_per_chunk):
        yield df.iloc[i:i + rows_per_chunk]


def format_elapsed(seconds: float) -> str:
    """Elapsed time the way Hasbrouck's etDuration.m prints it."""
    if seconds < 10:
        return f"{seconds:.4g} sec"
    if seconds < 3600:
        m, s = divmod(seconds, 60)
        return f"{int(m):02d}:{s:05.2f}"
    h, rem = divmod(int(seconds), 3600)
    m, s = divmod(rem, 60)
    return f"{h:02d}:{m:02d}:{s:02d}"


class Listing:
    """Print to stdout and append the same text to a listing file
    (the equivalent of Hasbrouck's MATLAB `diary`)."""

    def __init__(self, path: Optional[Path] = None):
        self.path = Path(path) if path else None
        if self.path:
            self.path.parent.mkdir(parents=True, exist_ok=True)

    def __call__(self, *parts, sep=" "):
        text = sep.join(str(p) for p in parts)
        print(text, flush=True)
        if self.path:
            with open(self.path, "a", encoding="utf-8") as f:
                f.write(text + "\n")

    def stamp(self) -> str:
        return datetime.now().strftime("%d-%b-%Y %H:%M:%S")


def ensure_project_on_path():
    """Make `import replication...` and `import wrds_creds` work from any script."""
    root = str(PROJECT_ROOT)
    if root not in sys.path:
        sys.path.insert(0, root)
