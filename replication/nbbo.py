"""
NBBO construction from exchange quotes, and price discretization.

build_nbbo  is a vectorized port of the %buildNBBO macro in Hasbrouck's
            taqExtract01.sas (the hash table of each exchange's current quote).
get_price   is a port of getPrice.m (single symbol / single date branch).

Both are shared by the PartSIP (participant vs SIP timestamps) and LexExl
(listing vs other exchanges) analyses.
"""
from typing import Tuple

import numpy as np
import pandas as pd

US_PER_SECOND = 1_000_000


def build_nbbo(quotes: pd.DataFrame, time_col: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """National best bid and offer from one symbol-date's exchange quotes.

    quotes   : one symbol on one date, in original TAQ order, with columns
               ex, bid, bidsiz, ask, asksiz and the timestamp column `time_col`
               (microseconds after midnight, may contain NaN).
    time_col : the timestamp that sequences the quotes, e.g. 'part_time_us'
               (participant) or 'time_m_us' (SIP).

    Returns (nbb, nbo): one row per record that changes the NBB (NBO), with
    columns time_us and price. This is SAS's temp.<dsOut>_nbb / _nbo.

    Replicates the SAS data step exactly, including:
      * records are stably sorted by `time_col`; missing times sort first;
      * an exchange enters the table only with a positive price AND size, but once
        in the table, any nonzero price replaces its quote (the SAS tests
        `bidsize`/`asksize`, which aren't cqm columns, so size is ignored there);
        a zero price removes that side;
      * if no exchange has a valid bid (ask), the NBB (NBO) keeps its last value;
      * the first NBB/NBO of the day is not output (no previous value to compare).
    """
    t = quotes[time_col].to_numpy(dtype="float64")
    order = np.argsort(np.where(np.isnan(t), -np.inf, t), kind="stable")
    t = t[order]
    ex_codes, _ = pd.factorize(quotes["ex"].to_numpy()[order])
    bid = quotes["bid"].to_numpy(dtype="float64")[order]
    bidsiz = quotes["bidsiz"].to_numpy(dtype="float64")[order]
    ask = quotes["ask"].to_numpy(dtype="float64")[order]
    asksiz = quotes["asksiz"].to_numpy(dtype="float64")[order]

    with np.errstate(invalid="ignore"):
        bid_pos, ask_pos = bid > 0, ask > 0
        new_bid_ok = bid_pos & (bidsiz > 0)
        new_ask_ok = ask_pos & (asksiz > 0)
    # Was the exchange already in the hash table before this record?
    added = (new_bid_ok | new_ask_ok).astype(np.int64)
    in_table = (pd.Series(added).groupby(ex_codes).cumsum().to_numpy() - added) > 0

    eff_bid = np.where(bid_pos & (in_table | new_bid_ok), bid, np.nan)
    eff_ask = np.where(ask_pos & (in_table | new_ask_ok), ask, np.nan)

    nbb_raw = _best_across_exchanges(ex_codes, eff_bid, np.fmax)
    nbo_raw = _best_across_exchanges(ex_codes, eff_ask, np.fmin)
    return _changes(t, nbb_raw), _changes(t, nbo_raw)


def _best_across_exchanges(ex_codes: np.ndarray, eff: np.ndarray, combine) -> np.ndarray:
    """At each record, combine (max for bids, min for asks) each exchange's current quote.
    An exchange's current quote is the value set by its most recent record (NaN = no quote)."""
    n = len(eff)
    best = np.full(n, np.nan)
    rows = np.arange(n)
    for e in np.unique(ex_codes):
        last = np.maximum.accumulate(np.where(ex_codes == e, rows, -1))
        cur = np.where(last >= 0, eff[np.maximum(last, 0)], np.nan)
        best = combine(best, cur)   # fmax/fmin ignore NaN
    return best


def _changes(t: np.ndarray, raw: np.ndarray) -> pd.DataFrame:
    """Carry the last valid value forward, then keep records where it changes."""
    level = pd.Series(raw).ffill().to_numpy()
    prev = np.concatenate(([np.nan], level[:-1]))
    keep = ~np.isnan(level) & ~np.isnan(prev) & (level != prev)
    return pd.DataFrame({"time_us": t[keep], "price": level[keep]})


def get_price(time_us: np.ndarray, price: np.ndarray, second_factor: int,
              start_s: float = 9 * 3600, end_s: float = 17 * 3600) -> Tuple[np.ndarray, np.ndarray]:
    """Discretize a price series to a grid of 1/second_factor seconds (getPrice.m).

    1. time2 = floor(time * second_factor): the grid index (intervals after midnight);
    2. keep the last record within each interval;
    3. keep start_s < time < end_s (9:00-17:00 by default);
    4. keep only records where the price changes.

    Returns (index, price) as int64 / float64 arrays."""
    if US_PER_SECOND % second_factor:
        raise ValueError(f"second_factor={second_factor} must divide 1,000,000")
    time_us = np.asarray(time_us, dtype="float64")
    price = np.asarray(price, dtype="float64")
    ok = ~np.isnan(time_us)
    time_us, price = time_us[ok].astype(np.int64), price[ok]
    if np.any(np.diff(time_us) < 0):
        raise ValueError("time_us must be sorted ascending")

    time2 = time_us // (US_PER_SECOND // second_factor)
    last_in_interval = np.append(np.diff(time2) > 0, True)
    time_us, time2, price = time_us[last_in_interval], time2[last_in_interval], price[last_in_interval]

    hours = (time_us > start_s * US_PER_SECOND) & (time_us < end_s * US_PER_SECOND)
    time2, price = time2[hours], price[hours]

    changed = np.append(True, np.diff(price) != 0) if len(price) else np.zeros(0, bool)
    return time2[changed], price[changed]
