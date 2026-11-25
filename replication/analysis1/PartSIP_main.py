#!/usr/bin/env python3
"""
Analysis 1 (PartSIP): price discovery between NBBOs built from participant
(exchange) timestamps and NBBOs built from SIP timestamps.

Replicates Hasbrouck, "Price Discovery in High Resolution", Table 5A and Figure 1:
    data       taqExtract01.sas  (%buildNBBOPartSIP: NBBOs from cqm quotes)
    estimation MVARiPartSIP01.m  (PDL VECM at 1s ... 10us, random-walk decomposition, IRFs)
    summary    summaryPartSIP01.m (sd of the random walk and grouped information-share bounds)

The four prices, in Hasbrouck's order, are nbbPart, nboPart, nbbSIP, nboSIP. Information-
share bounds are grouped as Part = {nbbPart, nboPart} and SIP = {nbbSIP, nboSIP}.

From your machine (repo root), as a WRDS Cloud batch job:
    python replication/analysis1/PartSIP_main.py setup            # one-time checks
    python replication/analysis1/PartSIP_main.py submit --wait    # Hasbrouck's sample, then fetch results
    python replication/analysis1/PartSIP_main.py submit -- --start 2016-10-03 --end 2016-11-11
    python replication/analysis1/PartSIP_main.py status | logs JOB_ID | fetch JOB_ID
    python replication/analysis1/PartSIP_main.py plot replication/analysis1/output/JOB_ID

Directly (on WRDS inside `qrsh`, or locally with a ~/.pgpass):
    python replication/analysis1/PartSIP_main.py run --symbols IBM --second-factors 1 10 --output-dir test_out
Re-estimate from NBBOs saved by an earlier run (no WRDS connection needed):
    python replication/analysis1/PartSIP_main.py run --prices-dir replication/analysis1/output/JOB_ID/prices

Output folder ($OUT_DIR on WRDS; fetched to replication/analysis1/output/<job_id>/):
    listing.txt                     run log (Hasbrouck's MATLAB diary)
    summary_table5A.csv             symbol, date, resolution, sdw, Part/SIP information-share bounds
    info_share_bounds.csv           bounds for every price and group, all resolutions
    prices/SYM_YYYYMMDD.csv.gz      the four NBBO change series (input to the estimation)
    results/SYM_YYYYMMDD_sfN.json   estimates, eCov, VMA sums, bounds and IRFs for one resolution
"""
import argparse
import json
import os
import sys
import time
import traceback
from datetime import date
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = ANALYSIS_DIR.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from replication.wrds_batch import BatchJob, add_batch_commands  # noqa: E402

# ---------------------------------------------------------------------------
# Analysis settings (MVARiPartSIP01.m). Dates/symbols/resolutions can be
# overridden on the command line.
# ---------------------------------------------------------------------------
START_DATE = date(2016, 10, 3)
END_DATE = date(2016, 10, 3)
SYMBOLS = ("IBM", "NVDA")
SECOND_FACTORS = (1, 10, 100, 1000, 10000, 100000)   # 1s, 100ms, 10ms, 1ms, 100us, 10us
LAG_PDL_MAX = 10 ** 6
FIRST_TIME_S = 9 * 3600 + 45 * 60                    # estimation range 9:45-16:00
LAST_TIME_S = 16 * 3600
PRICE_NAMES = ("nbbPart", "nboPart", "nbbSIP", "nboSIP")
IS_GROUPS = {"Part": ("nbbPart", "nboPart"), "SIP": ("nbbSIP", "nboSIP")}
IRF_SHOCK, IRF_RESPONSE = "nbbPart", "nbbSIP"        # Figure 1


def irf_horizon(second_factor: int) -> int:
    return min(500 * second_factor, 500 * 10 ** 4)


# ---------------------------------------------------------------------------
# WRDS batch job (shared machinery in replication/wrds_batch.py)
# ---------------------------------------------------------------------------
JOB = BatchJob(
    name="PartSIP",
    entry_script="replication/analysis1/PartSIP_main.py",
    local_output_dir=ANALYSIS_DIR / "output",
    extra_files=("blacklist_files/missing_days_cqm.txt",),
    cores=None,            # WRDS default (2 cores); the job is single-threaded
    mem_per_core=None,     # WRDS default (8G per core)
    notify_email=None,
    dropbox_dir=None,      # results are copied to your machine with `fetch`
)


# ---------------------------------------------------------------------------
# Data: participant and SIP NBBOs (taqExtract01.sas, %buildNBBOPartSIP)
# ---------------------------------------------------------------------------
def build_partsip_series(quotes):
    """The four NBBO change series for one symbol-date, as {name: DataFrame[time_us, price]}."""
    from replication.nbbo import build_nbbo
    nbb_part, nbo_part = build_nbbo(quotes, "part_time_us")
    nbb_sip, nbo_sip = build_nbbo(quotes, "time_m_us")
    return {"nbbPart": nbb_part, "nboPart": nbo_part, "nbbSIP": nbb_sip, "nboSIP": nbo_sip}


def save_series(series, path: Path):
    import pandas as pd
    frames = [df.assign(series=name) for name, df in series.items()]
    out = pd.concat(frames, ignore_index=True)[["series", "time_us", "price"]]
    out["time_us"] = out["time_us"].round().astype("Int64")
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, index=False, compression="gzip")


def load_series(path: Path):
    import pandas as pd
    df = pd.read_csv(path)
    return {name: df.loc[df["series"] == name, ["time_us", "price"]].astype("float64").reset_index(drop=True)
            for name in PRICE_NAMES}


def prices_file(prices_dir: Path, symbol: str, day: date) -> Path:
    return Path(prices_dir) / f"{symbol}_{day:%Y%m%d}.csv.gz"


# ---------------------------------------------------------------------------
# Estimation for one symbol, date and resolution (MVARiPartSIP01.m loop body)
# ---------------------------------------------------------------------------
def analyze(series, symbol: str, day: date, second_factor: int, log, with_irf: bool = True):
    import numpy as np
    from replication.helpers import format_elapsed
    from replication.mvar import PDLVECM, PriceSeries, pdl_polys
    from replication.nbbo import get_price
    from replication.random_walk import RandomWalkDecomp

    sf = second_factor
    log(f"\n{symbol} {day:%Y%m%d} secondFactor={sf} ({log.stamp()})")
    polys = pdl_polys(sf, LAG_PDL_MAX)
    for p in polys:
        log(f"name={p.name}; deg=0; n={p.n}; kOffset={p.k_offset}; vNames={p.var_name}.")

    prices = []
    for name in PRICE_NAMES:
        s = series[name]
        t2, v = get_price(s["time_us"].to_numpy(), s["price"].to_numpy(), sf)
        log(f"{day:%Y%m%d} {symbol} {name:8s} {len(s):7d} NBBO changes "
            f"({len(t2)} with nonzero price changes at this resolution and 9am-5pm)")
        prices.append(PriceSeries(name, t2, v))

    t0 = time.perf_counter()
    model = PDLVECM(prices, polys, sf, ecm=True, intercept=False)
    model.set_valid(FIRST_TIME_S, LAST_TIME_S, log)
    model.build_zpz(log=log)
    model.estimate()
    log(f"VAR/VECM estimates for secondFactor={sf}  ({model.n_segments:,} segments, "
        f"{model.last - model.first + 1:,} grid points)")
    log(model.estimates_table().to_string(float_format=lambda x: f"{x:.4g}"))
    log("eCov")
    log(_frame(model.ecov, model.y_names).to_string(float_format=lambda x: f"{x:.4g}"))
    log("eCorr")
    log(_frame(model.ecorr, model.y_names).to_string(float_format=lambda x: f"{x:.4g}"))
    log(f"Estimation finished for {symbol} {day:%Y%m%d}. secondFactor={sf}. "
        f"Elapsed time is {format_elapsed(time.perf_counter() - t0)}.")

    rw = RandomWalkDecomp(model.theta_sum()[0], model.ecov, model.price_names, sf)
    rw.grouped_bounds(IS_GROUPS)
    log(rw.report())

    if with_irf:
        t0 = time.perf_counter()
        n_ahead = irf_horizon(sf)
        model.irf_packed(n_ahead, log=log)
        log(f"IRF calculation finished ({n_ahead:,} steps). "
            f"Elapsed time is {format_elapsed(time.perf_counter() - t0)}.")

    result = {
        "symbol": symbol, "date": f"{day:%Y-%m-%d}", "second_factor": sf, "resolution_s": 1 / sf,
        "polys": [{"name": p.name, "n": p.n, "k_offset": p.k_offset} for p in polys],
        "price_names": model.price_names,
        "valid_range": [model.first, model.last], "n_obs": model.last - model.first + 1,
        "n_changes": {p.name: [int(len(p.t)), int(np.sum((p.t >= model.first) & (p.t <= model.last)))]
                      for p in prices},
        "ecm_means": model.ecm_means.tolist(),
        "x_names": model.x_names, "y_names": model.y_names,
        "b": model.b.tolist(), "t": model.tb.tolist(),
        "ecov": model.ecov.tolist(), "ecorr": model.ecorr.tolist(),
        "vma1": rw.vma1.tolist(), "var_w": rw.var_w, "var_w_annual": rw.var_w_annual,
        "sdw": float(np.sqrt(rw.var_w_annual)),
        "is_bounds": {k: [v["min"], v["max"]] for k, v in rw.is_bounds.iterrows()},
        "is_group_bounds": {k: [v["min"], v["max"]] for k, v in rw.is_group_bounds.iterrows()},
    }
    if with_irf:
        result["irft"] = model.irft.tolist()
        result["irf"] = model.irf.tolist()     # [price, horizon, shock]
    return result


def _frame(m, names):
    import pandas as pd
    return pd.DataFrame(m, index=names, columns=names)


def summary_rows(result):
    """One row of Table 5A, plus long-format bounds for every price and group."""
    g = result["is_group_bounds"]
    row = {"symbol": result["symbol"], "date": result["date"], "resolution_s": result["resolution_s"],
           "second_factor": result["second_factor"], "sdw": result["sdw"]}
    for label in IS_GROUPS:
        row[f"{label}_min"], row[f"{label}_max"] = g[label]
    long = []
    for kind, bounds in (("price", result["is_bounds"]), ("group", g)):
        for name, (lo, hi) in bounds.items():
            long.append({"symbol": result["symbol"], "date": result["date"],
                         "second_factor": result["second_factor"], "kind": kind, "name": name,
                         "min": lo, "max": hi})
    return row, long


# ---------------------------------------------------------------------------
# run: the job itself
# ---------------------------------------------------------------------------
def cmd_run(args):
    import pandas as pd
    from replication.helpers import Listing, format_elapsed, iter_days, load_blacklist

    out = args.output_dir.expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    log = Listing(out / "listing.txt")
    symbols = [s.upper() for s in args.symbols]
    sfs = sorted(set(args.second_factors))
    log(f"PartSIP run starting at {log.stamp()}")
    log(f"dates {args.start} -> {args.end} | symbols {' '.join(symbols)} | secondFactors {sfs} | "
        f"irf={'no' if args.no_irf else 'yes'} | prices from "
        f"{args.prices_dir if args.prices_dir else 'WRDS cqm'}")
    log(f"output: {out}")

    table, long_rows, failures = [], [], []
    t_all = time.perf_counter()
    for day in iter_days(args.start, args.end, load_blacklist("cqm")):
        # ---- data for all symbols on this day
        try:
            day_series = get_day_series(day, symbols, args.prices_dir, out, log)
        except Exception:
            log(f"FAILED to get data for {day}:\n{traceback.format_exc()}")
            failures.append((str(day), "*", "data"))
            continue
        # ---- estimation
        for symbol, series in day_series.items():
            for sf in sfs:
                try:
                    result = analyze(series, symbol, day, sf, log, with_irf=not args.no_irf)
                except Exception:
                    log(f"FAILED {symbol} {day} secondFactor={sf}:\n{traceback.format_exc()}")
                    failures.append((str(day), symbol, sf))
                    continue
                res_file = out / "results" / f"{symbol}_{day:%Y%m%d}_sf{sf}.json"
                res_file.parent.mkdir(exist_ok=True)
                res_file.write_text(json.dumps(result))
                row, long = summary_rows(result)
                table.append(row)
                long_rows += long
            # rewrite the summaries after every symbol-date so partial runs keep their results
            if table:
                pd.DataFrame(table).to_csv(out / "summary_table5A.csv", index=False)
                pd.DataFrame(long_rows).to_csv(out / "info_share_bounds.csv", index=False)

    if table:
        log("\nSummary (Table 5A layout: sd of random walk, annualized; information-share bounds)")
        log(pd.DataFrame(table).to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    log(f"\nPartSIP run finished at {log.stamp()}. Elapsed time is {format_elapsed(time.perf_counter() - t_all)}.")
    if failures:
        log(f"{len(failures)} failure(s): {failures}")
        sys.exit(1)
    if not table:
        log("No results were produced (no trading days in range?).")
        sys.exit(1)


def get_day_series(day, symbols, prices_dir, out, log):
    """{symbol: series} for one day, from saved NBBO files or from WRDS cqm quotes."""
    result = {}
    if prices_dir:
        for symbol in symbols:
            f = prices_file(prices_dir, symbol, day)
            if f.exists():
                result[symbol] = load_series(f)
            else:
                log(f"{day} {symbol}: no saved prices at {f}; skipping")
        return result

    from replication.taq_queries import cqm_exists, cqm_quotes
    if not cqm_exists(day):
        log(f"{day}: no cqm table on WRDS (add it to blacklist_files/missing_days_cqm.txt); skipping")
        return result
    t0 = time.perf_counter()
    quotes = cqm_quotes(day, symbols)
    log(f"\n{day}: {len(quotes):,} cqm quotes for {', '.join(symbols)} "
        f"({time.perf_counter() - t0:.0f}s)")
    for symbol in symbols:
        q = quotes[quotes["sym_root"] == symbol]
        if q.empty:
            log(f"{day} {symbol}: no quotes; skipping")
            continue
        n_missing_part = int(q["part_time_us"].isna().sum())
        series = build_partsip_series(q)
        log(f"{day} {symbol}: {len(q):,} quotes ({n_missing_part:,} without participant time); NBBO changes: "
            + ", ".join(f"{k}={len(v):,}" for k, v in series.items()))
        save_series(series, prices_file(out / "prices", symbol, day))
        result[symbol] = series
    return result


# ---------------------------------------------------------------------------
# plot: Figure 1 (summaryPartSIP01.m). Needs matplotlib.
# ---------------------------------------------------------------------------
def cmd_plot(args):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        sys.exit("plot needs matplotlib (pip install matplotlib)")
    results_dir = Path(args.results_dir) / "results"
    files = sorted(results_dir.glob("*.json"))
    if not files:
        sys.exit(f"no result files in {results_dir}")
    runs = {}
    for f in files:
        r = json.loads(f.read_text())
        if "irf" in r:
            runs.setdefault((r["symbol"], r["date"]), []).append(r)
    if not runs:
        sys.exit("no IRFs in the result files (was the run made with --no-irf?)")

    for (symbol, day), rs in sorted(runs.items()):
        fig, ax = plt.subplots(figsize=(6, 4))
        for r in sorted(rs, key=lambda r: r["second_factor"]):
            i = r["price_names"].index(args.response)
            s = r["price_names"].index(args.shock)
            x = [h / r["second_factor"] for h in r["irft"][1:]]
            y = [row[s] for row in r["irf"][i][1:]]
            ax.semilogx(x, y, label=_resolution_label(r["second_factor"]))
        ax.set_title(f"{symbol} {day}: initial impulse in {args.shock}; IRF for {args.response}", fontsize=9)
        ax.set_xlabel("Time")
        ax.set_ylabel("$/share")
        ticks = [10.0 ** k for k in range(-5, 3)]
        ax.set_xticks(ticks, ["10 µs", "100 µs", "1 ms", "10 ms", "100 ms", "1 sec", "10 sec", "100 sec"])
        ax.set_xlim(1e-6, 1300)
        ax.grid(True, which="major")
        ax.legend(title="resolution", fontsize=8)
        fig.tight_layout()
        path = Path(args.results_dir) / f"figure1_{symbol}_{day.replace('-', '')}.png"
        fig.savefig(path, dpi=300)
        plt.close(fig)
        print(f"wrote {path}")


def _resolution_label(sf: int) -> str:
    return {1: "1 sec", 10: "100 ms", 100: "10 ms", 1000: "1 ms", 10000: "100 µs", 100000: "10 µs"}.get(sf, f"1/{sf} s")


# ---------------------------------------------------------------------------
# command line
# ---------------------------------------------------------------------------
def main(argv=None):
    p = argparse.ArgumentParser(description="Hasbrouck HRVAR analysis 1: participant vs SIP NBBO information shares",
                                formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    sub = p.add_subparsers(dest="command", required=True)

    r = sub.add_parser("run", help="run the analysis here (this is what the WRDS job executes)")
    r.add_argument("--start", type=date.fromisoformat, default=START_DATE, help=f"first date (default {START_DATE})")
    r.add_argument("--end", type=date.fromisoformat, default=END_DATE, help=f"last date (default {END_DATE})")
    r.add_argument("--symbols", nargs="+", default=list(SYMBOLS), help=f"default {' '.join(SYMBOLS)}")
    r.add_argument("--second-factors", nargs="+", type=int, default=list(SECOND_FACTORS),
                   help="time units per second, e.g. 1 (1s) ... 100000 (10us); default all six")
    r.add_argument("--output-dir", type=Path, default=Path(os.environ.get("OUT_DIR", ANALYSIS_DIR / "output" / "local")),
                   help="where to write results (default $OUT_DIR on WRDS, else replication/analysis1/output/local)")
    r.add_argument("--prices-dir", type=Path, default=None,
                   help="re-use NBBO files saved by an earlier run instead of querying WRDS")
    r.add_argument("--no-irf", action="store_true", help="skip impulse responses (Figure 1); much faster at 100us/10us")
    r.set_defaults(func=cmd_run)

    pl = sub.add_parser("plot", help="draw Figure 1 from a fetched output folder (needs matplotlib)")
    pl.add_argument("results_dir", help="e.g. replication/analysis1/output/<job_id>")
    pl.add_argument("--shock", default=IRF_SHOCK)
    pl.add_argument("--response", default=IRF_RESPONSE)
    pl.set_defaults(func=cmd_plot)

    add_batch_commands(sub, JOB)

    args = p.parse_args(argv)
    if getattr(args, "end", None) and args.end < args.start:
        sys.exit(f"--end ({args.end}) is before --start ({args.start})")
    args.func(args)


if __name__ == "__main__":
    main()
