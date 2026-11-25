# Analysis 1: participant vs. SIP NBBOs

This analysis replicates the first analysis in Hasbrouck's *Price Discovery in High Resolution*:
**Table 5A** and **Figure 1**. Quotes in TAQ carry two timestamps: one from the exchange
(participant, `part_time`) and one from the SIP (`time_m`). Building the national best bid and
offer (NBBO) from each gives four prices: `nbbPart`, `nboPart`, `nbbSIP` and `nboSIP`. A VECM
at resolutions from 1 s to 10 µs shows how much price discovery happens in the participant
NBBO vs. the SIP NBBO.

For setup (Python environment, `~/.pgpass`) and how batch jobs work, see the
[project README](../../README.md) and the [replication README](../README.md).

## What it computes

For each symbol and date:

1. **Data** ([nbbo.py](../nbbo.py), [taq_queries.py](../taq_queries.py)). Pull the day's
   `cqm` quotes. Build the NBB and NBO twice: once ordering quotes by participant time and
   once by SIP time.
2. **Estimation** ([mvar.py](../mvar.py)), at each resolution: 1 s, 100 ms, 10 ms, 1 ms,
   100 µs and 10 µs.
   - The model is a VECM in the four price changes, with no intercept.
   - Its lags are grouped into PDL terms covering up to 10 seconds of lags.
   - It has three error-correction terms, `nbbPart − p_k`.
   - The estimation window is 9:45 to 16:00.
3. **Information shares** ([random_walk.py](../random_walk.py)). From the random-walk
   decomposition, compute the Cholesky information-share bounds over all 24 orderings of the
   four prices. They are reported by group: Part = {nbbPart, nboPart}, SIP = {nbbSIP, nboSIP}.
4. **Impulse responses.** Compute the response of `nbbSIP` to a shock in `nbbPart`, which is
   Figure 1.

Hasbrouck's corresponding code is `taqExtract01.sas`, `MVARiPartSIP01.m` and
`summaryPartSIP01.m` in [HRVAR_V01/](../../HRVAR_V01/).

## How to run

Run all commands from the **repo root**.

### First time: setup

```bash
python replication/analysis1/PartSIP_main.py setup
```

This asks for your WRDS password and a Duo approval. It then checks Python on WRDS and that
`~/.pgpass` exists there.

### Submit the batch job

The default run is Hasbrouck's sample: IBM and NVDA on 2016-10-03, at all six resolutions.

```bash
python replication/analysis1/PartSIP_main.py submit --wait
```

`--wait` waits for the job to finish, then copies the results to
`replication/analysis1/output/<job_id>/`. A first test can be smaller. Put the options for
`run` after `--`:

```bash
python replication/analysis1/PartSIP_main.py submit -- --symbols IBM --second-factors 1 10
```

Other samples:

```bash
python replication/analysis1/PartSIP_main.py submit -- --start 2016-10-03 --end 2016-10-31 --symbols IBM NVDA AAPL
python replication/analysis1/PartSIP_main.py submit -- --no-irf     # skip impulse responses (much faster)
```

### Check on a job and collect the results

Without `--wait`, `submit` prints a job ID. Use it with:

```bash
python replication/analysis1/PartSIP_main.py status            # queued / running jobs
python replication/analysis1/PartSIP_main.py logs JOB_ID       # the job's log
python replication/analysis1/PartSIP_main.py fetch JOB_ID      # copy results to your machine
python replication/analysis1/PartSIP_main.py fetch JOB_ID --force   # if some symbol-dates failed
```

### Draw Figure 1

This command needs matplotlib.

```bash
python replication/analysis1/PartSIP_main.py plot replication/analysis1/output/JOB_ID
```

It writes `figure1_*.png` to the same folder. `--shock` and `--response` choose a different
pair of prices.

### Re-estimate without WRDS

Each run saves the four NBBO series under `prices/`. To re-estimate from them on your own
computer, with no WRDS connection:

```bash
python replication/analysis1/PartSIP_main.py run --prices-dir replication/analysis1/output/JOB_ID/prices
```

### `run` options

| Option | Default | Meaning |
|---|---|---|
| `--start`, `--end` | 2016-10-03 | Date range (YYYY-MM-DD); blacklisted dates are skipped. |
| `--symbols` | IBM NVDA | Ticker symbols (`sym_root`). |
| `--second-factors` | 1 10 100 1000 10000 100000 | Resolutions as intervals per second (1 = 1 s, 100000 = 10 µs). |
| `--output-dir` | `$OUT_DIR` in a batch job, otherwise `output/local/` | Where results are written. |
| `--prices-dir` | none | Read saved NBBO series instead of querying TAQ. |
| `--no-irf` | off | Skip impulse responses. |

The defaults are constants at the top of [PartSIP_main.py](PartSIP_main.py).

## Output

`replication/analysis1/output/<job_id>/` contains:

| File | Contents |
|---|---|
| `summary_table5A.csv` | Table 5A: symbol, date, resolution, sd of the random walk (annualized), and the Part and SIP information-share bounds. |
| `info_share_bounds.csv` | Information-share bounds for every price and group at every resolution. |
| `listing.txt` | Full run log (the equivalent of Hasbrouck's MATLAB listings). |
| `prices/SYM_YYYYMMDD.csv.gz` | The four NBBO series for each symbol-date. |
| `results/SYM_YYYYMMDD_sfN.json` | Everything for one symbol-date and resolution: coefficients, residual covariance, VMA sums, bounds and impulse responses. |
| `figure1_*.png` | Figure 1, after running `plot`. |
| `_SUCCESS` | Present if the job finished without errors. |

### Checking against Hasbrouck

Hasbrouck's output for this sample is in
[HRVAR_V01/MatlabNovember2018/PartSIP01/listing/](../../HRVAR_V01/MatlabNovember2018/PartSIP01/listing/).
For example, NVDA at 1 s should give an sd of the random walk of about 12.25, with grouped
bounds of about 0.00028 to 0.99978.

## Run time

The NBBO construction and estimation take a few minutes per symbol-date. Most of the time
goes to the impulse responses at 100 µs and 10 µs (up to 5 million steps, about 2 minutes
at 10 µs). `--no-irf` skips them if you only need Table 5A.

## Differences from Hasbrouck's code

- **His SAS quirks are kept on purpose:**
  - His NBBO step tests `bidsize`/`asksize`, which aren't `cqm` columns. As a result, a
    zero-size quote still updates an exchange already in the NBBO.
  - If no exchange has a valid bid or offer, the NBB/NBO keeps its previous value.
  - The day's first NBBO is not output.
- **Suffixed share classes are dropped** (`sym_suffix IS NULL`). Hasbrouck doesn't filter
  them. This doesn't affect IBM or NVDA.
- **Timestamps are whole microseconds.** At 10 µs, a few intervals can differ by one from
  MATLAB, whose floating-point `fix(time*sf)` sometimes rounds down by one.
- **The regression data is built from intervals between price changes** rather than a full
  grid (see the [replication README](../README.md#a-note-on-the-time-grid)). The result is
  the same.

## Before the first full run

`taq_queries.cqm_quotes` converts `time_m` and `part_time` to microseconds with
`EXTRACT(EPOCH FROM ...)`, which assumes both are PostgreSQL `time` columns. If the first
job fails at the query, check those column types first. The small `--symbols IBM
--second-factors 1 10` test above shows this quickly.
