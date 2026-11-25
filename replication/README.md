# Replication of Hasbrouck, *Price Discovery in High Resolution*

This folder holds the Python replication of the analyses in Hasbrouck's paper. Hasbrouck built
his data in SAS and estimated his models in MATLAB; his code is in [../HRVAR_V01/](../HRVAR_V01/).
Here both steps happen in Python, run as batch jobs on the WRDS Cloud.

Set up your Python environment and `~/.pgpass` first, as described in the
[project README](../README.md#getting-started).

## Layout

```
replication/
├── helpers.py          shared utilities
├── taq_queries.py      TAQ queries on WRDS
├── nbbo.py             NBBO construction
├── mvar.py             PDL VECM estimation and impulse responses
├── random_walk.py      random-walk decomposition and information shares
├── wrds_batch.py       submit / monitor / fetch WRDS batch jobs
├── wrds_job.sh         Grid Engine wrapper that every batch job runs
├── analysis1/          participant vs. SIP NBBOs (Table 5A, Figure 1)
├── analysis2/          (not yet written)
└── analysis3/          (not yet written)
```

### Shared scripts

Every analysis can use these modules. Analysis-specific code goes in its `analysisN/` folder.

| Script | What it does | Hasbrouck's equivalent |
|---|---|---|
| [helpers.py](helpers.py) | Trading-day iteration that skips the dates in `blacklist_files/` (`load_blacklist`, `iter_days`), pulling data in chunks of days (`bulk_data_pull`), splitting large DataFrames (`chunk_dataframe`), and a run log that prints and writes to a file (`Listing`). | MATLAB `diary` |
| [taq_queries.py](taq_queries.py) | Queries TAQ on WRDS's PostgreSQL server. `cqm_quotes(day, symbols)` returns one day of quotes (`taqm_YYYY.cqm_YYYYMMDD`) with SIP (`time_m`) and participant (`part_time`) timestamps as integer microseconds. The connection opens on the first query and uses `~/.pgpass`. | SAS data steps in `taqExtract01.sas` |
| [nbbo.py](nbbo.py) | `build_nbbo(quotes, time_col)` builds the national best bid and offer from exchange quotes. It follows Hasbrouck's `%buildNBBO` SAS macro, including its quirks. `get_price(time_us, price, second_factor)` puts a price series on a time grid at a given resolution. | `%buildNBBO`, `getPrice.m` |
| [mvar.py](mvar.py) | `PDLVECM` estimates a vector error-correction model whose lags are grouped into polynomial distributed lag (PDL) terms, at any time resolution. It also computes the VMA sums and impulse responses. `pdl_polys(second_factor)` builds Hasbrouck's lag structure. | `mClasses/` (MVAR classes) |
| [random_walk.py](random_walk.py) | `RandomWalkDecomp` computes the random-walk variance and the Cholesky information-share bounds over all orderings of the prices, individually or by group. | `randomWalkDecomp` |
| [wrds_batch.py](wrds_batch.py) | `BatchJob` and `add_batch_commands`, which give each analysis's main script the `setup`, `submit`, `status`, `logs`, `fetch` and `disconnect` commands. See [Batch jobs](#batch-jobs-on-the-wrds-cloud). | |
| [wrds_job.sh](wrds_job.sh) | The Grid Engine script that `qsub` runs on a WRDS compute node. It runs the analysis's `run` command, writes output to scratch, and marks success with a `_SUCCESS` file. You never run it directly. | |

### A note on the time grid

At 10 µs a trading day has about 2.25 billion grid points, too many to build the regression
data in memory. Each regressor changes value only at a few shifted price-change times, so
`mvar.py` builds the moment matrix by summing over the intervals between those times. The
result is exactly what the full grid gives, with no approximation, in a fraction of the
memory. Hasbrouck's MATLAB code reaches the same goal with sparse matrices.

## Running an analysis

Each analysis has one main script in its folder, run from the **repo root**:

```bash
python replication/analysisN/<main script>.py <command> [options]
```

Every main script has these commands:

| Command | What it does |
|---|---|
| `setup` | One-time checks on WRDS: that Python, numpy, pandas and wrds work there, and that `~/.pgpass` exists. |
| `submit [--wait] [-- run options]` | Copies the code to WRDS and submits the batch job. `--wait` waits for the job to finish, then fetches the results. Options after `--` go to the analysis's `run` command. |
| `status` | Shows your queued and running jobs (`qstat`). |
| `logs JOB_ID` | Prints the job's log. |
| `fetch JOB_ID [--force]` | Copies the job's output from WRDS to `replication/analysisN/output/JOB_ID/`. `--force` fetches a job that did not finish cleanly. |
| `disconnect` | Closes the shared SSH connection to WRDS. |
| `run [options]` | Runs the analysis in the current session. A batch job runs this on WRDS; you can also run it yourself (see below). |

Run `python replication/analysisN/<main script>.py <command> --help` to list a command's options.

A typical session:

```bash
python replication/analysis1/PartSIP_main.py setup          # first time only
python replication/analysis1/PartSIP_main.py submit --wait  # submit, wait, fetch
```

If you don't use `--wait`, the job keeps running after you close the terminal. Check on it
and collect the results later. `submit` prints the job ID:

```bash
python replication/analysis1/PartSIP_main.py status
python replication/analysis1/PartSIP_main.py logs 12345678
python replication/analysis1/PartSIP_main.py fetch 12345678
```

## Batch jobs on the WRDS Cloud

The analyses read large amounts of quote data, so they run on WRDS's servers rather than
your computer. Each job has three parts:

1. **Your computer.** `submit` logs in to `wrds-cloud.wharton.upenn.edu` over SSH. It copies
   the shared scripts, `wrds_creds.py`, the analysis's main script and any files the analysis
   lists (such as its blacklist) to `~/wrds_batch_jobs/` on WRDS. Then it runs `qsub`.
2. **A WRDS compute node.** [wrds_job.sh](wrds_job.sh) runs `python3 <main script> run ...`.
   The job queries TAQ, estimates the models and writes to
   `/scratch/<group>/<user>/<job name>/<job_id>/`. When it finishes cleanly it writes `_SUCCESS`.
   The log goes to `~/wrds_batch_jobs/logs/`.
3. **Back on your computer.** `fetch` (or `submit --wait`) copies that folder to
   `replication/analysisN/output/<job_id>/`.

Things to know:

- **Logging in.** WRDS requires your password and a Duo approval. The first command of a
  session asks for both, then opens one shared SSH connection that later commands reuse for
  up to 8 hours of inactivity. You don't need SSH keys. `disconnect` closes the connection.
- **`.pgpass` on WRDS.** The job connects to the TAQ database from the compute node, so it
  needs `~/.pgpass` in your **WRDS** home directory as well as on your computer. `setup`
  checks for it. If it's missing, `setup` offers to copy your local file with permission 600,
  and it never overwrites an existing one.
- **Scratch is temporary.** WRDS deletes scratch files after about 7 days, so fetch results soon.
- **Results stay on your machine.** Nothing is uploaded elsewhere. `wrds_batch.py` can also
  copy results to Dropbox via rclone, but this is off: `dropbox_dir=None` in each analysis's
  `BatchJob`.
- **Resources.** Jobs use WRDS's defaults: 2 cores and 8 GB per core. Change `cores` or
  `mem_per_core` in the analysis's `BatchJob` (WRDS allows at most 8 cores and 48 GB per job).
- **Missing dates.** Dates in `blacklist_files/missing_days_*.txt` have no TAQ table on WRDS
  and are skipped.
- **Failures.** If one symbol-date fails, the error goes to the log and the job moves on to
  the next. The job then exits with an error, so use `fetch JOB_ID --force` to get the
  results that did finish.

### Running without a batch job

`run` works anywhere that can reach the WRDS database with a `~/.pgpass`. That includes an
interactive WRDS session (`ssh` to WRDS, then `qrsh`) or your own computer for small tests.
Its output goes to `--output-dir`, which defaults to `replication/analysisN/output/local/`.

### Adding a new analysis

Create `replication/analysisN/<name>_main.py` and:

1. define a `BatchJob` (name, `entry_script`, `local_output_dir`, any `extra_files`);
2. add a `run` subcommand that writes its results to `$OUT_DIR`;
3. call `add_batch_commands(subparsers, JOB)` to get the batch commands.

[analysis1/PartSIP_main.py](analysis1/PartSIP_main.py) is the template. Shared modules in
`replication/` are copied to WRDS automatically.

## Analyses

### Analysis 1: participant vs. SIP NBBOs

The information shares of NBBOs built from participant (exchange) timestamps vs. SIP
timestamps, at resolutions from 1 second to 10 microseconds (Table 5A and Figure 1).

- Main script: [analysis1/PartSIP_main.py](analysis1/PartSIP_main.py)
- Instructions: [analysis1/README.md](analysis1/README.md)

### Analysis 2

*Not yet written.* See [analysis2/README.md](analysis2/README.md).

### Analysis 3

*Not yet written.* See [analysis3/README.md](analysis3/README.md).
