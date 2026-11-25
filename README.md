# Hasbrouck_HRVAR

Replicating (and extending) Joel Hasbrouck's *Price Discovery in High Resolution* in Python,
using TAQ data on WRDS.

- [Paper](https://pages.stern.nyu.edu/~jhasbrou/Research/HRVAR/HRVAR08.pdf)
- [Hasbrouck's documentation](https://pages.stern.nyu.edu/~jhasbrou/Research/HRVAR/HRVARindex.html)

## Repository layout

| Path | Contents |
|---|---|
| [replication/](replication/) | Python replication of Hasbrouck's analyses: shared code plus one folder per analysis. See [replication/README.md](replication/README.md). |
| [HRVAR_V01/](HRVAR_V01/) | Hasbrouck's original SAS and MATLAB code and his output listings (reference only, not run). |
| [reference_materials/](reference_materials/) | The paper (`HRVAR08-2.pdf`) and Hasbrouck's code documentation. |
| [blacklist_files/](blacklist_files/) | Dates whose TAQ tables are missing on WRDS (`missing_days_cqm.txt`, etc.) and the script that generates them. These dates are skipped. |
| [wrds_creds.py](wrds_creds.py) | Reads your WRDS username and password from `~/.pgpass`. |
| [requirements.txt](requirements.txt) | Python packages for the local environment. |

## Getting started

### 1. Python environment

From the repo root:

```bash
python3 -m venv env
source env/bin/activate
pip install -r requirements.txt
```

The analyses need `numpy`, `pandas` and `wrds`. Plotting needs `matplotlib`.
Batch jobs also use the `ssh` and `scp` commands, which come with macOS and Linux.

### 2. WRDS account and `~/.pgpass`

You need a WRDS account with access to TAQ (millisecond and microsecond TAQ, the `taqm_YYYY`
schemas). The code never asks for your password interactively. It reads your credentials from
a `.pgpass` file in your home directory, the same file the `wrds` library uses.

Create `~/.pgpass` with one line, replacing the username and password with your own:

```
wrds-pgdata.wharton.upenn.edu:9737:wrds:your_wrds_username:your_wrds_password
```

Then make it readable only by you. PostgreSQL ignores the file otherwise:

```bash
chmod 600 ~/.pgpass
```

The format is `host:port:database:username:password`. If the file has several lines,
[wrds_creds.py](wrds_creds.py) uses the last one.

A batch job running on WRDS also needs a `.pgpass` in your WRDS home directory. The `setup`
command checks for it and offers to copy your local one (see
[replication/README.md](replication/README.md#batch-jobs-on-the-wrds-cloud)).

**Never commit `.pgpass` or put your password in the code.**

### 3. Run an analysis

Every command is run from the repo root. The general instructions are in
[replication/README.md](replication/README.md), and each analysis has its own README:

- [Analysis 1: participant vs. SIP NBBOs](replication/analysis1/README.md)
- [Analysis 2](replication/analysis2/README.md) (not yet written)
- [Analysis 3](replication/analysis3/README.md) (not yet written)

## Data and outputs

TAQ data is proprietary. Results that contain or are derived from it go to
`replication/**/output/` and job logs go to `logs/`. Both are in `.gitignore`, so keep them
out of the repository.
