"""
Run an analysis as a batch job on the WRDS Cloud from your local machine, and
copy the results back to your machine.

Each analysis describes its job with a BatchJob and calls add_batch_commands()
to give its main script these subcommands:

    setup        one-time checks: remote python/numpy/pandas/wrds, ~/.pgpass (and Dropbox if configured)
    submit       copy the code to WRDS and qsub the job (--wait: wait, then fetch the results)
    status       qstat
    logs JOB_ID  print the job's stdout/stderr
    fetch JOB_ID copy the job's output folder from WRDS scratch to your machine
    upload JOB_ID  copy the output to Dropbox via rclone (only if dropbox_dir is set)
    disconnect   close the shared SSH connection

The job itself runs `python3 <entry_script> run [args]` on a compute node through
the generic wrapper replication/wrds_job.sh, with its output folder in $OUT_DIR.

About SSH: WRDS needs your password plus Duo two-factor approval. One shared
("master") SSH connection is opened; you log in once and every later ssh/scp
call reuses it until control_persist has passed with no activity.
"""
import argparse
import re
import shlex
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Sequence

PROJECT_ROOT = Path(__file__).resolve().parent.parent
WRDS_HOST = "wrds-cloud.wharton.upenn.edu"
CONTROL_PATH = Path.home() / ".ssh" / "wrds-cm-%C"


@dataclass
class BatchJob:
    """Everything that differs between analyses' batch jobs."""
    name: str                          # Grid Engine job name (-N) and scratch subfolder
    entry_script: str                  # repo-relative path of the script the job runs
    local_output_dir: Path             # `fetch` copies results to <local_output_dir>/<job_id>/
    extra_files: Sequence[str] = ()    # repo-relative files besides the shared code (e.g. blacklists)
    job_command: Sequence[str] = ("run",)   # subcommand of entry_script the job executes
    cores: Optional[int] = None        # -pe onenode N (None = WRDS default, 2)
    mem_per_core: Optional[str] = None  # -l m_mem_free=8G (None = WRDS default, 8G)
    notify_email: Optional[str] = None  # Grid Engine begin/end/abort emails
    dropbox_dir: Optional[str] = None   # e.g. "/wrds_batch_jobs/partsip"; None = no Dropbox upload
    rclone_remote: str = "dropbox"      # name of the rclone remote on WRDS
    remote_dir: str = "wrds_batch_jobs"  # repo copy under your WRDS home directory
    wrds_host: str = WRDS_HOST
    job_script: str = "replication/wrds_job.sh"
    control_persist: str = "8h"
    poll_seconds: int = 120
    _username: Optional[str] = field(default=None, repr=False)

    # ------------------------------------------------------------------ files
    def files_to_sync(self) -> List[str]:
        """Shared code (replication/*.py, the wrapper, wrds_creds.py) + this analysis' files.
        Paths are repo-relative; the same layout is recreated under remote_dir."""
        shared = sorted(str(p.relative_to(PROJECT_ROOT)) for p in (PROJECT_ROOT / "replication").glob("*.py"))
        files = shared + [self.job_script, "wrds_creds.py", self.entry_script, *self.extra_files]
        return list(dict.fromkeys(files))   # de-duplicate, keep order

    # ------------------------------------------------------------------- SSH
    def username(self) -> str:
        if self._username is None:
            sys.path.insert(0, str(PROJECT_ROOT))
            from wrds_creds import read_pgpass
            self._username = read_pgpass()[0]
        return self._username

    def target(self) -> str:
        return f"{self.username()}@{self.wrds_host}"

    def ssh_opts(self) -> list:
        return ["-o", "ControlMaster=auto",
                "-o", f"ControlPath={CONTROL_PATH}",
                "-o", f"ControlPersist={self.control_persist}",
                "-o", "ServerAliveInterval=60"]

    def ensure_connected(self):
        """Open the shared SSH connection if it isn't already up (prompts for password + Duo)."""
        alive = subprocess.run(["ssh", *self.ssh_opts(), "-O", "check", self.target()],
                               capture_output=True).returncode == 0
        if alive:
            return
        CONTROL_PATH.parent.mkdir(mode=0o700, exist_ok=True)
        print(f"Connecting to {self.wrds_host}: enter your WRDS password and approve the Duo prompt.")
        if subprocess.run(["ssh", *self.ssh_opts(), "-N", "-f", self.target()]).returncode != 0:
            sys.exit("Couldn't open an SSH connection to WRDS.")

    def remote(self, cmd: str, *, capture=True, check=True, stdin=None, input_text=None):
        """Run `cmd` in a login shell on the WRDS head node (so qsub etc. are on PATH)."""
        self.ensure_connected()
        args = ["ssh", *self.ssh_opts(), self.target(), f"bash -lc {shlex.quote(cmd)}"]
        r = subprocess.run(args, text=True, capture_output=capture, stdin=stdin, input=input_text)
        if r.returncode == 255:
            sys.exit(f"SSH connection error while running: {cmd}\n{r.stderr or ''}")
        if check and r.returncode != 0:
            sys.exit(f"Remote command failed (exit {r.returncode}): {cmd}\n{r.stdout or ''}{r.stderr or ''}")
        return r

    def remote_test(self, path: str) -> bool:
        # Double quotes (not shlex.quote) so $HOME etc. expand on the remote side.
        return self.remote(f'test -e "{path}"', check=False).returncode == 0

    def scratch_root(self) -> str:
        """/scratch/<group>/<user>/<job name> on WRDS."""
        r = self.remote(f'echo "/scratch/$(id -gn)/$USER/{self.name}"')
        return _last_line(r.stdout)

    def job_out_dir(self, job_id: str) -> str:
        return f"{self.scratch_root()}/{job_id}"

    def script_cmd(self) -> str:
        return f"python {self.entry_script}"

    # -------------------------------------------------------------- commands
    def setup(self):
        self.ensure_connected()
        steps = 3 if self.dropbox_dir else 2

        print(f"\n[1/{steps}] Checking Python on WRDS...")
        r = self.remote('python3 -c "import sys, numpy, pandas, wrds; from importlib.metadata import version; '
                        'print(\'python\', sys.version.split()[0], \'| numpy\', numpy.__version__, '
                        '\'| pandas\', pandas.__version__, \'| wrds\', version(\'wrds\'))"', check=False)
        print("  " + (_last_line(r.stdout) if r.returncode == 0
                      else f"problem importing numpy/pandas/wrds:\n{r.stderr}"))

        print(f"\n[2/{steps}] Checking ~/.pgpass on WRDS (the job uses it to log in to the database)...")
        if self.remote_test("$HOME/.pgpass"):
            print("  found.")
        else:
            local_pgpass = Path.home() / ".pgpass"
            answer = input(f"  Not found. Copy your local {local_pgpass} to WRDS (readable only by you)? [y/N] ")
            if answer.strip().lower() == "y":
                with open(local_pgpass) as f:
                    # set -C (noclobber): never overwrite an existing .pgpass
                    self.remote("set -C && umask 077 && cat > ~/.pgpass && chmod 600 ~/.pgpass", stdin=f)
                print("  copied.")
            else:
                print("  skipped. The job can't connect to the database without it.")

        if self.dropbox_dir:
            print(f"\n[3/{steps}] Linking Dropbox (rclone remote '{self.rclone_remote}:' on WRDS)...")
            remotes = self.remote("rclone listremotes 2>/dev/null", check=False).stdout.split()
            if f"{self.rclone_remote}:" not in remotes:
                self._link_dropbox()
            dest = self._dropbox_dest()
            self.remote(f"rclone mkdir {shlex.quote(dest)}")
            self.remote(f"rclone lsd {shlex.quote(dest)}")
            print(f"  OK, results will also go to Dropbox folder {self.dropbox_dir}")
        print(f"\nResults are fetched to {self.local_output_dir}/<job_id>/")

    def sync_files(self):
        files = self.files_to_sync()
        dirs = sorted({str(Path(self.remote_dir, f).parent) for f in files} | {f"{self.remote_dir}/logs"})
        self.remote("mkdir -p " + " ".join(shlex.quote(d) for d in dirs))
        for rel in files:
            src = PROJECT_ROOT / rel
            if not src.exists():
                sys.exit(f"Missing local file: {src}")
            subprocess.run(["scp", "-q", *self.ssh_opts(), str(src),
                            f"{self.target()}:{self.remote_dir}/{rel}"], check=True)
        print(f"Copied {len(files)} files to ~/{self.remote_dir} on WRDS.")

    def submit(self, job_args: Sequence[str], wait: bool = False):
        self.sync_files()
        env = f"ENTRY_SCRIPT={self.entry_script},OUT_ROOT={self.scratch_root()}"
        qsub = ["qsub", "-N", self.name, "-v", env]
        if self.cores:
            qsub += ["-pe", "onenode", str(self.cores)]
        if self.mem_per_core:
            qsub += ["-l", f"m_mem_free={self.mem_per_core}"]
        if self.notify_email:
            qsub += ["-m", "abe", "-M", self.notify_email]
        qsub += [self.job_script, *self.job_command, *job_args]

        r = self.remote(f"cd {shlex.quote(self.remote_dir)} && {shlex.join(qsub)}")
        m = re.search(r"Your job (\d+)", r.stdout)
        if not m:
            sys.exit(f"Couldn't find the job ID in the qsub output:\n{r.stdout}{r.stderr}")
        job_id = m.group(1)
        print(f"Submitted job {job_id}.")
        print(f"  status: {self.script_cmd()} status")
        print(f"  logs:   {self.script_cmd()} logs {job_id}")
        print(f"  fetch:  {self.script_cmd()} fetch {job_id}")
        if wait:
            self.wait_for_job(job_id)
            self.fetch(job_id)
            if self.dropbox_dir:
                self.upload(job_id)
        return job_id

    def wait_for_job(self, job_id: str):
        print(f"Waiting for job {job_id} (checking every {self.poll_seconds}s). "
              "Ctrl-C stops waiting; the job keeps running.")
        while self.remote(f"qstat -j {job_id}", check=False).returncode == 0:
            time.sleep(self.poll_seconds)
        print("Job has left the queue.")

    def _wait_for_success(self, out_dir: str) -> bool:
        # WRDS says result files can take up to ~2 minutes to appear after a job ends.
        for attempt in range(6):
            if self.remote_test(f"{out_dir}/_SUCCESS"):
                return True
            if attempt < 5:
                time.sleep(30)
        return False

    def fetch(self, job_id: str, force: bool = False) -> Path:
        """Copy <scratch>/<job_id>/ to <local_output_dir>/<job_id>/."""
        out_dir = self.job_out_dir(job_id)
        if not self.remote_test(out_dir):
            sys.exit(f"{out_dir} doesn't exist on WRDS (scratch is wiped after 7 days).")
        if not self._wait_for_success(out_dir):
            msg = f"No _SUCCESS marker in {out_dir}; the job may have failed or still be running."
            if not force:
                sys.exit(f"{msg}\nCheck: {self.script_cmd()} logs {job_id}  (or fetch --force for partial output)")
            print(f"{msg} Fetching anyway (--force).")
        dest_parent = Path(self.local_output_dir)
        dest_parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(["scp", "-r", "-q", *self.ssh_opts(), f"{self.target()}:{out_dir}",
                        str(dest_parent)], check=True)
        dest = dest_parent / job_id
        print(f"Copied results of job {job_id} to {dest}")
        return dest

    def _dropbox_dest(self) -> str:
        return f"{self.rclone_remote}:{self.dropbox_dir.lstrip('/')}"

    def _link_dropbox(self):
        """Authorize Dropbox in a browser on this machine, then save the token in WRDS's rclone config."""
        if shutil.which("rclone") is None:
            sys.exit("  rclone isn't installed on this machine. Run `brew install rclone`, then rerun setup.")
        print("  Opening your browser: sign in to Dropbox and click Allow...")
        r = subprocess.run(["rclone", "authorize", "dropbox"], stdout=subprocess.PIPE, text=True)
        m = re.search(r"--->\s*(\{.*?\})\s*<---", r.stdout, re.DOTALL) or re.search(r"^(\{.*\})\s*$", r.stdout, re.M)
        if r.returncode != 0 or not m:
            sys.exit(f"  rclone authorize didn't return a token:\n{r.stdout}")
        # The token goes over stdin (never on a command line) into a config file only you can read.
        name = shlex.quote(self.rclone_remote)
        self.remote(f"""
set -e
umask 077
conf="$(rclone config file | tail -1)"
mkdir -p "$(dirname "$conf")"
if grep -q "^\\[{self.rclone_remote}\\]" "$conf" 2>/dev/null; then
    echo "a [{self.rclone_remote}] section already exists in $conf; not overwriting" >&2; exit 3
fi
{{ printf '\\n[%s]\\ntype = dropbox\\ntoken = ' {name}; cat; printf '\\n'; }} >> "$conf"
chmod 600 "$conf"
""", input_text=m.group(1).strip())
        print("  token saved on WRDS.")

    def upload(self, job_id: str, force: bool = False):
        """Copy a finished job's output to Dropbox (runs rclone on the login node)."""
        if not self.dropbox_dir:
            sys.exit("dropbox_dir is not set for this job; use `fetch` to copy results to your machine.")
        out_dir = self.job_out_dir(job_id)
        if not self._wait_for_success(out_dir):
            sys.exit(f"No _SUCCESS marker in {out_dir}; the job may have failed. "
                     f"Check: {self.script_cmd()} logs {job_id}")
        if not force and self.remote_test(f"{out_dir}/_UPLOADED"):
            print(f"Job {job_id} already uploaded to Dropbox:{self.dropbox_dir} (use --force to upload again).")
            return
        q_out = shlex.quote(out_dir)
        self.remote(f"rclone copy {q_out} {shlex.quote(self._dropbox_dest())}/{job_id} --exclude '_*' -P "
                    f"&& touch {q_out}/_UPLOADED", capture=False)
        print(f"Uploaded results of job {job_id} to Dropbox:{self.dropbox_dir}/{job_id}")

    def status(self):
        r = self.remote("qstat", check=False)
        print(r.stdout.strip() or "No jobs queued or running.")

    def logs(self, job_id: str):
        log = f"{self.remote_dir}/logs/{self.name}.o{job_id}"
        r = self.remote(f"cat {shlex.quote(log)}", check=False)
        print(r.stdout if r.returncode == 0 else f"No log yet at ~/{log}")

    def disconnect(self):
        subprocess.run(["ssh", *self.ssh_opts(), "-O", "exit", self.target()])


def add_batch_commands(subparsers, job: BatchJob):
    """Add setup/submit/status/logs/fetch/upload/disconnect subcommands to a script's parser.
    Each subcommand sets args.func; call args.func(args)."""
    run_cmd = " ".join(job.job_command)
    subparsers.add_parser("setup", help="one-time WRDS checks (python packages, ~/.pgpass)") \
        .set_defaults(func=lambda a: job.setup())

    s = subparsers.add_parser("submit", help="copy the code to WRDS and qsub the job",
                              description=f"Arguments after '--' are passed to `{run_cmd}` on WRDS.")
    s.add_argument("--wait", action="store_true", help="wait for the job to finish, then fetch the results")
    s.add_argument("job_args", nargs=argparse.REMAINDER, help=f"arguments for `{run_cmd}`, after '--'")
    s.set_defaults(func=lambda a: job.submit(a.job_args[1:] if a.job_args[:1] == ["--"] else a.job_args,
                                             wait=a.wait))

    subparsers.add_parser("status", help="show your queued and running jobs") \
        .set_defaults(func=lambda a: job.status())

    s = subparsers.add_parser("logs", help="print a job's log")
    s.add_argument("job_id")
    s.set_defaults(func=lambda a: job.logs(a.job_id))

    s = subparsers.add_parser("fetch", help=f"copy a job's output to {job.local_output_dir}/<job_id>/")
    s.add_argument("job_id")
    s.add_argument("--force", action="store_true", help="fetch even without the _SUCCESS marker")
    s.set_defaults(func=lambda a: job.fetch(a.job_id, force=a.force))

    s = subparsers.add_parser("upload", help="copy a job's output to Dropbox (only if dropbox_dir is set)")
    s.add_argument("job_id")
    s.add_argument("--force", action="store_true", help="upload even if already uploaded")
    s.set_defaults(func=lambda a: job.upload(a.job_id, force=a.force))

    subparsers.add_parser("disconnect", help="close the shared SSH connection") \
        .set_defaults(func=lambda a: job.disconnect())


def _last_line(text: str) -> str:
    lines = [ln for ln in (text or "").splitlines() if ln.strip()]
    return lines[-1].strip() if lines else ""
