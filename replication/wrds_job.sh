#!/bin/bash
#$ -cwd
#$ -j y
#$ -o logs/
#
# Generic WRDS Cloud (Grid Engine) wrapper shared by every analysis.
#
# replication/wrds_batch.py submits it from ~/wrds_batch_jobs (the synced repo copy) as
#   qsub -N <job name> -v ENTRY_SCRIPT=<script>,OUT_ROOT=<dir> [-pe onenode N] [-l m_mem_free=XG] \
#        replication/wrds_job.sh run [script args]
# and everything after the wrapper's name is passed to the entry script.
#
# The entry script finds its output folder in $OUT_DIR = $OUT_ROOT/$JOB_ID
# (WRDS scratch by default; scratch is wiped after 7 days). Markers written there:
#   _SUCCESS   the Python job finished cleanly
#   _UPLOADED  results were copied to Dropbox (only if a Dropbox folder is configured)
#
# Default resources are 2 cores and 8GB per core (WRDS limits: 8 cores, 48GB per job);
# set cores / mem_per_core in the analysis' BatchJob to change them.

set -uo pipefail

: "${ENTRY_SCRIPT:?ENTRY_SCRIPT is not set (pass it with qsub -v)}"

echo "Job ${JOB_ID} (${JOB_NAME}) starting on $(hostname) at $(date)"
echo "Entry script: ${ENTRY_SCRIPT} $*"

OUT_ROOT="${OUT_ROOT:-/scratch/$(id -gn)/${USER}/${JOB_NAME}}"
export OUT_DIR="${OUT_ROOT}/${JOB_ID}"
mkdir -p "$OUT_DIR"
echo "Output directory: $OUT_DIR"

# Keep numpy's BLAS within the cores Grid Engine gave us.
export OMP_NUM_THREADS="${NSLOTS:-1}" OPENBLAS_NUM_THREADS="${NSLOTS:-1}" MKL_NUM_THREADS="${NSLOTS:-1}"

python3 -u "$ENTRY_SCRIPT" "$@"
status=$?
if (( status != 0 )); then
    echo "${ENTRY_SCRIPT} failed with exit code $status at $(date)"
    exit "$status"
fi
touch "$OUT_DIR/_SUCCESS"

# Compute nodes have no internet access, so results are collected from your
# machine (scp through the login node) rather than pushed from here.
echo "To copy the results to your machine: python ${ENTRY_SCRIPT} fetch ${JOB_ID}"
echo "Job ${JOB_ID} finished at $(date)"
