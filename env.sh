#!/usr/bin/env bash
# CoralShot — shared environment preamble. SOURCE this first in every job/script.
#   source /scratch/sr7729/CORAL_P/env.sh
# Encodes the two cluster-specific fixes the design pass verified live on Torch:
#   (1) /home is at 100% inodes -> redirect ALL tool caches to /scratch.
#   (2) The /share/apps/bioinformatics bundle is a Singularity wrapper whose
#       bindpath excludes /tmp and /dev/shm -> temp MUST live on a BOUND fs
#       (node-local /state/partition1, 810 GB) or large samtools sorts corrupt.

export BASE=/scratch/sr7729/CORAL_P
export DB=$BASE/db

# --- SLURM accounting (user choice: general allocation) -------------------
# NOTE: #SBATCH lines cannot expand env vars; submit.sh passes these via the
# sbatch CLI (-A / --qos). Kept here for interactive use and reference.
export ACCT=torch_pr_113_general
export CPU_QOS=cpu48          # short/most stages; use cpu168 for humann/assembly
export CPU_QOS_LONG=cpu168    # >48 h stages (HUMAnN, assembly)
export GPU_QOS=gpu48          # genomic-LM jobs (verify exact GPU QoS at GNG-0)

# --- caches OFF /home (100% inodes) ---------------------------------------
export XDG_CACHE_HOME=$BASE/cache
export HF_HOME=$BASE/cache/hf
export CONDA_PKGS_DIRS=$BASE/cache/condapkgs
export PIP_CACHE_DIR=$BASE/cache/pip
export CONDA_ENVS_DIRS=/scratch/sr7729/conda_envs   # existing envs dir; make `-n` resolve
export MPLCONFIGDIR=$BASE/cache/mpl
export NUMBA_CACHE_DIR=$BASE/cache/numba
export R_LIBS_USER=$BASE/cache/Rlibs
export FONTCONFIG_PATH=$BASE/cache/fontconfig
export LMOD_CACHE_DIR=$BASE/cache/lmod
export NXF_HOME=$BASE/cache/nextflow
mkdir -p "$XDG_CACHE_HOME" "$HF_HOME" "$CONDA_PKGS_DIRS" "$PIP_CACHE_DIR" "$MPLCONFIGDIR" \
         "$NUMBA_CACHE_DIR" "$R_LIBS_USER" "$FONTCONFIG_PATH" "$LMOD_CACHE_DIR" "$NXF_HOME"

# --- node-local BOUND temp (verified: /state/partition1 = 810 GB, in bindpath) ---
# This cluster has no $SLURM_TMPDIR; /tmp and /dev/shm are NOT bound in the bundle container.
if [ -n "${SLURM_JOB_ID:-}" ] && [ -d /state/partition1 ]; then
  export TMPDIR=/state/partition1/$SLURM_JOB_ID
else
  export TMPDIR=$BASE/cache/tmp/${SLURM_JOB_ID:-$$}   # login-node fallback: PID-unique, never a shared 'login' dir
fi
mkdir -p "$TMPDIR"
# cleanup must survive scancel/OOM — but ONLY under SLURM. An unconditional EXIT trap on a login shell would
# wipe a shared temp dir under concurrency and silently clobber any pre-existing EXIT trap when env.sh is sourced.
if [ -n "${SLURM_JOB_ID:-}" ]; then
  trap 'rm -rf "$TMPDIR"' EXIT TERM INT HUP
fi

# --- compute nodes have no internet: force offline model/DB use ------------
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1

# --- conda --------------------------------------------------------------------
source /scratch/sr7729/miniforge3/etc/profile.d/conda.sh

# convenience: activate a coralshot env by short name
csenv() { conda activate /scratch/sr7729/conda_envs/"$1"; }
