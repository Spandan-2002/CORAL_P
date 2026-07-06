#!/usr/bin/env bash
# CORAL_P submission dispatcher.  Usage:  bash scripts/submit.sh <stage> [dependency_jobid]
# Submits ONE stage, sized to its list, with account/QoS via CLI (since #SBATCH can't expand vars).
# Stops between stages by design — you review the output, then submit the next.
#
# This repo is the read -> relative-abundance subset, so only these stages ship:
#   s0host                 # one-time bowtie2 host index build (after stage_dbs.sh)
#   s1 ; s2(dep s1) ; s3(dep s2)          # FastQC -> fastp trim -> host removal
#   s4(dep s3)                            # MetaPhlAn 4 per library
#   s4b(dep s4)                           # merge -> results/04_mpa/merged_sgb.tsv
# Then figures + phyloseq/stats layer are run directly (see README / RUN_AND_VERIFY.txt),
# not through this dispatcher.
set -euo pipefail
source /scratch/sr7729/CORAL_P/env.sh
cd "$BASE"
STAGE="${1:?usage: submit.sh <stage> [dep_jobid]}"; DEP="${2:-}"
dep_arg(){ [ -n "$DEP" ] && echo "--dependency=afterok:$DEP" || true; }
nmanifest=$(($(wc -l < manifest.tsv) - 1)); nfiles=$(wc -l < filelist.txt)
assert(){ [ "$1" -eq "$2" ] || { echo "ASSERT FAIL: $3 ($1 != $2)"; exit 1; }; }

case "$STAGE" in
  s0host) sbatch -A "$ACCT" --qos="$CPU_QOS" scripts/s0_build_host.sbatch ;;
  s1)  assert "$nfiles" $((2*nmanifest)) "filelist=2*manifest (paired)"; sbatch -A "$ACCT" --qos="$CPU_QOS" --array=1-${nfiles}%24 $(dep_arg) scripts/s1_qc.sbatch ;;
  s2)  assert "$nfiles" $((2*nmanifest)) "filelist=2*manifest (paired)"; sbatch -A "$ACCT" --qos="$CPU_QOS" --array=1-${nmanifest}%24 $(dep_arg) scripts/s2_fastp.sbatch ;;
  s3)  sbatch -A "$ACCT" --qos="$CPU_QOS" --array=1-${nmanifest}%12 $(dep_arg) scripts/s3_host.sbatch ;;
  s4)  sbatch -A "$ACCT" --qos="$CPU_QOS" --array=1-${nmanifest}%24 $(dep_arg) scripts/s4_metaphlan.sbatch ;;
  s4b) sbatch -A "$ACCT" --qos="$CPU_QOS" $(dep_arg) scripts/s4b_merge_metaphlan.sbatch ;;
  *) echo "unknown stage: $STAGE  (this repo ships: s0host s1 s2 s3 s4 s4b)"; exit 1 ;;
esac
