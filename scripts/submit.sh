#!/usr/bin/env bash
# CoralShot submission dispatcher.  Usage:  bash scripts/submit.sh <stage> [dependency_jobid]
# Submits ONE stage, sized to its list, with account/QoS via CLI (since #SBATCH can't expand vars).
# Stops at go/no-go gates by design — you review, then submit the next stage.
#
#   Recommended order (gates in CAPS):
#     s0host s0bracken           # one-time builds (after stage_dbs.sh)
#     -- GNG-0: bash scripts/gng0_smoketest.sbatch  (MUST pass) --
#     s1 ; s2(dep s1) ; s3(dep s2)
#     -- GNG-A: cat results/03_host/depth.d/*.row > results/03_host/microbial_depth.tsv ; review it;
#               set humann=yes for top pills; HARD-FAIL <10k --
#     s4(dep s3) ; s6(dep s3) ; s5(dep s4)        # profiling (s5 only after humann= set)
#     -- delete $DB/kraken2_pluspf ; bash scripts/stage_gtdbtk.sh --
#     s9squeegee(dep s6)                                  # optional: control-free contaminant cross-check
#     -- decontam: Rscript scripts/decontam.R + Rscript scripts/decontam_lda_carryover.R ; GNG-B --
#     s7(dep s3) ; s8(dep s7) ; s8b(dep s8) ; s8c(dep s8b)
#     -- GNG-C: count GTDB-unplaced MAGs -> fire gLMs or skip --
#     s9genomad(dep s7) ; s9annot(dep s8b)
#     stats: Rscript scripts/sap.R
set -euo pipefail
source /scratch/sr7729/CORAL_P/env.sh
cd "$BASE"
STAGE="${1:?usage: submit.sh <stage> [dep_jobid]}"; DEP="${2:-}"
dep_arg(){ [ -n "$DEP" ] && echo "--dependency=afterok:$DEP" || true; }
nmanifest=$(($(wc -l < manifest.tsv) - 1)); nfiles=$(wc -l < filelist.txt)
nunits=$( [ -f asm_units.tsv ] && echo $(($(wc -l < asm_units.tsv) - 1)) || echo 0 )   # assembly-only; not shipped in this MetaPhlAn-only repo
assert(){ [ "$1" -eq "$2" ] || { echo "ASSERT FAIL: $3 ($1 != $2)"; exit 1; }; }

case "$STAGE" in
  s0host)    sbatch -A "$ACCT" --qos="$CPU_QOS" scripts/s0_build_host.sbatch ;;
  s0bracken) sbatch -A "$ACCT" --qos="$CPU_QOS" scripts/s0_bracken.sbatch ;;
  s1)  assert "$nfiles" $((2*nmanifest)) "filelist=2*manifest (paired)"; sbatch -A "$ACCT" --qos="$CPU_QOS" --array=1-${nfiles}%24 $(dep_arg) scripts/s1_qc.sbatch ;;
  s2)  assert "$nfiles" $((2*nmanifest)) "filelist=2*manifest (paired)"; sbatch -A "$ACCT" --qos="$CPU_QOS" --array=1-${nmanifest}%24 $(dep_arg) scripts/s2_fastp.sbatch ;;
  s3)  sbatch -A "$ACCT" --qos="$CPU_QOS" --array=1-${nmanifest}%12 $(dep_arg) scripts/s3_host.sbatch ;;
  s4)  sbatch -A "$ACCT" --qos="$CPU_QOS" --array=1-${nmanifest}%24 $(dep_arg) scripts/s4_metaphlan.sbatch ;;
  s5)  sbatch -A "$ACCT" --qos="$CPU_QOS" --array=1-${nmanifest}%6  $(dep_arg) scripts/s5_humann.sbatch ;;
  s6)  sbatch -A "$ACCT" --qos="$CPU_QOS" --array=1-${nmanifest}%6  $(dep_arg) scripts/s6_kraken_bracken.sbatch ;;
  s7)  sbatch -A "$ACCT" --qos="$CPU_QOS" --array=1-${nunits}%8 $(dep_arg) scripts/s7_assembly.sbatch ;;
  s8)  sbatch -A "$ACCT" --qos="$CPU_QOS" --array=1-${nunits}%6 $(dep_arg) scripts/s8_binning.sbatch ;;
  s8b) sbatch -A "$ACCT" --qos="$CPU_QOS" $(dep_arg) scripts/s8b_mags.sbatch ;;
  s8c) sbatch -A "$ACCT" --qos="$CPU_QOS" $(dep_arg) scripts/s8c_mag_control_screen.sbatch ;;
  s9squeegee) sbatch -A "$ACCT" --qos="$CPU_QOS" $(dep_arg) scripts/s9b_squeegee.sbatch ;;   # optional control-free contaminant cross-check
  s9genomad) sbatch -A "$ACCT" --qos="$CPU_QOS" --array=1-${nunits}%8 $(dep_arg) scripts/s9_genomad.sbatch ;;
  s9annot)
      MAGN=$(ls results/08_mags_drep/dereplicated_genomes/*.fa 2>/dev/null | wc -l)
      [ "$MAGN" -gt 0 ] || { echo "no dRep MAGs yet"; exit 1; }
      sbatch -A "$ACCT" --qos="$CPU_QOS" --array=1-"$MAGN"%8 $(dep_arg) scripts/s9_annot_mags.sbatch ;;
  dark_orfs)  sbatch -A "$ACCT" --qos="$CPU_QOS" $(dep_arg) scripts/lm/make_dark_orfs.sbatch ;;   # ESM-2 input: homology baseline -> dark ORFs
  evo2|genomeocean|esm2)
      sbatch -A "$ACCT" --qos="$GPU_QOS" $(dep_arg) scripts/lm/${STAGE}.sbatch ;;
  *) echo "unknown stage: $STAGE"; exit 1 ;;
esac
