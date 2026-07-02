#!/usr/bin/env bash
#SBATCH -p cs -c 16 --mem=24G -t 1:00:00 --array=1-53%24 -J cs_fastp -o logs/%x_%A_%a.out
# Adapter + quality + poly-G/X trimming. Keeps orphan (singleton) reads — every read counts at low biomass.
# Array size = 53 libraries; authoritatively (re)set by scripts/submit.sh (1-$nmanifest).
set -euo pipefail
source /scratch/sr7729/CORAL_P/env.sh
module load bioinformatics/20260224
IFS=$'\t' read -r SID SUBJ ARM PART TYPE I7 I5 HOPF HU R1 R2 \
  < <(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$BASE/manifest.tsv")
[ -z "${SID:-}" ] && { echo "FATAL: empty manifest row"; exit 1; }
[[ "$SID" =~ ^[PC][0-9]{2}[MB]_[CSF]$ ]] || { echo "FATAL: bad SID '$SID'"; exit 1; }
mkdir -p "$BASE/work/trim" "$BASE/results/02_fastp"
fastp -i "$R1" -I "$R2" \
  -o "$BASE/work/trim/${SID}_R1.fq.gz" -O "$BASE/work/trim/${SID}_R2.fq.gz" \
  --unpaired1 "$BASE/work/trim/${SID}_U.fq.gz" --unpaired2 "$BASE/work/trim/${SID}_U.fq.gz" \
  --detect_adapter_for_pe \
  --adapter_sequence       CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
  --adapter_sequence_r2    CTGTCTCTTATACACATCTGACGCTGCCGACGA \
  --trim_poly_g --poly_g_min_len 10 --trim_poly_x \
  --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 \
  --length_required 60 --thread 16 \
  -j "$BASE/results/02_fastp/${SID}.json" -h "$BASE/results/02_fastp/${SID}.html"
