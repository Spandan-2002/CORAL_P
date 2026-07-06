#!/usr/bin/env bash
#SBATCH -p cs -c 24 --mem=48G -t 3:00:00 --array=1-53%24 -J cs_mpa -o logs/%x_%A_%a.out
# PRIMARY taxonomy: MetaPhlAn 4.1 (SGB relative abundance) -> one profile per library.
# Array size is authoritatively (re)set by scripts/submit.sh (1-$nmanifest); this header is for direct `sbatch`.
set -euo pipefail
source /scratch/sr7729/CORAL_P/env.sh
conda activate /scratch/sr7729/conda_envs/cs_profile
IFS=$'\t' read -r SID REST < <(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$BASE/manifest.tsv")
[ -z "${SID:-}" ] && { echo "FATAL empty row"; exit 1; }
zcat "$BASE/work/clean/${SID}_R1.fq.gz" "$BASE/work/clean/${SID}_R2.fq.gz" > "$TMPDIR/${SID}.fq"
if [ -s "$BASE/work/clean/${SID}_U.fq.gz" ]; then zcat "$BASE/work/clean/${SID}_U.fq.gz" >> "$TMPDIR/${SID}.fq"; fi
metaphlan "$TMPDIR/${SID}.fq" --input_type fastq --nproc 24 \
  --bowtie2db "$DB/metaphlan" --index mpa_vJun23_CHOCOPhlAnSGB_202403 --unclassified_estimation \
  --bowtie2out "$BASE/results/04_mpa/${SID}.bt2.bz2" \
  -o "$BASE/results/04_mpa/${SID}.profile.tsv"
# After the array finishes, merge with scripts/s4b_merge_metaphlan.sbatch:
#   merge_metaphlan_tables.py results/04_mpa/*.profile.tsv > results/04_mpa/merged_sgb.tsv
