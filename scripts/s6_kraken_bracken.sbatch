#!/usr/bin/env bash
#SBATCH -p cs -c 64 --mem=220G -t 3:00:00 --array=1-53%6 -J cs_kraken -o logs/%x_%A_%a.out
# Kraken2 PlusPF (DECONTAM SUBSTRATE) + Bracken. On the 3 TB cl nodes (DB is mem-mapped).
# Array size is authoritatively (re)set by scripts/submit.sh (1-$nmanifest); this header is for direct `sbatch`.
set -euo pipefail
source /scratch/sr7729/CORAL_P/env.sh
module load bioinformatics/20260224
IFS=$'\t' read -r SID SUBJ ARM PART TYPE I7 I5 HOPF HU R1 R2 \
  < <(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$BASE/manifest.tsv")
[ -z "${SID:-}" ] && { echo "FATAL empty row"; exit 1; }
mkdir -p "$BASE/results/06_kraken" "$BASE/results/06_bracken"
kraken2 --db "$DB/kraken2_pluspf" --threads 64 --paired --gzip-compressed \
  --confidence 0.15 --memory-mapping \
  --report "$BASE/results/06_kraken/${SID}.kreport" --output /dev/null \
  "$BASE/work/clean/${SID}_R1.fq.gz" "$BASE/work/clean/${SID}_R2.fq.gz"
conda activate /scratch/sr7729/conda_envs/cs_profile
[ "$TYPE" = "feces" ] && THR=10 || THR=20    # higher threshold for low-biomass pills
bracken -d "$DB/kraken2_pluspf" -i "$BASE/results/06_kraken/${SID}.kreport" -r 150 -l S -t "$THR" \
  -o "$BASE/results/06_bracken/${SID}.S.bracken"
bracken -d "$DB/kraken2_pluspf" -i "$BASE/results/06_kraken/${SID}.kreport" -r 150 -l G -t "$THR" \
  -o "$BASE/results/06_bracken/${SID}.G.bracken"
# After array (run FROM the bracken dir so columns are basenamed SIDs, not path-prefixed):
#   ( cd results/06_bracken && combine_bracken_outputs.py --files *.S.bracken -o merged.S.tsv )
# THEN delete $DB/kraken2_pluspf (160 GB) before staging GTDB-Tk (scripts/stage_gtdbtk.sh).
