#!/usr/bin/env bash
#SBATCH -p cs -c 32 --mem=64G -t 4:00:00 --array=1-53%12 -J cs_host -o logs/%x_%A_%a.out
# Human+decoy removal (paired + orphan), no SAM on disk, bound node-local temp.
# Array size = 53 libraries; authoritatively (re)set by scripts/submit.sh (1-$nmanifest).
# %12: each task spins up ~4 bundle containers; bound overlay mounts must stay modest.
set -euo pipefail
source /scratch/sr7729/CORAL_P/env.sh
module load bioinformatics/20260224
IFS=$'\t' read -r SID SUBJ ARM PART TYPE I7 I5 HOPF HU R1 R2 \
  < <(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$BASE/manifest.tsv")
[ -z "${SID:-}" ] && { echo "FATAL: empty manifest row"; exit 1; }
ST="$TMPDIR/$SID"; mkdir -p "$ST"          # BOUND node-local (/state/partition1), NOT /tmp
mkdir -p "$BASE/work/clean" "$BASE/results/03_host/depth.d"
IDX="$DB/host/host_decoy"

# paired: keep read PAIRS where BOTH mates are unmapped (-f 12), drop secondary (-F 256)
bowtie2 -p 32 --very-sensitive -x "$IDX" \
    -1 "$BASE/work/trim/${SID}_R1.fq.gz" -2 "$BASE/work/trim/${SID}_R2.fq.gz" \
    2> "$BASE/results/03_host/${SID}.bt2.log" \
  | samtools view -b -f 12 -F 256 -@ 8 - \
  | samtools sort -n -@ 8 -m 2G -T "$ST/sortP" - \
  | samtools fastq -@ 8 -1 "$BASE/work/clean/${SID}_R1.fq.gz" -2 "$BASE/work/clean/${SID}_R2.fq.gz" \
      -0 /dev/null -s /dev/null -n -

# orphans (single-end): keep unmapped (-f 4)
if [ -s "$BASE/work/trim/${SID}_U.fq.gz" ]; then
  bowtie2 -p 32 --very-sensitive -x "$IDX" -U "$BASE/work/trim/${SID}_U.fq.gz" \
      2>> "$BASE/results/03_host/${SID}.bt2.log" \
    | samtools view -b -f 4 -F 256 -@ 8 - \
    | samtools fastq -@ 8 - | gzip > "$BASE/work/clean/${SID}_U.fq.gz"
fi

# OPERATIONAL biomass = post-host microbial read PAIRS (not file size; pills are host-heavy)
PAIRS=$(( $(zcat "$BASE/work/clean/${SID}_R1.fq.gz" | wc -l) / 4 ))
# per-task row (NFS append is not atomic at %12 concurrency); GNG-A merges depth.d/*.row -> microbial_depth.tsv
echo -e "${SID}\t${TYPE}\t${PART}\t${PAIRS}" > "$BASE/results/03_host/depth.d/${SID}.row"
rm -f "$BASE/work/trim/${SID}_R1.fq.gz" "$BASE/work/trim/${SID}_R2.fq.gz" "$BASE/work/trim/${SID}_U.fq.gz"
