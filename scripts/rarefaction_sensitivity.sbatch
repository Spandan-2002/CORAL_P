#!/usr/bin/env bash
#SBATCH -p cs -c 16 --mem=48G -t 2:00:00 --array=1-30%12 -J cs_rare -o logs/%x_%A_%a.out
# READ-LEVEL RAREFACTION sensitivity check for alpha diversity.
# METHOD (literature-backed): rarefy the INPUT READS, not the MetaPhlAn profile table -- the
# MetaPhlAn developer's own guidance (N. Segata, metaphlan-users forum: "you should rarefy the
# input metagenomes to the desired depth ... use seqtk or other tools"). MetaPhlAn reports relative
# abundance, so there is no count table to rarefy; rarefaction is only needed for the depth-correlated
# statistic (richness; Liu et al. 2022, Genome, PMID 35939836).
# Subsample each clean-10 library's host-removed reads (R1+R2+U, as MetaPhlAn consumes them)
# to a COMMON depth, then re-profile with MetaPhlAn (same DB/params as s4). Downstream we
# recompute Observed/Shannon/InvSimpson on the rarefied profiles and compare to native depth.
# N is set below every clean-10 sample's total combined reads (min pairs 6.90M -> ~14.06M combined).
set -euo pipefail
source /scratch/sr7729/CORAL_P/env.sh
conda activate /scratch/sr7729/conda_envs/cs_profile
N=13000000; SEED=42
SK="$(dirname "$(command -v metaphlan)")/seqkit"
SID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${LIST:-$BASE/rare_samples.txt}")   # LIST overridable for retries
[ -z "${SID:-}" ] && { echo "FATAL empty SID for task $SLURM_ARRAY_TASK_ID"; exit 1; }
mkdir -p "$BASE/results/rarefaction"

comb="$TMPDIR/${SID}.comb.fq"
zcat "$BASE/work/clean/${SID}_R1.fq.gz" "$BASE/work/clean/${SID}_R2.fq.gz" > "$comb"
[ -s "$BASE/work/clean/${SID}_U.fq.gz" ] && zcat "$BASE/work/clean/${SID}_U.fq.gz" >> "$comb"
total=$(( $(wc -l < "$comb") / 4 ))
[ "$total" -ge "$N" ] || { echo "FATAL $SID total=$total < N=$N (raise floor)"; exit 1; }
# Fixed-depth subsampling to N reads per library (fixed seed). seqkit '-p' does Bernoulli subsampling
# (STREAMS, low memory; '-n' loads ALL reads into RAM and OOMs on deep libraries); oversample to
# ~1.05N in a FILE, then truncate to the first N. Every library ends at EXACTLY N reads = equal depth.
# (This is a fixed-depth subsample, not a uniform exact-size SRS, but is adequate for depth equalization.)
# File input -> seqkit head exits cleanly; no pipe -> no SIGPIPE.
p=$(awk -v n=$N -v t=$total 'BEGIN{v=(n*1.05)/t; if(v>1)v=1; printf "%.6f", v}')
"$SK" sample -s $SEED -p "$p" "$comb" -o "$TMPDIR/${SID}.over.fq"
over=$(( $(wc -l < "$TMPDIR/${SID}.over.fq") / 4 ))
[ "$over" -ge "$N" ] || { echo "FATAL $SID oversample=$over < N=$N (raise p margin)"; exit 1; }
"$SK" head -n $N "$TMPDIR/${SID}.over.fq" -o "$TMPDIR/${SID}.rare.fq"
got=$(( $(wc -l < "$TMPDIR/${SID}.rare.fq") / 4 ))
echo "[$SID] total_combined=$total  p=$p  oversample=$over  rarefied=$got  (target N=$N)"
[ "$got" -eq "$N" ] || { echo "FATAL $SID rarefied=$got != $N"; exit 1; }

metaphlan "$TMPDIR/${SID}.rare.fq" --input_type fastq --nproc 16 \
  --bowtie2db "$DB/metaphlan" --index mpa_vJun23_CHOCOPhlAnSGB_202403 --unclassified_estimation \
  --bowtie2out "$BASE/results/rarefaction/${SID}.bt2.bz2" \
  -o "$BASE/results/rarefaction/${SID}.profile.tsv"
echo "[$SID] done -> results/rarefaction/${SID}.profile.tsv"
# After the array: merge_metaphlan_tables.py results/rarefaction/*.profile.tsv > results/rarefaction/merged_rare13M.tsv
