#!/usr/bin/env bash
#SBATCH -p cs -c 4 --mem=16G -t 0:30:00 -J cs_mpamerge -o logs/%x_%j.out
# Merge ALL per-sample MetaPhlAn profiles into results/04_mpa/merged_sgb.tsv
# (the MetaPhlAn SGB relative-abundance table). Mirror of s6b_merge_bracken.sbatch.
# Submit with: --dependency=afterok:<s4_array_jobid> so it runs after the MetaPhlAn array.
set -euo pipefail
source /scratch/sr7729/CORAL_P/env.sh
module load bioinformatics/20260224
conda activate /scratch/sr7729/conda_envs/cs_profile
cd "$BASE/results/04_mpa"
n=$(ls -1 *.profile.tsv 2>/dev/null | wc -l)
echo "[$(date +%T)] merging $n MetaPhlAn per-sample profiles"
[ "$n" -ge 1 ] && merge_metaphlan_tables.py *.profile.tsv > merged_sgb.tsv
echo "[$(date +%T)] done -> results/04_mpa/merged_sgb.tsv  (cols: $(head -2 merged_sgb.tsv | tail -1 | tr '\t' '\n' | wc -l))"
