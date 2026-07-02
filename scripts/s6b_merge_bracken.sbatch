#!/usr/bin/env bash
#SBATCH -p cs -c 4 --mem=16G -t 0:30:00 -J cs_brkmerge -o logs/%x_%j.out
# Re-merge ALL per-sample Bracken outputs into merged.S.tsv / merged.G.tsv (now n=53).
# Submit with: --dependency=afterok:<s6_array_jobid> so it runs only after the new-5
# Kraken/Bracken array completes. MUST run FROM the bracken dir so combine_bracken_outputs
# names columns by basename SID (not path-prefixed) — sap.R basename()s anyway, but keep it clean.
set -euo pipefail
source /scratch/sr7729/CORAL_P/env.sh
module load bioinformatics/20260224
conda activate /scratch/sr7729/conda_envs/cs_profile
cd "$BASE/results/06_bracken"
nS=$(ls -1 *.S.bracken 2>/dev/null | wc -l); nG=$(ls -1 *.G.bracken 2>/dev/null | wc -l)
echo "[$(date +%T)] merging $nS species-level + $nG genus-level Bracken files"
[ "$nS" -ge 1 ] && combine_bracken_outputs.py --files *.S.bracken -o merged.S.tsv
[ "$nG" -ge 1 ] && combine_bracken_outputs.py --files *.G.bracken -o merged.G.tsv
echo "[$(date +%T)] merged.S.tsv samples (num cols): $(head -1 merged.S.tsv | tr '\t' '\n' | grep -c '_num')"
echo "[$(date +%T)] done -> results/06_bracken/merged.S.tsv  merged.G.tsv"
