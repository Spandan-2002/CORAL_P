#!/usr/bin/env bash
#SBATCH -p cs -c 64 --mem=220G -t 8:00:00 -J cs_bracken_build -o logs/%x_%j.out
# Build the Bracken k-mer distribution for the Kraken2 PlusPF DB (read length 150). One-time.
set -euo pipefail
source /scratch/sr7729/CORAL_P/env.sh
module load bioinformatics/20260224
conda activate /scratch/sr7729/conda_envs/cs_profile
[ -e "$DB/kraken2_pluspf/database150mers.kmer_distrib" ] && { echo "bracken db exists"; exit 0; }
bracken-build -d "$DB/kraken2_pluspf" -k 35 -l 150 -t 64
echo "bracken db built"
