#!/usr/bin/env bash
#SBATCH -p cs -c 8 --mem=12G -t 1:00:00 --array=1-106%24 -J cs_qc_raw -o logs/%x_%A_%a.out
# Raw read QC on all 106 fastq files (53 libraries x 2; one per array task).
# Array size is authoritatively (re)set by scripts/submit.sh (1-$nfiles); this header is for direct `sbatch`.
set -euo pipefail
source /scratch/sr7729/CORAL_P/env.sh
module load bioinformatics/20260224
mkdir -p "$BASE/results/01_fastqc_raw"
F=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BASE/filelist.txt")
[ -z "$F" ] && { echo "FATAL: empty filelist row $SLURM_ARRAY_TASK_ID"; exit 1; }
fastqc -t 8 -o "$BASE/results/01_fastqc_raw" "$F"
# After the array finishes, run once (see submit.sh post-hooks):
#   multiqc -f -o results/01_fastqc_raw results/01_fastqc_raw
#   seqkit stats -a -T -j 8 Shotgun/*.fastq.gz > results/01_fastqc_raw/seqkit_raw.tsv
