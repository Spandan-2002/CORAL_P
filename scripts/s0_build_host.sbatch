#!/usr/bin/env bash
#SBATCH -p cs -c 32 --mem=200G -t 6:00:00 -J cs_hostidx -o logs/%x_%j.out
# Build the combined human+decoy bowtie2 index (GRCh38 + T2T-CHM13 + PhiX).
set -euo pipefail
source /scratch/sr7729/CORAL_P/env.sh
module load bioinformatics/20260224
cd "$DB/host"
[ -e host_decoy.1.bt2l ] && { echo "host index exists"; exit 0; }
zcat GRCh38.fa.gz chm13v2.0.fa.gz phix.fa.gz > "$TMPDIR/host_decoy.fa"
bowtie2-build --large-index --threads 32 "$TMPDIR/host_decoy.fa" host_decoy
echo "built $DB/host/host_decoy.*"
