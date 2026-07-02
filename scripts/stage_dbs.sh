#!/usr/bin/env bash
# Stage the reference databases for the reads -> relative-abundance pipeline.
# RUN ON THE LOGIN NODE (compute nodes have no internet).
#   bash scripts/stage_dbs.sh 2>&1 | tee logs/stage_dbs.log
# Idempotent-ish: skips a target if its sentinel exists. Logs sizes to manifest_dbs.tsv.
# Stages ONLY the DBs this subset needs: host (GRCh38 + T2T-CHM13 + PhiX), MetaPhlAn 4 SGB,
# and Kraken2 PlusPF. (The full CoralShot project also stages HUMAnN / CheckM2 / GTDB-Tk /
# geNomad / CheckV / AMRFinder / antiSMASH — those are out of scope here and are omitted.)
set -euo pipefail
source /scratch/sr7729/CORAL_P/env.sh
module load bioinformatics/20260224 || true   # for datasets/bowtie2-build if needed
mkdir -p "$DB"/{host,metaphlan,kraken2_pluspf} "$BASE/logs"
MAN=$BASE/manifest_dbs.tsv
[ -f "$MAN" ] || echo -e "name\tpath\tbytes\tmd5_or_note\tdate" > "$MAN"
log(){ local n="$1" p="$2"; local b; b=$(du -sb "$p" 2>/dev/null | cut -f1 || echo 0); \
       echo -e "${n}\t${p}\t${b}\t$(date +%F)\t$(date +%F)" >> "$MAN"; echo ">> staged $n ($((b/1024/1024/1024)) GB)"; }
need(){ local need_gb="$1"; local avail; avail=$(df -BG --output=avail "$DB" | tail -1 | tr -dc '0-9'); \
        [ "$avail" -ge "$need_gb" ] || { echo "FATAL: need ${need_gb}G, have ${avail}G free"; exit 1; }; }

# 1. GRCh38 (no-alt analysis set) -------------------------------------------------
if [ ! -e "$DB/host/GRCh38.fa.gz" ]; then need 10
  wget -q -O "$DB/host/GRCh38.fa.gz" \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
  log GRCh38 "$DB/host/GRCh38.fa.gz"; fi
# 2. T2T-CHM13v2.0 ---------------------------------------------------------------
if [ ! -e "$DB/host/chm13v2.0.fa.gz" ]; then need 5
  wget -q -O "$DB/host/chm13v2.0.fa.gz" \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz"
  log T2T-CHM13 "$DB/host/chm13v2.0.fa.gz"; fi
# 2b. PhiX (NC_001422) -----------------------------------------------------------
if [ ! -e "$DB/host/phix.fa.gz" ]; then
  wget -q -O "$DB/host/phix.fa" \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=NC_001422.1&rettype=fasta&retmode=text"
  gzip -f "$DB/host/phix.fa"   # efetch returns PLAINTEXT fasta -> compress to a valid .gz so zcat in s0_build_host works
  log PhiX "$DB/host/phix.fa.gz"; fi
# NOTE: the combined bowtie2 host index is BUILT on a compute node -> scripts/s0_build_host.sbatch

# 3. MetaPhlAn 4.1 SGB DB (download + index build; ~25 GB) -----------------------
if [ ! -e "$DB/metaphlan/mpa_vJun23_CHOCOPhlAnSGB_202403.pkl" ]; then need 30
  conda activate /scratch/sr7729/conda_envs/cs_profile
  metaphlan --install --index mpa_vJun23_CHOCOPhlAnSGB_202403 --bowtie2db "$DB/metaphlan"
  conda deactivate; log MetaPhlAn "$DB/metaphlan"; fi

# 4. Kraken2 PlusPF (~160 GB; pin exact dated build) -----------------------------
# Pick the latest k2_pluspf_YYYYMMDD from https://benlangmead.github.io/aws-indexes/k2 and set K2URL.
K2URL="${K2URL:-https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20240904.tar.gz}"
if [ ! -e "$DB/kraken2_pluspf/hash.k2d" ]; then need 200
  wget -q -O "$DB/kraken2_pluspf.tar.gz" "$K2URL"
  tar xzf "$DB/kraken2_pluspf.tar.gz" -C "$DB/kraken2_pluspf" && rm -f "$DB/kraken2_pluspf.tar.gz"
  echo "$K2URL" > "$DB/kraken2_pluspf/SOURCE_URL.txt"; log Kraken2_PlusPF "$DB/kraken2_pluspf"; fi
# NOTE: the Bracken k-mer distribution is BUILT on a compute node -> scripts/s0_bracken.sbatch

echo "DB staging complete (host + MetaPhlAn + Kraken2 PlusPF). See $MAN"
