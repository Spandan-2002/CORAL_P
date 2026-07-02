# CORAL_P

**Shotgun-metagenomic reads → taxonomic relative abundance** for the CORAL ingestible
core/shell sampling device (Ramadi Lab, NYU). This is the minimal profiling pipeline: it takes raw
shotgun FASTQs for the device **core**, **shell**, **stool**, and **control** samples and produces
per-sample and merged **relative-abundance tables** (MetaPhlAn SGB + Kraken2/Bracken).

> Extracted from the full CoralShot analysis project — this repo contains **only** the
> read → relative-abundance stages (QC → host removal → taxonomic profiling). Assembly/MAGs, strain
> analysis, function (HUMAnN), virome, and statistics are **not** included here.

## What it produces

| Output | Tool | Level |
|---|---|---|
| `results/04_mpa/merged_sgb.tsv` | MetaPhlAn 4 | species / SGB relative abundance (%) |
| `results/06_bracken/merged.S.tsv` | Kraken2 + Bracken | species (counts + fraction) |
| `results/06_bracken/merged.G.tsv` | Kraken2 + Bracken | genus (counts + fraction) |
| `figures/stacked_composition.png` | `make_stacked_composition.py` | stacked genus composition — controls vs stool / core / shell |

Each column is one library; library naming encodes compartment and group (below).

## Sample naming — control / core / stool / shell

Library IDs are `<subject><arm>_<compartment>`:

- **compartment suffix:** `_C` = **core**, `_S` = **shell**, `_F` = **stool (feces)**
- **`group` column:** `control` (un-ingested control units — core+shell only) vs `participant` (core+shell+stool)
- examples: `P01M_C` = participant P01, arm M, **core**; `C01B_S` = control C01, arm B, **shell**

See `sample_manifest.tsv` for the full sheet (subject, arm, compartment, i7/i5, read counts).

## Pipeline

```
raw FASTQ  (Shotgun/)
  s1_qc            FastQC on raw reads
  s2_fastp         adapter / quality trim (fastp)
  s3_host          host-read removal (bowtie2 vs T2T-CHM13 + GRCh38) -> microbial reads
  ├─ s4_metaphlan      MetaPhlAn 4 per sample
  │    s4b_merge_metaphlan   merge_metaphlan_tables.py -> results/04_mpa/merged_sgb.tsv
  └─ s6_kraken_bracken Kraken2 (PlusPF) + Bracken per sample
        s6b_merge_bracken  combine_bracken_outputs.py -> results/06_bracken/merged.{S,G}.tsv

  merged_sgb.tsv -> make_stacked_composition.py -> figures/stacked_composition.png
```

Setup stages: `s0_build_host` (bowtie2 host index) · `s0_bracken` (Bracken k-mer DB) ·
`stage_dbs.sh` (stage MetaPhlAn + Kraken2 PlusPF DBs) · `build_manifest.py` (sample sheet).

## Quick start (NYU Torch HPC / SLURM)

```bash
# 0. paths + caches + offline fixes  (edit BASE in env.sh for another system)
source env.sh                       # sets BASE, ACCT, CPU_QOS
mkdir -p logs                       # SLURM opens -o logs/... at job start (also kept via logs/.gitkeep)

# 1. conda env for profiling (MetaPhlAn / Kraken2 / Bracken)
conda env create -p /scratch/sr7729/conda_envs/cs_profile -f envs/cs_profile.yml

# 2. stage reference DBs, then build the host index + Bracken DB
bash   scripts/stage_dbs.sh 2>&1 | tee logs/stage_dbs.log   # host + MetaPhlAn + Kraken2 PlusPF only
sbatch -A "$ACCT" --qos="$CPU_QOS" scripts/s0_build_host.sbatch
sbatch -A "$ACCT" --qos="$CPU_QOS" scripts/s0_bracken.sbatch

# 3. place raw FASTQs in Shotgun/  (named <SID>_S##_L002_R{1,2}_001.fastq.gz).
#    manifest.tsv + filelist.txt for the original 53-library run are COMMITTED, so you can run as-is
#    once the FASTQs are present. To REGENERATE the sheet for a NEW run, build_manifest.py also needs
#    the bcl2fastq run-report HTML at $BASE/<run>.html (barcode/QC metadata are merged from it):
python scripts/build_manifest.py    # -> manifest.tsv, filelist.txt, asm_units.tsv, barcode_audit.txt

# 4. run: QC -> trim -> host removal -> profiling  (submit.sh chains SLURM dependencies)
J1=$(bash scripts/submit.sh s1     | awk '{print $NF}')
J2=$(bash scripts/submit.sh s2 $J1 | awk '{print $NF}')
J3=$(bash scripts/submit.sh s3 $J2 | awk '{print $NF}')
J4=$(bash scripts/submit.sh s4 $J3 | awk '{print $NF}')     # MetaPhlAn relative abundance
J6=$(bash scripts/submit.sh s6 $J3 | awk '{print $NF}')     # Kraken2 / Bracken

# 5. merge per-sample profiles into the two relative-abundance tables
sbatch -A "$ACCT" --qos="$CPU_QOS" --dependency=afterok:$J4 scripts/s4b_merge_metaphlan.sbatch  # -> merged_sgb.tsv
sbatch -A "$ACCT" --qos="$CPU_QOS" --dependency=afterok:$J6 scripts/s6b_merge_bracken.sbatch    # -> merged.S.tsv / merged.G.tsv

# 6. stacked genus-composition figure (login node; needs cs_viz = pandas + numpy + matplotlib)
conda env create -p /scratch/sr7729/conda_envs/cs_viz -f envs/cs_viz.yml
conda run -p /scratch/sr7729/conda_envs/cs_viz python scripts/make_stacked_composition.py  # -> figures/stacked_composition.png
```

> `submit.sh` also lists downstream stages (`s5`, `s7`, `s8`, …) inherited from the parent project —
> **only `s0`–`s6` apply here**; the other branches reference scripts not included in this repo.

## Requirements

- **Cluster:** NYU Torch HPC (SLURM). `env.sh` encodes Torch-specific fixes (tool caches redirected
  off an inode-full `/home`; node-local bound `$TMPDIR`; offline model/DB caches). For another
  system, edit `BASE`, the SLURM account/QOS, and the temp/cache paths in `env.sh`.
- **conda:** `envs/cs_profile.yml` — MetaPhlAn 4, Bracken, seqkit. **Kraken2**, plus the QC/host
  tools (FastQC, fastp, bowtie2, samtools), come from the Torch `module load bioinformatics/20260224`
  bundle — install or module-load equivalents on another system.
- **conda (figure):** `envs/cs_viz.yml` — pandas + numpy + matplotlib for `make_stacked_composition.py`.
- **Databases (not in git):** MetaPhlAn `mpa_vJun23_CHOCOPhlAnSGB_202403`; Kraken2 PlusPF; Bracken
  k-mer distribution; bowtie2 host index (T2T-CHM13 + GRCh38). Staged by `stage_dbs.sh` / `s0_*`.

## Not included

Raw FASTQs, reference DBs, and pipeline outputs are **gitignored** (regenerated). This repository is
code only. Downstream analysis (assembly, MAGs, strain sharing, function, virome, statistics) lives
in the parent CoralShot project.

## License

See [`LICENSE`](LICENSE). Default is MIT — **confirm or replace** per Ramadi Lab / NYU policy before
any public release.
