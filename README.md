# CORAL_P

**Shotgun reads → taxonomic relative abundance → composition figure** for the CORAL ingestible
core/shell sampling device (Ramadi Lab, NYU). Minimal, reproducible pipeline:
**fastp → bowtie2 host removal → MetaPhlAn 4 → merged relative-abundance table → stacked genus
composition** (drawn both with matplotlib and with phyloseq).

> 📋 **Copy-paste runbook:** [`RUN_AND_VERIFY.txt`](RUN_AND_VERIFY.txt) — every command, grouped and ready.
> Extracted from the parent **CoralShot** project — this repo is the read→relative-abundance subset only
> (no assembly/MAGs, strain, HUMAnN, virome, or stats).

---

## ▶ Run it (NYU Torch — fastest path)

The MetaPhlAn DB, prebuilt host index, and raw reads are already staged in the parent project, so on
Torch you reuse them and go straight to running — **no database download, no index build:**

```bash
cd /scratch/sr7729/CORAL_P
source env.sh                                     # sets BASE, ACCT, CPU_QOS, caches, offline fixes
rmdir db Shotgun 2>/dev/null                      # clear empty placeholders (else ln nests the link inside)
ln -s /scratch/sr7729/coralshot/db      db        # MetaPhlAn DB + prebuilt bowtie2 host index
ln -s /scratch/sr7729/coralshot/Shotgun Shotgun   # raw FASTQs (manifest.tsv + filelist.txt are committed)

# pipeline: QC -> trim -> host removal -> MetaPhlAn -> merge   (SLURM; ~a few hours on 53 libraries)
J1=$(bash scripts/submit.sh s1     | awk '{print $NF}')
J2=$(bash scripts/submit.sh s2 $J1 | awk '{print $NF}')
J3=$(bash scripts/submit.sh s3 $J2 | awk '{print $NF}')
J4=$(bash scripts/submit.sh s4 $J3 | awk '{print $NF}')
sbatch -A "$ACCT" --qos="$CPU_QOS" --dependency=afterok:$J4 scripts/s4b_merge_metaphlan.sbatch
squeue -u $USER

# figures (after the merge finishes) -> results/04_mpa/merged_sgb.tsv must exist
conda run -p /scratch/sr7729/conda_envs/cs_viz python scripts/make_stacked_composition.py         # matplotlib
BASE=$PWD /scratch/sr7729/conda_envs/cs_profile/bin/Rscript scripts/make_phyloseq.R                # phyloseq object + figure
```

The `cs_profile` and `cs_viz` conda envs already exist on Torch. Create them from `envs/*.yml` if missing
(see *Requirements*). **Already have `merged_sgb.tsv`?** Skip the pipeline — just run the two figure
commands (seconds).

---

## What it produces

| Output | Made by | What it is |
|---|---|---|
| `results/04_mpa/merged_sgb.tsv` | MetaPhlAn 4 (`s4` + `s4b`) | species/SGB **relative abundance (%)**, one column per library |
| `figures/stacked_composition.png` | `make_stacked_composition.py` | stacked genus composition — controls \| stool \| core \| shell |
| `figures/stacked_composition_phyloseq.png` | `make_phyloseq.R` | same, drawn with **phyloseq** (faceted CONTROL \| core \| shell \| stool) |
| `results/phyloseq/coralp_ps.rds` | `make_phyloseq.R` | **phyloseq object**: 817 species × 53 samples, full taxonomy + sample_data |

## Verify the relative abundance independently (phyloseq)

phyloseq is a *container*, not a profiler — load the same MetaPhlAn table and it reports the same
numbers. Two gotchas that otherwise cause a false mismatch: **filter to one rank** (genus = `|g__`,
not `|s__`) and **close each sample to sum = 1**. `phyloseq_relabund.R` does this and reproduces the
plotted numbers to **5×10⁻⁹** (verified):

```bash
BASE=$PWD /scratch/sr7729/conda_envs/cs_profile/bin/Rscript scripts/phyloseq_relabund.R
#   -> results/04_mpa/phyloseq_genus_relabund.tsv   (should diff to ~0 against ours)
```
Load the object in R: `ps <- readRDS("results/phyloseq/coralp_ps.rds")` — then `tax_glom(ps,"Genus")`
relative abundance equals the `|g__` numbers exactly (Bifidobacterium P01M_C = 0.08849).

## Sample naming — control / core / stool / shell

Library IDs are `<subject><arm>_<compartment>`:
`_C` = **core** · `_S` = **shell** · `_F` = **stool (feces)**.
Controls = `C0x` (core+shell only); participants = `P0x` (core+shell+stool).
Full sheet with subject/arm/i7/i5/read-counts: `sample_manifest.tsv`.

---

## Pipeline (from scratch — non-Torch or a new run)

```
raw FASTQ  (Shotgun/)
  s1_qc            FastQC on raw reads
  s2_fastp         adapter / quality trim (fastp)
  s3_host          host-read removal (bowtie2 vs T2T-CHM13 + GRCh38)
  s4_metaphlan     MetaPhlAn 4 per sample
  s4b_merge        merge_metaphlan_tables.py -> results/04_mpa/merged_sgb.tsv
  figures          make_stacked_composition.py (matplotlib) + make_phyloseq.R (phyloseq)
```

Full setup when the DBs are **not** already staged:

```bash
source env.sh
mkdir -p logs
conda env create -p /scratch/sr7729/conda_envs/cs_profile -f envs/cs_profile.yml   # MetaPhlAn (+ R/phyloseq on Torch)
conda env create -p /scratch/sr7729/conda_envs/cs_viz     -f envs/cs_viz.yml       # matplotlib figure
bash   scripts/stage_dbs.sh 2>&1 | tee logs/stage_dbs.log                          # host + MetaPhlAn DBs
sbatch -A "$ACCT" --qos="$CPU_QOS" scripts/s0_build_host.sbatch                    # build bowtie2 host index
# put FASTQs in Shotgun/ (named <SID>_S##_L002_R{1,2}_001.fastq.gz); manifest.tsv + filelist.txt are committed.
# To regenerate the sheet for a NEW run, build_manifest.py additionally needs the bcl2fastq run-report HTML:
python scripts/build_manifest.py     # -> manifest.tsv, filelist.txt, asm_units.tsv, barcode_audit.txt
# then run the s1->s4->s4b chain + figures as in "Run it" above.
```

> `submit.sh` still lists parent-project stages (`s5`, `s6`, `s7`, …) — **only `s1`–`s4` apply here**;
> the other cases reference scripts not shipped in this subset.

## Requirements

- **Cluster:** NYU Torch HPC (SLURM). `env.sh` encodes Torch-specific fixes (caches off an inode-full
  `/home`; node-local bound `$TMPDIR`; offline caches). For another system, edit `BASE`, the SLURM
  account/QOS, and the cache/temp paths in `env.sh`.
- **conda `cs_profile`** (`envs/cs_profile.yml`): MetaPhlAn 4 (+ seqkit). QC/host tools (FastQC, fastp,
  bowtie2, samtools) come from the Torch `module load bioinformatics/20260224` bundle.
- **conda `cs_viz`** (`envs/cs_viz.yml`): pandas + numpy + matplotlib (for `make_stacked_composition.py`).
- **R + phyloseq** (for `make_phyloseq.R` / `phyloseq_relabund.R`): on Torch these run via
  `cs_profile`'s `Rscript` with phyloseq on the `coral_reef/R_libs` path (set inside the scripts);
  elsewhere install `phyloseq` + `ggplot2`.
- **Databases (not in git):** MetaPhlAn `mpa_vJun23_CHOCOPhlAnSGB_202403`; bowtie2 host index
  (T2T-CHM13 + GRCh38). Staged by `stage_dbs.sh` / `s0_build_host` (or reused via the symlinks above).

## Not included

Raw FASTQs, reference DBs, and pipeline outputs are **gitignored** (regenerated) — except the two
figures and the phyloseq `.rds`, committed as ready-to-use artifacts. This repo is otherwise code
only; downstream analysis (assembly, MAGs, strain, function, virome, stats) lives in parent CoralShot.

## License

See [`LICENSE`](LICENSE) — MIT.
