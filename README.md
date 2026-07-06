# CORAL_P

**Shotgun reads → taxonomic relative abundance → composition figure** for the CORAL ingestible
core/shell sampling device (Ramadi Lab, NYU). Minimal, reproducible pipeline:
**fastp → bowtie2 host removal → MetaPhlAn 4 → merged relative-abundance table → stacked genus
composition** (drawn both with matplotlib and with phyloseq).

> 📋 **Copy-paste runbook:** [`RUN_AND_VERIFY.txt`](RUN_AND_VERIFY.txt) — every command, grouped and ready.
> Extracted from the parent **CoralShot** project — the read→relative-abundance subset **plus a
> phyloseq/vegan verification layer** (alpha / beta / participant-blocked PERMANOVA + MaAsLin2
> differential abundance). The heavy downstream (assembly/MAGs, strain, HUMAnN, virome, AMR/BGC)
> stays in CoralShot.

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
conda run -p /scratch/sr7729/conda_envs/cs_viz python scripts/make_stacked_composition.py            # matplotlib (cs_viz)
BASE=$PWD /scratch/sr7729/conda_envs/cs_stats/bin/Rscript scripts/make_phyloseq.R                    # phyloseq object + figure (cs_stats)
```

The `cs_profile` (MetaPhlAn), `cs_viz` (matplotlib figure), and `cs_stats` (the R/phyloseq layer) conda
envs already exist on Torch. Create any that are missing from `envs/*.yml` (see *Requirements*).
**Already have `merged_sgb.tsv`?** Skip the pipeline — just run the two figure commands (seconds).

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
BASE=$PWD /scratch/sr7729/conda_envs/cs_stats/bin/Rscript scripts/phyloseq_relabund.R
#   -> results/04_mpa/phyloseq_genus_relabund.tsv   (should diff to ~0 against ours)
```
Load the object in R: `ps <- readRDS("results/phyloseq/coralp_ps.rds")` — then `tax_glom(ps,"Genus")`
relative abundance equals the `|g__` numbers exactly (Bifidobacterium P01M_C = 0.08849).

## Statistics — phyloseq + vegan verification layer

An independent phyloseq + vegan layer reproduces the CoralShot compartment analysis and prints each
number next to its CoralShot reference. **Environment:** create `cs_stats` from
[`envs/cs_stats.yml`](envs/cs_stats.yml) — one env with phyloseq, vegan, ggplot2, and MaAsLin2 runs
**all four** R scripts. No site-specific R-library path is required: the scripts prepend the Torch
`R_libs` path only if it exists, otherwise they use the conda env's own libraries.

```bash
conda env create -p /scratch/sr7729/conda_envs/cs_stats -f envs/cs_stats.yml   # if missing
RS=/scratch/sr7729/conda_envs/cs_stats/bin/Rscript
BASE=$PWD $RS scripts/make_phyloseq.R        # -> results/phyloseq/coralp_ps.rds (817 sp × 53) + figures/stacked_composition_phyloseq.png
BASE=$PWD $RS scripts/phyloseq_relabund.R    # genus relabund cross-check (sample sums = 1)
BASE=$PWD $RS scripts/downstream_analysis.R  # alpha (Friedman) + Bray-Curtis + participant-blocked PERMANOVA + PCoA -> figures/pcoa_braycurtis_phyloseq.png
BASE=$PWD $RS scripts/maaslin_da.R           # REAL MaAsLin2 (v1.18.0): all-13 + clean-10
```

| Result (clean-10 unless noted) | Value |
|---|---|
| alpha Shannon — core / shell / stool | 3.31 / 3.77 / 3.25 (Friedman p ≈ 0.02–0.05) |
| Bray-Curtis medians — core-stool / core-shell / shell-stool | 0.360 / 0.301 / 0.334 |
| PERMANOVA omnibus (adonis2, blocked within participant) | R² 0.016, F 0.22, p ≈ 0.21 |
| Differential abundance (MaAsLin2) | *Streptococcus mitis* ↑ core — all-13 q = 0.0022; clean-10 q = 0.026 |

`Observed` alpha = **# detected species/SGBs** (MetaPhlAn presence), *not* rarefied count richness
(MetaPhlAn emits no read counts to rarefy) — read it as detected taxa.

**Depth-sensitivity (rarefaction) check.** Since MetaPhlAn has no count table to rarefy, depth-robustness
is validated at the **read level** — the MetaPhlAn developer's recommended approach (N. Segata: *"rarefy
the input metagenomes … use seqtk"* — this repo uses **seqkit**, an equivalent subsampler).
`scripts/rarefaction_sensitivity.sbatch` applies **fixed-depth
subsampling** of each clean-10 library's reads to **13M** (fixed seed) and re-profiles;
`scripts/rarefaction_alpha.R` recomputes alpha and compares to native depth. Result:
Shannon / InvSimpson essentially unchanged (per-sample Spearman ρ = 0.999; Friedman p = 0.025 in both),
and Observed richness keeps its median ranking (shell highest; **core~stool is n.s.** — paired Wilcoxon
core-stool p≈0.08–1.0 depending on metric) and its **omnibus** significance (Friedman p = 0.045 native and
rarefied) — so the diversity findings are **not explained by sequencing depth** (they persist at equal
depth; this supports, but does not by itself prove, a biological rather than technical origin). The one
compartment that stands out is the **shell**; core and stool are statistically indistinguishable. The result
text is committed (`results/rarefaction/RAREFACTION_RESULT.txt`), as are the plots from
`scripts/rarefaction_plot.R` (`figures/rarefaction_alpha_box.png` = per-compartment boxplots native
vs rarefied; `figures/rarefaction_alpha_scatter.png` = per-sample native-vs-rarefied scatter with ρ;
`figures/rarefaction_richness.png` = paired Observed-richness native→rarefied per library; `figures/rarefaction_richness_box.png` = dedicated richness boxplot by compartment, native vs rarefied).
The intermediate `merged_rare13M.tsv` and per-sample profiles are gitignored and regenerated by
re-running the scripts. Refs: MetaPhlAn 4
(Blanco-Míguez 2023, *Nat Biotechnol*); depth/richness (Liu 2022, *Genome*, PMID 35939836).

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
  s3_host          host/decoy removal (bowtie2 vs GRCh38 + T2T-CHM13 + PhiX)
  s4_metaphlan     MetaPhlAn 4 per sample
  s4b_merge        merge_metaphlan_tables.py -> results/04_mpa/merged_sgb.tsv
  figures          make_stacked_composition.py (matplotlib) + make_phyloseq.R (phyloseq)
```

Full setup when the DBs are **not** already staged:

```bash
source env.sh
mkdir -p logs
conda env create -p /scratch/sr7729/conda_envs/cs_profile -f envs/cs_profile.yml   # MetaPhlAn
conda env create -p /scratch/sr7729/conda_envs/cs_viz     -f envs/cs_viz.yml       # matplotlib figure
conda env create -p /scratch/sr7729/conda_envs/cs_stats   -f envs/cs_stats.yml     # R: phyloseq + vegan + MaAsLin2
bash   scripts/stage_dbs.sh 2>&1 | tee logs/stage_dbs.log                          # host + MetaPhlAn DBs
sbatch -A "$ACCT" --qos="$CPU_QOS" scripts/s0_build_host.sbatch                    # build bowtie2 host index
# put FASTQs in Shotgun/ (named <SID>_S##_L002_R{1,2}_001.fastq.gz); manifest.tsv + filelist.txt are committed.
# To regenerate the sheet for a NEW run, build_manifest.py additionally needs the bcl2fastq run-report HTML:
python scripts/build_manifest.py     # -> manifest.tsv, filelist.txt, asm_units.tsv, barcode_audit.txt
# then run the s1->s4->s4b chain + figures as in "Run it" above.
```

> `submit.sh` dispatches only the stages shipped here: **`s0host s1 s2 s3 s4 s4b`** (any other name
> is rejected). The phyloseq/stats layer is run directly, not through the dispatcher.

## Requirements

- **Cluster:** NYU Torch HPC (SLURM). `env.sh` encodes Torch-specific fixes (caches off an inode-full
  `/home`; node-local bound `$TMPDIR`; offline caches). For another system, edit `BASE`, the SLURM
  account/QOS, and the cache/temp paths in `env.sh`.
- **conda `cs_profile`** (`envs/cs_profile.yml`): MetaPhlAn 4 (+ seqkit). QC/host tools (FastQC, fastp,
  bowtie2, samtools) come from the Torch `module load bioinformatics/20260224` bundle.
- **conda `cs_viz`** (`envs/cs_viz.yml`): pandas + numpy + matplotlib (for `make_stacked_composition.py`).
- **conda `cs_stats`** (`envs/cs_stats.yml`): the R/phyloseq layer — phyloseq + vegan + ggplot2 +
  MaAsLin2. Runs all four R scripts (`make_phyloseq.R`, `phyloseq_relabund.R`,
  `downstream_analysis.R`, `maaslin_da.R`). No site-specific R-library path is required: the scripts
  prepend the Torch `coral_reef/R_libs` path only if phyloseq is not otherwise available.
- **Databases (not in git):** MetaPhlAn `mpa_vJun23_CHOCOPhlAnSGB_202403`; bowtie2 host+decoy index
  (GRCh38 + T2T-CHM13 + PhiX/NC_001422, built as `host_decoy`). Staged by `stage_dbs.sh` /
  `s0_build_host` (or reused via the symlinks above).

## Not included

Raw FASTQs, reference DBs, and pipeline outputs are **gitignored** (regenerated) — except these
committed ready-to-use artifacts: the two figures, the phyloseq `.rds`, the two MaAsLin2
`significant_results.tsv` tables, and the rarefaction result (`results/rarefaction/RAREFACTION_RESULT.txt`).
This repo is otherwise code only; the heavy downstream (assembly, MAGs, strain, function, virome,
AMR/BGC) lives in parent CoralShot.

## License

See [`LICENSE`](LICENSE) — MIT.
