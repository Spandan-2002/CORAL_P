#!/usr/bin/env Rscript
# ============================================================================
# REAL MaAsLin2 differential abundance off the CORAL_P phyloseq object, using
# CoralShot's exact parameters (Maaslin2 1.18.0, the env CoralShot itself used).
#   fixed_effects="type"(=compartment) · random_effects="subject" · reference stool
#   normalization=TSS · transform=LOG · analysis_method=LM · min_prevalence=0.10
# Run BOTH cohorts so the comparison to CoralShot is honest:
#   (a) all-13 participants  -> matches CoralShot's sample set (N=39)
#   (b) clean-10 participants -> consistent with the rest of the downstream
# CoralShot's committed run (maaslin_metaphlan_scrub) used SCRuB-decontaminated,
# estimated-count MetaPhlAn on all 13; this uses the RAW MetaPhlAn from the phyloseq
# object -> any difference is attributable to SCRuB (which CoralShot found over-corrects).
# Run: /scratch/sr7729/conda_envs/cs_stats/bin/Rscript scripts/maaslin_da.R
# ============================================================================
suppressMessages({ library(phyloseq); library(Maaslin2) })
BASE  <- Sys.getenv("BASE", "/scratch/sr7729/CORAL_P")
ps    <- readRDS(file.path(BASE, "results/phyloseq/coralp_ps.rds"))
sd    <- as.data.frame(as(sample_data(ps), "matrix"))
otu   <- as(otu_table(ps), "matrix")                       # taxa x samples (MetaPhlAn relab)
CLEAN <- c("P01","P02","P03","P04","P05","P08","P09","P10","P12","P13")
ALL13 <- sprintf("P%02d", 1:13)

run_maaslin <- function(subjects, label) {
  keep <- rownames(sd)[sd$subject %in% subjects]
  inp  <- as.data.frame(t(otu[, keep]))                    # samples x taxa
  meta <- data.frame(type = sd[keep, "compartment"],
                     subject = sd[keep, "subject"], row.names = keep)
  outdir <- file.path(BASE, "results/phyloseq", paste0("maaslin_", label))
  fit <- Maaslin2(input_data = inp, input_metadata = meta, output = outdir,
                  fixed_effects = "type", random_effects = "subject",
                  reference = "type,stool", normalization = "TSS", transform = "LOG",
                  analysis_method = "LM", min_prevalence = 0.10, max_significance = 0.05,
                  plot_heatmap = FALSE, plot_scatter = FALSE, standardize = FALSE)
  sig <- tryCatch(read.delim(file.path(outdir, "significant_results.tsv")), error = function(e) data.frame())
  cat(sprintf("\n=== MaAsLin2 [%s]  N=%d samples ===\n", label, length(keep)))
  cat(sprintf("  significant (q<0.05): %d\n", nrow(sig)))
  if (nrow(sig)) print(head(sig[, c("feature","metadata","value","coef","qval")], 10), row.names = FALSE)
}

run_maaslin(ALL13, "all13_raw")     # CoralShot sample set, RAW MetaPhlAn (no SCRuB)
run_maaslin(CLEAN, "clean10_raw")   # clean-10 cohort, consistent with alpha/beta/PERMANOVA
cat("\nCoralShot committed ref (maaslin_metaphlan_scrub, SCRuB'd est-counts, N=39): 1 sig = Streptococcus mitis (core, q=0.00219)\n")
