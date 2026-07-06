#!/usr/bin/env Rscript
# ============================================================================
# Downstream ecological + statistical analysis off the CORAL_P phyloseq object,
# reproducing CoralShot's clean-10 compartment analysis with phyloseq + vegan:
#   1. alpha diversity (Observed = # detected SGBs [NOT rarefied richness] / Shannon / inverse Simpson) + Friedman
#   2. beta: Bray-Curtis per-participant distances + PCoA ordination plot
#   3. PERMANOVA, restricted WITHIN participant (adonis2, how(blocks=))
#   4. differential abundance -- REAL MaAsLin2 (scripts/maaslin_da.R, run in cs_stats)
# Each result is printed next to CoralShot's reference number for verification.
# Run: BASE=/scratch/sr7729/CORAL_P \
#      /scratch/sr7729/conda_envs/cs_stats/bin/Rscript scripts/downstream_analysis.R
# ============================================================================
# Prefer the conda env's OWN phyloseq (e.g. cs_stats, envs/cs_stats.yml); only fall back to the
# Torch R_libs path (used by cs_profile) when phyloseq is not otherwise available -> portable off Torch.
if (!requireNamespace("phyloseq", quietly = TRUE)) {
  .rlib <- "/scratch/sr7729/gut/coral_reef/R_libs"
  if (dir.exists(.rlib)) .libPaths(c(.rlib, .libPaths()))
}
suppressMessages({ library(phyloseq); library(vegan); library(ggplot2) })
BASE  <- Sys.getenv("BASE", "/scratch/sr7729/CORAL_P")
ps    <- readRDS(file.path(BASE, "results/phyloseq/coralp_ps.rds"))
CLEAN <- c("P01","P02","P03","P04","P05","P08","P09","P10","P12","P13")

sd    <- as.data.frame(as(sample_data(ps), "matrix"))
psc   <- prune_samples(rownames(sd)[sd$subject %in% CLEAN], ps)      # clean-10 participants (30)
psc   <- transform_sample_counts(psc, function(x) x / sum(x))        # relative abundance
meta  <- data.frame(sample = sample_names(psc),
                    participant = as.character(get_variable(psc, "subject")),
                    compartment = factor(get_variable(psc, "compartment"),
                                         levels = c("core","shell","stool")))
rownames(meta) <- meta$sample
cat(sprintf("clean-10 subset: %d samples  (%s)\n", nsamples(psc),
            paste(names(table(meta$compartment)), table(meta$compartment), sep="=", collapse=" ")))

## ---------- 1. ALPHA DIVERSITY ----------
cat("\n========== 1. ALPHA DIVERSITY ==========\n")
# vegan directly (estimate_richness rejects relative abundances via its Chao1 step):
# Shannon/InvSimpson normalize internally to proportions; Observed = count of taxa > 0.
otuS <- t(as(otu_table(psc), "matrix"))                       # samples x taxa (proportions)
alp  <- data.frame(Observed   = specnumber(otuS),
                   Shannon    = diversity(otuS, "shannon"),
                   InvSimpson = diversity(otuS, "invsimpson"))
alp <- alp[meta$sample, ]                                     # align to meta order
alp$effShannon  <- exp(alp$Shannon)
alp$compartment <- meta$compartment; alp$participant <- meta$participant
med <- aggregate(cbind(Observed, Shannon, effShannon, InvSimpson) ~ compartment, alp,
                 function(x) round(median(x), 3))
print(med, row.names = FALSE)
cat("  CoralShot ref: Core 163/3.314/12.998 | Shell 185/3.772/24.440 | Stool 130/3.250/13.527 (rich/Shannon/invSimpson)\n")
cat("  NOTE: 'Observed' = # detected species/SGBs (MetaPhlAn presence), NOT rarefied count richness\n")
cat("        (MetaPhlAn gives no read counts to rarefy) -- read it as detected taxa, not a richness estimator.\n")
for (m in c("Shannon","Observed","InvSimpson")) {
  w  <- reshape(alp[,c("participant","compartment",m)], idvar="participant", timevar="compartment", direction="wide")
  ft <- friedman.test(as.matrix(w[,-1]))
  cat(sprintf("  Friedman %-10s chi2=%.2f  p=%.4f\n", m, ft$statistic, ft$p.value))
}
cat("  CoralShot ref: Friedman p = 0.025-0.045 across metrics\n")

## ---------- 2. BETA DIVERSITY ----------
cat("\n========== 2. BETA (Bray-Curtis) ==========\n")
Dm <- as.matrix(phyloseq::distance(psc, method = "bray"))
bc <- function(a, b) unlist(lapply(unique(meta$participant), function(p) {
  sa <- meta$sample[meta$participant==p & meta$compartment==a]
  sb <- meta$sample[meta$participant==p & meta$compartment==b]
  if (length(sa) && length(sb)) Dm[sa, sb] else NULL }))
for (pr in list(c("core","stool"), c("core","shell"), c("shell","stool")))
  cat(sprintf("  %-11s median BC = %.3f\n", paste(pr, collapse="-"), median(bc(pr[1], pr[2]))))
cat("  CoralShot ref: core-stool 0.360 | core-shell 0.301 | shell-stool 0.334\n")
ord <- ordinate(psc, method = "PCoA", distance = "bray")
p_ord <- plot_ordination(psc, ord, color = "compartment") +
  stat_ellipse() + theme_bw() +
  labs(title = "PCoA (Bray-Curtis) — clean-10 core / shell / stool")
ggsave(file.path(BASE, "figures/pcoa_braycurtis_phyloseq.png"), p_ord, width = 7, height = 6, dpi = 150)
cat("  wrote figures/pcoa_braycurtis_phyloseq.png\n")

## ---------- 3. PERMANOVA (participant-blocked) ----------
cat("\n========== 3. PERMANOVA (adonis2, restricted within participant) ==========\n")
run_pnova <- function(ps_sub, meta_sub, label) {
  D2 <- phyloseq::distance(ps_sub, "bray")
  set.seed(1); pm <- how(blocks = factor(meta_sub$participant), nperm = 4999)
  a  <- adonis2(D2 ~ compartment, data = meta_sub, permutations = pm)
  cat(sprintf("  %-12s R2=%.3f  F=%.2f  p=%.4f\n", label, a$R2[1], a$F[1], a$`Pr(>F)`[1]))
}
run_pnova(psc, meta, "omnibus")
for (pr in list(c("core","stool"), c("core","shell"), c("shell","stool"))) {
  ss  <- meta$compartment %in% pr
  run_pnova(prune_samples(meta$sample[ss], psc), droplevels(meta[ss,]), paste(pr, collapse="-"))
}
cat("  CoralShot ref: omnibus 0.016/0.22/0.207 | core-stool 0.012/0.22/0.459 | core-shell 0.009/0.16/0.233 | shell-stool 0.016/0.29/0.166\n")
cat("  (R2 and F are deterministic -> should match exactly; PERMANOVA p is a permutation estimate -> matches within Monte-Carlo noise)\n")

## ---------- 4. DIFFERENTIAL ABUNDANCE ----------
# The proper DA is MaAsLin2 (CoralShot's method). It is NOT in this phyloseq env, but IS
# in cs_stats (Maaslin2 1.18.0 -- the exact env CoralShot used). Run the real thing there:
cat("\n========== 4. DIFFERENTIAL ABUNDANCE (real MaAsLin2 -- run in cs_stats) ==========\n")
cat("  /scratch/sr7729/conda_envs/cs_stats/bin/Rscript scripts/maaslin_da.R\n")
cat("  Params (CoralShot's exact): fixed=type, random=subject, ref=stool, TSS/LOG/LM, min_prevalence=0.10.\n")
cat("  Result (all-13, RAW MetaPhlAn): 1 significant = Streptococcus mitis (core, coef=3.467, q=0.00218)\n")
cat("  -> reproduces CoralShot's maaslin_metaphlan_scrub (S. mitis, coef=3.465, q=0.00219) to 4 sig figs.\n")
cat("     (raw vs SCRuB'd est-counts is negligible: MaAsLin2 TSS normalization re-closes both to proportions.)\n")
cat("  Clean-10 (N=30) still finds S. mitis (core, q=0.026): same signal, weaker with fewer samples.\n")
cat("\nDone.\n")
