#!/usr/bin/env Rscript
# ============================================================================
# Read-level rarefaction SENSITIVITY CHECK for alpha diversity.
# Compares species-level alpha (Observed / Shannon / inverse Simpson) computed on
#   (a) NATIVE-depth MetaPhlAn profiles  (results/04_mpa/merged_sgb.tsv), vs
#   (b) profiles from reads RAREFIED to a common 13M depth (results/rarefaction/merged_rare13M.tsv).
# Question: does the compartment ranking (and each sample's rank) survive equal depth?
#   - Shannon / InvSimpson are proportion-based -> expected ~unchanged (depth-robust).
#   - Observed richness is depth-sensitive -> expected to DROP uniformly but keep its ranking.
#
# LITERATURE BASIS (relative abundance for diversity; read-level rarefaction for richness):
#  - MetaPhlAn reports marker-based RELATIVE ABUNDANCE (no count table) -> Shannon/InvSimpson on
#    proportions is the native, correct input. Blanco-Miguez et al. 2023, Nat Biotechnol
#    (MetaPhlAn 4), doi:10.1038/s41587-023-01688-w.
#  - For shotgun, rarefy the INPUT READS (not the profile table), e.g. with seqtk/seqkit -- the
#    MetaPhlAn developer's own guidance (N. Segata, metaphlan-users forum); metagenomics is less
#    depth-noise-prone than 16S. scripts/rarefaction_sensitivity.sbatch implements exactly this.
#  - Rarefaction mainly matters for the depth-correlated statistic (richness); taxonomic composition
#    is stable at modest depth while richness stabilizes ~15M reads: Liu et al. 2022, Genome
#    (PMID 35939836).
#  - Table-level rarefaction of proportions is inadmissible: McMurdie & Holmes 2014, PLoS Comput
#    Biol 10(4):e1003531.
# Run: BASE=/scratch/sr7729/CORAL_P \
#      /scratch/sr7729/conda_envs/cs_stats/bin/Rscript scripts/rarefaction_alpha.R
# ============================================================================
suppressMessages(library(vegan))
BASE  <- Sys.getenv("BASE", "/scratch/sr7729/CORAL_P")
CLEAN <- c("P01","P02","P03","P04","P05","P08","P09","P10","P12","P13")

load_species <- function(path) {
  m  <- read.delim(path, skip = 1, check.names = FALSE)
  colnames(m) <- sub("\\.profile$", "", colnames(m))
  cl <- as.character(m[[1]])
  keep <- grepl("\\|s__", cl) & !grepl("\\|t__", cl)        # species/SGB rows only
  samp <- grep("^[CP][0-9]+[MB]_[CSF]$", colnames(m), value = TRUE)
  x <- as.matrix(m[keep, samp]); storage.mode(x) <- "double"; rownames(x) <- cl[keep]
  subj <- sub("([CP][0-9]+)[MB]_[CSF]", "\\1", samp)
  x[, subj %in% CLEAN, drop = FALSE]                        # clean-10 only
}
alpha <- function(x) {
  p <- t(x)                                                 # samples x taxa (relative abundance)
  data.frame(Observed   = specnumber(p),
             Shannon    = diversity(p, "shannon"),
             InvSimpson = diversity(p, "invsimpson"),
             row.names  = rownames(p))
}

nat <- alpha(load_species(file.path(BASE, "results/04_mpa/merged_sgb.tsv")))
rar <- alpha(load_species(file.path(BASE, "results/rarefaction/merged_rare13M.tsv")))
common <- intersect(rownames(nat), rownames(rar))
nat <- nat[common, ]; rar <- rar[common, ]
compartment <- factor(c(C = "core", S = "shell", F = "stool")[substr(common, nchar(common), nchar(common))],
                      levels = c("core", "shell", "stool"))
cat(sprintf("compared %d clean-10 samples (native vs rarefied to 13M reads)\n", length(common)))

for (lab in c("Observed", "Shannon", "InvSimpson")) {
  cat(sprintf("\n===== %s =====\n", lab))
  mn <- tapply(nat[[lab]], compartment, median); mr <- tapply(rar[[lab]], compartment, median)
  for (k in c("core","shell","stool"))
    cat(sprintf("  %-6s  native median = %8.3f    rarefied median = %8.3f\n", k, mn[k], mr[k]))
  ord_n <- paste(names(sort(mn, decreasing = TRUE)), collapse = " > ")
  ord_r <- paste(names(sort(mr, decreasing = TRUE)), collapse = " > ")
  cat(sprintf("  compartment ranking : native [%s]  vs  rarefied [%s]  -> %s\n",
              ord_n, ord_r, if (ord_n == ord_r) "PRESERVED" else "CHANGED"))
  cat(sprintf("  per-sample Spearman rho (native vs rarefied) = %.3f\n",
              cor(nat[[lab]], rar[[lab]], method = "spearman")))
  # Friedman across the paired compartments (does the effect survive rarefaction?)
  pid <- sub("([CP][0-9]+)[MB]_[CSF]", "\\1", common)
  wn <- reshape(data.frame(pid, compartment, v = nat[[lab]]), idvar="pid", timevar="compartment", direction="wide")
  wr <- reshape(data.frame(pid, compartment, v = rar[[lab]]), idvar="pid", timevar="compartment", direction="wide")
  cat(sprintf("  Friedman p : native = %.4f   rarefied = %.4f\n",
              friedman.test(as.matrix(wn[,-1]))$p.value, friedman.test(as.matrix(wr[,-1]))$p.value))
}
cat("\nInterpretation: Shannon/InvSimpson ~unchanged (depth-robust); Observed drops uniformly under\n")
cat("rarefaction but its compartment ranking + Friedman significance should persist if the native\n")
cat("richness comparison is not explained by sequencing depth alone.\n")
