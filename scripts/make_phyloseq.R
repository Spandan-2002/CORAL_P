#!/usr/bin/env Rscript
# ============================================================================
# Build a phyloseq object from the CoralShot/CORAL_P MetaPhlAn table and make a
# phyloseq-native stacked composition figure.
#   object : species/SGB-level otu_table + full tax_table (Kingdom..Species) +
#            sample_data (group/subject/arm/compartment)  -> results/phyloseq/coralp_ps.rds
#   figure : tax_glom to Genus -> relative abundance -> plot_bar (top-14 + Other)
#            -> figures/stacked_composition_phyloseq.png
# Run: BASE=/scratch/sr7729/CORAL_P \
#      /scratch/sr7729/conda_envs/cs_profile/bin/Rscript scripts/make_phyloseq.R
# ============================================================================
.libPaths(c("/scratch/sr7729/gut/coral_reef/R_libs", .libPaths()))
suppressMessages({ library(phyloseq); library(ggplot2) })
BASE  <- Sys.getenv("BASE", "/scratch/sr7729/CORAL_P")
mpa   <- read.delim(file.path(BASE, "results/04_mpa/merged_sgb.tsv"), skip = 1, check.names = FALSE)
colnames(mpa) <- sub("\\.profile$", "", colnames(mpa))
clade <- as.character(mpa[[1]])
samp  <- grep("^[CP][0-9]+[MB]_[CSF]$", colnames(mpa), value = TRUE)

## --- species/SGB rows (|s__, not |t__) -> otu_table + tax_table ---------------
is_s  <- grepl("\\|s__", clade) & !grepl("\\|t__", clade)
otu   <- as.matrix(mpa[is_s, samp]); storage.mode(otu) <- "double"
ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
pref  <- c("k__","p__","c__","o__","f__","g__","s__")
parse_lineage <- function(cl){
  parts <- strsplit(cl, "\\|")[[1]]; out <- setNames(rep(NA_character_, 7), ranks)
  for (i in seq_along(pref)) { h <- grep(pref[i], parts, value = TRUE); if (length(h)) out[i] <- sub(pref[i], "", h[1]) }
  out
}
tax   <- t(sapply(clade[is_s], parse_lineage)); rownames(tax) <- NULL
taxid <- make.unique(tax[, "Species"]); rownames(otu) <- taxid; rownames(tax) <- taxid

## --- sample_data from the committed manifest ---------------------------------
sm <- read.delim(file.path(BASE, "sample_manifest.tsv"), check.names = FALSE)
sm <- sm[match(samp, sm$sample_id), ]; rownames(sm) <- samp
sm$compartment <- c(C = "core", S = "shell", F = "stool")[substr(samp, nchar(samp), nchar(samp))]

ps <- phyloseq(otu_table(otu, taxa_are_rows = TRUE),
               tax_table(as.matrix(tax)),
               sample_data(sm))
dir.create(file.path(BASE, "results/phyloseq"), showWarnings = FALSE, recursive = TRUE)
saveRDS(ps, file.path(BASE, "results/phyloseq/coralp_ps.rds"))
cat(sprintf("saved: %d species x %d samples -> results/phyloseq/coralp_ps.rds\n", ntaxa(ps), nsamples(ps)))

## --- phyloseq-native genus composition figure --------------------------------
psg  <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
psgr <- transform_sample_counts(psg, function(x) x / sum(x))               # relative abundance
top  <- names(sort(taxa_sums(psgr), decreasing = TRUE))[1:14]
tt   <- as.data.frame(as(tax_table(psgr), "matrix"))
tt$Genus[!(taxa_names(psgr) %in% top)] <- "Other"
tax_table(psgr) <- tax_table(as.matrix(tt))
psgr <- tax_glom(psgr, taxrank = "Genus", NArm = FALSE)                     # lump the "Other"

p <- plot_bar(psgr, fill = "Genus") +
  facet_grid(~ compartment, scales = "free_x", space = "free") +
  labs(title = "Genus composition (phyloseq) — relative abundance",
       x = NULL, y = "relative abundance") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 4, hjust = 1),
                     legend.text = element_text(size = 7))
ggsave(file.path(BASE, "figures/stacked_composition_phyloseq.png"), p, width = 16, height = 7, dpi = 150)
cat("wrote figures/stacked_composition_phyloseq.png\n")
