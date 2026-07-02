#!/usr/bin/env Rscript
# ============================================================================
# Independent replication of the CoralShot / CORAL_P GENUS relative abundance
# via phyloseq (the same numbers make_stacked_composition.py plots).
#
# phyloseq does NOT recompute abundance from reads -- it loads the MetaPhlAn
# table and reports it. So it reproduces our relative abundance EXACTLY, as
# long as two things are handled (both are common sources of a false mismatch):
#
#   (1) merged_sgb.tsv mixes ALL taxonomic ranks as rows -> filter to ONE rank.
#       Here: genus = rows whose clade_name has |g__ and NOT |s__.
#       (Loading all rows and normalizing double-counts every read at every rank.)
#   (2) MetaPhlAn values are ALREADY relative abundance (%). We close each sample
#       to sum = 1 over the chosen rank  ==  the figure's  g / colSums(g).
#
# Run (phyloseq env on Torch):
#   BASE=/scratch/sr7729/CORAL_P \
#   /scratch/sr7729/conda_envs/cs_profile/bin/Rscript scripts/phyloseq_relabund.R
# ============================================================================
.libPaths(c("/scratch/sr7729/gut/coral_reef/R_libs", .libPaths()))
suppressMessages(library(phyloseq))

BASE  <- Sys.getenv("BASE", "/scratch/sr7729/CORAL_P")
mpa   <- read.delim(file.path(BASE, "results/04_mpa/merged_sgb.tsv"), skip = 1, check.names = FALSE)
colnames(mpa) <- sub("\\.profile$", "", colnames(mpa))
clade <- as.character(mpa[[1]])

# (1) GENUS-level rows only
is_g  <- grepl("\\|g__", clade) & !grepl("\\|s__", clade)
gname <- sub(".*g__", "", clade[is_g])
samp  <- grep("^[CP][0-9]+[MB]_[CSF]$", colnames(mpa), value = TRUE)
otu   <- as.matrix(mpa[is_g, samp]); rownames(otu) <- gname
storage.mode(otu) <- "double"

# build the phyloseq object; (2) close each sample to relative abundance (sum = 1)
ps    <- phyloseq(otu_table(otu, taxa_are_rows = TRUE))
psrel <- transform_sample_counts(ps, function(x) x / sum(x))

relab <- as.data.frame(as(otu_table(psrel), "matrix"))
out   <- file.path(BASE, "results/04_mpa/phyloseq_genus_relabund.tsv")
write.table(cbind(genus = rownames(relab), round(relab, 8)), out,
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("phyloseq: %d genera x %d samples -> %s\n", nrow(relab), ncol(relab), out))
cat(sprintf("sample column sums (all should be 1.0): %s ...\n",
            paste(round(head(colSums(relab), 6), 6), collapse = " ")))
