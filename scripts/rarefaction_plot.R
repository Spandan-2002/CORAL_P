#!/usr/bin/env Rscript
# ============================================================================
# Visualize the read-level rarefaction sensitivity check for alpha diversity:
#   (a) figures/rarefaction_alpha_box.png     — per-compartment boxplots, native
#       depth vs rarefied-to-13M, one panel per metric (shows the ranking holds
#       and how little Shannon/InvSimpson move while Observed drops uniformly).
#   (b) figures/rarefaction_alpha_scatter.png — per-sample native vs rarefied
#       scatter (diagonal = no change), one panel per metric, Spearman rho annotated.
# Reads the same native + rarefied merged tables as rarefaction_alpha.R.
# Run: BASE=/scratch/sr7729/CORAL_P \
#      /scratch/sr7729/conda_envs/cs_stats/bin/Rscript scripts/rarefaction_plot.R
# ============================================================================
suppressMessages({ library(vegan); library(ggplot2) })
BASE  <- Sys.getenv("BASE", "/scratch/sr7729/CORAL_P")
CLEAN <- c("P01","P02","P03","P04","P05","P08","P09","P10","P12","P13")

load_species <- function(path) {
  m  <- read.delim(path, skip = 1, check.names = FALSE)
  colnames(m) <- sub("\\.profile$", "", colnames(m))
  cl <- as.character(m[[1]])
  keep <- grepl("\\|s__", cl) & !grepl("\\|t__", cl)
  samp <- grep("^[CP][0-9]+[MB]_[CSF]$", colnames(m), value = TRUE)
  x <- as.matrix(m[keep, samp]); storage.mode(x) <- "double"
  subj <- sub("([CP][0-9]+)[MB]_[CSF]", "\\1", samp)
  x[, subj %in% CLEAN, drop = FALSE]
}
alpha <- function(x) {
  p <- t(x)
  data.frame(sample = rownames(p),
             Observed   = specnumber(p),
             Shannon    = diversity(p, "shannon"),
             InvSimpson = diversity(p, "invsimpson"))
}
nat <- alpha(load_species(file.path(BASE, "results/04_mpa/merged_sgb.tsv")))
rar <- alpha(load_species(file.path(BASE, "results/rarefaction/merged_rare13M.tsv")))
common <- intersect(nat$sample, rar$sample)
nat <- nat[match(common, nat$sample), ]; rar <- rar[match(common, rar$sample), ]
comp <- factor(c(C="core", S="shell", F="stool")[substr(common, nchar(common), nchar(common))],
               levels = c("core","shell","stool"))
mets <- c("Observed","Shannon","InvSimpson")

# ---- long frame for boxplots (native vs rarefied) --------------------------
long <- do.call(rbind, lapply(mets, function(m) rbind(
  data.frame(compartment=comp, metric=m, depth="native (full)",   value=nat[[m]]),
  data.frame(compartment=comp, metric=m, depth="rarefied (13M)",  value=rar[[m]]))))
long$metric <- factor(long$metric, levels = mets)

pbox <- ggplot(long, aes(compartment, value, fill = depth)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.8), width = 0.7, alpha = 0.6) +
  geom_point(aes(color = depth), position = position_dodge(0.8), size = 1.1) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_manual(values = c("native (full)"="#4C72B0","rarefied (13M)"="#DD8452")) +
  scale_color_manual(values = c("native (full)"="#31517d","rarefied (13M)"="#a85c34")) +
  labs(title = "Alpha diversity is depth-robust — native depth vs read-level rarefaction to 13M reads (clean-10)",
       subtitle = "Each metric's compartment ranking is preserved under rarefaction (shell highest throughout; core~stool, n.s.); Shannon/InvSimpson barely move, Observed drops uniformly",
       x = NULL, y = "alpha diversity", fill = NULL, color = NULL) +
  theme_bw() + theme(legend.position = "top",
                     plot.subtitle = element_text(size = 9),
                     strip.text = element_text(face = "bold"))
ggsave(file.path(BASE, "figures/rarefaction_alpha_box.png"), pbox, width = 11, height = 4.5, dpi = 150)
cat("wrote figures/rarefaction_alpha_box.png\n")

# ---- per-sample scatter native vs rarefied (with Spearman rho) --------------
sc <- do.call(rbind, lapply(mets, function(m)
  data.frame(compartment=comp, metric=m, native=nat[[m]], rarefied=rar[[m]])))
sc$metric <- factor(sc$metric, levels = mets)
rho <- sapply(mets, function(m) cor(nat[[m]], rar[[m]], method = "spearman"))
labs_rho <- data.frame(metric = factor(mets, levels = mets),
                       lab = sprintf("Spearman rho = %.3f", rho))

psc <- ggplot(sc, aes(native, rarefied, color = compartment)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
  geom_point(size = 1.6, alpha = 0.85) +
  facet_wrap(~ metric, scales = "free") +
  geom_text(data = labs_rho, aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE,
            hjust = -0.08, vjust = 1.6, size = 3.2) +
  scale_color_manual(values = c(core="#55A868", shell="#C44E52", stool="#8172B3")) +
  labs(title = "Per-sample alpha: native depth vs rarefied to 13M (points on the dashed line = no change)",
       x = "native (full depth)", y = "rarefied (13M reads)", color = NULL) +
  theme_bw() + theme(legend.position = "top", strip.text = element_text(face = "bold"))
ggsave(file.path(BASE, "figures/rarefaction_alpha_scatter.png"), psc, width = 11, height = 4.5, dpi = 150)
cat("wrote figures/rarefaction_alpha_scatter.png\n")

# ---- dedicated Observed-RICHNESS view: paired native -> rarefied per library --
rich <- rbind(
  data.frame(sid = common, compartment = comp, depth = "native (full)",  richness = nat$Observed),
  data.frame(sid = common, compartment = comp, depth = "rarefied (13M)", richness = rar$Observed))
rich$depth <- factor(rich$depth, levels = c("native (full)","rarefied (13M)"))
medlab <- aggregate(richness ~ compartment + depth, rich, function(v) round(median(v)))
prich <- ggplot(rich, aes(depth, richness, group = sid)) +
  geom_line(color = "grey72", linewidth = 0.4) +
  geom_point(aes(color = compartment), size = 2) +
  geom_text(data = medlab, aes(depth, richness, label = richness), inherit.aes = FALSE,
            vjust = -0.9, size = 3, fontface = "bold") +
  facet_wrap(~ compartment) +
  scale_color_manual(values = c(core="#55A868", shell="#C44E52", stool="#8172B3")) +
  labs(title = "Observed species richness — native depth vs rarefied to 13M (each line = one library; bold = median)",
       subtitle = "Richness falls when depth is equalized (expected: fewer reads detect fewer rare taxa), but the drop is ~uniform so the compartment ordering is preserved (Friedman p=0.045 native & rarefied; core~stool n.s.)",
       x = NULL, y = "observed species / SGBs", color = NULL) +
  theme_bw() + theme(legend.position = "none", plot.subtitle = element_text(size = 8),
                     strip.text = element_text(face = "bold"))
ggsave(file.path(BASE, "figures/rarefaction_richness.png"), prich, width = 10, height = 4.2, dpi = 150)
cat("wrote figures/rarefaction_richness.png\n")

# ---- richness BOXPLOT: Observed by compartment, native vs rarefied ----
pbox_rich <- ggplot(rich, aes(compartment, richness, fill = depth)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.8), width = 0.7, alpha = 0.6) +
  geom_point(aes(color = depth), position = position_dodge(0.8), size = 1.4) +
  scale_fill_manual(values = c("native (full)"="#4C72B0","rarefied (13M)"="#DD8452")) +
  scale_color_manual(values = c("native (full)"="#31517d","rarefied (13M)"="#a85c34")) +
  labs(title = "Observed species richness by compartment — native depth vs rarefied to 13M",
       subtitle = "Boxplots across the 10 clean-10 libraries (points = libraries); richness drops ~uniformly under rarefaction; shell highest, core~stool (n.s.)",
       x = NULL, y = "observed species / SGBs", fill = NULL, color = NULL) +
  theme_bw() + theme(legend.position = "top", plot.subtitle = element_text(size = 8.5))
ggsave(file.path(BASE, "figures/rarefaction_richness_box.png"), pbox_rich, width = 8, height = 4.6, dpi = 150)
cat("wrote figures/rarefaction_richness_box.png\n")
