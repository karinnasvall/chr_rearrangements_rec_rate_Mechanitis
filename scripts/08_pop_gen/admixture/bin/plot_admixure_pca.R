#!/usr/bin/env Rscript
# Usage: Rscript plot_pca_admixture.R [prefix] [K] [outdir]
#   prefix : PLINK/ADMIXTURE file prefix (default "pruned")
#   K      : number of ancestral populations to plot (default: lowest CV error)
#   outdir : output directory (default "plots")
 
suppressPackageStartupMessages(library(tidyverse))
 
args   <- commandArgs(trailingOnly = TRUE)
prefix <- if (length(args) >= 1) args[1] else "output/"
K      <- if (length(args) >= 2) as.integer(args[2]) else NA_integer_
outdir <- if (length(args) >= 3) args[3] else "plots"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
 
fam <- read_table(paste0(prefix, "pruned.fam"), col_names = FALSE,
                  show_col_types = FALSE)$X2
 
## ---- PCA -------------------------------------------------------------
 
# Metadata: tab-delimited, columns "sample" and "population"
meta <- read_tsv("metadata.tsv", show_col_types = FALSE, col_names = c("sample", "population"))

# PLINK 2 eigenvec has a header; "#FID"/"#IID" are cleaned to "FID"/"IID"
pca <- read_tsv(paste0(prefix, "pca.eigenvec"), show_col_types = FALSE) %>%
  rename_with(~ sub("^#", "", .x)) %>%
  rename(sample = IID) %>%
  left_join(meta, by = "sample")

# Check that every sample was matched
if (any(is.na(pca$population))) {
  warning("Samples without population: ",
          paste(pca$sample[is.na(pca$population)], collapse = ", "))
}

eig <- read_table(paste0(prefix, "pca.eigenval"), col_names = "ev", show_col_types = FALSE)$ev
pve <- round(100 * eig / sum(eig), 1)

p_pca <- ggplot(pca, aes(PC1, PC2, colour = population)) +
  geom_point(size = 2) +
  scale_colour_brewer(palette = "Set1", na.value = "grey60") +
  labs(x = paste0("PC1 (", pve[1], "%)"), y = paste0("PC2 (", pve[2], "%)"),
       colour = "Population") +
  theme_classic()
ggsave(file.path(outdir, "pca.pdf"), p_pca, width = 6, height = 5)

## ---- Cross-validation error -------------------------------------------
cv <- read_table(paste0(prefix, "cv.txt"), col_names = c("K", "cv"), show_col_types = FALSE)
 
p_cv <- ggplot(cv, aes(K, cv)) +
  geom_line() + geom_point() +
  scale_x_continuous(breaks = cv$K) +
  labs(x = "K", y = "CV error") +
  theme_classic()
ggsave(file.path(outdir, "cv_error.pdf"), p_cv, width = 5, height = 4)
 
if (is.na(K)) K <- cv$K[which.min(cv$cv)]
message("Plotting K = ", K)
 
## ---- ADMIXTURE barplot ------------------------------------------------
q_wide <- read_table(paste0(prefix, "pruned.", K, ".Q"),
                     col_names = paste0("anc", seq_len(K)),
                     show_col_types = FALSE) %>%
  mutate(sample = fam)
 
# Order samples by dominant ancestry, then by its proportion
q_wide <- q_wide %>%
  mutate(top = max.col(select(., starts_with("anc")), ties.method = "first"),
         topprop = apply(select(., starts_with("anc")), 1, max)) %>%
  arrange(top, desc(topprop)) %>%
  mutate(sample = factor(sample, levels = sample))
 
q <- q_wide %>%
  pivot_longer(starts_with("anc"), names_to = "ancestry", values_to = "prop")
 
p_adm <- ggplot(q, aes(sample, prop, fill = ancestry)) +
  geom_col(width = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Ancestry proportion", fill = NULL,
       title = paste0("K = ", K)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))
ggsave(file.path(outdir, paste0("admixture_K", K, ".pdf")), p_adm,
       width = max(6, 0.15 * nrow(q_wide)), height = 3.5, limitsize = FALSE)
 
