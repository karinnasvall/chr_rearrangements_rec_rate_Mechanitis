library(dplyr)
library(tidyr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

map_directory <- args[1]
prefix <- args[2]

files <- list.files(map_directory,
                    pattern = "^ordered.*\\.fitted$",
                    full.names = TRUE)

cat("=== Input Structure Inspection ===\n")
cat("Files found:", length(files), "\n\n")

# Inspect first file
if (length(files) > 0) {
  cat("Inspecting first file:", basename(files[1]), "\n")
  
  # Read with header to see structure
  first_file <- read.csv(files[1], header = FALSE, sep = "\t", nrows = 5)
  cat("\nFirst 5 rows:\n")
  print(first_file)
  
  cat("\nDimensions:", nrow(first_file), "rows x", ncol(first_file), "columns\n")
  cat("Column names/classes:\n")
  for (i in 1:ncol(first_file)) {
    cat(sprintf("  V%d: %s\n", i, class(first_file[[i]])))
  }
}

cat("\n=== Reading all fitted files ===\n")
df <- purrr::map_dfr(files, ~ read.csv(.x, header = FALSE, skip = 1, sep = "\t"))

cat("Loaded data: ", nrow(df), "rows x", ncol(df), "columns\n\n")

cat("=== Column Selection ===\n")
cat("Selecting columns V1 (chromosome), V6, V7 (genotypes)\n\n")

df <- df %>% select(V1, V6, V7) %>% distinct() %>% mutate(V6=sub(" .*", "", V6), 
                                                          V7=sub(" .*", "", V7)) 

cat("After selection and cleaning:\n")
print(head(df, 10))
cat("\nDimensions:", nrow(df), "rows x", ncol(df), "columns\n")
cat("Unique chromosomes:", n_distinct(df$V1), "\n")
cat("V6 (first genotype string) sample lengths:\n")
print(table(nchar(df$V6[!is.na(df$V6)])))
cat("V7 (second genotype string) sample lengths:\n")
print(table(nchar(df$V7[!is.na(df$V7)])))

df2 <- df %>%
  transmute(
    chrom = V1,
    geno = V6,
    geno2 = V7
  )

cat("\nGenotype string structure confirmed:\n")
cat("  V6 (Family A):", nchar(df2$geno[1]), "individuals\n")
cat("  V7 (Family B):", nchar(df2$geno2[1]), "individuals\n\n")

# Split genotype string to one column per individual (PATERNAL HAPLOTYPE ONLY = first half)
cat("\nExtracting paternal haplotype (first half of each string)...\n")

# Family A: first 48 of 96
geno_strings_A <- df2$geno
geno_strings_A <- substr(geno_strings_A, 1, nchar(df2$geno[1])/2)
geno_mat <- stringr::str_split(geno_strings_A, "", simplify = TRUE)

cat("Family A (paternal): ", nrow(geno_mat), "x", ncol(geno_mat), "\n")
n_A <- ncol(geno_mat)
colnames(geno_mat) <- sprintf("A_%02d", seq_len(n_A))

# Family B: first 49 of 98
geno_strings_B <- df2$geno2
geno_strings_B <- substr(geno_strings_B, 1, nchar(df2$geno2[1])/2)
geno_mat2 <- stringr::str_split(geno_strings_B, "", simplify = TRUE)

cat("Family B (paternal): ", nrow(geno_mat2), "x", ncol(geno_mat2), "\n")
n_B <- ncol(geno_mat2)
colnames(geno_mat2) <- sprintf("B_%02d", seq_len(n_B))

long_gtA <- as_tibble(geno_mat) %>%
  mutate(chrom = df2$chrom, pos_index = row_number()) %>%
  pivot_longer(-c(chrom, pos_index), names_to = "individual", values_to = "gt") %>%
  mutate(
    gt = case_when(
      gt == "0" ~ 0L,
      gt == "1" ~ 1L,
      TRUE ~ NA_integer_
    ),
    family = "A"
  ) %>%
  select(family, individual, chrom, pos_index, gt)


long_gtB <- as_tibble(geno_mat2) %>%
  mutate(chrom = df2$chrom, pos_index = row_number()) %>%
  pivot_longer(-c(chrom, pos_index), names_to = "individual", values_to = "gt") %>%
  mutate(
    gt = case_when(
      gt == "0" ~ 0L,
      gt == "1" ~ 1L,
      TRUE ~ NA_integer_
    ),
    family = "B"
  ) %>%
  select(family, individual, chrom, pos_index, gt)

long_gt <- rbind(long_gtA, long_gtB)

# Handle missing data (-) represented as NA in integer conversion
long_gt <- long_gt %>%
  mutate(gt = na_if(gt, NA))  # "-" already converts to NA during as.integer()

cat("\n=== Counting Recombination Events ===\n")

count_recombinations <- function(gt_vec) {
  last_known <- NA_integer_
  recombination <- 0L

  for (gt in gt_vec) {
    if (is.na(gt)) {
      next
    }

    if (!is.na(last_known) && gt != last_known) {
      recombination <- recombination + 1L
    }

    last_known <- gt
  }

  recombination
}

# Count recombination events per chromosome and individual.
# Missing data are skipped, and each known genotype is compared with the last
# known genotype within the chromosome/individual series.
rec_counts <- long_gt %>%
  arrange(family, individual, chrom, pos_index) %>%
  group_by(family, individual, chrom) %>%
  summarise(
    recombination = count_recombinations(gt),
    n_markers = sum(!is.na(gt)),
    .groups = "drop"
  )

# function to clean chromosome names (remove prefixes, and zeropad numbers but not letters etc.)
clean_chr <- function(x) {
  # Remove everything up to and including the last underscore
  stripped <- sub(".*_", "", x)
  # Zero-pad 1-digit + letter (e.g. "2A" -> "02A"), leave others as-is
  ifelse(grepl("^[0-9][A-Za-z]$|^[0-9]$", stripped),
         sub("^([0-9])", "0\\1\\2", stripped),
         stripped)
}

rec_counts$chrom <- clean_chr(rec_counts$chrom)

cat("Calculated recombination counts\n")
print(head(rec_counts, 20))


wide_rec <- rec_counts %>%
  select(family, individual, chrom, recombination) %>%
  pivot_wider(names_from = chrom, values_from = recombination, values_fill = 0)

cat("\n=== Wide Format (Recombinations per Individual per Chromosome) ===\n")
print(head(wide_rec, 20))

cat("\nSummary statistics per individual:\n")
summary_stats <- wide_rec %>%
  mutate(
    total = rowSums(select(., -family, -individual)),
    mean = rowMeans(select(., -family, -individual))
  ) %>%
  select(family, individual, total, mean)
print(summary_stats)

write.table(wide_rec, paste0("tables/recombination_summary_wide_", prefix, ".tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(paste0("\nWrote recombination summary to: tables/recombination_summary_wide_", prefix, ".tsv\n"))


cat("\n=== Detailed Summary Statistics ===\n")

sink(paste0("tables/recombination_summary_statistics_", prefix, ".txt"))

cat("\n=== Summary by Chromosome ===\n")
chr_summary <- rec_counts %>%
  group_by(chrom) %>%
  summarise(
    mean_rec = mean(recombination),
    median_rec = median(recombination),
    sd_rec = sd(recombination),
    n_individuals = n(),
    count_0 = sum(recombination == 0),
    count_1 = sum(recombination == 1),
    count_2 = sum(recombination == 2),
    count_3 = sum(recombination == 3),
    count_4_plus = sum(recombination >= 4),
    .groups = "drop"
  ) %>%
  arrange(chrom)

print(chr_summary, n = 50)

cat("\n=== Summary by Family and Individual ===\n")
indiv_summary <- rec_counts %>%
  group_by(family, individual) %>%
  summarise(
    total_recomb = sum(recombination),
    mean_rec_per_chrom = mean(recombination),
    n_chrom = n(),
    .groups = "drop"
  ) %>%
  arrange(family, individual)

print(indiv_summary, n = 200)

cat("\n=== Distribution: Count of individuals with N total recombinations ===\n")
distribution <- indiv_summary %>%
  group_by(family, total_recomb) %>%
  summarise(n_individuals = n(), .groups = "drop") %>%
  arrange(family, total_recomb)

print(distribution, n = 200)

sink()

cat(paste0("Wrote summary statistics to: tables/recombination_summary_statistics_", prefix, ".txt\n"))

# ============================================================================
# Optional: Create plots if ggplot2 is available
# ============================================================================

if (require(ggplot2, quietly = TRUE)) {
  cat("\nGenerating plots...\n")
  
  # Histogram per individual (per chromosome)
  p_indiv <- ggplot(rec_counts, aes(recombination, fill = family)) +
    geom_histogram(binwidth = 0.5, color = "black", alpha = 0.7) +
    facet_wrap(~family + individual, ncol = 12) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 7),
      strip.text = element_text(size = 7),
      legend.position = "top"
    ) +
    labs(title = "Recombination count distribution per individual across chromosomes",
         x = "Number of recombination events",
         y = "Number of chromosomes")
  
  ggsave(paste0("plots/recombination_histogram_per_individual_", prefix, ".pdf"), p_indiv, width = 20, height = 24)
  cat(paste0("Wrote plot to: plots/recombination_histogram_per_individual_", prefix, ".pdf\n"))
  
  # Distribution of recombination counts per chromosome
  p1 <- ggplot(rec_counts, aes(recombination, fill = family)) +
    geom_histogram(binwidth = 1, color = "black", alpha = 0.7) +
    facet_wrap(~chrom) +
    theme_bw() +
    labs(title = "Distribution of recombination counts per chromosome",
         x = "Number of recombination events",
         y = "Count")
  
  ggsave(paste0("plots/recombination_histogram_per_chromosome_", prefix, ".pdf"), p1, width = 14, height = 10)
  cat(paste0("Wrote plot to: plots/recombination_histogram_per_chromosome_", prefix, ".pdf\n"))
  
  # Total recombinations per chromosome
  p2 <- rec_counts %>%
    group_by(chrom, family) %>%
    summarise(total_recomb = sum(recombination), .groups = "drop") %>%
    ggplot(aes(chrom, total_recomb, fill = family)) +
    geom_bar(stat = "identity", position = "dodge", color = "black") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Total recombination events per chromosome",
         x = "Chromosome",
         y = "Total recombination count")
  
  ggsave(paste0("plots/recombination_total_per_chromosome_", prefix, ".pdf"), p2, width = 10, height = 6)
  cat(paste0("Wrote plot to: plots/recombination_total_per_chromosome_", prefix, ".pdf\n"))
}
