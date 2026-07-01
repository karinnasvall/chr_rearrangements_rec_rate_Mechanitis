#!/usr/bin/env Rscript 

#Rscript to eveluate map order agaist physical order from lapmap3 maps

library(dplyr)
library(ggplot2)
library(ggpubr)

args <- commandArgs(trailingOnly = TRUE)

mapdirectory <- args[1]
prefix <- args[2]
#for chrom lengths
fai_file <- args[3]


# or a file with a list of files to read in, one per line
if (file.exists(mapdirectory) && file.info(mapdirectory)$isdir) {
  print(paste("Reading maps from directory:", mapdirectory))
  files <- list.files(mapdirectory, pattern = paste0("ordered.*.fitted$"), full.names = TRUE)  
} else if (file.exists(mapdirectory)) {
  print(paste("Reading maps from file list:", mapdirectory))
  files <- readLines(mapdirectory)
} else {
  stop("Provided mapdirectory is neither a directory nor a file.")
}

clean_chr <- function(x) {
  # Remove everything up to and including the last underscore
  stripped <- sub(".*_", "", x)
  # Zero-pad numbers (e.g. "1" -> "01"), leave letters as-is
  ifelse(grepl("^\\d+$", stripped),
         sprintf("%02d", as.integer(stripped)),
         stripped)
}


print(paste("Using prefix for plots and tables:", prefix))

print(files)

# readin the files and combine to one df give the number in the filename as a column, 
# keep only the first 4 columns and rename them to marker, position, genetic_position_male, genetic_position_female and add LG column with the number in the filename
maps <- lapply(files, function(file) {
  map <- read.table(file, header = TRUE, stringsAsFactors = FALSE, skip = 3)
  map <- map[, 1:4]
  colnames(map) <- c("chrom", "position", "genetic_position_male", "genetic_position_female")
  map$LG <- as.numeric(gsub(".*ordered\\.([0-9]+)\\.fitted", "\\1", file))
  map$position <- as.numeric(gsub("\\*", "", map$position))
  return(map)
})  %>% bind_rows()

str(maps)


maps$chrom <- clean_chr(maps$chrom)
maps$chrom <- factor(maps$chrom)

maps$LG <- as.factor(maps$LG)

str(maps)

### get chr lengths
print(paste("Reading .fai file:", fai_file))
chr_lengths <- read.table(fai_file, header = FALSE, stringsAsFactors = FALSE)[, 1:2]
colnames(chr_lengths) <- c("chrom", "length")


# clean chromosme names to be just the number and zeropadded if needed
chr_lengths$chrom <- clean_chr(chr_lengths$chrom)
chr_lengths$chrom <- factor(chr_lengths$chrom)

print("Chromosome lengths:")
print(chr_lengths)

# add to df
maps <- maps %>% left_join(chr_lengths, by = "chrom")


# add supprted markers

files_sup <- list.files(mapdirectory, pattern = paste0("marker_support_*"), full.names = TRUE)

# read in files and add to maps df
mark_sup <- lapply(files_sup, function(file) {
  map <- read.table(file, header = FALSE, stringsAsFactors = FALSE, skip = 1)
  colnames(map) <- c("chrom", "position", "support")
  map$LG <- as.numeric(gsub(".*marker_support_([0-9]+)\\.txt", "\\1", file))
  map$position <- as.numeric(gsub("\\*", "", map$position))
  return(map)
})  %>% bind_rows()

head(mark_sup)

mark_sup$chrom <- clean_chr(mark_sup$chrom)
mark_sup$chrom <- factor(mark_sup$chrom)
mark_sup$LG <- as.factor(mark_sup$LG)

sup_summary <- mark_sup %>% filter(support > 0.5) %>% group_by(LG, chrom) %>% summarise(supported_markers = n())
head(sup_summary)

plot_title=prefix
# print information about the maps df
print(paste("Result for:", plot_title))

#Number of markers in map
print(paste("Number of markers in LG:", length(maps$position)))   
apply(maps[c("LG")], 2, table)
apply(maps[c("chrom")], 2, table)

#maplength
print("Map length per LG min and max")
#cbind(aggregate(map_ordered$distance_min, list(map_ordered$lg), max), aggregate(map_ordered$distance_max, list(map_ordered$lg), max))
aggregate(maps$genetic_position_male, list(maps$LG), max)
paste("Total map length (cM): ",sum(aggregate(maps$genetic_position_male, list(maps$LG), max)$x))

maps_auto <- maps %>% filter(!chrom %in%  c("Z", "Z1", "Z2","Z3","W", "W1", "W2", "W3", "W4", "W5")) 

paste("Total autosomal map length (cM): ", sum(aggregate(maps_auto$genetic_position_male, list(maps_auto$LG), max)$x))

paste("Appr autosomal recombination rate (cM/Mb):", sum(aggregate(maps_auto$genetic_position_male, list(maps_auto$LG), max)$x)/sum(aggregate(maps_auto$length, list(maps_auto$chrom), max)$x)*1000000)

# print as table
summary_maps <- maps %>% group_by(LG, chrom) %>% summarise(map_length = max(genetic_position_male), physical_length = max(length), markers = n(), rec_rate = map_length / physical_length * 1000000) %>% filter(markers > 100) %>% arrange(desc(physical_length))

summary_maps <- summary_maps %>% left_join(sup_summary, by = c("LG", "chrom")) 

write.table(summary_maps, file = paste0("tables/summary_ordered_fitted_", prefix, ".txt"), quote = FALSE, row.names = FALSE, sep = "\t")
str(maps)

# plot the number of markers per chrom and per LG
ggplot(maps, aes(x = chrom, fill = LG)) +
  geom_bar(stat = "count", position = "stack") +
  labs(title = plot_title) +
  theme(panel.background = element_blank(),
        axis.line = element_line(linewidth = 1),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10))

ggsave(paste0("plots/marker_distribution_scaffolds_fitted_", prefix, ".png"), width = 10, height = 5)

ggplot(maps, aes(LG, fill=chrom)) +
  geom_bar(stat = "count", position = "stack") +
  labs(title = plot_title) +
  theme(panel.background = element_blank(),
        axis.line = element_line(linewidth = 1),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10))

ggsave(paste0("plots/marker_distribution_lg_fitted", prefix, ".png"), width = 10, height = 5)


#in one plot by lg
ggplot(maps, aes(position, genetic_position_male)) +
  geom_point(aes(colour=LG)) +
  #geom_smooth(se = FALSE) +
  facet_wrap(~chrom) +
  labs(title = plot_title) +
  xlab("Position (bp)") +
  ylab("Distance (cM)") +
  theme_classic()

ggsave(paste0("plots/physical_vs_genetic_position_fitted_", prefix, ".png"), width = 15, height = 10)


#plot rec_rate per chrom vs chrom length
ggplot(summary_maps, aes(x = physical_length/1000000, y = rec_rate)) +
  geom_point(size = 5, shape = 21) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_line(aes(x = physical_length/1000000, y = (50/physical_length)*1000000), linetype = "dashed", color = "tan1") +
  geom_line(aes(x = physical_length/1000000, y = (100/physical_length)*1000000), linetype = "dashed", color = "tan4") +
  stat_cor(method = "spearman", cor.coef.name = "rho", label.x = 3, label.y = max(summary_maps$rec_rate) * 0.9) +
  geom_text(aes(label = chrom), size = 3) +
  labs(title = paste("Recombination rate vs chromosome length for", plot_title),
       x = "Chromosome length (Mb)",
       y = "Recombination rate (cM/Mb)") +
  ylim(0, max(summary_maps$rec_rate) * 1.3) +
  theme_classic()

ggsave(paste0("plots/rec_rate_vs_chrom_length_fitted_", prefix, ".png"), width = 10, height = 5)