#!/usr/bin/env Rscript 

#Rscript to eveluate map order agaist physical order from lapmap3 maps
# the ordered maps needs to be mapped with physical position (map_physical_pos_to_ordered_maps.sh)
# Usage: Rscript lepmap_plot_order.R <map_directory_or_file_list> <prefix_for_plots_and_tables>

library(dplyr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

mapdirectory <- args[1]
# or a file with a list of files to read in, one per line
if (file.exists(mapdirectory) && file.info(mapdirectory)$isdir) {
  print(paste("Reading maps from directory:", mapdirectory))
  files <- list.files(mapdirectory, pattern = paste0("ordered.*_mapped$"), full.names = TRUE)  
} else if (file.exists(mapdirectory)) {
  print(paste("Reading maps from file list:", mapdirectory))
  files <- readLines(mapdirectory)
} else {
  stop("Provided mapdirectory is neither a directory nor a file.")
}

prefix <- args[2]

print(paste("Using prefix for plots and tables:", prefix))

print(files)

# readin the files and combine to one df give the number in the filename as a column, 
# keep only the first 4 columns and rename them to marker, position, genetic_position_male, genetic_position_female and add LG column with the number in the filename
maps <- lapply(files, function(file) {
  map <- read.table(file, header = TRUE, stringsAsFactors = FALSE, skip = 2)
  map <- map[, 1:4]
  colnames(map) <- c("chrom", "position", "genetic_position_male", "genetic_position_female")
  map$LG <- as.numeric(gsub(".*ordered\\.([0-9]+).*", "\\1", file))
  map$position <- as.numeric(gsub("\\*", "", map$position))
  return(map)
})  %>% bind_rows()

str(maps)
# clean chromosme names to be just the number and zeropadded if needed
# substitute everything before _

clean_chr <- function(x) {
  # Remove everything up to and including the last underscore
  stripped <- sub(".*_", "", x)
  # Zero-pad numbers (e.g. "1" -> "01"), leave letters as-is
  ifelse(grepl("^\\d+$", stripped),
         sprintf("%02d", as.integer(stripped)),
         stripped)
}

maps$chrom <- clean_chr(maps$chrom)
unique(maps$chrom)

maps$chrom <- factor(maps$chrom)

maps$LG <- as.factor(maps$LG)

str(maps)

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

paste("Appr recombination rate (cM/Mb):", sum(aggregate(maps$genetic_position_male, list(maps$LG), max)$x)/sum(aggregate(maps$position, list(maps$chrom), max)$x)*1000000)

# print as table, note this is by lg and chrom so not the total size if markers are on multiple chromsomes per lg
summary_maps <- maps %>% group_by(LG, chrom) %>% summarise(map_length = max(genetic_position_male), physical_length = max(position), markers = n()) %>% filter(markers > 100) %>% arrange(desc(physical_length))

write.table(summary_maps, file = paste0("tables/summary_ordered_", prefix, ".txt"), quote = FALSE, row.names = FALSE, sep = "\t")
str(maps)

# plot the number of markers per chrom and per LG
ggplot(maps, aes(x = chrom, fill = LG)) +
  geom_bar(stat = "count", position = "stack") +
  labs(title = plot_title) +
  theme(panel.background = element_blank(),
        axis.line = element_line(size = 1),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10))

ggsave(paste0("plots/marker_distribution_scaffolds_", prefix, ".png"), width = 10, height = 5)

ggplot(maps, aes(LG, fill=chrom)) +
  geom_bar(stat = "count", position = "stack") +
  labs(title = plot_title) +
  theme(panel.background = element_blank(),
        axis.line = element_line(size = 1),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10))

ggsave(paste0("plots/marker_distribution_lg_", prefix, ".png"), width = 10, height = 5)


#in one plot by lg
ggplot(maps, aes(position, genetic_position_male)) +
  geom_point(aes(colour=chrom)) +
  #geom_smooth(se = FALSE) +
  facet_wrap(~LG) +
  guides(colour=FALSE) +
  labs(title = plot_title) +
  xlab("Position (bp)") +
  ylab("Distance (cM)") +
  theme_classic()

ggsave(paste0("plots/physical_vs_genetic_position_", prefix, ".png"), width = 15, height = 10)
