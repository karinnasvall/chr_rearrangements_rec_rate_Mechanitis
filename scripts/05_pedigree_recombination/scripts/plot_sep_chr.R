#!/usr/bin/env Rscript   

# script to plot the number of markers assigned to each LG with different LodLimit values in the separate chromosomes step, 
# summarise the ouput from SeparateChromosomes first with create_lod_summary.sh

library(ggplot2)
library(viridis)
library("tidyr")
library(dplyr)

# arguments

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: plot_sep_chrom.R <maps_sum_file> <output_file>", call. = FALSE)
}

maps_sum_file <- args[1]
output_file <- args[2]

maps_sum <- read.csv(maps_sum_file, sep = "\t", header = T, fill=T)[,-1]

# Pivot to long and parse the "markers LG chrname" cell values
maps_long <- maps_sum %>%
  mutate(rank = row_number()) %>%
  pivot_longer(-rank, names_to = "lod_lim", values_to = "cell") %>%
  mutate(lod_lim = gsub("lodLimit_", "", lod_lim)) %>%
  filter(!is.na(cell) & cell != "") %>%  
  mutate(cell = trimws(cell)) %>%   
  separate(cell, into = c("nr_markers", "lg", "chr_name"), sep = "\\s+", extra = "merge") %>%
  mutate(nr_markers = as.numeric(nr_markers),
         lg         = as.numeric(lg))

# Tile label: markers on top, LG + chr below
maps_long <- maps_long %>%
  mutate(label = paste0(nr_markers, "\nLG", lg, " ", chr_name))


ggplot(maps_long, aes(x = lod_lim, y = rank)) +
  geom_tile(aes(fill = lg), colour = "grey80") +
  geom_text(aes(label = label), size = 2, lineheight = 0.85) +
  scale_fill_gradient(low = "yellow", high = "blue", name = "LG",
                       na.value = "grey90") +
  scale_x_discrete(position = "top") +
  scale_y_reverse() +
  labs(x = "LOD limit", y = "Rank") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 15),
  axis.text.y = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_blank())

ggsave(output_file, width = 10, height = 12)

