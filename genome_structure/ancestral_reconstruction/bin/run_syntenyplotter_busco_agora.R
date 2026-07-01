#!/usr/bin/env Rscript --vanilla

#script by Karin Näsvall using syntenyPlotter
#https://github.com/Farre-lab/syntenyPlotteR
#This script uses the output from busco
#Creates chain files for syntenyPlotter.

#install.packages("devtools")
#devtools::install_github("marta-fb/syntenyPlotteR")
#library(syntenyPlotteR)
#install.packages("pals")
#library(pals)
library(tidyverse)
library(ggplot2)

rm(list = ls())

#dir organisation: intermediate/, plots/


######
####If using command line and argument uncomment these (and comment away the variable lines below):
cmd_args <- commandArgs(trailingOnly = TRUE)
table_file <- cmd_args[1]
align_length=cmd_args[2]
#taxa1=cmd_args[3]
#taxa2=cmd_args[4]

# table file must be in the format dir/taxa1_taxa2.paf
taxa1=str_split(sub(".paf", "", sub(".*/","", table_file)), pattern = "_")[[1]][1]
taxa2=str_split(sub(".paf", "", sub(".*/","", table_file)), pattern = "_")[[1]][2]

#################
# Variables if using hardcoding, change these variables:
#output from minimap2 paf-format

#################
output_path="/lustre/scratch125/tol/teams/meier/users/kn9/recombination/breakpoints/anc_reconstruction_pipeline_260609/analysis_break/"
table_file <- paste(output_path, "/tables/agora_link_parsed_flipped.csv", sep = "")

busco_table <- read.csv2(table_file, sep = ",", header = T)

colnames(busco_table)
 [1] "chr"         "start"       "end"         "strand"      "agora_node" 
 [6] "geneid"      "length"      "anc_element" "n"           "Node"       
[11] "node"        "length_max"  "start_init"  "end_init"  

# trying out a first ordering file chr_size file for syntenyplotter
busco_table %>% select(Node, chr, start, geneid, length) %>% filter(Node %in% c("n1", "n4", "ilMecPoly1", "ilMecLysi212")) %>% 
  pivot_wider(names_from = Node, values_from = c(chr, start, length)) %>% arrange(chr_ilMecPoly1, start_ilMecPoly1) %>% 
  pivot_longer(cols = -geneid, names_to = c(".value","Node"), names_sep = "_") %>% 
  select(chr, length, Node) %>% distinct() %>% arrange(factor(Node, levels = c("ilMecPoly1", "ilMecLysi212", "n4", "n1"))) %>%
  mutate(heterozygous = "N") %>% filter(!is.na(chr)) %>% 
  mutate(length = case_when(length < 10000 ~ length*50000, TRUE ~ length)) %>% 
   write.table(paste(output_path, "/intermediate/chr_length_all_auto_ordered.txt", sep = ""),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)




unique(busco_table$Node)


# ─────────────────────────────────────────────────────────────────────────────
# loop over taxa pairs
# ─────────────────────────────────────────────────────────────────────────────
taxa_pairs <- list(
  c("n1", "n4"),
  c("n4", "ilMecPoly1"),
  c("n4", "ilMecLysi212")
)

source("~/scripts/recombination/syntenyplotter_function_busco.R")

anc_colour_file = "/lustre/scratch125/tol/teams/meier/users/kn9/recombination/breakpoints/odp_ms_mec-hap1/list_anc_colour_rainbow.txt"
for (taxa_pair in taxa_pairs) {
  prep_and_plot_synteny(
    taxa1       = taxa_pair[1],
    taxa2       = taxa_pair[2],
    busco_table = busco_table,
    output_path = output_path,
    anc_colour_file = anc_colour_file
  )
}

# run one multiplot
# concatenate size files and remove duplicates (since the same chromosome sizes are used for both comparisons)
sizefile4 <- paste(output_path, "/intermediate/chr_length_n4n1.txt", sep = "")
sizefile1 <- paste(output_path, "/intermediate/chr_length_ilMecLysi212n4.txt", sep = "")
sizefile2 <- paste(output_path, "/intermediate/chr_length_ilMecPoly1n4.txt", sep = "")

list_sizefiles <- c(sizefile1, sizefile2, sizefile4)
# read all size files to one dataframe, remove duplicates, and write to new file
read_sizefiles <- function(sizefiles) {
  do.call(rbind, lapply(sizefiles, function(f) {
    utils::read.delim(f, header = FALSE, colClasses = "character")
  }))
}

all_sizes <- read_sizefiles(list_sizefiles) %>%
  distinct() # remove duplicates

sizefile_out <- paste(output_path, "/intermediate/chr_length_all.txt", sep = "")
write.table(all_sizes, file = sizefile_out, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


  anc_colours <- read.table(anc_colour_file, header = FALSE, sep = "\t",
                             col.names = c("anc_element", "anc_colour"),
                             comment.char = "")
  anc_palette <- setNames(anc_colours$anc_colour, anc_colours$anc_element)


size_file  <- paste(output_path, "/intermediate/chr_length_all_auto_ordered.txt", sep = "")
chain_files <- file.path(output_path, "intermediate", paste0("chain_*.txt"))
fileformat <- "pdf"
w          <- 12
h          <- 8
show_chromosome_names <- TRUE



chain_file_list <- as.list(Sys.glob(file.path(output_path, "intermediate", "chain_*.txt")))
size_file  <- paste(output_path, "/intermediate/chr_length_all_ordered.txt", sep = "")

source("~/scripts/recombination/syntenyplotter_agora_multiplot.R")
do.call(draw.linear, c(
  list(
    directory             = file.path(output_path, "plots"),
    output                = paste0("synt_all_", Sys.Date()),
    sizefile              = size_file,
    fileformat            = fileformat,
    w                     = w,
    h                     = h,
    show_chromosome_names = show_chromosome_names,
    anc_colours           = anc_palette
  ),
  chain_file_list
))

# make palette for the ancestral node
# I neeed to map the chromosmes to the colours



# The first 8 (black through reddish purple) are Okabe-Ito, 
# # The remaining 7 are from Paul Tol's muted palette, which is also colourblind-safe. 
# colour_pal <-
# c(
#   "reddish_purple" = "#CC79A7",
#   "vermillion"     = "#D55E00",
#   "orange"         = "#E69F00",
#   "sand"           = "#DDCC77",
#   "yellow"         = "#F0E442",
#   "teal"           = "#44AA99",
#   "bluish_green"   = "#009E73",
#   "sky_blue"       = "#56B4E9",
#   "black"          = "#000000",
#   "blue"           = "#0072B2",
#   "indigo"         = "#332288",
#   "wine"           = "#882255",
#   "purple"         = "#AA4499",  
#   "cyan"           = "#88CCEE",
#   "green"          = "#117733"

# )

colour_pal <- c("#e75632ff", "#ff9d00" ,"#be8d3dff" ,  "#ffd700" ,  "#4caf50" ,  "#20c997",  "#00cfe8" , "#42a5f5" ,  "#3f51b5" ,"#051886ff", "#8968c1ff", "#3e1585ff" , "#6a1b78ff" ,  "#ab47bc" ,  "#d63384" ,"grey20",  "#f06292",   "#f48fb1", "black")

colour_pal <- c("#42a5f5", "lightblue", "cyan2", "#117733", "lightgreen", "#FFD700", "orange", "#F3E081", "tan", "tan4", "red3","grey40", "grey80", "pink","#ab47bc", "black")

colour_pal <- c(
  "royalblue4", "#1F77B4",  "#17BECF", "#9EDAE5", "#2CA02C", 
  "#98DF8A", "#BCBD22", "#DBDB8D",  "#FFD700", "#FF7F0E","#8C564B", 
  "#C49C94",  "#D62728", "pink","#ab47bc", "#202020"
)

size_file  <- paste(output_path, "/intermediate/chr_length_all_auto_ordered_and_manual.txt", sep = "")

#get the chromosome names for anc node, in the order of the plot
n1_chr <- read.table(size_file, header = FALSE, sep = "\t") %>%
  distinct() %>% filter(V3 == "n1") %>% select(V1) # remove duplicates

palette_n1 <- setNames(colour_pal[1:length(n1_chr$V1)], n1_chr$V1)
palette_n1

write.table(
  data.frame(name = names(palette_n1), value = palette_n1),
  file = paste(output_path, "tables/palette_n1.txt", sep = ""),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

chain_file_list <- as.list(Sys.glob(file.path(output_path, "intermediate", "chain_*.txt")))
fileformat <- "pdf"
w          <- 18
h          <- 6
show_chromosome_names <- TRUE
chr_colours <- palette_n1

source("~/scripts/recombination/syntenyplotter_agora_multiplot.R")
do.call(draw.linear, c(
  list(
    directory             = file.path(output_path, "plots"),
    output                = paste0("synt_all_N2_colour", Sys.Date()),
    sizefile              = size_file,
    fileformat            = fileformat,
    w                     = w,
    h                     = h,
    show_chromosome_names = show_chromosome_names,
    chr_colours           = chr_colours
  ),
  chain_file_list
))

