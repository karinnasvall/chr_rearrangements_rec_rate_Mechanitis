#!/usr/bin/env Rscript --vanilla

#script by Karin Näsvall using syntenyPlotter
#https://github.com/Farre-lab/syntenyPlotteR
#This script uses the output from minimap2 whole genome alignments.
#Creates chain files for syntenyPlotter.

#install.packages("devtools")
#devtools::install_github("marta-fb/syntenyPlotteR")
#library(syntenyPlotteR)
#install.packages("pals")
library(pals)
library(tidyverse)
library(ggplot2)


rm(list = ls())

#dir organisation: intermediate/, plots/


######
####If using command line and argument uncomment these (and comment away the variable lines below):
cmd_args <- commandArgs(trailingOnly = TRUE)
table_file <- cmd_args[1]
align_length=as.numeric(cmd_args[2])
# taxa1=cmd_args[3]
# taxa2=cmd_args[4]
# table file must be in the format dir/taxa1_taxa2.paf

#################
# Variables if using hardcoding, change these variables:
#output from minimap2 paf-format
#table_file <- "ilMecLysi317.hap1_ilMecLysi317.hap2.paf"
# # #reference (the taxa in the sixth column of the output in minimap2)
# taxa1="ilTitTari1"
# # #query (the taxa in the first column of the output in minimap2)
# taxa2="ilTitHarm1"
#align_length=100000
#################

taxa1=str_split(sub(".paf", "", sub(".*/","", table_file)), pattern = "_")[[1]][1]
taxa2=str_split(sub(".paf", "", sub(".*/","", table_file)), pattern = "_")[[1]][2]


paf_table <- read.csv2(table_file, sep = "\t", header = F)[,1:10]

colnames(paf_table) <- c("query", "Qseq_length", "Qstart", "Qend", "strand", "reference", "Rseq_length", "Rstart", "Rend",  "matches")

head(paf_table)

#align_length=100000
length(paf_table[paf_table$matches > align_length,c("matches")])

#filter , used 100-500 kb alignments since closely related taxa, needs to be adjusted
paf_table <- paf_table[paf_table$matches > align_length,]
length(paf_table$matches)

#remove small scaffolds
paf_table <- filter(paf_table, !grepl("unloc|SCAFFOLD", query))
paf_table <- filter(paf_table, !grepl("unloc|SCAFFOLD", reference))

paf_table$refID <- taxa1
paf_table$queryID <- taxa2

#correct column order for the chainfile
chain_table <- paf_table[, c("reference", "Rstart", "Rend", "query", "Qstart", "Qend", "strand", "refID", "queryID", "Rseq_length", "Qseq_length")]

# check nr of scffolds longer than min_scaffold_length
min_scaffold_length=500000
length(chain_table[chain_table$Rseq_length > min_scaffold_length, "Rseq_length"])
length(chain_table[chain_table$Qseq_length > min_scaffold_length, "Qseq_length"])
# remove small scaffolds
chain_table <- chain_table[chain_table$Qseq_length > min_scaffold_length,]
chain_table <- chain_table[chain_table$Rseq_length > min_scaffold_length,]

#change names of chr OBS customise!!!
unique(chain_table$reference)
unique(chain_table$query)

#zeropadding
pad_chr <- function(x, width = 2) {
  ifelse(
    grepl("^[ZW]", x),
    x,
    sprintf(paste0("%0", width, "d"), as.integer(x))
  )
}
#if ending in _HAP1 or _HAP2, remove that part of the name
chain_table$reference <- sub("_HAP1|_HAP2", "", chain_table$reference)
chain_table$query <- sub("_HAP1|_HAP2", "", chain_table$query)

chain_table$reference <- as.character(pad_chr(sub("^[^_]*_", "", chain_table$reference)))
chain_table$query <- as.character(pad_chr(sub("^[^_]*_", "", chain_table$query)))


head(chain_table)

#order df after chr length in ref taxa, will also order the other taxa after first appearance in the chain file
chain_table <- chain_table %>% arrange(desc(Rseq_length))


#change orientation if majority is on negative strand 
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


chain_1 <- chain_table
chain_1$Rstart_init <- chain_1$Rstart
chain_1$Rend_init <- chain_1$Rend

chain_1$Qstart_init <- chain_1$Qstart
chain_1$Qend_init <- chain_1$Qend


for (i in unique(chain_table$query)) {
  if (Mode(chain_table[chain_table$query==i, "strand"])=="-") {
    print(i)
    
    chain_1[chain_1$query==i, "Qend"] <- chain_1[chain_1$query==i, "Qseq_length"] - chain_1[chain_1$query==i, "Qstart"]

    chain_1[chain_1$query==i, "Qstart"] <- chain_1[chain_1$query==i, "Qseq_length"] - chain_1[chain_1$query==i, "Qend_init"]
    
    chain_1[chain_1$query==i & chain_1$strand=="-", "strand"] <- "temp"
    chain_1[chain_1$query==i & chain_1$strand=="+", "strand"] <- "-"
    chain_1[chain_1$query==i & chain_1$strand=="temp", "strand"] <- "+"
  }
}

chain_table <- chain_1

chain_table <- chain_table %>%
  mutate(
    aln_len = abs(Rend - Rstart)
  )

chr_pair_lengths <- chain_table %>%
  group_by(reference, query) %>%
  summarise(
    total_aligned = sum(aln_len, na.rm = TRUE),
    n_blocks = n(),
    .groups = "drop"
  ) %>%
  # zeropadding reference and query names to ensure correct ordering
  arrange(desc(total_aligned)) %>%
  arrange(reference)

  best_chr_by_alignment <- function(df, chr_col, score_col = total_aligned) {
  chr_col <- enquo(chr_col)
  score_col <- enquo(score_col)
  
  df %>%
    mutate(.row_id = row_number()) %>%
    group_by(!!chr_col) %>%
    slice_max(!!score_col, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(.row_id) %>%
    select(-.row_id)
}


#make size_file for input (and order of chr in the plot)
chr_size_q <- unique(chain_table[,c("query", "Qseq_length", "queryID")])
chr_size_q <- chr_size_q %>% arrange(match(query, best_chr_by_alignment(chr_pair_lengths, query)$query))
colnames(chr_size_q) <- c("chr", "seq_length", "ID")
chr_size_r <- unique(chain_table[,c("reference", "Rseq_length", "refID")])
chr_size_r <- chr_size_r %>% arrange(match(reference, best_chr_by_alignment(chr_pair_lengths, reference)$reference))
colnames(chr_size_r) <- c("chr", "seq_length", "ID")
chr_size <- rbind(chr_size_q, chr_size_r)
chr_size$polymorphic <- "N"


#writing input files
names(chr_size) <- NULL
write.table(chr_size, file = paste("intermediate/chr_length_", taxa2, taxa1, ".txt", sep = ""), sep = "\t", row.names = F)

names(chain_table) <- NULL
write.table(chain_table, file = paste("intermediate/chain_", taxa1, taxa2, ".txt", sep = ""), sep = "\t", row.names = F)

source("/nfs/users/nfs_k/kn9/scripts/comp_gen_scripts_working/04_genome_alignment/syntenyplotter_mod_v2.R")
#produce the synteny figures
draw.linear(directory = "plots", output = paste("synt_", taxa1, taxa2, Sys.Date(), sep = ""), 
            paste("intermediate/chr_length_", taxa2, taxa1, ".txt", sep = ""), 
            paste("intermediate/chain_", taxa1, taxa2, ".txt", sep = ""),  
            fileformat = "pdf", w=13,h=4, show_chromosome_names = F)

draw.linear(directory = "plots", output = paste("synt_chrnames", taxa1, taxa2, Sys.Date(), sep = ""), 
            paste("intermediate/chr_length_", taxa2, taxa1, ".txt", sep = ""), 
            paste("intermediate/chain_", taxa1, taxa2, ".txt", sep = ""),  
            fileformat = "pdf", w=13,h=4, show_chromosome_names = T)

# ggsave(filename = paste("plots/synt_", taxa1, taxa2, Sys.Date(), ".pdf", sep = ""), 
#        device = "png", width = 13, height = 5)

