---
  title: "combined_analyses"
author: "KNasvall"
date: "2026-06-12"
output: html_document 
editor_options: 
  chunk_output_type: console
---
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r setup2, echo=FALSE}
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(data.table)
library(ggpattern)

source("bin/combined_analyses_functions.R")
```

# Combined analyses of recombination rate, genomic features and breakpoints


```{r set_variables, echo=FALSE}

GENOME_SIZE="../breakpoints/anc_reconstruction_pipeline_260609/results/01_genomes/%s_hap1.fa.fai"
GC="../breakpoints/composition/output/%s_gc_content_100000.txt"
REPEATS="/lustre/scratch125/tol/teams/meier/projects/ithomiini/comp_gen/09_repeats/pantera/RepeatMasker/%s_2/%s_repeat_occupancy.csv"
GENE="/lustre/scratch125/tol/teams/meier/projects/ithomiini/comp_gen/03_gene_annotation/helixer/Mec_CDS_density_100kb/%s_1_CDS_100000bp_density.tsv"

BREAKPOINTS="../breakpoints/anc_reconstruction_pipeline_260609/analysis_break/tables/anc_breakpoints_formatted_after_correction.tsv"
ANCESTRAL="../breakpoints/anc_reconstruction_pipeline_260609/analysis_break/tables/agora_link_parsed_flipped.csv"

COMPARTMENTS="input/%s_all_chroms_eigenvectors_corrected_100kb.tsv"

DIVERSITY="input/%s_filtered_concat.vcf.out.csv.gz"

REC_RATE="../pedigree/lepmap/lysimnia_test/analysis_rec_rate/tables/scaled_pop_rec_rate_100kb_windows.csv"

WINDOW=100000

```




```{r read_in_data, echo=FALSE}
ids <- c("ilMecPoly1", "ilMecLysi212")

######### Size ##############
size_df <- do.call(rbind, lapply(ids, function(id) {
  path <- sprintf(
    GENOME_SIZE,
    id
  )
  d <- read.csv(path, header = FALSE, sep = "\t")[,1:2]
  d$genome <- id
  d
}))

# use only the first three columns, add start for each chr
size_df$V3 <- 1
colnames(size_df) <- c("chr", "end","genome","start")
size_df <- size_df[,c("chr", "start", "end","genome")]
str(size_df)

######### Breakpoints ##############
break_df <- read.csv(BREAKPOINTS, header = T, sep = "\t")

break_df <- break_df %>% filter(type != "inversion")


# merge event_category and parent node and change parent node to syngraph child (N1 -> n4)
break_df <-  break_df %>% mutate(event_timing = case_when(event_parent_node == "N1" ~ paste(event_category, "_n4", sep = ""),
                                                          TRUE ~ event_category)) %>% select(!c(event_parent_node, event_category, all_brp_ids))


break_df %>% filter(event_timing != "anc_proxy" & type != "Complex") %>% 
  ggplot() +
  geom_bar(aes(x = event_timing, fill = type), position = "dodge") +
  facet_wrap(~species) +
  theme_bw() +
  labs(title = "Counts of breakpoint events by species, type, and parent node",
       x = "Parent node",
       y = "Count of events",
       fill = "Event type")


break_df <-  break_df %>% rename(genome = species)

break_df<- break_df[,c("chr", "start", "end", "event_timing", "type", "genome")]
str(break_df)

########## Ancestral ############

anc_df <- read.csv(ANCESTRAL, sep = ",")

# node for synteny comparison and chromosome painting
anc_node <- "n1"

unique(anc_df$Node)
anc_chr.df <- anc_df %>% filter(Node == anc_node) %>% select(chr, geneid) %>% rename(anc_chr = chr)

anc_df <- anc_df %>% select(chr, start_init, end_init, Node, geneid) %>% left_join(., anc_chr.df) %>% rename(start = start_init,
                                                                                                             end = end_init,
                                                                                                             genome = Node)

anc_df <-  anc_df %>% filter(!is.na(anc_df$anc_chr))

str(anc_df)

# #anc_col_palette <- cols <- c(
#   "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
#   "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
#   "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
#   "#E5C494", "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3",
#   "#E7298A", "#66A61E"
# )

```


```{r read_in_data, echo=FALSE}
######### Repeats ##############
rpt_df <- do.call(rbind, lapply(ids, function(id) {
  path <- sprintf(
    REPEATS,
    id, id
  )
  d <- read.csv(path)
  d$genome <- id
  d
}))

# summarise repeat content
sum_rpt_df <- rpt_df %>%
  group_by(genome, chrom) %>% summarise(occupancy = sum(overlap_bases))



# Add a row for total repeats
totals <- rpt_df %>%
  group_by(chrom, start, end, genome) %>%
  summarise(overlap_bases = sum(overlap_bases),
            fraction_covered = sum(fraction_covered),
            .groups = "drop") %>%
  mutate(repeat_family = "Total_repeats")

rpt_df <- bind_rows(rpt_df, totals)
rpt_df <- rpt_df %>% select(!overlap_bases)
rpt_df <- rpt_df %>% rename(value = fraction_covered)
rpt_df <- rpt_df %>% rename(chr = chrom)
rpt_df <- rpt_df %>% rename(feature = repeat_family)

# complete the df with missing values as 0 in this case
rpt_df <- rpt_df %>%
  complete(
    nesting(chr, start, end, genome),   # keep valid window combinations only
    feature,                             # expand across all feature types
    fill = list(value = 0)
  )

str(rpt_df)



######### GC ##############
gc_df <- do.call(rbind, lapply(ids, function(id) {
  path <- sprintf(
    GC,
    id
  )
  d <- read.csv(path, skip = 1, header = FALSE, sep = "\t")[,c(1:3,5,7,8)]
  d$genome <- id
  d
}))

gc_df$occupancy <- (gc_df$V7 + gc_df$V8)

gc_df$feature <- "GC"
colnames(gc_df) <- c("chr", "start", "end", "fraction_covered", "V5","V6", "genome","occupancy", "feature")
gc_df <- gc_df[,c("chr", "start", "end","feature", "occupancy", "fraction_covered", "genome")]

sum_gc_df <- gc_df %>%
  group_by(genome, chr, feature) %>% summarise(occupancy = sum(occupancy))

gc_df <- gc_df %>% rename(value = fraction_covered) %>% select(!occupancy)

str(gc_df)

######### Genes ##############
gene_df <- do.call(rbind, lapply(ids, function(id) {
  path <- sprintf(
    GENE,
    id
  )
  d <- read.csv(path, sep = "\t", header = FALSE)
  d$genome <- id
  d
}))

colnames(gene_df) <- c("chr", "start", "end", "count","occupancy", "genome")
sum_gene_df <- gene_df %>% group_by(genome, chr) %>% summarise(occupancy = sum(occupancy))

gene_df$fraction_covered <- (gene_df$occupancy/WINDOW)
gene_df$feature <- "Genes"
gene_df <- gene_df[,c("chr", "start", "end", "feature", "occupancy", "fraction_covered", "genome")]
gene_df <- gene_df %>% rename(value = fraction_covered) %>% select(!occupancy)

str(gene_df)
```


```{r read_in_data, echo=FALSE}
######### Compartments ##############

comp_df_raw <- do.call(rbind, lapply(ids, function(id) {
  path <- sprintf(
    COMPARTMENTS,
    id
  )
  d <- read.csv(path, sep = "\t", header = TRUE)
  d$genome <- id
  d
}))

comp_df_raw <- comp_df_raw %>% rename(chr=chrom)

spec="ilMecPoly1"
comp_df_raw %>% select(genome, chr, eig_used) %>% unique()
ggarrange(comp_df_raw %>% filter(genome == spec) %>%
            ggplot() + 
            geom_line(aes(start, E1)) + 
            facet_wrap(~chr, scales = "free_x", nrow = 1),
          comp_df_raw %>% filter(genome == spec) %>%
            ggplot() + 
            geom_line(aes(start, E2)) + 
            facet_wrap(~chr, scales = "free_x", nrow = 1),
          comp_df_raw %>% filter(genome == spec) %>%
            ggplot() + 
            geom_line(aes(start, E3)) + 
            facet_wrap(~chr, scales = "free_x", nrow = 1),
          nrow = 3)


#correct pattern and PC
comp_df_raw_corr <- comp_df_raw %>% mutate(E2_corr=case_when(genome == "ilMecPoly1" & chr %in% c("SUPER_2", "SUPER_8") ~ E1,
                                                             genome == "ilMecPoly1" & chr %in% c("SUPER_14") ~ E3,
                                                             TRUE ~ E2)) %>%
  mutate(E2_corr=case_when( genome == "ilMecPoly1" & chr %in% c("SUPER_8", "SUPER_9", "SUPER_10", "SUPER_11", "SUPER_Z") ~ E2_corr*-1,
                            TRUE ~ E2_corr)) %>%
  mutate(E2_corr=case_when( genome == "ilMecLysi212" & chr %in% c("SUPER_9","SUPER_Z2") ~ E3,
                            TRUE ~ E2_corr))

#correct E3
comp_df_raw_corr <- comp_df_raw_corr %>% mutate(E3=case_when(genome == "ilMecLysi212" & chr %in% c("SUPER_9","SUPER_Z2") ~ E2,
                                                             genome == "ilMecPoly1" & chr %in% c("SUPER_14") ~ E2,
                                                             TRUE ~ E3)) %>%
  #correct PC1
  mutate(E1=case_when(genome == "ilMecPoly1" & chr %in% c("SUPER_2", "SUPER_8") ~ E2,
                      TRUE ~ E1))

spec="ilMecPoly1"

ggarrange(comp_df_raw_corr %>% filter(genome == spec) %>%
            ggplot() + 
            geom_line(aes(start, E1)) + 
            facet_wrap(~chr, scales = "free_x", nrow = 1),
          comp_df_raw_corr %>% filter(genome == spec) %>%
            ggplot() + 
            geom_line(aes(start, E2_corr)) + 
            facet_wrap(~chr, scales = "free_x", nrow = 1),
          comp_df_raw_corr %>% filter(genome == spec) %>%
            ggplot() + 
            geom_line(aes(start, E3)) + 
            facet_wrap(~chr, scales = "free_x", nrow = 1),
          nrow = 3)

ggsave(paste("plots/pc_corrected_", spec, Sys.Date(), ".pdf"),
       height = 12,
       width = 24)                 
comp_df <- comp_df_raw_corr %>% select(!c(eig_used, compartment, E2)) %>% rename(eig_scores = eig_value) %>% pivot_longer(cols = c("eig_scores", "E2_corr", "E3"), names_to = "feature", values_to = "value")
comp_df <- comp_df[,c("chr", "start", "end", "feature", "value", "genome")]

############ Rec rate ##############

rec_df <- read.csv(REC_RATE, sep = ",", header = TRUE)

rec_df <- rec_df %>% rename(value = weighted_rec_rate)
rec_df$feature <- "Recombination_rate"
rec_df$end <- rec_df$win_start + WINDOW
rec_df$start <- rec_df$win_start

# remove Z 
rec_df <- rec_df %>% filter(!chr %in% c("Z", "Z1", "Z2"))

rec_df <- rec_df[,c("chr", "start", "end", "feature", "value", "genome")]
str(rec_df)
```


```{r}



```


```{r diversity}
######### Diversity ##############

library(GGally)


div_df <- do.call(rbind, lapply(ids, function(id) {
  path <- sprintf(
    DIVERSITY,
    id, id
  )
  d <- read.csv(path, sep = ",", header = TRUE)
  d$genome <- id
  d
}))


div_df <- div_df %>% rename(chr=scaffold,
                            pi_Mpol=pi_Mechanitis_polymnia,
                            pi_Mlys=pi_Mechanitis_lysimnia,
                            dxy=dxy_Mechanitis_polymnia_Mechanitis_lysimnia,
                            Fst=Fst_Mechanitis_polymnia_Mechanitis_lysimnia)
summary(div_df)

colnames(div_df)



my_scatter <- function(data, mapping) {
  ggplot(data, mapping) +
    geom_point(size = 0.5, alpha = 0.2, colour = "grey60") +
    geom_smooth(method = "loess", se = FALSE, colour = "grey20", linewidth = 0.8) +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(linewidth = 0.1, colour = "grey20"))
}

my_density <- function(data, mapping) {
  ggplot(data, mapping) +
    geom_density(linewidth = 0.5, alpha = 0.4, colour = "grey20", fill = "grey60") +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(linewidth = 0.1, colour = "grey20"))
}

ggpairs(title = "ilMecPoly1", 
        div_df %>% filter(genome=="ilMecPoly1"), 
        columns = c("sites", "pi_Mpol", "pi_Mlys",  "dxy", "Fst"),
        upper = list(continuous = "cor"),
        lower = list(continuous = my_scatter),
        diag  = list(continuous = my_density))

ggsave(filename = paste("plots/div_dxy_fst_cor_", "ref_ilMecPoly1_", Sys.Date(), ".png", sep = ""),
       height = 10,
       width = 10)


ggpairs(title = "ilMecLysi212", 
        div_df %>% filter(genome=="ilMecLysi212"), 
        columns = c("sites", "pi_Mpol", "pi_Mlys",  "dxy", "Fst"),
        upper = list(continuous = "cor"),
        lower = list(continuous = my_scatter),
        diag  = list(continuous = my_density))

ggsave(filename = paste("plots/div_dxy_fst_cor_", "ref_ilMecLysi212_", Sys.Date(), ".png", sep = ""),
       height = 10,
       width = 10)

div_df <- div_df %>% select("chr", "start", "end","sites", "pi_Mpol", "pi_Mlys",  "dxy", "Fst", "genome") %>%
  pivot_longer(cols = c("sites", "pi_Mpol", "pi_Mlys",  "dxy", "Fst"), names_to = "feature", values_to = "value")

div_df <- div_df[,c("chr", "start", "end", "feature", "value", "genome")]

#convert to bedformat
div_df$start = div_df$start -1

str(div_df)


```


```{r read_in_data, echo=FALSE}
######### Combine ##############

dfs <- list(
  break_df   = break_df,
  size_df    = size_df,
  gc_df      = gc_df,
  rpt_df     = rpt_df,
  gene_df    = gene_df,
  anc_df     = anc_df,
  comp_df    = comp_df,
  div_df     = div_df,
  rec_df     = rec_df
)

lapply(dfs, function(d) names(d))
dfs <- lapply(dfs, function(d) d %>% mutate(chr = clean_chr(chr)))

# Unpack back to individual variables
list2env(dfs, envir = environment())

combined_df <- rbind(gc_df, gene_df, rpt_df, 
                     div_df, 
                     comp_df, rec_df) 

combined_df <- combined_df %>% filter(!str_detect(pattern = "SCAFFOLD|scaffold", chr))

write_csv(combined_df, 
          file = paste("tables/combined_df_all_features_", Sys.Date(), ".csv", sep = "")
)

str(combined_df)
```




```{r brp_df}

#df from data_preparation_summary.R
# combined_df <- read.csv(combined_df, 
#           file = paste("tables/combined_df_all_features_", Sys.Date(), ".csv", sep = "")
#           )

# break_df

str(combined_df)
#make a df

combined_dt <- combined_df
break_dt <- break_df
setDT(combined_dt)
setDT(break_dt)

setkey(combined_dt, genome, chr, start, end)
setkey(break_dt, genome, chr, start, end)
combined_dt[, window_id := .I]


ov <- foverlaps(break_dt, combined_dt, type = "any", nomatch = 0L)

combined_dt <- ov[combined_dt, on = c("window_id","chr", "genome", "start", "end", "feature", "value")]

combined_dt <- combined_dt[is.na(type), type := "colinear"]
combined_dt <- combined_dt[is.na(event_timing), event_timing := "colinear"]
str(combined_dt)

combined_df_brp <- as.data.frame(combined_dt)

#which are colinear and which are rearranged?
df_colinear_only <- combined_df_brp %>% filter(event_timing %in% c("focal", "colinear")) %>% select(genome, chr, type) %>% unique() %>%
  group_by(genome, chr) %>%   # adjust column names as needed
  filter(all(type == "colinear")) %>%
  ungroup()

chr_status <- combined_df_brp %>% select(genome, chr) %>% unique() %>% 
  left_join(., df_colinear_only) %>% 
  mutate(type = case_when(is.na(type) ~"Rearranged",
                          chr %in% c("W", "W1", "Z2") ~ NA,
                          TRUE ~ type))

#write.csv(chr_status, file = paste("tables/chr_status_", Sys.Date(), ".csv", sep = ""), row.names = F, quote = F)
#add rel_pos

#genome_df <- size_df %>% select(chr, end, genome)

combined_df_brp <- left_join(combined_df_brp, size_df %>% select(-start) %>% dplyr::rename(chr_length = end))


combined_df_brp <- combined_df_brp %>% mutate(window_mid = (start + end) / 2,
                                              dist_telo = (chr_length/2) - abs((chr_length/2) - window_mid),
                                              rel_position = dist_telo /chr_length)

# add distance to breakpoint

combined_df_brp <- combined_df_brp %>%
  group_by(chr) %>%
  group_modify(~{
    bp_mid <- .x$window_mid[.x$type != "colinear"]
    
    .x$dist_bp_kb <- sapply(
      .x$window_mid,
      function(x) min(abs(x - bp_mid))
    )
    
    .x
  }) %>%
  ungroup()

combined_df_brp_test <- combined_df_brp %>%
  rowwise() %>%
  mutate(
    # Find the index of the closest breakpoint using absolute difference
    closest_idx = which.min(abs(start - dist_bp_kb)),
    # Extract the distance and type
    closest_bp = dist_bp_kb[closest_idx])
,
    breakpoint_type = type[closest_idx]
  ) %>%
  ungroup()

library(dplyr)
library(purrr)

library(dplyr)
library(data.table)

# 1. Extract true breakpoint locations based on break_status
# (Change "normal" to whatever your background/non-break label is if needed)
true_breakpoints <- combined_df_brp %>%
  filter(!is.na(break_status) & break_status != "colinear") %>% 
  mutate(bp_mid = (i.start + i.end) / 2) %>%
  select(genome, chr, bp_mid, closest_break_type = break_status) %>%
  distinct()

# 2. Convert to data.table for an ultra-fast rolling join
dt_windows <- as.data.table(combined_df_brp)
dt_breaks  <- as.data.table(true_breakpoints)

# Set keys for matching groups
setkey(dt_windows, genome, chr, window_mid)
setkey(dt_breaks, genome, chr, bp_mid)

# 3. Perform a rolling join to the nearest 'bp_mid' value
# 'roll = "nearest"' pairs each window_mid to its mathematically closest bp_mid
result_dt <- dt_breaks[dt_windows, roll = "nearest", on = .(genome, chr, bp_mid = window_mid)]

result_dt$window_mid <- (result_dt$start + result_dt$end)/2 
# 4. Clean up columns and convert back to a tibble
final_df <- result_dt %>%
  as_tibble() %>%
  # Put columns back in their original order plus your new variable
  select(all_of(names(combined_df_brp)), closest_break_type)


final_df[is.na(final_df$closest_break_type),]$closest_break_type <- "colinear"

combined_df_brp <- final_df
# add ends

ends=200000
combined_df_brp <- combined_df_brp %>%
  mutate(event_timing = case_when(window_mid < ends & type == "colinear" ~ "end", 
                                  (chr_length - window_mid) < ends & type == "colinear" ~ "end",
                                  TRUE ~ event_timing))

write_csv(combined_df_brp %>% select(-window_id), file = paste("tables/combined_df_incl_brp_rel_pos", Sys.Date(), ".csv", sep = ""))

str(combined_df_brp)

```



```{r summary}

sum_rpt_df$feature <- "Repeats"

sum_gene_df <- sum_gene_df %>% rename(chrom = chr)
sum_gene_df$feature <- "CDS"

sum_gc_df <- sum_gc_df[,c("genome", "chr", "occupancy", "feature")] %>% rename(chrom = chr)

sum_feat <- rbind(sum_gene_df, sum_rpt_df, sum_gc_df)
sum_feat$chrom <- clean_chr(sum_feat$chrom)

genome_size_df <- size_df %>% group_by(genome) %>% summarise(genome_len = sum(end))



rpt_gene_size_df <- left_join(sum_feat, size_df, by = c("chrom" = "chr", "genome" = "genome")) %>% filter(!grepl("SCAFFOLD", chrom) & !grepl("scaffold", chrom)) %>% mutate(percent=occupancy/end)
# ilMecLysi212  279425840
# ilMecPoly1    282253049


# summary per genome

sink(file = "tables/sum_rpt_gene_content.txt")
rpt_gene_size_df %>% group_by(genome, feature) %>% summarise(total=sum(occupancy)) %>% left_join(., genome_size_df) %>% mutate(percent=total/genome_len)
print(rpt_gene_size_df, n=100)
sink()


rpt_gene_size_df %>% 
  ggplot(aes(end/1000000, percent*100)) +
  geom_point() +
  geom_smooth(method = "lm", linetype = "dashed", se = F) + 
  facet_grid(feature~genome, scales = "free_y") +
  stat_cor(method = "spearman", cor.coef.name = "rho") +
  theme_pubr() +
  theme(panel.background = element_rect(fill = "white", colour = "grey30")) +
  scale_y_continuous("Fraction covered bases (%)") +
  scale_x_continuous("Chromosome length (Mb)")

ggsave(paste("plots/gene_rpt_vs_chr_length_", Sys.Date(), ".pdf", sep = ""),
       height = 8,
       width = 8)

