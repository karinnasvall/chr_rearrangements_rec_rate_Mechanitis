library(ggplot2)

library(cobs)
library(dplyr)
library(ggtext)

args <- commandArgs(trailingOnly = TRUE)

mapdirectory="fitted_maps_for_rec_rate/"
prefix="M_polymnia"

mapdirectory <- args[1]
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

prefix <- args[2]
print(paste("Using prefix for plots and tables:", prefix))

print(files)

# readin the files and combine to one df give the number in the filename as a column, 
# keep only the first 4 columns and rename them to marker, position, genetic_position_male, genetic_position_female and add LG column with the number in the filename
maps <- lapply(files, function(file) {
  map <- read.table(file, header = TRUE, stringsAsFactors = FALSE, skip = 3)
  map <- map[, 1:4]
  colnames(map) <- c("chrom", "position", "genetic_position_male", "genetic_position_female")
  map$position <- as.numeric(gsub("\\*", "", map$position))
  return(map)
})  %>% bind_rows()

str(maps)
# clean chromosome names to be just the number and zeropadded if needed
# substitute everything before _


clean_chr <- function(x) {
  # Remove everything up to and including the last underscore
  stripped <- sub(".*_", "", x)
  # Zero-pad 1-digit + letter (e.g. "2A" -> "02A"), leave others as-is
  ifelse(grepl("^[0-9][A-Za-z]$|^[0-9]$", stripped),
         sub("^([0-9])", "0\\1\\2", stripped),
         stripped)
}


maps$chrom <- clean_chr(maps$chrom)
maps <- maps %>% arrange(chrom)

unique(maps$chrom)

maps$chrom <- factor(maps$chrom)

maps <- maps %>% group_by(chrom) %>% mutate(phys_length=max(position)-min(position)) %>% ungroup()




str(maps)
unique(maps$chrom)

plot_title=prefix


# print information about the maps df
sink(file = paste("tables/test_summary_map_", prefix, Sys.Date(), ".txt", sep = ""))

print(paste("Result for:", plot_title))

#Number of markers in map
print(paste("Number of markers in LG:", length(maps$position)))   
apply(maps[c("chrom")], 2, table)

#maplength
print("Map length per chr min and max")
#cbind(aggregate(map_ordered$distance_min, list(map_ordered$lg), max), aggregate(map_ordered$distance_max, list(map_ordered$lg), max))
aggregate(maps$genetic_position_male, list(maps$chrom), max)
paste("Total map length (cM): ",sum(aggregate(maps$genetic_position_male, list(maps$chrom), max)$x))

paste("Appr recombination rate (cM/Mb):", sum(aggregate(maps$genetic_position_male, list(maps$chrom), max)$x)/sum(aggregate(maps$phys_length, list(maps$chrom), max)$x)*1000000)
sink()

# print as table
summary_maps <- maps %>% group_by(chrom) %>% summarise(map_length = max(genetic_position_male), physical_length = max(phys_length), markers = n()) %>% filter(markers > 100) %>% arrange(chrom) %>% 
  mutate(rec_rate=1000000*map_length/physical_length)

write.table(summary_maps, file = paste0("tables/summary_ordered_", prefix, Sys.Date(), ".txt"), quote = FALSE, row.names = FALSE, sep = "\t")
str(maps)



map.df <- maps

text_size = 14
my_theme <- theme(plot.title = element_markdown(),
                  panel.background = element_blank(),
                  strip.background = element_rect(linewidth = 0.2, fill = "grey90", colour = "black"),
                  axis.title = element_text(size = text_size),
                  axis.text = element_text(size = text_size),
                  strip.text = element_text(size = text_size),
                  panel.border = element_rect(linewidth = 0.2),
                  axis.line = element_line(linewidth = 0.2), 
                  axis.ticks = element_line(linewidth = 0.2),
                  panel.grid.major = element_line(linewidth = 0.1, colour = "grey90")
)

# plot marey map
ggplot(map.df, aes(x=position/1000000, y=genetic_position_male)) +
  geom_point(colour = "grey40", alpha=0.5, size = 1) +
  facet_wrap(~chrom) + 
  labs(title = paste("Marey map <i>M. polymnia<i>"),
       y = "Genetic position (cM)",
       x = "Physical position (Mb)") +
  my_theme
  
ggsave("plots_ms/marey_map_polymnia.pdf",
       width = 14,
       height = 10)


map.df %>% filter(is.na(position))

#map with cobs

smooth_map<-data.frame()

for(i in unique(maps$chrom)){

  map_chr<-maps[maps$chrom == i,]

  print(i)

  if (i == "14"){
    fit = cobs(map_chr$position, map_chr$genetic_position_male,
               method="uniform",
               constraint= "increase", 
               lambda=0.5, 
               degree=1, # for L1 roughness
               knots=seq(min(map_chr$position),max(map_chr$position),length.out=10), # desired nr of knots 
               tau=0.5) # to predict median   
  }else{
fit = cobs(map_chr$position, map_chr$genetic_position_male,
           method="uniform",
           constraint= "increase", 
           lambda=0.5, 
           degree=1, # for L1 roughness
          #knots=seq(min(map_chr$position),max(map_chr$position),length.out=30), # desired nr of knots 
           tau=0.5) # to predict median
}

X = predict(fit,interval="none",z=map_chr$position)[,1]
predY = predict(fit,interval="none",z=map_chr$position)[,2]

#ggplot(data=NULL, aes(x=X, y=predY))+geom_point()

testDF<-as.data.frame(cbind(X=X, Y=predY))
testDF<-testDF[order(testDF$X),]

all(round(testDF$Y,5) == cummax(round(testDF$Y,5)))

# set values smaller than 1e-10 (also includes negative values!) to 0 and round to 10 digits
testDF$Y2 <- ifelse(testDF$Y<1e-10, 0, testDF$Y)
testDF$Y2 <- round(testDF$Y2, 10)

print(all(testDF$Y2== cummax(testDF$Y2)))

testDF$chrom <- i

smooth_map <- rbind(smooth_map, testDF)
}

ggplot(smooth_map, aes(x=X/1000000, y=Y2))+
  geom_point(data=maps, aes(x=position/1000000, y=genetic_position_male), alpha=0.4, colour = "grey70")+
  geom_line(col="blue")+
  facet_wrap(~chrom) +
  labs(x="Physical position (Mb)", y="Genetic position (cM)") +
  my_theme


ggsave(paste("plots_ms/smoothed_map_", prefix, ".pdf", sep = ""),
       width = 14,
       height = 10 )


# Fit cobs spline (example: to marey map data for one LG/sex)

# Predict derivative (deriv=1) over a fine grid
maps_deriv <- data.frame()
for (chr in unique(maps$chrom)){
  
  map_chr <- maps[maps$chrom == chr, ]
  
  if (chr == "14"){
    fit = cobs(map_chr$position, map_chr$genetic_position_male,
               method="uniform",
               constraint= "increase", 
               lambda=0.5, 
               degree=1, # for L1 roughness
               knots=seq(min(map_chr$position),max(map_chr$position),length.out=10), # desired nr of knots 
               tau=0.5) # to predict median   
  }else{
    fit = cobs(map_chr$position, map_chr$genetic_position_male,
               method="uniform",
               constraint= "increase", 
               lambda=0.5, 
               degree=1, # for L1 roughness
               #knots=seq(min(map_chr$position),max(map_chr$position),length.out=30), # desired nr of knots 
               tau=0.5) # to predict median
  }
  
x_grid <- seq(min(map_chr$position), max(map_chr$position),100000)
deriv_pred <- predict(fit, x_grid, deriv = 1)

deriv_df <- data.frame(
    chrom = chr,
  pos_bp = deriv_pred[, 1],
  rate   = deriv_pred[, 2]   # cM per bp
)
maps_deriv <- rbind(maps_deriv, deriv_df)
}

ggplot(maps_deriv, aes(x=pos_bp/1e6, y=rate*1e6)) +
  geom_point(colour = "grey", alpha=0.7, size = 1) +
  geom_smooth(span = 0.2, se = F, colour = "royalblue", linewidth = 1) +
  facet_wrap(~chrom) +
  scale_y_continuous(limits = c(0,max(maps_deriv$rate*1e6, na.rm = T))) +
  labs(title = paste("Recombiation rate in <i>M. polymnia<i>"), x = "Position (Mb)", y = "Recombination rate (cM/Mb)") +
  my_theme

ggsave(paste0("plots_ms/Recombination_rate_", prefix, ".pdf"), 
       width=14, 
       height=10)


####### rec rate vs chr length

fai_dir <- "../../rawdata_to_vcf/03_variant_call_polymnia_rptmask/ref/"
### get chr lengths
print(paste("Reading .fai from directory:", fai_dir))
files <- list.files(fai_dir, pattern = paste0(".*\\.fai$"), full.names = TRUE)  

chr_lengths <- lapply(files, function(file) {
  df <- read.table(file, header = FALSE, stringsAsFactors = FALSE)[, 1:2]
  colnames(df) <- c("chrom", "length")
  df
}) %>% bind_rows()

# clean chromosme names to be just the number and zeropadded if needed
chr_lengths$chrom <- clean_chr(chr_lengths$chrom)
chr_lengths$chrom <- factor(chr_lengths$chrom)

# add chromosome length for chr 2A and 2B, that do not exist in the fai files

chr_lengths <- chr_lengths %>%
  add_row(chrom = "02A", length = 13051474) %>%
  add_row(chrom = "02B", length = chr_lengths[chr_lengths$chrom=="02", "length"] - 13051474) %>%
  filter(chrom != "02")

# add to df
summary_maps <-  summary_maps %>% left_join(., chr_lengths, by = "chrom")

# test nonlinear model
m_fit <- nls(rec_rate ~ a * (length/1000000)^b, data = summary_maps, start = list(a = 1, b = 1))
sink(file = paste("tables/test_summary_map_", prefix, Sys.Date(), ".txt", sep = ""), append = T)
summary(m_fit)
sink()

pred_model <- data.frame(x = summary_maps$length, y = predict(m_fit))

eq_txt <- sprintf("y = %.3f + %.3f x", coef(m_fit)[1], coef(m_fit)[2])
p_txt  <- sprintf("p-value = %.3g", coef(summary(m_fit))[2, 4])


ggplot(summary_maps, aes(x = length/1000000, y = rec_rate)) +
  geom_point(size = 8, shape = 21, fill= "grey80", alpha = 0.5) +
  #geom_smooth(se = FALSE, method = "gam", linetype = "dashed", linewidth = 1) +
  geom_line(data = pred_model, aes(x/1000000, y)) +
  geom_line(aes(x = length/1000000, y = (50/length)*1000000), linetype = "dotted", linewidth = 1, color = "tan1") +
  geom_line(aes(x = length/1000000, y = (100/length)*1000000), linetype = "dotted", linewidth = 1, color = "tan4") +
  geom_text(aes(label = chrom), size = 4) +
  annotate("text", x = 20, y = Inf, hjust = -0.0, vjust = 2, label = eq_txt, size = text_size/2.8 ) +
  annotate("text", x = 20, y = Inf, hjust = -0.0, vjust = 4, label = p_txt, size = text_size/2.8) +
  labs(title = paste("Recombination rate vs chromosome length, pedigree <i>M. polymnia<i>"),
       x = "Chromosome length (Mb)",
       y = "Recombination rate (cM/Mb)") +
  scale_y_continuous(limits = c(1, 7)) +
  my_theme

ggsave(paste0("plots_ms/rec_rate_vs_chrom_length_", Sys.Date(), ".pdf"), width = 8, height = 6)


##### get 100kb window

# 2) Build all 100 kb windows from 0 for each chromosome
all_windows <- chr_lengths %>%
  rowwise() %>%
  mutate(start = list(seq(0, length, by = 100000))) %>%
  unnest(start) %>%
  mutate(end = pmin(start + 100000, length)) %>%
  ungroup() %>%
  select(chrom, start, end) %>% filter(chrom != "W")

# 3) Bin recombination positions to 100 kb windows
# maps_deriv columns: chrom, pos_bp, rate
map_windows <- maps_deriv %>%
  mutate(start = (pos_bp %/% 100000) * 100000) %>%
  group_by(chrom, start) %>%
  summarise(rate = mean(rate, na.rm = TRUE), .groups = "drop")

# 4) Join to keep all genome windows; uncovered windows become NA
maps_100kb_full <- all_windows %>%
  left_join(map_windows, by = c("chrom", "start")) %>%
  arrange(chrom, start) 

summary(maps_100kb_full)
ggplot(maps_100kb_full, aes(x = start, y = rate)) +
  geom_point(colour= "grey80", alpha = 0.5) +
  facet_wrap(~chrom)

write.csv(maps_100kb_full, file = "tables/rec_rate_pedigree_polymnia_100kb.csv", row.names = F, quote = F)

##########
#compare with population data

COMBINED="../../combined_analysis/tables/combined_df_all_features_2026-05-26.csv"

rec_df <- read.csv(COMBINED) %>% filter(feature == "Recombination_rate" & genome=="ilMecPoly1") 

rec_df <- rec_df %>%
  mutate(chr=case_when(chr == "02" & start < 13051474 ~ "02A",
                       chr == "02" & start > 13051474 ~ "02B",
                       TRUE ~ chr)) %>% filter(!chr %in% c("02", "W", "Z"))


chr_status <- read.csv(file = paste("tables/chr_status_2026-06-03.csv"))


left_join(maps_100kb_full %>% filter(chrom != "Z"), rec_df %>% filter(chr != "Z"), by = c("chrom"="chr", "start"="start")) %>%
  ggplot(aes(value, rate*1000000)) +
  geom_point(colour = "grey30", alpha = 0.7) +
  geom_smooth(method ="lm", ) +
  stat_cor(method = "spearman") +
  facet_wrap(~chrom) +
  labs(title = "Population vs pedigree <i>M. polymnia<i>",
       x = "Population recombination rate",
       y = "Pedigree recombination rate") +
  my_theme

ggsave(filename = paste("plots_ms/corr_pop_ped_per_chr", Sys.Date(), ".pdf", sep = ""),
       height = 12,
       width = 16)

left_join(maps_100kb_full %>% filter(chrom != "Z"), rec_df, by = c("chrom"="chr", "start"="start")) %>%
  ggplot(aes(value, rate*1000000)) +
  geom_point(colour = "grey30", alpha = 0.7) +
  geom_smooth(method ="lm", ) +
  stat_cor(method = "spearman") +
  #facet_wrap(~chrom) +
  labs(title = "Population vs pedigree <i>M. polymnia<i>",
       x = "Population recombination rate",
       y = "Pedigree recombination rate") +
  my_theme

ggsave(filename = paste("plots_ms/corr_pop_ped_", Sys.Date(), ".pdf", sep = ""),
       height = 12,
       width = 16)

rec_df %>% select(genome, chr, start, value) %>% left_join(., chr_status) %>%
  group_by(chr, type) %>% summarise(mean_value=mean(value, na.rm = T)) %>%
  left_join(., summary_maps, by = c("chr"="chrom")) %>% 
  mutate(type=case_when(chr %in% c("02A", "02B") ~ "Rearranged",
                        TRUE ~ type)) %>%
  ggplot(aes(mean_value, rec_rate)) +
  geom_point(aes(fill= type), shape = 21, size = 8, stroke = 1, colour = "#2166AC") +
  geom_smooth(method ="lm", ) +
  stat_cor(method = "spearman") +
  geom_text(aes(label = chr), size = 4) +
  #facet_wrap(~chrom) +
  scale_fill_manual(name = "",
                    values = c("#2166AC", "white", "white"),
                    labels = c("Colinear", "Rearranged", "NA")) +
  labs(title = "Population vs pedigree <i>M. polymnia<i>",
       x = "Population recombination rate",
       y = "Pedigree recombination rate") +
  my_theme  


ggsave(paste("plots_ms/pop_vs_ped_chr_mean_", Sys.Date(), ".pdf", sep = ""),
       height = 6,
       width = 8)

