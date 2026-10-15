library(dplyr)
library(stringr)
library(ggpubr)

#get population data

lys_dir="/lustre/scratch125/tol/teams/meier/users/ac85/pyrho_project/maps_miss0.1_ecuador/"

files <- list.files(lys_dir, pattern = paste0("\\.rmap$"), full.names = TRUE)  

rec_rate_pop_lys <- lapply(files, function(file) {
  df <- read.table(file, header = FALSE, stringsAsFactors = FALSE)
  colnames(df) <- c("start", "end", "rec_rate")
  df$chr = sub(".*SUPER_(.*?)\\.rmap$", "\\1", basename(file))
  df$genome = "ilMecLysi212"
  df
}) %>% bind_rows()

pol_dir="/lustre/scratch125/tol/teams/meier/projects/ithomiini/mechanitis/recombination_mutation/results/polymnia/population/05_pyrho/output_optimize_ecu_col/"

files_pol <- list.files(pol_dir, pattern = paste0("\\.rmap$"), full.names = TRUE)  

rec_rate_pop_pol <- lapply(files_pol, function(file) {
  df <- read.table(file, header = FALSE, stringsAsFactors = FALSE)
  colnames(df) <- c("start", "end", "rec_rate")
  df$chr =sub(".*SUPER_(.*?)\\.rmap$", "\\1", basename(file))
  df$genome = "ilMecPoly1"
  df
}) %>% bind_rows()


rec_rate_pop <- rbind(rec_rate_pop_lys, rec_rate_pop_pol)


#mean rec rate excluding sex chr
mean_rec_auto <- rec_rate_pop %>% filter(!chr %in% c("Z", "Z1", "Z2", "W", "W1", "W2")) %>%
  mutate(distance = end - start) %>%
  group_by(genome) %>%
  summarise(
    weighted_rec_rate =
      sum(rec_rate * distance, na.rm = TRUE) /
      sum(distance, na.rm = TRUE)
  )


# scale  
ped_mean_pol=3.386
ped_mean_lys=4.118

rec_rate_pop <- rec_rate_pop %>% mutate(rec_rate_sc=case_when(genome == "ilMecPoly1" ~ rec_rate*ped_mean_pol/mean_rec_auto[mean_rec_auto$genome=="ilMecPoly1",]$weighted_rec_rate,
                                              genome == "ilMecLysi212" ~ rec_rate*ped_mean_lys/mean_rec_auto[mean_rec_auto$genome=="ilMecLysi212",]$weighted_rec_rate))


rec_rate_pop$chr <- clean_chr(rec_rate_pop$chr)



write.csv(rec_rate_pop, paste("tables/scaled_pop_rec_rate_", Sys.Date(), ".csv",sep = ""),
          quote = F, row.names = F)

chr_rec <- rec_rate_pop %>% 
  mutate(distance = end - start) %>%
  group_by(genome, chr) %>%
  summarise(
    weighted_rec_rate =
      sum(rec_rate_sc * distance, na.rm = TRUE) /
      sum(distance, na.rm = TRUE)
  )

chr_rec %>%
  ggplot() +
  geom_point(aes(chr, weighted_rec_rate)) 

table(rec_rate_pop$chr, rec_rate_pop$genome)

rec_pop_window <- rec_rate_pop %>% 
  mutate(distance = end - start) %>%
  mutate(win_start = (start %/% 100000) * 100000) %>%
  group_by(genome, chr, win_start) %>%
  summarise(
    weighted_rec_rate =
      sum(rec_rate_sc * distance, na.rm = TRUE) /
      sum(distance, na.rm = TRUE)
)

write.csv(rec_pop_window, paste("tables/scaled_pop_rec_rate_100kb_windows_", Sys.Date(), ".csv", sep = ""),
          quote = F, row.names = F)

rec_pop_window %>%
  ggplot() +
  geom_point(aes(win_start, weighted_rec_rate)) +
  facet_wrap(~genome+chr)

# read in polymnia lm


rec_rate_pol <- read.csv("../polymnia/tables/rec_rate_pedigree_polymnia_100kb.csv")
head(rec_rate_pol)
rec_rate_pol$genome <- "ilMecPoly1"

rec_rate_lys <- read.csv("tables/rec_rate_pedigree_M_lysimnia_100kb.csv")
head(rec_rate_lys)
rec_rate_lys$genome <- "ilMecLysi212"

map_all <- rbind(rec_rate_pol, rec_rate_lys)

# merge dfs
head(map_all)
head(rec_pop_window)

map_ped_pop <- left_join(map_all, rec_pop_window, by = c("genome" = "genome", "chrom" = "chr", "start" = "win_start"))

map_ped_pop %>% filter(!chrom %in% c("Z", "Z1", "Z2", "W", "W1", "W2")) %>%
  ggplot(aes(weighted_rec_rate, rate*1000000)) +
  geom_point(aes(colour = genome), size = 1, alpha = 0.5) +
  geom_smooth(aes(colour = genome), method ="lm") +
  stat_cor(method = "spearman", cor.coef.name = "Rho", aes(colour = genome), label.x = 7, show.legend = F) +
  scale_colour_manual(name = "Species",
                      values = c("orange", "#2166AC"),
                      labels = c(expression(italic("Mechanitis lysimnia")), expression(italic("Mechanitis polymnia")))) +
  labs(title = "Population vs pedigree",
       x = "Population recombination rate",
       y = "Pedigree recombination rate") +
  my_theme

ggsave(paste("plots/pop_vs_ped_", Sys.Date(), ".pdf", sep = ""),
       height = 6,
       width = 8)

####
# Chromosome level

summary_maps_pol <- read.csv("../polymnia/tables/summary_ordered_M_polymnia2026-06-03.txt", sep = "\t")
summary_maps_pol$genome <- "ilMecPoly1"
summary_maps_lys <- read.csv("tables/summary_ordered_M_lysimnia2026-06-30.txt", sep = "\t")
summary_maps_lys$genome <- "ilMecLysi212"


summary_maps_pop_ped <- rbind(summary_maps_pol, summary_maps_lys) %>% left_join(., chr_rec, by = c("genome"="genome", "chrom"="chr"))

chr_status <- read.csv(file = paste("../polymnia/tables/chr_status_2026-06-03.csv"))

summary_maps_pop_ped <- summary_maps_pop_ped %>% left_join(., chr_status, by = c("genome"="genome", "chrom"="chr"))
View(summary_maps_pop_ped)

summary_maps_pop_ped %>% filter(!chrom %in% c("Z", "Z1", "Z2", "W", "W1", "W2") & !is.na(type)) %>%
  ggplot(aes(rec_rate, weighted_rec_rate  )) +
  geom_point(aes(shape = type, colour = genome), size = 8, stroke = 1) +
  geom_smooth(method ="lm", ) +
  stat_cor(method = "spearman", cor.coef.name = "Rho", label.y = 6) +
  geom_text(aes(label = chrom), size = 4) +
  scale_shape_manual(name = "Chromosome type",
                     values = c(19, 21, 21),
                     labels = c("Colinear", "Rearranged", "NA")) +
  scale_colour_manual(name = "Species",
                      values = c("orange", "#2166AC"),
                      labels = c(expression(italic("Mechanitis lysimnia")), expression(italic("Mechanitis polymnia")))) +
  labs(title = "Population vs pedigree",
       x = "Population recombination rate",
       y = "Pedigree recombination rate") +
  my_theme

ggsave("plots/rec_pop_ped_comp.pdf",
       height = 6,
       width = 8)

summary_maps_pop_ped %>% filter(!chrom %in% c("Z", "Z1", "Z2", "W", "W1", "W2") & !is.na(type)) %>%
  ggplot(aes(rec_rate, weighted_rec_rate  )) +
  geom_point(aes(shape = type, colour = genome), size = 8, stroke = 1) +
  geom_smooth(method ="lm", ) +
  stat_cor(method = "spearman", cor.coef.name = "Rho", label.y = 6) +
  geom_text(aes(label = chrom), size = 4) +
  facet_wrap(~type+genome) +
  scale_shape_manual(name = "Chromosome type",
                     values = c(19, 21, 21),
                     labels = c("Colinear", "Rearranged", "NA")) +
  scale_colour_manual(name = "Species",
                      values = c("orange", "#2166AC"),
                      labels = c(expression(italic("Mechanitis lysimnia")), expression(italic("Mechanitis polymnia")))) +
  labs(title = "Population vs pedigree",
       x = "Population recombination rate",
       y = "Pedigree recombination rate") +
  my_theme

ggsave("plots/rec_pop_ped_comp_type_genome.pdf",
       height = 6,
       width = 8)
