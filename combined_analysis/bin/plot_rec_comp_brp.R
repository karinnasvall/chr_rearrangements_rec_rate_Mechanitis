
# script to format and plot compartments agains recombination rate and breakpoints

source("bin/combined_analyses_functions.R")


ids <- c("ilMecPoly1", "ilMecLysi212")

rec_vs_comp.df <- read.csv("tables/combined_df_all_features_2026-05-26.csv") %>% filter(feature %in% c("Recombination_rate", "eig_scores")) %>% select(-window_id) %>% pivot_wider(names_from = feature, values_from = value)

COMPARTMENTS="input/%s_all_chroms_eigenvectors_corrected_100kb.tsv"

comp_df_raw <- do.call(rbind, lapply(ids, function(id) {
  path <- sprintf(
    COMPARTMENTS,
    id, id
  )
  d <- read.csv(path, sep = "\t", header = TRUE)
  d$genome <- id
  d
}))

comp_df_raw <- comp_df_raw %>% rename(chr=chrom)

comp_df_raw$chr <- clean_chr(comp_df_raw$chr)

rec_vs_comp.df <- left_join(rec_vs_comp.df, comp_df_raw)

#########
# Plot genome tracks


size_q1 <- size_df %>% filter(genome=="ilMecPoly1") %>% arrange(chr)

# track_df: one row per window
p_rec <- plot_genome_track(rec_vs_comp.df %>% filter(genome == "ilMecPoly1") %>% select(chr, start, end, Recombination_rate), size_q1,
                           y_col      = "Recombination_rate",
                           #colour_col = "compartment",
                           line       = FALSE,
                           point      = TRUE,
                           loess = TRUE,
                           title = "ilMecPoly1",
                           y_label    = "Recombination rate (cM/Mb)"
)

p_comp <- plot_genome_track(rec_vs_comp.df %>% filter(genome == "ilMecPoly1") %>% select(chr, start, end, compartment), size_q1,
                            y_col      = "compartment",
                            colour_col = "compartment",
                            line       = FALSE,
                            point      = FALSE,
                            vline       = TRUE,
                            shape = 15,
                            loess = FALSE,
                            y_label    = "Compartment"
) + scale_y_discrete(breaks = NULL, labels = NULL) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank())

size_q2 <- size_df %>% filter(genome=="ilMecLysi212") %>% arrange(chr)

# track_df: one row per window
p_rec_lys <- plot_genome_track(rec_vs_comp.df %>% filter(genome == "ilMecLysi212") %>% select(chr, start, end, Recombination_rate), size_q2,
                               y_col      = "Recombination_rate",
                               #colour_col = "compartment",
                               line       = FALSE,
                               point      = TRUE,
                               loess = TRUE,
                               title = "ilMecLysi212",
                               y_label    = "Recombination rate (cM/Mb)"
)

p_comp_lys <- plot_genome_track(rec_vs_comp.df %>% filter(genome == "ilMecLysi212") %>% select(chr, start, end, compartment), size_q2,
                                y_col      = "compartment",
                                colour_col = "compartment",
                                line       = FALSE,
                                point      = FALSE,
                                vline = TRUE,
                                shape = 15,
                                loess = FALSE,
                                y_label    = "Compartment"
) + scale_y_discrete(breaks = NULL, labels = NULL) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

ggarrange(p_rec, p_comp, p_rec_lys, p_comp_lys, ncol = 1, heights = c(2,1, 2, 1), 
          align = "v", common.legend = T, legend = "bottom")

ggsave(paste("plots/genome_track_rec_vs_comp_", Sys.Date(), ".png", sep = ""),
       height = 10,
       width = 20)


#######
# compartment lengths


comp_length <- rec_vs_comp.df %>% select(genome, chr, start, end, compartment, Recombination_rate) %>% arrange(genome, chr, start) %>%
  group_by(genome, chr) %>%
  mutate(run_id = cumsum(row_number() == 1 | compartment != lag(compartment))) %>%
  group_by(genome, chr, run_id, compartment) %>%
  summarise(start = min(start),
            end = max(end),
            .groups = "drop") %>%
  arrange(genome, chr, start)

comp_length$comp_length <- comp_length$end - comp_length$start

comp_length <- comp_length[, c("genome", "chr", "start", "end", "compartment", "comp_length")]

write_tsv(comp_length, file = "tables/compartment_lengths.bed", col_names = F)

comp_length %>% 
  ggplot(aes(y = compartment)) +
  geom_segment(
    aes(x = start, xend = end, yend = compartment, colour = compartment),
    linewidth = 3,
    lineend = "butt"
  ) +
  facet_wrap(~ genome + chr, scales = "free_x", nrow = 2) +
  labs(x = "Genomic position", fill = "Compartment") +
  theme_bw()


# excluding the W
comp_len_distr <- comp_length %>% filter(!is.na(compartment) & !chr %in% c("W", "W1", "W2")) %>% 
  ggplot(aes(compartment, comp_length/1000000)) +
  geom_boxplot(aes(fill=compartment), outliers = TRUE, alpha = 0.7, outlier.color = "grey70") +
  stat_compare_means() +
  stat_summary(
    geom = "text",
    fun.data = function(x) data.frame(y = -Inf, label = paste0("n=", length(x))),
    size = 5, vjust = -0.5
  ) +
  facet_grid(~genome) +
  scale_y_continuous(#transform = "log10",
    name = "Compartment length (Mb)") +
  scale_fill_manual(name = "Compartment", 
                    values = c("red2", "royalblue")) +
  theme_classic() +
  theme(axis.title.x = element_blank())

comp_length %>% filter(!is.na(compartment)) %>% 
  ggplot(aes(comp_length/1000000)) +
  geom_histogram(aes(fill = compartment), alpha = 0.8, bins = 100) +
  xlim(0,5) +
  facet_wrap(~genome)

comp_length %>% group_by(genome, chr, compartment) %>% count()
comp_length %>% group_by(genome, compartment) %>% count()

ggsave(comp_len_distr, filename = paste("plots/comp_comb_length_box_", Sys.Date(), ".png", sep = ""),
       height = 4,
       width = 6)

```

```{r rec_comp}

rec_vs_comp.df %>%
  ggplot() +
  geom_density(aes(Recombination_rate, fill = compartment), alpha = 0.5) +
  geom_density(aes(Recombination_rate), fill = "grey", alpha = 0.5) +
  facet_wrap(~genome) +
  theme_bw()

rec_comp_plot <- rec_vs_comp.df %>%
  ggplot(aes(compartment, Recombination_rate)) +
  geom_violin(aes(compartment, Recombination_rate, fill = compartment), alpha = 0.5) +
  geom_boxplot( outliers = FALSE, width = 0.3, fill = "grey90") +
  stat_compare_means() +
  facet_grid(~genome, labeller = my_labeller_comp) +
  scale_y_continuous(name = "Recombination rate") +
  scale_x_discrete(name = "Compartment") +
  scale_fill_manual(values = c("red2", "royalblue")) +
  theme_classic() + theme(legend.position = "none",
                          axis.title = element_blank(),
                          axis.text = element_text(size = textsize),
                          strip.text = element_text(size = textsize))



legend_comp <- get_legend(brp_comp_plot)

ggarrange(comp_len_distr + theme(legend.position = "none"), 
          rec_comp_plot  + theme(legend.position = "none"), 
          brp_comp_plot  + theme(legend.position = "none"), 
          legend_comp)

ggsave(paste("plots/comp_comb_", Sys.Date(), ".png", sep = ""),
       height = 8,
       width = 12)


```

```{r}


library(data.table)

# Expected columns:
# roi_df:  chr, start, end, type
# feat_df: chr, start, end, feature  (values: "A", "B")

roi_df <- break_df %>% select(chr, start, end, type, event_timing, genome)
feat_df <- comp_length

setDT(roi_df);  setDT(feat_df)
roi_df[, roi_id := .I]  # unique ID per ROI

hits <- roi_df[feat_df,
               on = .(genome, chr, start < end, end > start),
               nomatch = NULL,
               .(roi_id = x.roi_id, type = x.type, event_timing = x.event_timing, genome = x.genome, feature = i.compartment)]

roi_features <- hits[, .(features = list(sort(unique(feature)))), by = .(roi_id, type, event_timing, genome)]
roi_features[, overlap_class := fcase(
  sapply(features, function(f) identical(f, "A")),        "A",
  sapply(features, function(f) identical(f, "B")),        "B",
  sapply(features, function(f) all(c("A", "B") %in% f)), "A or B",
  default = "other"
)]

result <- roi_features[, .N, by = .(genome, type, event_timing, overlap_class)]
dcast(result, type ~ overlap_class, value.var = "N", fill = 0)

brp_comp_plot_all <- ggplot(result, 
                            aes(x = event_timing, y = N, fill = overlap_class)) +
  geom_col(alpha = 0.7) +
  scale_fill_manual(
    values = c(
      "A"       = "red2",
      "B"       = "royalblue",
      "A or B" = "grey70"
    ),
    name = "Overlaps"
  ) +
  labs(x = "ROI type", y = "Count") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(brp_comp_plot_all, filename = paste("plots/comp_comb_brp_all_", Sys.Date(), ".png", sep = ""),
       height = 8,
       width = 12)


brp_comp_plot <- ggplot(result %>% filter(event_timing %in% c( "focal")), 
                        aes(x = type, y = N, fill = overlap_class)) +
  geom_col(alpha = 0.5) +
  facet_wrap(~event_timing) +
  scale_fill_manual(
    values = c(
      "A"       = "red2",
      "B"       = "royalblue",
      "A or B" = "grey70"
    ),
    name = "Compartment"
  ) +
  labs(x = "", y = "Count") +
  scale_x_discrete(labels = c("Fission", "Fission-fusion", "Fusion")) +
  theme_classic() +
  theme(axis.text.x = element_text(size = textsize),
        strip.background = element_blank(),
        strip.text = element_blank())

```

```{r comp_vs_breakpoint}


comp_length

brp_comp_plot_anc <- anc_only %>% filter(!is.na(compartment) & event_timing == "focal", ) %>%
  ggplot(aes(type)) + 
  geom_bar(aes(fill = compartment)) +
  facet_grid(~event_timing) +
  scale_fill_manual(values = c("red2", "royalblue"),
                    name = "Compartment") +
  scale_y_continuous(name = "Count 100kb windows") +
  theme_classic() +
  theme(axis.title.x = element_blank())


anc_only <- rec_vs_comp.df %>% 
  filter(event_timing == "anc_proxy") %>%
  anti_join(
    rec_vs_comp.df %>% filter(event_timing == "focal"),
    by = c("chr", "start", "end", "type", "genome")
  ) %>% select(genome, chr, i.start,  i.end, event_timing, type, compartment) 

View(anc_only)

ilMecPoly1
11
5360096
5365971
anc_proxy
fusion
A

ilMecLysi2… 13    7.1 e6 7.20e6 7192328 7.39e6 focal        fusi…   12790513    7150000   5640513
2 ilMecLysi2… 13    7.20e6 7.3 e6 7192328 7.39e6 focal        fusi…   12790513    7250000   5540513
3 ilMecLysi2… 13    7.3 e6 7.4 e6 7192328 7.39e6 focal        fusi…   12790513    7350000   5440513
rec_vs_comp.df %>% 
  +     filter(event_timing == "focal" & chr == "13") %>% select(compartment)
# A tibble: 3 × 1
compartment
<chr>      
  1 A          
2 B          
3 B          


ilMecLysi212
04
11195219
11223413
anc_proxy
fission-fusion
A

ilMecPoly1  11    5.3 e6 5.40e6 5360096 5.37e6 focal        fiss…   16717337    5350000   5350000 A ilMecPoly1  03 all B
View(anc_only)


```

```{r brp-comp}

brp_comp_plot <- ggplot(result %>% filter(event_timing %in% c( "focal")), 
                        aes(x = type, y = N, fill = overlap_class)) +
  geom_col(alpha = 0.5) +
  facet_wrap(~event_timing) +
  scale_fill_manual(
    values = c(
      "A"       = "red2",
      "B"       = "royalblue",
      "A or B" = "grey70"
    ),
    name = "Compartment"
  ) +
  labs(x = "", y = "Count") +
  scale_x_discrete(labels = c("Fission", "Fission-fusion", "Fusion")) +
  theme_classic() +
  theme(axis.text.x = element_text(size = textsize),
        strip.background = element_blank(),
        strip.text = element_blank())

legend_comp <- get_legend(brp_comp_plot)

```

```{r}

my_labeller_comp1 <- labeller(
  genome = c(
    ilMecPoly1   = "italic('Mechanitis polymnia')",
    ilMecLysi212 = "italic('Mechanitis lysimnia')"
  ),
  .default = label_parsed
)
rec_comp_plot <- rec_vs_comp.df %>%
  ggplot(aes(compartment, Recombination_rate)) +
  geom_violin(aes(compartment, Recombination_rate, fill = compartment), alpha = 0.5) +
  geom_boxplot( outliers = FALSE, width = 0.3, fill = "grey90") +
  stat_compare_means() +
  facet_grid(~genome, labeller = my_labeller_comp1) +
  scale_y_continuous(name = "Recombination rate") +
  scale_x_discrete(name = "Compartment") +
  scale_fill_manual(values = c("red2", "royalblue")) +
  theme_classic() + theme(legend.position = "none",
                          axis.title = element_blank(),
                          axis.text = element_text(size = textsize),
                          strip.text = element_text(size = textsize))




```


```{r}
taxa1 = "Mechanitis polymnia"
taxa2 = "Mechanitis lysimnia"

my_labeller_comp2 <- labeller(
  species = c(
    ilMecPoly1   = "italic('Mechanitis polymnia')",
    ilMecLysi212 = "italic('Mechanitis lysimnia')"
  ),
  .default = label_parsed
)

bin_data <- read.csv(paste("overlapping_compartments", REF, QUERY, Sys.Date(), "_table.csv", sep = ""))

pol_df <- rec_vs_comp.df %>% 
  filter(genome == "ilMecPoly1" & event_timing == "focal") %>% add_column(species = "ilMecPoly1") %>% dplyr::rename(chrom=chr) 

lys_pol_ref_df <- rec_vs_comp.df %>% 
  filter(genome == "ilMecPoly1" & event_timing == "anc_proxy") %>% add_column(species = "ilMecLysi212") %>% dplyr::rename(chrom=chr) %>% rbind(pol_df)

lys_pol_ref_df$species <- factor(lys_pol_ref_df$species, levels = c("ilMecPoly1", "ilMecLysi212"))

p <- ggplot(bin_data) +
  geom_rect(
    aes(xmin = x_start, xmax = x_end,
        ymin = 0,        ymax = 1,
        fill = comp_clean),
    colour = NA, alpha = 0.7
  ) +
  geom_vline(data = lys_pol_ref_df,
             aes(xintercept = (start + end)/2000000),
             linewidth = 1, colour = "black") +
  scale_fill_manual(
    values   = c("red2", "royalblue", "#B0C4DE"),
    #labels   = TRACK_LABELS,
    name     = "Compartment",
    na.value = "grey92"
  ) +
  scale_x_continuous(
    #labels = function(x) paste0(x, " Mb"),
    #expand = c(0.005, 0)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_grid(species ~ chrom, scales = "free_x", space = "free_x", switch = "both", labeller = my_labeller_comp2) +
  labs(
    # title    = paste("Compartment conservation:", REF, "vs", QUERY, sep = " "),
    # subtitle = paste0(
    #   "Compartments lifted via synteny blocks | bin size = ",
    #   BIN_SIZE / 1000, " kb | grey = no data"
    # ),
    x =  bquote("Position in" ~ italic(.(taxa1)) ~ "genome (Mb)"),
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    panel.border    = element_rect(colour = "grey40", fill = NA, linewidth = 0.5),
    axis.text.y     = element_blank(),
    axis.ticks.y    = element_blank(),
    axis.text.x     = element_blank(),
    axis.title.x     = element_text(size = textsize),
    strip.text.x    = element_text(size = textsize),
    strip.text.y    = element_text(size = textsize, angle = 0),
    legend.position = "bottom",
    legend.key.size = unit(0.5, "cm"),
    legend.text     = element_text(size = 8),
    #plot.title      = element_text(face = "bold", size = 13),
    #plot.subtitle   = element_text(size = 9, colour = "grey40"),
    plot.margin     = margin(10, 10, 10, 10)
  )

p
ggarrange(
  ggarrange(brp_comp_plot  + theme(legend.position = "none"),
            rec_comp_plot  + theme(legend.position = "none"),
            legend_comp, 
            nrow = 1, labels = c("A", "B"), widths = c(1,1,0.3)),
  p + theme(legend.position = "none", plot.margin = margin(t=20,r=1,b=1,l=20)), nrow = 2, labels = c("", "C"))

ggsave(paste("plots/comp_comb_", Sys.Date(), ".png", sep = ""),
       height = 12,
       width = 16)

ggsave(paste("plots/comp_comb_", Sys.Date(), ".png", sep = ""),
       height = 12,
       width = 16)



