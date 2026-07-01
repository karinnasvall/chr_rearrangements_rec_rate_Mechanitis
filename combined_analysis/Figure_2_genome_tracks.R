
library(tidyverse)
library(patchwork)

# use df:s from data_preparation_summary.R

# # colour based on tab20, 09 is Z, 15 CARs
# anc_col_palette <- c(
#   "royalblue4", "#1F77B4",  "#17BECF", "#9EDAE5", "#2CA02C", 
#  "#98DF8A", "#BCBD22", "#DBDB8D", "#8C564B", 
#   "#C49C94", "#FFBB78", "#FF7F0E", "#D62728", "#FF9896",   "#202020"
# 
# #  "#F7B6D2", "#E377C2", "#9467BD", "#C5B0D5", "#AEC7E8", "#7F7F7F", "#C7C7C7"
# )

anc_col_palette <- c(
  "royalblue4", "#1F77B4",  "#17BECF", "#9EDAE5", "#2CA02C", 
  "#98DF8A", "#BCBD22", "#DBDB8D",  "#FFD700", "#FF7F0E","#8C564B", 
  "#C49C94",  "#D62728", "pink","#ab47bc", "#202020"
)

unique(anc_df[, c("anc_chr", "chr", "genome")])
unique(anc_df$anc_chr)
anc_colors <- setNames(anc_col_palette, unique(anc_df$anc_chr))

anc_colors <- read.table("/lustre/scratch125/tol/teams/meier/users/kn9/recombination/breakpoints/anc_reconstruction_pipeline_260609/analysis_break/tables/palette_n1.txt", header = TRUE, sep = "\t", comment.char = "", colClasses = c(name = "character"))


anc_colors <- setNames(anc_colors$value, anc_colors$name)




# ── Compute cumulative offsets from size_df ──────────────────────────────────
gap <- 1000000  # gap between chromosomes

make_offsets <- function(size_df) {
  size_df <- size_df %>% arrange(chr)
  size_df %>%
    mutate(chr_len = end - start,
           offset = cumsum(lag(chr_len + gap, default = 0)),
           chr_mid = offset + chr_len / 2)
}

# ── Plot one genome ─────────────────────────────────────────────────────────
plot_one_genome <- function(species_name, genome_name, size_df, gene_df, combined_df, break_df, show_legend = TRUE, feature_order = NULL, smooth = FALSE, smooth_span = 0.3, area_shade=0, plot_paint=TRUE, text_size=12) {
  
  offsets <- make_offsets(size_df)
  offset_map <- setNames(offsets$offset, offsets$chr)
  
  # Add offsets to all data
  size_off <- offsets %>%
    mutate(x_start = offset, x_end = offset + chr_len)
  
  gene_off <- gene_df %>%
    filter(genome == genome_name) %>%
    mutate(x_start = start + offset_map[chr],
           x_end   = end   + offset_map[chr])
  
  feat_off <- combined_df %>%
    filter(genome == genome_name) %>%
    mutate(x_mid = (start + end) / 2 + offset_map[chr])
  
  break_off <- break_df %>%
    filter(genome == genome_name, 
           !str_detect(pattern = "inversion", type)) %>%
    mutate(x_start = start + offset_map[chr],
           x_end   = end   + offset_map[chr])
  
  # Chr label positions
  chr_labels <- size_off %>% select(chr, chr_mid)
  
  # Create alternating background rectangles from size_off
  chr_bg <- size_off %>%
    mutate(fill = ifelse(row_number() %% 2 == 0, "grey85", "white"))
  
  
  # ── Chromosome painting ────────────────────────────────────────────────
  p_paint <- ggplot() +
    # Add as first layer in both p_paint and track plots:
    geom_rect(data = size_off,
              aes(xmin = x_start, xmax = x_end, ymin = -0.1, ymax = 0.1),
              fill = "white", color = "black", linewidth = 0.3) +
    # Ancestral chr painting
    geom_rect(data = gene_off,
              aes(xmin = x_start-100000, xmax = x_end+100000, ymin = -0.08, ymax = 0.08,
                  fill = anc_chr)) +
    geom_point(data = break_off %>% filter(event_timing == "focal"),
               aes(x=x_start, y=0.12, shape=type), color="black", fill = "black", alpha=0.7, size = 2, stroke = 0.5,
               show.legend = c(colour = FALSE, fill = FALSE)) +
    scale_fill_manual(values = anc_colors, name = "Inferred ancestral LG") +
    scale_x_continuous(breaks = chr_labels$chr_mid,
                       labels = chr_labels$chr,
                       expand = c(0.01, 0)) +
    scale_shape_manual(name = "Breakpoint event", 
                       labels = c("Fission", "Fusion", "Fission and fusion"),
                       limits = c("fission", "fusion", "fission-fusion"),
                       values = c(25, 24, 23, 6, 2),
                       drop = FALSE) +
    labs(title = species_name, x = NULL, y = NULL) +
    coord_cartesian(clip = "off") +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(face = "bold", size = text_size),
          axis.ticks.y = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(face = "italic", size = text_size+2), 
          legend.text = element_text(size = text_size),
          legend.title = element_text(size = text_size))
  
  if (!show_legend) {
    p_paint <- p_paint + theme(legend.position = "none")
  }
  
  # ── Feature tracks ────────────────────────────────────────────────────
  
  # feat_colors <- c(
  #   "GC"            = "orange3",
  #   "Genes"         = "goldenrod1",
  #   "Total_repeats" = "tan4",
  #   "Diversity"     = "tan",
  #   "E1"   = "grey30"
  # )
  
  mono_feat_col="grey30"
  feat_colors <- c(
    "GC"            = mono_feat_col,
    "Genes"         = mono_feat_col,
    "Total_repeats" = mono_feat_col,
    "pi_Mpol" = mono_feat_col,
    "pi_Mlys" = mono_feat_col,  
    "dxy" = mono_feat_col, 
    "Fst" = mono_feat_col,
    "eig_scores"   = mono_feat_col,
    "sites" = mono_feat_col,
    "Recombination_rate" = mono_feat_col
  )
  
  make_track <- function(fdf, feat, size_off, break_off, chr_labels, area_shade, smooth_span) {
    
    chr_bg <- size_off %>%
      mutate(fill = ifelse(row_number() %% 2 == 0, "grey85", "white"))
    
    base <- ggplot(fdf, aes(x = x_mid, group = chr)) +
      # geom_vline(data = size_off, aes(xintercept = x_start),
      #            colour = "grey80", linewidth = 0.2) +
      geom_rect(data = chr_bg,
                aes(xmin = x_start, xmax = x_end, ymin = -Inf, ymax = Inf),
                fill = chr_bg$fill, alpha = 0.5, inherit.aes = FALSE) +
      geom_vline(data = break_off,
                 aes(xintercept = x_start, linetype = event_timing),
                 colour = "black", linewidth = 0.4) +
      scale_x_continuous(breaks = chr_labels$chr_mid,
                         labels = chr_labels$chr,
                         expand = c(0.01, 0)) +
      scale_linetype_manual(name = "Breakpoint", 
                            values = c("dotted","dashed"),
                            labels = c("Origin (ancestral proxy)", "Focal")
                            ) +
      theme_minimal() +
      theme(axis.text.x = element_blank(),
            axis.title.y = element_text(size = text_size),
            axis.text.y = element_text(size = text_size),
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            legend.text = element_text(size = text_size),
            legend.title = element_text(size = text_size))
    
    if (feat == "eig_scores") {
      # Positive red, negative blue
      # base +
      #   geom_area(data = fdf %>% mutate(pos = pmax(value, 0)),
      #             aes(y = pos), fill = "firebrick", alpha = 0.4) +
      #   geom_area(data = fdf %>% mutate(neg = pmin(value, 0)),
      #             aes(y = neg), fill = "dodgerblue3", alpha = 0.4) +
      #   #geom_line(aes(y = value), linewidth = 0.25, color = "grey30") +
      #   geom_hline(yintercept = 0, linewidth = 0.3, color = "black") +
      #   scale_y_continuous(limits = c(-1.5,1.5)) +
      #   labs(x = NULL, y = "Eigenvector")
      
      p <- base +
        geom_col(aes(y = value, fill = value > 0), width = 40000, alpha = 0.7, show.legend = FALSE) +
        scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "dodgerblue3"),
                          guide = "none") +
        geom_hline(yintercept = 0, linewidth = 0.3, color = "black") +
        scale_y_continuous(limits = c(-1.5, 1.5)) +
        labs(x = NULL, y = "Eigenvector") +
        theme(legend.position = "none")
      
    } else if (feat %in% c("pi_Mpol", "pi_Mlys")) {
      p <- base +
        geom_point(aes(x = x_mid, y = value),
                   colour = feat_colors[[feat]], alpha = 0.2, size = 0.5) +
        # geom_line(aes(y = value),
        #           linewidth = 0.25, color = feat_colors[[feat]]) +
        labs(x = NULL, y = bquote(pi)) +
        ylim(c(0, 0.04))
      
      
    } else if (feat == "Recombination_rate") {
      p <- base +
        geom_point(aes(x = x_mid, y = value),
                   colour = feat_colors[[feat]], alpha = 0.2, size = 0.5) +
        # geom_area(aes(y = value),
        #           fill = feat_colors[[feat]], alpha = area_shade) +
        # geom_line(aes(y = value),
        #           linewidth = 0.25, color = feat_colors[[feat]]) +
        labs(x = NULL, y = "Recombination rate") +
        ylim(c(0, 12))
      
    } else if  (feat == "GC") {
      p <- base +
        geom_point(aes(x = x_mid, y = value),
                   colour = feat_colors[[feat]], alpha = 0.2, size = 0.5) +
        # geom_area(aes(y = value*100),
        #           fill = feat_colors[[feat]], alpha = area_shade) +
        # geom_line(aes(y = value*100),
        #           linewidth = 0.25, color = feat_colors[[feat]]) +
        labs(x = NULL, y = paste0(feat, " (%)")) +
        scale_y_continuous(limits = c(min(fdf[fdf$feature=="GC", "value"], na.rm = T)*0.8*100, max(fdf[fdf$feature=="GC", "value"], na.rm = T)*1.2*100))
      
    } else if (feat == "Total_repeats"){
      # Genes, Repeats — percentage tracks
      p <- base +
        geom_point(aes(x = x_mid, y = value*100),
                   colour = feat_colors[[feat]], alpha = 0.2, size = 0.5) +
        # geom_area(aes(y = value*100),
        #           fill = feat_colors[[feat]], alpha = area_shade) +
        #geom_line(aes(y = value*100),
        #         linewidth = 0.25, color = feat_colors[[feat]]) +
        labs(x = NULL, y = paste0(feat, " (%)")) +
        ylim(c(0,50))
      
    } else {
      # Genes, Repeats — percentage tracks
      p <- base +
        geom_point(aes(x = x_mid, y = value*100),
                   colour = feat_colors[[feat]], alpha = 0.2, size = 0.5) +
        # geom_area(aes(y = value*100),
        #           fill = feat_colors[[feat]], alpha = area_shade) +
        #geom_line(aes(y = value*100),
        #         linewidth = 0.25, color = feat_colors[[feat]]) +
        labs(x = NULL, y = paste0(feat, " (%)"))
    }
    
    if (feat %in% c("Fst", "dxy")) {
      
      base_fst <- ggplot(fdf, aes(x = x_mid, group = chr)) +
        # geom_vline(data = size_off, aes(xintercept = x_start),
        #            colour = "grey80", linewidth = 0.2) +
        geom_rect(data = chr_bg,
                  aes(xmin = x_start, xmax = x_end, ymin = -Inf, ymax = Inf),
                  fill = chr_bg$fill, alpha = 0.5, inherit.aes = FALSE) +
        geom_vline(data = break_off,
                   aes(xintercept = x_start, colour = event_timing, linetype = event_timing),
                   linewidth = 0.6) +
        scale_x_continuous(breaks = chr_labels$chr_mid,
                           labels = chr_labels$chr,
                           expand = c(0.01, 0)) +
        scale_colour_manual(name = "Breakpoint", 
                            values = c("#EC7014", "#2166AC", "black"), 
                            labels = c(expression(italic("M. lysimnia")), expression(italic("M. polymnia")), "Ancestral")) +
        scale_linetype_manual(name = "Breakpoint", 
                              labels = c(expression(italic("M. lysimnia")), expression(italic("M. polymnia")), "Ancestral"),
                              values = c("dashed", "dashed", "dotted")) +
        theme_minimal() +
        theme(axis.text.x = element_blank(),
              axis.title.y = element_text(size = text_size),
              axis.text.y = element_text(size = text_size),
              panel.grid.minor = element_blank(),
              panel.grid.major.x = element_blank(),
              legend.text = element_text(size = text_size),
              legend.title = element_text(size = text_size))
      
      
    }
    
    
    if (feat %in% c("Fst")) {
      p <- base_fst +
        geom_point(aes(x = x_mid, y = value),
                   colour = feat_colors[[feat]], alpha = 0.2, size = 0.5) +
        # geom_line(aes(y = value),
        #           linewidth = 0.25, color = feat_colors[[feat]]) +
        labs(#title = bquote("Reference " * italic(.(species_name))), 
          x = NULL, y = expression(F[ST])) +
        theme(legend.position = "none")
    } else if (feat %in% c("dxy")) {
      
      p <- base_fst +
        geom_point(aes(x = x_mid, y = value),
                   colour = feat_colors[[feat]], alpha = 0.2, size = 0.5) +
        # geom_line(aes(y = value),
        #           linewidth = 0.25, color = feat_colors[[feat]]) +
        labs(x = NULL, y = expression(D[XY]))
    }
    
    
    # Add per-chromosome smooth if requested
    if (smooth && feat != "eig_scores") {
      y_expr <- if (feat %in% c("Genes", "Total_repeats", "GC")) quote(value * 100) else quote(value)
      p <- p + geom_smooth(aes(y = !!y_expr),
                           method = "loess", span = smooth_span, se = FALSE,
                           linewidth = 0.4, color = "black", alpha = 1)
    }
    
    p
  }
  
  # Build tracks
  features <- unique(feat_off$feature)
  if (!is.null(feature_order)) {
    features <- feature_order[feature_order %in% features]
  }
  
  track_plots <- map(features, function(feat) {
    fdf <- feat_off %>% filter(feature == feat)
    make_track(fdf, feat, size_off, break_off, chr_labels, area_shade, smooth_span)
  })
  
  # Last track gets x-axis labels
  track_plots[[length(track_plots)]] <- track_plots[[length(track_plots)]] +
    theme(axis.text.x = element_text(size = text_size))
  
  # Stack all
  if (plot_paint){
    Reduce(`/`, c(list(p_paint), track_plots)) +
      plot_layout(heights = c(1, rep(1, length(track_plots))))
    
  } else {  
    Reduce(`/`, c(track_plots)) +
      plot_layout(heights = c(rep(1, length(track_plots))))  
  }  
  
}


# ── Build both genomes side by side ──────────────────────────────────────────
# Filter size_df per genome if needed, or pass the same one

feature_order_q1 <- c("Recombination_rate", "Genes", "Total_repeats", "pi_Mpol")
combined_df_filt_q1 <- combined_df %>% filter(feature %in% feature_order_q1)
size_q1 <- size_df %>% filter(genome=="ilMecPoly1")
p_q1 <- plot_one_genome("Mechanitis polymnia", "ilMecPoly1", size_q1, anc_df, 
                        combined_df_filt_q1, break_df %>% filter(event_timing != "shared_event_n4" & type != "Complex"),
                        feature_order = feature_order_q1,
                        smooth = TRUE,
                        smooth_span=0.2,
                        area_shade=0.2,
                        text_size = 12)

feature_order_q2 <- c("Recombination_rate", "Genes", "Total_repeats", "pi_Mlys")
combined_df_filt_q2 <- combined_df %>% filter(feature %in% feature_order_q2)
size_q2 <- size_df %>% filter(genome=="ilMecLysi212")
p_q2 <- plot_one_genome("Mechanitis lysimnia","ilMecLysi212", size_q2, anc_df, 
                        combined_df_filt_q2, break_df %>% filter(event_timing != "shared_event_n4" & type != "Complex"),
                        feature_order = feature_order_q2,
                        smooth = TRUE,
                        smooth_span=0.2,
                        area_shade=0.2,
                        text_size = 12)



feature_order_q3 <- c("Fst", "dxy")
combined_df_filt_q3 <- combined_df %>% filter(feature %in% feature_order_q3)
p_q3 <- plot_one_genome("Mechanitis polymnia","ilMecPoly1", size_q1, anc_df, 
                        combined_df_filt_q3, break_df %>% filter(type != "Complex" & event_timing != "shared_event_n4"),
                        feature_order = feature_order_q3,
                        smooth = TRUE,
                        smooth_span=0.1,
                        area_shade=0.2, plot_paint=FALSE,
                        text_size = 12)



plots <- lapply(p_q2, function(p) p + theme(axis.title.y = element_blank()))
p_q2_corr <- wrap_plots(plots, ncol = 1, guides = "collect")


genome_tracks <- p_q1 | p_q2_corr 

genome_tracks <- genome_tracks + plot_layout(guides = "collect") +
  theme(legend.position = "right",
        legend.direction = "vertical") + 
  guides(fill = guide_legend(ncol = 2, title.position = "top",
                             override.aes = list(shape = NA, size = 0)))

genome_tracks_legend <- get_legend(p_q2_corr)

legend_pop <- get_legend(p_q3 + theme(legend.direction = "horizontal",
                                      legend.title.position = "top"))


fig_2 <- ggarrange(genome_tracks, 
                   ggarrange(
  ggarrange(legend_pop, 
            (p_q3 + theme(legend.position = "none",
                          plot.margin = margin(30,1,1,1))), 
            nrow = 2, heights = c(1,3)),
  genome_tracks_legend, 
  nrow = 1, widths = c(1, 1.25)
), 
nrow = 2, heights = c(1.5,1), labels = c("A", "C"))

ggsave(plot = fig_2, paste("plots/Fig_2_", Sys.Date(), ".pdf", sep = ""),
       width = 18, height = 12)

```
