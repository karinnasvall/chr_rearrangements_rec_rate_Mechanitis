# Functions to evaluate and plot the output from syngraph
# Karin Näsvall 2025-07-05

get_cluster_all <- function(prefix = "output_folder", 
                            file_prefix = "infer_m", 
                            param_range = 5:60) {
  cluster_all <- NULL
  for (i in stringr::str_pad(param_range, 2, pad = "0")) {
    file_path <- paste0(prefix, i, "/", file_prefix, i, ".clusters.tsv")
    clusters <- read.table(file = file_path)
    clusters$param_nr <- i
    if (i == min(stringr::str_pad(param_range, 2, pad = "0"))){
      cluster_all <- clusters
    } else {
      cluster_all <- rbind(cluster_all, clusters)
    }
  }
  return(cluster_all)
}

# Example usage:
#param_prefix="output_2025-07-05/2025-07-05_"
#output_folder=paste(work_dir, param_prefix, sep="/")
#param_range= 5:60

#cluster_all <- get_cluster_all(param_prefix, "test_", param_range)

plot_linkage_clusters <- function(cluster_data, outfile_prefix = "plots/linkage_cluster", height= 8, width = 16) {
  p <- ggplot(cluster_data) +
    geom_bar(aes(as.numeric(param_nr), V3), stat = "identity") +
    facet_wrap(~V2, ncol = 6) +
    theme_pubr()
  ggsave(
    plot = p,
    filename = paste(outfile_prefix, Sys.Date(), "pdf", sep = "."),
    height = height, width = width
  )
}

# Example usage:
#plot_linkage_clusters(cluster_all)


#concatenates the rearrangement files from syngraph output to a table for plotting and writes a summary table
get_rearrangement_table <- function(param_prefix, param_range = 5:60, file_prefix = "infer_m", out_sum_table = "tables/summary_rearrangements") {
  rearr_all <- NULL
  padded_range <- stringr::str_pad(param_range, 2, pad = "0")
  for (i in padded_range) {
    file_path <- paste0(param_prefix, i, "/", file_prefix, i, ".rearrangements.tsv")
    rearr <- read.csv2(file_path, sep = "\t", header = TRUE)
    rearr_sum <- rearr %>%
      group_by(X.parent, child, event) %>%
      summarise(rearr = sum(multiplicity), .groups = "drop")
    rearr_sum$param_nr <- i
    if (is.null(rearr_all)) {
      rearr_all <- rearr_sum
    } else {
      rearr_all <- rbind(rearr_all, rearr_sum)
    }
  }
    
  # Get all unique values for each variable for full grid (to include 0 rearrangements)
  # Only observed (X.parent, child) pairs, but all events and param_nr
  parent_child <- unique(rearr_all[, c("X.parent", "child")])
  all_event <- unique(rearr_all$event)
  all_param <- unique(rearr_all$param_nr)
  full_grid <- merge(parent_child, expand.grid(event = all_event, param_nr = all_param, stringsAsFactors = FALSE))
  rearr_all_full <- dplyr::left_join(full_grid, rearr_all, by = c("X.parent", "child", "event", "param_nr"))
  rearr_all_full$rearr[is.na(rearr_all_full$rearr)] <- 0


  # Write summary table
  summary_tbl <- rearr_all_full %>%
    group_by(child, event) %>%
    summarise(
      mean_rearr = mean(rearr),
      median_rearr = median(rearr),
      sd_rearr = sd(rearr),
      max_rearr = max(rearr),
      min_rearr = min(rearr),
      .groups = "drop"
    )
  write.csv(summary_tbl, file = paste(out_sum_table, Sys.Date(), "tsv", sep="."), row.names = FALSE)
  return(rearr_all_full)
}



plot_rearrangements_param <- function(rearr_all, 
                                      outfile_prefix = "plots/rearrangements_param_test", 
                                      height = 8, width = 16) {
  p <- ggplot(rearr_all) +
    geom_bar(aes(as.numeric(param_nr), rearr, fill = event), stat = "identity") +
    scale_fill_manual(values = c("goldenrod", "tan4"),
                      labels = c("Fission", "Fusion")) +
    facet_wrap(~child) +
    theme_pubr() +
    theme(
      panel.grid.major.x = element_line(linewidth = 0.2),
      panel.grid.minor.x = element_line(linewidth = 0.2),
      legend.position = "top",
      legend.text = element_text(size = 14),
      legend.title = element_blank()
    )
  ggsave(
    plot = p,
    filename = paste(outfile_prefix, Sys.Date(), "pdf", sep = "."),
    height = height, width = width
  )
  return(p)
}

get_nr_markers <- function(prefix = "output_folder", 
                            file_prefix = "infer_", 
                            param_range = 5:60) {
  nr_markers_all <- NULL
  for (i in stringr::str_pad(param_range, 2, pad = "0")) {
    file_path <- paste0(prefix, i, "/", file_prefix, i, ".table.tsv")
    print(i)
    nr_markers <- read.table(file = file_path, header = T, comment.char = "")
    nr_markers <- nr_markers %>% 
      dplyr::rename(marker = `X.marker`) %>% 
      select(marker, ends_with("_seq"))
    nr_markers$param_nr <- i
    nr_markers <- nr_markers %>% pivot_longer(ends_with("_seq"), values_to = "chr", names_to = "node")
    nr_markers_summary <- nr_markers%>%   filter(!is.na(chr)) %>% 
      group_by(node, chr, param_nr) %>% count()
    nr_markers_summary$node <- sub("_seq", "", nr_markers_summary$node)
    
    if (i == min(stringr::str_pad(param_range, 2, pad = "0"))){
      nr_markers_all <- nr_markers_summary
    } else {
      nr_markers_all <- rbind(nr_markers_all, nr_markers_summary)
    }
  }
  write.csv(nr_markers_all, file = paste(prefix, "nr_markers_chr_all.csv", sep = ""),
          quote = F, row.names = F)
  return(nr_markers_all)
}

plot_syngraph_tree_features <- function(tree, rearr_all, cluster_all, param_values, 
                                        out_prefix = "plots/Syngraph_clado_param", 
                                        rearr_ylim = c(0, 60), 
                                        bar_cols = c("goldenrod", "tan4"),
                                        legend_plot = NULL) {
  require(ggtree)
  require(ggimage)
  require(ggpubr)
  require(dplyr)
  require(tidyr)
  
  # Prepare wide rearrangement summary for nodebar
  rearr_sum_wide <- data.frame(tidyr::pivot_wider(rearr_all, names_from = event, values_from = rearr, values_fill = 0))
  rearr_sum_wide <- dplyr::inner_join(
    rearr_sum_wide, 
    data.frame(child = c(tree$tip.label, tree$node.label), node = seq(1, length(c(tree$tip.label, tree$node.label))))
  )
  
  # Prepare chromosome cluster info
  anc_chr <- cluster_all[, 2:4]
  colnames(anc_chr) <- c("Node", "Chr_nr", "param_nr")
  
  # If legend_plot is not provided, create a dummy for the legend
  if (is.null(legend_plot)) {
    legend_plot <- ggplot(rearr_all) +
      geom_bar(aes(as.numeric(param_nr), rearr, fill = event), stat = "identity") +
      scale_fill_manual(values = bar_cols, labels = c("Fission", "Fusion")) +
      theme_pubr()
    fill_legend <- get_legend(legend_plot)
  } else {
    fill_legend <- get_legend(legend_plot)
  }
  
  # Loop over parameter values
  for (i in stringr::str_pad(param_values, 2, pad = "0")) {
    # Node bars for current param
    bars <- nodebar(rearr_sum_wide[rearr_sum_wide$param_nr == i, 4:6], cols = 1:2, position = "dodge")
    bars <- lapply(bars, function(g) 
      g + scale_fill_manual(values = bar_cols) +
        geom_text(aes(label = .data$value), colour = "black", size = 5, vjust = -0.15, position = position_dodge(.9)) +
        scale_y_continuous(limits = rearr_ylim)
    )
    
    # Cladogram with chromosome number labels
    p_clado <- ggtree(tree, ladderize = TRUE, branch.length = "none", layout = "rectangular", right = TRUE) %<+% 
      anc_chr[anc_chr$param_nr == i, ] +
      geom_tiplab(size = 5) +
      #geom_nodelab(nudge_x = 0.15, size = 5) +
      geom_label(aes(label = Chr_nr), size = 5, nudge_x = 0.3, nudge_y = 0, fill = "yellow", alpha = 0.9) +
      xlim(0, 14)
    
    p_plot_cladogram <- p_clado + 
      geom_inset(bars, width = .05, height = .12, x = "branch", vjust = -1.2, hjust = 0.2) +
      theme(legend.position = "top")
    
    ggsave(
      ggarrange(p_plot_cladogram, fill_legend, nrow = 2, heights = c(7, 1)),
      filename = paste0(out_prefix, i, ".pdf"),
      height = 18, width = 12
    )
  }
}

clean_chr <- function(x) {
  # Remove everything up to and including the last underscore
  stripped <- sub(".*_", "", x)
  # Zero-pad numbers (e.g. "1" -> "01"), leave letters as-is
  ifelse(grepl("^\\d+$", stripped),
         sprintf("%02d", as.integer(stripped)),
         stripped)
}

# chromosome painting by ancestral elements

get_link_long <- function(i, anc_busco, param_prefix = param_prefix, file_prefix = "infer_m", tree = tree) {

  link_table <- read.table(paste0(param_prefix, i, "/", file_prefix, i, ".table.tsv"),
                           header = FALSE, sep = "\t")

  node_cols <- names(link_table)[sapply(link_table, function(x)
    any(grepl("^n\\d+_\\d+$", x), na.rm = TRUE))]

  tip_cols <- names(link_table)[sapply(link_table, function(x)
    any(grepl("^il.*SUPER", x), na.rm = TRUE))]

  link_long <- link_table %>%
    select(V1, all_of(c(tip_cols, node_cols))) %>%
    pivot_longer(-V1, names_to = "col", values_to = "chr_assignment") %>%
    filter(!is.na(chr_assignment)) %>%
    mutate(
      Node = gsub("_.*", "", chr_assignment),
      chr  = clean_chr(gsub(".*_", "", chr_assignment))
    ) %>%
    left_join(anc_busco, by = c("V1" = "geneid")) %>%
    select(-col, -chr_assignment) %>% dplyr::left_join(., 
    data.frame(Node = c(tree$tip.label, tree$node.label), node = seq(1, length(c(tree$tip.label, tree$node.label))))
  ) %>% filter(!is.na(anc_element)) 

  chr_order <- link_long %>%
    count(Node, chr) %>%
    arrange(Node, desc(n))

  link_long %>% left_join(chr_order, by = c("Node", "chr"))
}

plot_syngraph_tree_features_chr_painting <- function(tree, rearr_all, chr_per_param, param_values, 
                                        out_prefix = "plots/Syngraph_clado_param", 
                                        rearr_ylim = c(0, 30),
                                        plot_xlim =  c(0,10),
                                        plot_ylim =  c(0,10),
                                        bar_cols = c("goldenrod", "tan4"),
                                        legend_plot = NULL,
                                        link_df = NULL,
                                        anc_busco = NULL,
                                        param_prefix = NULL,
                                        anc_colours = NULL,
                                        agora = FALSE) {
  require(ggtree)
  require(ggimage)
  require(ggpubr)
  require(dplyr)
  require(tidyr)
  
  # Prepare wide rearrangement summary for nodebar
  rearr_sum_wide <- data.frame(tidyr::pivot_wider(rearr_all, names_from = event, values_from = rearr, values_fill = 0))
  rearr_sum_wide <- dplyr::inner_join(
    rearr_sum_wide, 
    data.frame(child = c(tree$tip.label, tree$node.label), node = seq(1, length(c(tree$tip.label, tree$node.label))))
  )
  
  # Prepare chromosome cluster info
  anc_chr <- chr_per_param %>%
    rename(Node = node)
  
  #   anc_chr <- cluster_all[, 2:4]
  # colnames(anc_chr) <- c("Node", "Chr_nr", "param_nr")
  
  # If legend_plot is not provided, create a dummy for the legend
  if (is.null(legend_plot)) {
    legend_plot <- ggplot(rearr_all) +
      geom_bar(aes(as.numeric(param_nr), rearr, fill = event), stat = "identity") +
      scale_fill_manual(values = bar_cols, labels = c("Fission", "Fusion")) +
      theme_pubr()
    fill_legend <- get_legend(legend_plot)
  } else {
    fill_legend <- get_legend(legend_plot)
  }

  # Loop over parameter values
  for (i in stringr::str_pad(param_values, 2, pad = "0")) {
    # Node bars for current param
    bars <- nodebar(rearr_sum_wide[rearr_sum_wide$param_nr == i, 4:6], cols = 1:2, position = "dodge")
    bars <- lapply(bars, function(g) 
      g + scale_fill_manual(values = bar_cols) +
        geom_text(aes(label = .data$value), colour = "black", size = 5, vjust = -0.15, position = position_dodge(.9)) +
        scale_y_continuous(limits = rearr_ylim)
    )

    # Cladogram with chromosome number labels
    p_clado <- ggtree(tree, ladderize = TRUE, branch.length = "none", layout = "rectangular", right = TRUE) %<+% 
      anc_chr[anc_chr$param_nr == i, ] +
      geom_tiplab(size = 5, #offset = 0.1,
                  nudge_x = 0.1, nudge_y = -0.2) +
      geom_label(aes(label = Chr_nr), size = 5, 
                nudge_x = -0.15, nudge_y = 0.15
                ) +
      xlim(plot_xlim) #+
      # ylim(plot_ylim)
      
    
    p_plot_cladogram <- p_clado + 
      geom_inset(bars, 
                 width = .08, height = .12, 
                 x = "branch",
                 vjust = -0.52, hjust = 0.1
                 ) +
      theme(legend.position = "top")

  
    if (agora == TRUE) {

      link_long <- link_df


      bars_chr <- lapply(sort(unique(link_long$node)), function(nd) {
      ggplot(link_long[link_long$node == nd, ]) +
        geom_rect(aes(x = chr, width = 0.8, ymin = start, ymax = end, fill = anc_element, colour = anc_element)) +
        scale_fill_manual(values = anc_palette) +
        scale_colour_manual(values = anc_palette) +
        theme_void() +
        theme(legend.position = "none")
    })

    names(bars_chr) <- sort(unique(link_long$node))




    } else {
    # if a link file is given as a parameter, use that. But if not, prepare the long format table for chromosome painting from the table file in the syngraph output for the current parameter value, and the anc_busco table for ancestral element info.
    if (is.null(link_df)) {
      link_long <- get_link_long(i = i, anc_busco = anc_busco, param_prefix = param_prefix, tree = tree)
    } else {
      #or if a link_file is already prepared:
      link_long <- link_df
    }

    anc_palette <- anc_colours$anc_colour
    names(anc_palette) <- anc_colours$anc_element

    bars_chr <- lapply(sort(unique(link_long$node)), function(nd) {
      ggplot(link_long[link_long$node == nd, ],
             aes(x = reorder(chr, -n), fill = anc_element)) +
        geom_bar(alpha = 0.7) +
        scale_fill_manual(values = anc_palette) +
        theme_void() +
      theme(legend.position = "none")
    })
    names(bars_chr) <- sort(unique(link_long$node))
    } 

    max_width <- 0.15  # inset width for the node with most chromosomes
    node_gap  <- 0.02 # fixed distance from node to near edge of inset

    n_chr <- sapply(names(bars_chr), function(nd) {
      length(unique(link_long$chr[link_long$node == nd]))
    })

    widths <- setNames(n_chr / max(n_chr) * max_width, names(bars_chr))

    for (nd in names(bars_chr)) {
      # keep near edge at constant distance: hjust = -(node_gap + width/2)
      hjust_nd <- -(node_gap + widths[[nd]]*4)

      p_plot_cladogram <- p_plot_cladogram +
        geom_inset(bars_chr[nd],
                   x = "node",
                   vjust = -0.1, hjust = hjust_nd,
                   width = widths[[nd]], height = 0.05)
    }


  
    legend_plot <- ggplot(data.frame(anc_element = names(anc_palette)),
                      aes(x = 1, y = anc_element, fill = anc_element)) +
    geom_tile() +
    scale_fill_manual(values = anc_palette) +
    theme_void() +
    guides(fill = guide_legend(title = "Ancestral element", direction = "horizontal", nrow = 2))

    anc_legend <- get_legend(legend_plot)

    ggsave(
      ggarrange(anc_legend, p_plot_cladogram, fill_legend, nrow = 3, heights = c(1, 7, 1)),
      filename = paste0(out_prefix, i, "_painting.pdf"),
      height = 20, width = 20
    )
  }
}


plot_syngraph_tree_features_chr_painting_one_plot <- function(tree, rearr_all, chr_per_param, param_values, 
                                        out_prefix = "plots/Syngraph_clado_param", 
                                        rearr_ylim = c(0, 30),
                                        plot_xlim =  c(0,10),
                                        plot_ylim =  c(0,10),
                                        bar_cols = c("goldenrod", "tan4"),
                                        legend_plot = NULL,
                                        link_df = NULL,
                                        anc_busco = NULL,
                                        param_prefix = NULL,
                                        anc_colours = NULL,
                                        agora = FALSE, text_size = 5) {
  require(ggtree)
  require(ggimage)
  require(ggpubr)
  require(dplyr)
  require(tidyr)
  
  # Prepare wide rearrangement summary for nodebar
  rearr_sum_wide <- data.frame(tidyr::pivot_wider(rearr_all, names_from = event, values_from = rearr, values_fill = 0))
  rearr_sum_wide <- dplyr::inner_join(
    rearr_sum_wide, 
    data.frame(child = c(tree$tip.label, tree$node.label), node = seq(1, length(c(tree$tip.label, tree$node.label))))
  )

  # Event column names (everything between param_nr and node in rearr_sum_wide)
  event_cols <- setdiff(names(rearr_sum_wide)[4:ncol(rearr_sum_wide)], "node")

  event_cols_map <- c(
  fission = "goldenrod",
  fusion  = "tan4"
)

  # Range across all param values per node/event — used for error bars
  rearr_range <- rearr_all %>%
    dplyr::group_by(child, event) %>%
    dplyr::summarise(min_rearr = min(rearr), max_rearr = max(rearr), .groups = "drop") %>%
    dplyr::inner_join(
      data.frame(child = c(tree$tip.label, tree$node.label),
                 node  = seq_len(length(c(tree$tip.label, tree$node.label)))),
      by = "child"
    ) %>%
    dplyr::mutate(event = factor(event, levels = event_cols))
  
  # Prepare chromosome cluster info
  anc_chr <- chr_per_param %>%
    rename(Node = node) 

  anc_chr_range <- anc_chr %>% 
        group_by(Node) %>%
        summarise(Chr_nr_max = max(Chr_nr),
                  Chr_nr_min = min(Chr_nr))

  anc_chr <- left_join(anc_chr, anc_chr_range, by = "Node") %>%
    mutate(Chr_nr_label = ifelse(Chr_nr_max == Chr_nr_min, as.character(Chr_nr_max), paste0(Chr_nr_min, "-", Chr_nr_max)))

                
  
  # If legend_plot is not provided, create a dummy for the legend
  if (is.null(legend_plot)) {
    legend_plot <- ggplot(rearr_all) +
      geom_bar(aes(as.numeric(param_nr), rearr, fill = event), stat = "identity") +
      scale_fill_manual(name = "Event", values = event_cols_map, labels = c(fission = "Fission", fusion = "Fusion")) +
      theme_pubr() +
      theme(legend.text = element_text(size = text_size*2.8), 
            legend.title = element_text(size = text_size*2.8),
            legend.title.position = "top",
            legend.direction = "vertical")
    fill_legend <- get_legend(legend_plot) 
  } else {
    fill_legend <- get_legend(legend_plot)
  }

  # Loop over parameter values
  for (i in stringr::str_pad(param_values, 2, pad = "0")) {
    # Build bar + errorbar plots per node (avoids scale conflicts from nodebar + geom_errorbar)
    node_ids <- rearr_sum_wide[rearr_sum_wide$param_nr == i, "node"]
    bars <- setNames(lapply(node_ids, function(nd) {
      nd_data <- rearr_sum_wide[rearr_sum_wide$param_nr == i & rearr_sum_wide$node == nd, ] %>%
        tidyr::pivot_longer(cols = dplyr::all_of(event_cols), names_to = "event", values_to = "rearr") %>%
        dplyr::mutate(event = factor(event, levels = event_cols))
      nd_range <- rearr_range[rearr_range$node == nd, ]

      ggplot(nd_data, aes(x = event, y = rearr, fill = event)) +
        geom_col(position = position_dodge(0.9)) +
        geom_errorbar(data = nd_range,
                      aes(x = event, ymin = min_rearr, ymax = max_rearr),
                      width = 0.5, linewidth = 0.8, inherit.aes = FALSE, colour = "grey20") +
        geom_text(aes(label = rearr), colour = "black", size = text_size, vjust = -1, hjust = 0.5,
                  position = position_dodge(0.9)) +
        scale_fill_manual(values = event_cols_map) +
        scale_y_continuous(limits = rearr_ylim) +
        ggtree::theme_inset()
    }), as.character(node_ids))

    # Cladogram with chromosome number labels
    p_clado <- ggtree(tree, ladderize = TRUE, branch.length = "none", layout = "rectangular", right = TRUE) %<+% 
      anc_chr[anc_chr$param_nr == i, ] +
      geom_tiplab(size = text_size, #offset = 0.1,
                  nudge_x = 0.1, nudge_y = -0.2) +
      geom_label(aes(label = Chr_nr_label), size = text_size, 
                nudge_x = -0.15, nudge_y = 0.15
                ) +
      xlim(plot_xlim) #+
      # ylim(plot_ylim)
      
    
    p_plot_cladogram <- p_clado + 
      geom_inset(bars, 
                 width = .08, height = .12, 
                 x = "branch",
                 vjust = -0.52, hjust = 0.1
                 ) +
      theme(legend.position = "top")

    if (agora == TRUE) {

      link_long <- link_df


      bars_chr <- lapply(sort(unique(link_long$node)), function(nd) {
      ggplot(link_long[link_long$node == nd, ]) +
        geom_rect(aes(x = chr, width = 0.8, ymin = start, ymax = end, fill = anc_element, colour = anc_element)) +
        scale_fill_manual(values = anc_palette) +
        scale_colour_manual(values = anc_palette) +
        theme_void() +
        theme(legend.position = "none")
    })

    names(bars_chr) <- sort(unique(link_long$node))




    } else {
    # if a link file is given as a parameter, use that. But if not, prepare the long format table for chromosome painting from the table file in the syngraph output for the current parameter value, and the anc_busco table for ancestral element info.
    if (is.null(link_df)) {
      link_long <- get_link_long(i = i, anc_busco = anc_busco, param_prefix = param_prefix, tree = tree)
    } else {
      #or if a link_file is already prepared:
      link_long <- link_df
    }

    anc_palette <- anc_colours$anc_colour
    names(anc_palette) <- anc_colours$anc_element

    bars_chr <- lapply(sort(unique(link_long$node)), function(nd) {
      ggplot(link_long[link_long$node == nd, ],
             aes(x = reorder(chr, -n), fill = anc_element)) +
        geom_bar(alpha = 0.7) +
        scale_fill_manual(values = anc_palette) +
        theme_void() +
      theme(legend.position = "none")
    })
    names(bars_chr) <- sort(unique(link_long$node))
    } 

    max_width <- 0.15  # inset width for the node with most chromosomes
    node_gap  <- 0.02 # fixed distance from node to near edge of inset

    n_chr <- sapply(names(bars_chr), function(nd) {
      length(unique(link_long$chr[link_long$node == nd]))
    })

    widths <- setNames(n_chr / max(n_chr) * max_width, names(bars_chr))

    for (nd in names(bars_chr)) {
      # keep near edge at constant distance: hjust = -(node_gap + width/2)
      hjust_nd <- -(node_gap + widths[[nd]]*4)

      p_plot_cladogram <- p_plot_cladogram +
        geom_inset(bars_chr[nd],
                   x = "node",
                   vjust = -0.1, hjust = hjust_nd,
                   width = widths[[nd]], height = 0.05)
    }


  
    legend_plot <- ggplot(data.frame(anc_element = names(anc_palette)),
                      aes(x = 1, y = anc_element, fill = anc_element)) +
    geom_tile() +
    scale_fill_manual(values = anc_palette) +
    theme_void() +
          theme(legend.text = element_text(size = text_size*2.8), legend.title = element_text(size = text_size*2.8)) +
          guides(fill = guide_legend(title = "Ancestral element", title.position = "top", label.position = "bottom", direction = "vertical", nrow = 2))

    anc_legend <- get_legend(legend_plot)

    ggsave(
      ggarrange(
        ggarrange(fill_legend, anc_legend, ncol = 2, widths = c(1,2), align = "h"),
         p_plot_cladogram, nrow = 2, heights = c(1,7)),
      filename = paste0(out_prefix, i, "_painting.pdf"),
      height = 12, width = 14
    )
  }
}


