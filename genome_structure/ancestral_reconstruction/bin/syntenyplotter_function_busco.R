# ─────────────────────────────────────────────────────────────────────────────
# prep_and_plot_synteny()
#
# Args:
#   taxa1           : reference taxon name (must match Node column in busco_table)
#   taxa2           : query taxon name
#   busco_table     : data.frame from agora_link_parsed_flipped.csv
#   output_path     : root directory; must contain intermediate/ and plots/ subdirs
#   anc_colour_file : path to two-column TSV with anc_element / hex colour
#   fileformat      : passed to draw.linear (default "pdf")
#   w, h            : plot dimensions (default 13 x 4)
#   show_chromosome_names : logical (default TRUE)
# ─────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(ggplot2)

prep_and_plot_synteny <- function(
    taxa1,
    taxa2,
    busco_table,
    output_path,
    anc_colour_file,
    fileformat = "pdf",
    w = 13,
    h = 4,
    show_chromosome_names = TRUE
) {

  print("prep_and_plot_synteny() called")
  print(paste("taxa1:", taxa1, "taxa2:", taxa2))


  best_chr_by_alignment <- function(df, chr_col, score_col = total_aligned) {
    chr_col   <- enquo(chr_col)
    score_col <- enquo(score_col)
    df %>%
      mutate(.row_id = row_number()) %>%
      group_by(!!chr_col) %>%
      slice_max(!!score_col, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      arrange(.row_id) %>%
      select(-.row_id)
  }

  # ── build chain table ──────────────────────────────────────────────────────
  print(paste("Building chain table for:", taxa1, "and", taxa2, "..."))
  
  chain_table <- busco_table %>%
    filter(Node %in% taxa1) %>%    
    inner_join(
      busco_table %>% filter(Node %in% taxa2),
      by = c("geneid", "anc_element")
    ) %>%
    mutate(strand = if_else(strand.x == strand.y, "+", "-")) %>%
    rename(
      reference   = chr.x,
      Rstart      = start.x,
      Rend        = end.x,
      query       = chr.y,
      Qstart      = start.y,
      Qend        = end.y,
      refID       = Node.x,
      queryID     = Node.y,
      Rseq_length = length.x,
      Qseq_length = length.y
    ) %>%
     select(!ends_with(".x") & !ends_with(".y")) %>%
     select(reference, Rstart, Rend, query, Qstart, Qend,
            strand, refID, queryID, Rseq_length, Qseq_length, anc_element, geneid) %>%
     arrange(desc(Rseq_length))

    print(head(chain_table))

  # ── scale reconstructed (gene-step) coordinates ───────────────────────────
  # if any sequences are shorter than 10 kb, we assume they are in gene coordinates and scale them up to 50 kb (arbitrary choice to make the plot look nice)
  print("Scaling reference coordinates...")

  if (any(chain_table$Rseq_length < 10000)) {
    chain_table$Rseq_length <- chain_table$Rseq_length * 50000
    chain_table$Rstart      <- chain_table$Rstart      * 50000
    chain_table$Rend        <- chain_table$Rend        * 50000
  } else {
    chain_table$Qend <- chain_table$Qend + 200000
  }

print("Scaling query coordinates...")
  if (any(chain_table$Qseq_length < 10000)) {
    chain_table$Qseq_length <- chain_table$Qseq_length * 50000
    chain_table$Qstart      <- chain_table$Qstart      * 50000
    chain_table$Qend        <- chain_table$Qend        * 50000
  } else {
    chain_table$Qend <- chain_table$Qend + 200000
  }

  # ── chromosome ordering by alignment coverage ─────────────────────────────
  print("Ordering chromosomes by alignment coverage...")

  chain_table <- chain_table %>%
    mutate(aln_len = abs(Rend - Rstart))

  chr_pair_lengths <- chain_table %>%
    group_by(reference, query) %>%
    summarise(
      total_aligned = sum(aln_len, na.rm = TRUE),
      n_blocks      = n(),
      .groups       = "drop"
    ) %>%
    arrange(reference, desc(total_aligned))

    print(head(chr_pair_lengths))

  chr_size_q <- unique(chain_table[, c("query", "Qseq_length", "queryID")])
  chr_size_q <- chr_size_q %>%
    arrange(match(query, best_chr_by_alignment(chr_pair_lengths, query)$query))
  colnames(chr_size_q) <- c("chr", "seq_length", "ID")

  chr_size_r <- unique(chain_table[, c("reference", "Rseq_length", "refID")])
  chr_size_r <- chr_size_r %>%
    arrange(match(reference, best_chr_by_alignment(chr_pair_lengths, reference)$reference))
  colnames(chr_size_r) <- c("chr", "seq_length", "ID")

  chr_size             <- rbind(chr_size_q, chr_size_r)
  chr_size$polymorphic <- "N"

  print(head(chr_size))


  # ── write intermediate files ───────────────────────────────────────────────
  size_file  <- file.path(output_path, "intermediate", paste0("chr_length_", taxa2, taxa1, ".txt"))
  chain_file <- file.path(output_path, "intermediate", paste0("chain_",      taxa1, taxa2, ".txt"))

  write.table(chr_size,   file = size_file,  sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(chain_table, file = chain_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

  print(paste("Intermediate files written to:", output_path, "intermediate/"))

  # ── colour palette ─────────────────────────────────────────────────────────
  anc_colours <- read.table(anc_colour_file, header = FALSE, sep = "\t",
                             col.names = c("anc_element", "anc_colour"),
                             comment.char = "")
  anc_palette <- setNames(anc_colours$anc_colour, anc_colours$anc_element)

  print(paste("Colour palette:", paste(names(anc_palette), anc_palette, sep = "=", collapse = ", ")))

  # ── plot ───────────────────────────────────────────────────────────────────
  source("~/scripts/recombination/syntenyplotter_agora.R")
  draw.linear(
    directory             = file.path(output_path, "plots"),
    output                = paste0("synt_", taxa1, taxa2, "_", Sys.Date()),
    size_file,
    chain_file,
    fileformat            = fileformat,
    w                     = w,
    h                     = h,
    show_chromosome_names = show_chromosome_names,
    anc_colours           = anc_palette
  )

  invisible(NULL)

  print(paste("Synteny plot generated for:", taxa1, "and", taxa2, "in", fileformat, "format."))
}
