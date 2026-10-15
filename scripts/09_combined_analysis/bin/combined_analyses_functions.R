# helper functions for combined analyses

library(dplyr)
library(ggplot2)


clean_chr <- function(x) {
  # Remove everything up to and including the last underscore
  stripped <- sub(".*_", "", x)
  # Zero-pad numbers (e.g. "1" -> "01"), leave letters as-is
  ifelse(grepl("^\\d+$", stripped),
         sprintf("%02d", as.integer(stripped)),
         stripped)
}




library(ggplot2)
library(dplyr)

# size_df: data.frame with columns chr, chr_len (in bp)
# track_df: data.frame with columns chr, start, end, value (+ any extra aesthetics)
# value_col: name of the column to plot on y axis (as string)

make_offsets <- function(size_df, gap = 0) {
  size_df %>%
    mutate(chr_len = end - start + 1) %>%
    arrange(chr) %>%
    mutate(
      offset  = cumsum(dplyr::lag(chr_len + gap, default = 0)),
      chr_mid = offset + chr_len / 2
    )
}

plot_genome_track <- function(
    track_df,
    size_df,
    y_col      = "value",
    colour_col = NULL,
    point      = FALSE,
    shape      = NULL,
    line       = TRUE,
    loess      = FALSE,
    vline       = FALSE,
    gap        = 1000000,        # gap between chromosomes in bp
    y_label    = y_col,
    title      = NULL,
    alpha      = 0.6,
    text_size  = 12
) {
  
  offsets <- make_offsets(size_df, gap = gap)
  
  # Build background rectangles (alternating)
  chr_bg <- offsets %>%
    mutate(
      x_start  = offset,
      x_end    = offset + chr_len,
      bg_group = factor(row_number() %% 2)
    )
  
  # Add genome-wide x coordinate to track data
  offset_map <- setNames(offsets$offset, offsets$chr)
  track <- track_df %>%
    mutate(
      x_mid = unname(offset_map[as.character(chr)] + (start + end) / 2)
    )
  
  # Chr label positions
  chr_labels <- offsets %>% select(chr, chr_mid)
  
  # Base plot
  p <- ggplot() +
    # alternating chr background
    geom_rect(
      data = chr_bg,
      aes(xmin = x_start, xmax = x_end, ymin = -Inf, ymax = Inf,
          fill = bg_group),
      alpha = 0.4, inherit.aes = FALSE
    ) +
    scale_fill_manual(
      values = c("0" = "grey88", "1" = "white"),
      guide = "none"
    )
  
  # Track layer
  aes_base <- aes(x = x_mid, y = .data[[y_col]])
  
  if (!is.null(colour_col)) {
    aes_base <- modifyList(aes_base, aes(colour = .data[[colour_col]]))
  }
  
  if (line)  p <- p + geom_line(data = track, mapping = aes_base, alpha = alpha)
  if (point) p <- p + geom_point(data = track, mapping = aes_base, alpha = alpha, size = 0.4, shape = shape)
  if (loess) p <- p + geom_smooth(
    data    = track,
    mapping = modifyList(aes_base, aes(group = chr)),
    method  = "loess",
    se      = FALSE,
    span    = 0.15
  )

  if (vline) {
    # Train the y scale when vertical markers are the only visible data layer.
    if (!line && !point && !loess) {
      p <- p + geom_point(
        data = track,
        mapping = aes_base,
        alpha = 0,
        size = 0,
        show.legend = FALSE,
        inherit.aes = FALSE
      )
    }

    vline_map <- aes(x = x_mid, xend = x_mid, y = -Inf, yend = Inf)

    if (!is.null(colour_col)) {
      vline_map <- modifyList(vline_map, aes(colour = .data[[colour_col]]))
    }

    p <- p + geom_segment(
      data = track,
      mapping = vline_map,
      alpha = alpha,
      linewidth = 0.4,
      inherit.aes = FALSE
    )
  }
  
  p +
    scale_x_continuous(
      breaks = chr_labels$chr_mid,
      labels = chr_labels$chr,
      expand = c(0.01, 0)
    ) +
    labs(x = NULL, y = y_label, title = title) +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.x        = element_text(size = text_size),
      axis.text.y = element_text(size = text_size)
    )
}
