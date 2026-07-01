library(dplyr)
library(tidyverse)

draw.linear <- function(output, sizefile, ..., directory = NULL, fileformat = "png", anc_colours = NULL, chr_colours = NULL, w = 13, h = 5, opacity = .5, show_chromosome_names = TRUE, threshold=0.80, colour_by_genome = NULL) {
  
  if (is.null(directory)) {
    directory <- tempdir()
  }
  
  synteny.data.reframing <- function(data, tar.y, ref.y, compiled.size) {
    synteny <- data.frame()
    for (i in c(1:nrow(data))) {
      reference <- data[i, "ref.species"]
      target <- data[i, "tar.species"]
      tar_chr <- data[i, "tarchr"]
      ref_chr <- data[i, "refchr"]
      dir <- data[i, "dir"]
      fill_group <- data[i, "fill_group"]
      gene_id <- data[i, "gene_id"]
      tar_sizes <- compiled.size[compiled.size$species == target, ]
      names(tar_sizes) <- c("tarchr", "size", "species", "heterozygous", "xstart", "xend")
      ref_sizes <- compiled.size[compiled.size$species == reference, ]
      names(ref_sizes) <- c("refchr", "size", "species", "heterozygous", "xstart", "xend")
      tar_add <- tar_sizes[as.character(tar_sizes$tarchr) == as.character(tar_chr), ]$xstart
      ref_add <- ref_sizes[as.character(ref_sizes$refchr) == as.character(ref_chr), ]$xstart      
      tar_y <- tar.y
      ref_y <- ref.y
      tar_xstart <- data[i, "tarstart"] + tar_add
      tar_xend <- data[i, "tarend"] + tar_add
      ref_xstart <- data[i, "refstart"] + ref_add
      ref_xend <- data[i, "refend"] + ref_add
      anc_element <- data[i, "anc_element"]
      
      inverted <- grepl("-", dir, fixed = TRUE)
      if (inverted == TRUE) {
        df <- data.frame(
          x = c(tar_xstart, tar_xend, ref_xstart, ref_xend), y = c(tar_y, tar_y, ref_y, ref_y),
          fill = ref_chr, tar_chr = tar_chr, group = paste0("s", i), ref = reference, tar = target, dir = dir, fill_group = fill_group, anc_element = anc_element, gene_id = gene_id
        )
      } else {
        df <- data.frame(
          x = c(tar_xstart, ref_xstart, ref_xend, tar_xend), y = c(tar_y, ref_y, ref_y, tar_y),
          fill = ref_chr, tar_chr = tar_chr, group = paste0("s", i), ref = reference, tar = target, dir = dir, fill_group = fill_group, anc_element = anc_element, gene_id = gene_id
        )
      }
      synteny <- rbind(synteny, df)
    }
    return(synteny)
  }
  
  bezier <- function(p0, p1, p2, p3, n = n_points) {
    t <- seq(0, 1, length.out = n)
    data.frame(
      x = (1 - t)^3*p0[1] + 3*(1 - t)^2*t*p1[1] + 3*(1 - t)*t^2*p2[1] + t^3*p3[1],
      y = (1 - t)^3*p0[2] + 3*(1 - t)^2*t*p1[2] + 3*(1 - t)*t^2*p2[2] + t^3*p3[2]
    )
  }
  
  sum_nonoverlap <- function(starts, ends) {
    intervals <- data.frame(start = starts, end = ends)
    intervals <- intervals[order(intervals$start), ]
    
    total <- 0
    current_start <- -Inf
    current_end <- -Inf
    
    for (i in 1:nrow(intervals)) {
      s <- intervals$start[i]
      e <- intervals$end[i]
      
      if (s > current_end) {
        total <- total + (e - s + 1)
        current_start <- s
        current_end <- e
      } else if (e > current_end) {
        total <- total + (e - current_end)
        current_end <- e
      }
    }
    return(total)
  }
  
  classify_reciprocal_conservation <- function(data, threshold = 0.80) {
    
    # reverse inverted segments
    df <- data %>%
      mutate(
        refstart_norm = pmin(refstart, refend),
        refend_norm   = pmax(refstart, refend),
        tarstart_norm = pmin(tarstart, tarend),
        tarend_norm   = pmax(tarstart, tarend)
      )
    
    # Reference chromosomes
    ref_totals <- df %>%
      group_by(refchr) %>%
      summarise(
        ref_total = sum_nonoverlap(refstart_norm, refend_norm),
        .groups = "drop"
      )
    
    # Target chromosomes
    tar_totals <- df %>%
      group_by(tarchr) %>%
      summarise(
        tar_total = sum_nonoverlap(tarstart_norm, tarend_norm),
        .groups = "drop"
      )
    
    pair_totals <- df %>%
          group_by(refchr, tarchr) %>%
          summarise(
            pair_ref_aln = sum_nonoverlap(refstart_norm, refend_norm),
            pair_tar_aln = sum_nonoverlap(tarstart_norm, tarend_norm),
            .groups = "drop"
          ) %>%
      left_join(ref_totals, by = "refchr") %>%
      left_join(tar_totals, by = "tarchr") %>%
      mutate(
        percent_ref = pair_ref_aln / ref_total,
        percent_tar = pair_tar_aln / tar_total)
        
    conserved.df <- pair_totals %>% filter(
           percent_ref >= threshold | percent_tar >= 0.20
      ) %>%
      add_count(refchr, name = "ref_count") %>%   # adds a column with how many times refchr appears
      add_count(tarchr, name = "tar_count") %>%  # adds a column with how many times tarchr appears
      filter(ref_count == 1 & tar_count == 1) %>%
      filter(percent_ref >= threshold & percent_tar >= threshold)
    
    conserved.df$fill_group <- "Conserved"

    data <- left_join(data, conserved.df %>% select(refchr, fill_group)) %>%
      mutate(fill_group = if_else(is.na(fill_group), "Rearranged", fill_group))
    
    data
  }
  
  
  colours.default <- c(
    "01" = "#BFD73B", "02" = "#39ACE2", "03" = "#F16E8A",
    "04" = "#2DB995", "05" = "#CC99CC", "06" = "#A085BD",
    "07" = "#2EB560", "08" = "#D79128", "09" = "#FDBB63",
    "10" = "#AFDFE5", "11" = "#BF1E2D", "12" = "purple4",
    "13" = "#B59F31", "14" = "#F68B1F", "15" = "#EF374B",
    "16" = "#D376FF", "17" = "#9590FF", "18" = "#CE4699",
    "19" = "#7C9ACD", "20" = "#84C441", "21" = "#404F23",
    "22" = "#607F4B", "23" = "#EBB4A9", "24" = "#F6EB83",
    "25" = "#915F6D", "26" = "#602F92", "27" = "#81CEC6",
    "28" = "#F8DA04", "29" = "peachpuff2", "30" = "gray85", "33" = "peachpuff3",
    "W" = "brown","W1" = "brown","W2" = "#855823","W3" = "#855823","W4" = "#855823","W5" = "#855823", 
    "Z" = "black", "Z1" = "black", "Z2" = "grey70", "Z3" = "grey70", "Z4" = "grey70", "Z5" = "grey70", 
    "Y" = "#9590FF", "X" = "#666666",
    "LGE22" = "grey", "LGE64" = "gray64",
    "1A" = "pink", "1B" = "dark blue", "4A" = "light green",
    "Gap" = "white", "LG2" = "black", "LG5" = "#CC99CC"
  )
  
  xstart <- xend <- refchr <- tarchr <- x <- y <- group <- fill <- chromosome <- species <- heterozygous <- NULL
  sizes <- utils::read.delim(sizefile, header = FALSE, colClasses = "character") # to be consistent with naming in EH
  names(sizes) <- c("chromosome", "size", "species", "heterozygous")
  sizes$size <- as.numeric(gsub(",", "", sizes$size))
  
  count <- 0
  compiled.size <- data.frame()
  for (i in unique(sizes$species)) {
    size.intermediate <- sizes[sizes$species == i, ]
    for (x in c(1:nrow(size.intermediate))) {
      if (x == 1) {
        total_start <- 1
        total_end <- size.intermediate[x, "size"]
      } else {
        total_start <- total_end + 6000000
        total_end <- total_start + size.intermediate[x, "size"]
      }
      size.intermediate[x, "xstart"] <- total_start
      size.intermediate[x, "xend"] <- total_end
    }
    compiled.size <- rbind(compiled.size, size.intermediate)
  }
  
  for (z in unique(compiled.size$species)) {
    compiled.size$y[compiled.size$species == z] <- count
    count <- count + 2
  }

  list.of.files <- list()
  for (i in list(...)) {
    list.of.files[[i]] <- i
  }
  
  listsynt <- list()
  for (i in 1:length(list.of.files)) {
    num <- i
    file <- list.of.files[[num]]
    dataTMP <- utils::read.delim(file, header = FALSE, colClasses = "character")
    data2 <- dataTMP[, c(4, 5, 6, 1, 2, 3, 7, 8, 9, 12, 13)]
    colnames(data2) <- c("tarchr", "tarstart", "tarend", "refchr", "refstart", "refend", "dir", "ref.species", "tar.species", "anc_element", "gene_id")
    data2$tarstart <- as.numeric(gsub(",", "", data2$tarstart))
    data2$tarend <- as.numeric(gsub(",", "", data2$tarend))
    data2$refstart <- as.numeric(gsub(",", "", data2$refstart))
    data2$refend <- as.numeric(gsub(",", "", data2$refend))
    reference <- data2[1, "ref.species"]
    target <- data2[1, "tar.species"]
    ref_y <- compiled.size[compiled.size$species == reference, "y"]
    tar_y <- compiled.size[compiled.size$species == target, "y"]
    if (tar_y[1] > ref_y[1]){
      ref_y <- ref_y[1] + 0.1
      tar_y <- tar_y[1]
    } else{
      ref_y <- ref_y[1]
      tar_y <- tar_y[1] + 0.1
    }
    
    data2 <- classify_reciprocal_conservation(as.data.frame(data2), threshold)
    data2 <- as.data.frame(data2)
    #write.csv(data2, file = paste("data", i, ".csv", sep = ""))
    x <- synteny.data.reframing(data2, tar_y, ref_y, compiled.size)
    x$fill <- as.factor(x$fill)
    x$fill_group <- as.factor(x$fill_group)
    listsynt[[i]] <- x
  }
  
  compiled.size$chromosome <- as.factor(compiled.size$chromosome)
  compiled.size <- left_join(compiled.size, sizes[, c("species","chromosome", "heterozygous")])

  # Build geneid -> chr lookup from the colour reference genome
  # Default: last species in size file (= target in first chain file)
  colour_ref_species <- if (is.null(colour_by_genome)) tail(unique(sizes$species), 1) else colour_by_genome
  geneid_chr_map <- NULL
  for (.i in seq_along(listsynt)) {
    .d <- listsynt[[.i]]
    if (.d$ref[1] == colour_ref_species) {
      geneid_chr_map <- distinct(.d, gene_id, chr_colour = as.character(fill))
      break
    } else if (.d$tar[1] == colour_ref_species) {
      geneid_chr_map <- distinct(.d, gene_id, chr_colour = as.character(tar_chr))
      break
    }
  }
  if (is.null(geneid_chr_map)) stop(paste0("colour_by_genome species '", colour_ref_species, "' not found in any chain file."))

  p <- ggplot2::ggplot()
  
  for (i in 1:length(listsynt)) {
    data <- listsynt[[i]]
    reference <- data[1, "ref"]
    target <- data[1, "tar"]
    ref_sizes <- compiled.size[compiled.size$species == reference, ]
    tar_sizes <- compiled.size[compiled.size$species == target, ]
    #ref_sizes <- left_join(ref_sizes, unique(data2[,c("refchr", "fill_group")]), by = c("chromosome"="refchr"))
    ref_sizes <- ref_sizes %>% mutate(
      fill_chr = case_when(
        str_detect(chromosome, "^Z") ~ "Z chromosome",
        str_detect(chromosome, "^W") ~ "W chromosome",
        heterozygous=="Y" ~ "Heterozygous",
        TRUE ~ "Homozygous"
      )
    )
    ref_sizes$fill_chr <- as.factor(ref_sizes$fill_chr)
    
    #tar_sizes <- left_join(tar_sizes, unique(data2[,c("tarchr", "fill_group")]), by = c("chromosome"="tarchr"))
    tar_sizes <- tar_sizes %>% mutate(
      fill_chr = case_when(
        str_detect(chromosome, "^Z") ~ "Z chromosome",
        str_detect(chromosome, "^W") ~ "W chromosome",
        heterozygous=="Y" ~ "Heterozygous",
        TRUE ~ "Homozygous"
      )
    )
    tar_sizes$fill_chr <- as.factor(tar_sizes$fill_chr)

    # Convert rect bounds to polygon corners for geom_shape (which needs x/y, not xmin/xmax)
    rect_to_shape <- function(df, chr_height = 0.10, chr_width_extra = 0.1) {
      do.call(rbind, lapply(seq_len(nrow(df)), function(k) {
        r <- df[k, ]
        ext <- (r$xend - r$xstart) * chr_width_extra
        data.frame(
          x         = c(r$xstart - ext, r$xend + ext, r$xend + ext, r$xstart - ext),
          y         = c(r$y,      r$y,    r$y + chr_height, r$y + chr_height),
          fill_chr  = r$fill_chr,
          species   = r$species,
          chromosome = r$chromosome,
          xstart    = r$xstart,
          xend      = r$xend,
          ybase     = r$y,
          group_id  = k
        )
      }))
    }

    ref_poly <- rect_to_shape(ref_sizes)
    tar_poly <- rect_to_shape(tar_sizes)

    p <- p +
      ggforce::geom_shape(
        data = ref_poly,
        mapping = ggplot2::aes(x = x, y = y, group = group_id, fill = fill_chr),
        color = "grey20", linewidth = 0.4, radius = unit(2, "pt")
      ) +
      ggplot2::geom_text(data = ref_sizes, mapping = ggplot2::aes(x = 2, y = y + 0.05, label = species), size = 3, hjust = 1) +
      ggforce::geom_shape(
        data = tar_poly,
        mapping = ggplot2::aes(x = x, y = y, group = group_id, fill = fill_chr),
        color = "grey20", linewidth = 0.4, radius = unit(2, "pt")
      ) +
      ggplot2::geom_text(data = tar_sizes, mapping = ggplot2::aes(x = 2, y = y + 0.05, label = species), size = 3, hjust = 1) +
      scale_fill_manual(name = "Chromosomes",
                        limits = c("Homozygous", "Heterozygous", "Z chromosome", "W chromosome"),
                        values = c( "grey95", "grey30","goldenrod1","goldenrod3")) +
      ggnewscale::new_scale_fill()
    
    # Function to generate Bezier ribbon polygons from synteny data
    make_synteny_ribbons <- function(data, n_points = 50, h = 0.01) {
      ribbons <- data %>%
        group_by(group, fill, fill_group, chr_colour) %>%
        group_modify(~{
          d <- .
          
          # Separate top and bottom x/y
          x_top <- d$x[d$y == max(d$y)]
          x_bottom <- d$x[d$y == min(d$y)]
          y_top <- max(d$y)
          y_bottom <- min(d$y)
          
          # Control y for Bezier curve
          ctrl_y <- (y_top + y_bottom)/2
          
          if (unique(d$dir) == "+"){
            # Top curve: left to right
            left_curve <- bezier(
              p0 = c(min(x_top), y_top),
              p1 = c(min(x_top), ctrl_y),
              p2 = c(min(x_bottom), ctrl_y),
              p3 = c(min(x_bottom), y_bottom + h),
              n = n_points
            )
            
            # Bottom curve: right to left
            right_curve <- bezier(
              p0 = c(max(x_top), y_top),
              p1 = c(max(x_top), ctrl_y),
              p2 = c(max(x_bottom), ctrl_y),
              p3 = c(max(x_bottom), y_bottom + h),
              n = n_points
            )
          } else {
            left_curve <- bezier(
              p0 = c(min(x_top), y_top),
              p1 = c(min(x_top), ctrl_y),
              p2 = c(max(x_bottom), ctrl_y),
              p3 = c(max(x_bottom), y_bottom + h),
              n = n_points
            )
            
            # Bottom curve: right to left
            right_curve <- bezier(
              p0 = c(max(x_top), y_top),
              p1 = c(max(x_top), ctrl_y),
              p2 = c(min(x_bottom), ctrl_y),
              p3 = c(min(x_bottom), y_bottom + h),
              n = n_points
            )
            
          }
          
          poly <- rbind(left_curve, right_curve[nrow(right_curve):1, ])
          poly
        }) %>%
        ungroup()
      
      return(ribbons)
    }
    
    if (!is.null(anc_colours)) {
      # User-supplied anc_element -> colour: colour ribbons directly by anc_element
      data$chr_colour <- as.character(data$anc_element)
      link_palette <- anc_colours
      fill_legend <- "Ancestral linkage group"
    } else {
      # Default: colour by chromosome of the colour reference genome via gene ID
      data <- left_join(data, geneid_chr_map, by = "gene_id")
      data$chr_colour[is.na(data$chr_colour)] <- "unknown"
      link_palette <- if (!is.null(chr_colours)) chr_colours else colours.default
      fill_legend <- paste0("Chromosome (", colour_ref_species, ")")
    }

    ribbon_poly <- make_synteny_ribbons(data)

    p <-  p + 
      geom_polygon(data = ribbon_poly, 
                   aes(x = x, y = y, group = group, fill = chr_colour), colour = NA) + 
      scale_fill_manual(name = fill_legend,
                      values = link_palette) +
      ggnewscale::new_scale_fill()
    
    if (show_chromosome_names) {
      p <- p + 
        ggplot2::geom_text(data = ref_sizes, ggplot2::aes(x = (xstart + xend) / 2, y = y + 0.05, label = chromosome), size = 3, angle = 0, fontface = "bold") +
        ggplot2::geom_text(data = tar_sizes, ggplot2::aes(x = (xstart + xend) / 2, y = y + 0.05, label = chromosome), size = 3, angle = 0, fontface = "bold")
    }
  }
  
  p <- p +
    #ggplot2::scale_fill_manual(values = colours) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      legend.position = "right",
      legend.title = element_text()
    )
  
  message(paste0("Saving linear image to ", directory))

  ggplot2::ggsave(paste0(directory,"/",output, ".", fileformat), p, device = fileformat, width = w, height = h)
}
