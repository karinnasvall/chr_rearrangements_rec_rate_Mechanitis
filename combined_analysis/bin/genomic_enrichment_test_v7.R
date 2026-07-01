#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(optparse)
})

# ============================================================
# Debugging helpers
# ============================================================

DEBUG <- interactive()

msg <- function(...) {
  cat(
    sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")),
    sprintf(...),
    "\n"
  )
}

dbg <- function(...) {
  if (DEBUG)
    message(sprintf(...))
}

check_cols <- function(dt, cols, name) {
  miss <- setdiff(cols, names(dt))
  if (length(miss))
    stop(name, " missing columns: ",
         paste(miss, collapse = ", "))
}

# ============================================================
# Background preprocessing
# ============================================================


precompute_runs <- function(bg) {

  setorder(bg, genome, chr, start)

  ws <- median(bg$end - bg$start)

  bg[, gap :=
       c(TRUE,
         genome[-1] != genome[-.N] |
         chr[-1] != chr[-.N] |
         diff(start) != ws)]

  bg[, run_id := cumsum(gap)]
  bg[, gap := NULL]

  bg
}

precompute_starts <- function(bg, block_sizes) {

  out <- vector("list", length(block_sizes))
  names(out) <- as.character(block_sizes)

  run_starts <- split(seq_len(nrow(bg)), bg$run_id)

  for (k in block_sizes) {

    starts <- integer()

    for (idx in run_starts) {

      n <- length(idx)

      if (n >= k) {
        starts <- c(
          starts,
          idx[1] + seq_len(n - k + 1) - 1
        )
      }
    }

    out[[as.character(k)]] <- starts
  }

  out
}

sample_block <- function(bg, k, valid_starts) {

  starts <- valid_starts[[as.character(k)]]

  if (is.null(starts) || length(starts) == 0)
    return(NULL)

  s <- sample(starts, 1)

  bg[s:(s + k - 1)]
}

sample_values <- function(bg,
                          block_sizes,
                          valid_starts) {

  vals <- vector("list", length(block_sizes))

  for (i in seq_along(block_sizes)) {

    k <- block_sizes[i]

    b <- sample_block(
      bg,
      k,
      valid_starts
    )

    if (!is.null(b))
      vals[[i]] <- b$value
  }

  unlist(vals, use.names = FALSE)
}

# ============================================================
# Permutation test
# ============================================================

permutation_test <- function(
    roi,
    bg,
    valid_starts,
    n_perm = 1000,
    seed = 42) {

  set.seed(seed)

  roi[, breakpoint_id :=
        paste(genome,
              chr,
              roi_start,
              roi_end,
              sep = ":")]

  block_sizes <- roi[
    ,
    .N,
    by = breakpoint_id
  ]$N

  stopifnot(sum(block_sizes) == nrow(roi))

  obs_med <- median(roi$value,
                    na.rm = TRUE)

  stat_fun <- if (obs_med == 0)
    mean else median

  stat_name <- if (obs_med == 0)
    "mean" else "median"

  obs <- stat_fun(
    roi$value,
    na.rm = TRUE
  )

  dbg("Running %d permutations", n_perm)

  null <- vapply(
    seq_len(n_perm),
    function(i) {

      vals <- sample_values(
        bg,
        block_sizes,
        valid_starts
      )

      if (length(vals) == 0)
        return(NA_real_)

      stat_fun(vals,
               na.rm = TRUE)

    },
    numeric(1)
  )

  null_center <- stat_fun(
    null,
    na.rm = TRUE
  )

  p <- (sum(abs(null - null_center) >= abs(obs - null_center), 
  na.rm = TRUE) + 1) / 
     (sum(is.finite(null)) + 1)
  

  list(
    observed = obs,
    null_dist = null,
    null_stat = null_center,
    background_stat =
      stat_fun(bg$value,
               na.rm = TRUE),
    p_two_tailed = p,
    n_roi_windows = nrow(roi),
    n_background_windows = nrow(bg),
    stat_function = stat_name
  )
}

# ============================================================
# Plotting
# ============================================================

plot_null <- function(store, out) {

  pd <- rbindlist(
    lapply(names(store), function(n) {

      data.table(
        panel = n,
        val = store[[n]]$null_dist,
        obs = store[[n]]$observed
      )
    })
  )

  # 1. Subset and clean on the fly using data.table syntax
  plot_data <- pd[!is.na(val) & !is.na(obs)]

  # 2. Force variables to drop unused factor levels 
  # (Crucial if var1 or var2 were factors before filtering)
  plot_data[, panel := factor(panel)]


  # 3. Safe check: Only plot if rows actually exist
  if (nrow(plot_data) > 0) {
     

  p <- ggplot(plot_data,
              aes(val)) +
    geom_histogram(
      bins = 50
    ) +
    geom_vline(
      aes(xintercept = obs),
      colour = "blue"
    ) +
    facet_wrap(
      ~panel,
      scales = "free"
    ) +
    theme_bw()

  ggsave(
    out,
    p,
    width = 12,
    height = 8
  )

    } else {
    message("Skipping plot: No valid combinations available.")
  }
}

# ============================================================
# Main
# ============================================================

main <- function(opt) {

  msg("Loading data")

  combined <- if (is.character(opt$features)) {
    fread(opt$features)
  } else {
    as.data.table(opt$features)
  }

  bg <- if (is.character(opt$bg_bed)) {
    fread(opt$bg_bed)
  } else {
    as.data.table(opt$bg_bed)
  }

  setDT(combined)
  setDT(bg)

  check_cols(
    combined,
    c(
      "chr",
      "genome",
      "start",
      "end",
      "feature",
      "value",
      "window_id",
      "roi_start",
      "roi_end",
      "event_timing",
      "type"
    ),
    "combined"
  )

  check_cols(
    bg,
    c(
      "chr",
      "genome",
      "start",
      "end"
    ),
    "bg"
  )

  # Define sampling universe by overlap with bg
  setkey(combined, genome, chr, start, end)
  setkey(bg, genome, chr, start, end)

  bg <- foverlaps(
      combined,
      bg,
      by.x = c("genome", "chr", "start", "end"),
      by.y = c("genome", "chr", "start", "end"),
      nomatch = 0L
  )

  bg[, `:=`(
      start = i.start,
      end   = i.end
  )]

  bg[, c("i.start", "i.end") := NULL]


  results <- list()
  store <- list()

  print(unique(combined$feature))

  for (feat in unique(combined$feature)) {

    msg("Feature: %s", feat)


    feat_df <- combined[
      feature == feat
    ]

    for (ev in unique(
      feat_df$event_timing)) {

      for (tp in unique(
        feat_df$type)) {

        roi <- feat_df[
          event_timing == ev &
            type == tp
        ]

        print(summary(roi))

        if (!nrow(roi))
          next

        # process background for this feature/timing/type combination
        # filter roi from bg
        # Replace your current bg2 assignment inside the nested loop with this:
        bg2 <- bg[feature == feat]

        # Create a unique text fingerprint of your ROI windows to ensure drop accuracy
        roi_keys <- roi[, paste(genome, chr, start, end, sep = "_")]
        bg2_keys  <- bg2[, paste(genome, chr, start, end, sep = "_")]

        # Explicitly exclude them
        bg2 <- bg2[!bg2_keys %in% roi_keys]


        print(head(bg2))

        msg("Precomputing runs")
        bg2 <- precompute_runs(bg2)

        max_block_sizes <- roi[
          ,
          .N,
          by = .(
            genome,
            chr,
            roi_start,
            roi_end
          )
        ]$N

        valid_starts <- precompute_starts(
          bg2,
          unique(max_block_sizes)
        )


        if (!nrow(bg2))
          next

        
        res <- permutation_test(
          roi = roi,
          bg = bg2,
          valid_starts =
            valid_starts,
          n_perm =
            opt$n_perm,
          seed =
            opt$seed
        )

        results[
           length(results)+1] <- list(
          data.table(
            feature = feat,
            event_timing = ev,
            type = tp,
            observed_stat =
              res$observed,
            background_stat =
              res$background_stat,
            null_stat =
              res$null_stat,
            p_two_tailed =
              res$p_two_tailed,
            n_roi_windows =
              res$n_roi_windows,
            n_background_windows =
              res$n_background_windows,
            stat_function =
              res$stat_function
          )
        )

        store[paste(
            feat,
            ev,
            tp,
            sep = "|"
          )] <- list(res)
      }
    }
  }

  head(store[[1]]$null_dist)

  out <- rbindlist(
    results,
    fill = TRUE
  )

  fwrite(
    out,
    opt$output,
    sep = "\t"
  )

  print(head(store[[1]]))

 

  if (!is.null(opt$plot))
    plot_null(
      store,
      opt$plot
    )

  msg("Done")
}

# ============================================================
# Interactive vs CLI
# ============================================================

if (!interactive()) {

    opt <- parse_args(
        OptionParser(
            option_list = list(
                make_option("--features", type = "character"),
                make_option("--bg-bed", type = "character", dest = "bg_bed"),
                make_option("--output", type = "character"),
                make_option("--plot", type = "character", default = NULL),
                make_option("--n-perm", type = "integer",
                            default = 1000, dest = "n_perm"),
                make_option("--seed", type = "integer", default = 42)
            )
        )
    )

    main(opt)
}


