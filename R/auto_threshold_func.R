#!/usr/bin/env R

########### ~~~~ Automated threshold ~~~~ ###########  

# Load required packages 
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(LaplacesDemon))


### Computes kernel density and trims ends (+ and - %1)
# v: Numeric vector of values to estimate density from
# ...: Additional arguments passed to stats::density()
# Returns a list containing trimmed density info 
density_with_trim <- function(v, ...) {
  dens <- stats::density(v, ...)
  x <- dens$x
  y <- dens$y
  
  if (length(x) < 3) {
    return(list(
      dens = dens,
      x_trim = x,
      y_trim = y,
      keep_mask = rep(TRUE, length(x))
    ))
  }
  
  min_x <- min(x)
  max_x <- max(x)
  margin <- 0.01 * (max_x - min_x)
  
  keep_mask <- x > (min_x + margin) & x < (max_x - margin)
  
  list(
    dens = dens,
    x_trim = x[keep_mask],
    y_trim = y[keep_mask],
    keep_mask = keep_mask,
    min_x = min_x,
    max_x = max_x,
    margin = margin,
    x = x,
    y = y
  )
}

### Computes automated marker expression thresholds across control samples
# df: A dataframe containing marker expression columns and metadata
# markers: Character vector of marker column names to threshold
# control_list: Character vector of control sample identifiers
# batch_list: Character vector of batch identifiers
# batch_colm: Column name in df containing batch identifiers
# subsample: Can be either "yes" or "no". Whether to subsample per batch prior to thresholding
# seed: Integer seed for reproducibility of random sampling
# num: Number of cells to subsample per batch (only done if subsample = "yes")
# output_dir: Directory where the CSV will be saved
# Returns a dataframe with one row per marker and the median cutoff value
automated_threshold_DF <- function(df,
                                    markers,
                                    control_list,
                                    batch_list,
                                    batch_colm,
                                    subsample = "yes",
                                    seed = 450,
                                    num = 20000,
                                    output_dir) {
  
  # Allowed option values
  valid_options <- c("yes", "no")

  subsample    <- tolower(subsample)

  # Input validation
  if (!subsample %in% valid_options) {
    stop("Invalid input: please specify 'yes' or 'no'")
  }
  
  # Ensure all controls exist
  if (!all(control_list %in% unique(df$control))) {
    stop("Check if dataframe has all unique control samples")
  }
  
  auto_cutoff <- data.frame()
  
  # Select thresholding strategy by dataset type
  for (i in seq_along(control_list)) {
    
    tmp <- df[df$control == control_list[i], ]
    
    # Optional per-batch subsampling
    if (subsample == "yes") {
      tmp <- data.frame(
        subsample_batch(
          df = tmp,
          batch_list = batch_list,
          batch_colm = batch_colm,
          num = num,
          seed = seed
        )
      )
    }
    
    for (m in seq_along(markers)) {
      
      cutoff_val <- auto_threshold_func(
        tmp[[markers[m]]],
        marker_name = markers[m]
      )
      
      auto_cutoff <- rbind(
        auto_cutoff,
        data.frame(
          marker  = markers[m],
          cutoff  = cutoff_val,
          control = as.character(control_list[i]),
          stringsAsFactors = FALSE
        )
      )
    }
  }
  
  # Aggregate median cutoff per marker
  auto_cutoff <- auto_cutoff %>%
    dplyr::group_by(marker) %>%
    dplyr::summarise(cutoff = median(cutoff), .groups = "drop")
  
  # Round and order by input marker list
  auto_cutoff$cutoff <- round(auto_cutoff$cutoff, 3)
  auto_cutoff <- auto_cutoff[match(markers, auto_cutoff$marker), ]
  
  # Output CSV
  file_path <- file.path(output_dir, "automated_marker_thresholds.csv")
  readr::write_csv(auto_cutoff, file = file_path)
  
  auto_cutoff
}

### Determines a cutoff for unimodal marker distributions using density shoulders
# vec: Numeric vector of marker expression values
# marker_name: Character string identifying the marker
# shoulder_frac: Fraction of peak density height used to define the shoulder
# Returns a numeric cutoff rounded to three decimal places
unimodal_threshold <- function(vec,
                               marker_name,
                               shoulder_frac = 0.05) {
  
  # Kernel density estimation with trimmed tails
  dpack  <- density_with_trim(vec)
  x_trim <- dpack$x_trim
  y_trim <- dpack$y_trim
  
  # Fallback if density is too sparse (should not happen)
  if (length(x_trim) < 5) {
    return(round(stats::quantile(vec, 0.995, na.rm = TRUE), 3))
  }
  
  # Identify main density peak
  peak_idx <- which.max(y_trim)
  peak_y   <- y_trim[peak_idx]
  peak_x   <- x_trim[peak_idx]
  
  # Approximate integration step
  dx <- diff(x_trim)[1]
  
  # Mass on either side of the peak
  mass_left  <- sum(y_trim[1:peak_idx]) * dx
  mass_right <- sum(y_trim[peak_idx:length(y_trim)]) * dx
  mass_ratio <- mass_right / (mass_left + 1e-9)
  
  # Quantile-based dynamic range
  q_low  <- stats::quantile(vec, 0.02, na.rm = TRUE)
  q_high <- stats::quantile(vec, 0.98, na.rm = TRUE)
  dyn    <- q_high - q_low
  
  # Fallback if dynamic range collapses (should not happen)
  if (!is.finite(dyn) || dyn <= 0) {
    return(round(stats::quantile(vec, 0.995, na.rm = TRUE), 3))
  }
  
  # Normalize peak position on 0-1 scale
  peak_position <- (peak_x - q_low) / dyn
  peak_position <- max(0, min(1, peak_position))
  
  # Define shoulder height
  threshold_height <- peak_y * shoulder_frac
  
  # Initial unimodal classification
  # Determine if unimodal positive or unimodal negative
  if (peak_position < 0.4) {
    mode_type <- "negative"
  } else if (peak_position > 0.8) {
    mode_type <- "positive"
  } else {
    if (mass_ratio < 0.8) {
      mode_type <- "negative"
    } else if (mass_ratio > 1.3) {
      mode_type <- "positive"
    } else {
      if (peak_position < 0.6) {
        mode_type <- "negative"
      } else {
        mode_type <- "positive"
      }
    }
  }
  
  # Override if the density is centered and wide 
  # In this case, unimodal distributions are likely positive
  if (mode_type == "negative" &&
      peak_position > 0.4 &&
      mass_ratio < 1.3 &&
      dyn > 1.2 &&
      q_low > 0) {
    
    mode_type <- "positive"
  }
  
  # Get cutoff for negative or positive unimodal
  if (mode_type == "negative") {
    
    # Search right shoulder
    search_idx <- seq(peak_idx, length(y_trim))
    below_idx  <- search_idx[y_trim[search_idx] <= threshold_height]
    
    if (length(below_idx) > 0) {
      cutoff <- x_trim[min(below_idx)]
    } else {
      # Fallback
      cutoff <- stats::quantile(vec, 0.995, na.rm = TRUE)
    }
    
  } else {
    
    # Search left shoulder
    search_idx <- seq_len(peak_idx)
    below_idx  <- search_idx[y_trim[search_idx] <= threshold_height]
    
    if (length(below_idx) > 0) {
      cutoff <- x_trim[max(below_idx)]
    } else {
      # Fallback
      cutoff <- stats::quantile(vec, 0.005, na.rm = TRUE)
    }
  }
  
  # Return cutoff, rounded
  round(cutoff, 3)
}

### Computes an automated cutoff for bimodal or multimodal marker distributions
# exp_vector: Numeric vector of marker expression values
# marker_name: Character string identifying the marker
# width_cutoff: Minimum allowable peak width (set to exclude spurious narrow peaks)
# sep_peaks: Minimum separation required between peaks
# dip_cutoff: Minimum dip strength (valley) required between peaks 
# min_peak_height_frac: Minimum relative peak height as a fraction of max peak
# Returns a numeric cutoff rounded to three decimal places
auto_threshold_func <- function(exp_vector, 
                                 marker_name = "UNKNOWN_MARKER", 
                                 width_cutoff = 0.4,
                                 sep_peaks = 0.5,
                                 dip_cutoff = 0.15,
                                 min_peak_height_frac = 0.008) {
  
  # Density on cleaned data
  dpack <- density_with_trim(exp_vector)
  x <- dpack$x
  y <- dpack$y
  x_trim <- dpack$x_trim
  y_trim <- dpack$y_trim
  
  if (length(y) < 5) {
    return(unimodal_threshold(exp_vector, marker_name))
  }
  
  # Identify local maxima
  peak_idx <- which(diff(sign(diff(y))) == -2) + 1
  
  # Estimate peak width by descending slopes
  estimate_peak_width <- function(x, y, peak_idx) {
    
    left <- peak_idx
    while (left > 2 && y[left - 1] < y[left]) {
      left <- left - 1
    }
    
    right <- peak_idx
    while (right < (length(y) - 1) && y[right + 1] < y[right]) {
      right <- right + 1
    }
    
    x[right] - x[left]
  }
  
  if (length(peak_idx) == 0) {
    return(unimodal_threshold(exp_vector, marker_name))
  }
  
  peak_positions_all <- x[peak_idx]
  peak_heights_all   <- y[peak_idx]
  
  # Remove peaks near density edges (+ or - 1%)
  min_x  <- dpack$min_x
  max_x  <- dpack$max_x
  margin <- dpack$margin
  
  keep_mask <- peak_positions_all > (min_x + margin) &
    peak_positions_all < (max_x - margin)
  
  peak_positions_all <- peak_positions_all[keep_mask]
  peak_heights_all   <- peak_heights_all[keep_mask]
  peak_idx           <- peak_idx[keep_mask]
  
  # Peak width filtering
  # Peaks must be wider then this cutoff
  peak_widths <- sapply(peak_idx, function(i) estimate_peak_width(x, y, i))
  
  keep_width <- peak_widths > width_cutoff
  
  peak_positions_all <- peak_positions_all[keep_width]
  peak_heights_all   <- peak_heights_all[keep_width]
  peak_idx           <- peak_idx[keep_width]
  
  # If there are less then 2 peaks left, revert to unimodal
  if (length(peak_positions_all) < 2) {
    return(unimodal_threshold(exp_vector, marker_name))
  }
  
  # Relative peak height filtering
  # Next tallest peak must be taller then this cutoff
  height_cut <- min_peak_height_frac * max(peak_heights_all)
  
  keep_rel <- peak_heights_all >= height_cut
  
  peak_positions <- peak_positions_all[keep_rel]
  peak_heights   <- peak_heights_all[keep_rel]
  
  # If there are less then 2 peaks left, revert to unimodal
  if (length(peak_positions) < 2) {
    return(unimodal_threshold(exp_vector, marker_name))
  }
  
  # Order peaks by descending height
  ord <- order(peak_heights, decreasing = TRUE)
  
  for (i in 2:length(ord)) {
    
    p1 <- peak_positions[ord[1]]
    p2 <- peak_positions[ord[i]]
    
    h1 <- peak_heights[ord[1]]
    h2 <- peak_heights[ord[i]]
    
    sep <- abs(p1 - p2)
    
    if (sep < sep_peaks) next
    
    lower <- min(p1, p2)
    upper <- max(p1, p2)
    
    window_mask <- x_trim > lower & x_trim < upper
    if (!any(window_mask)) next
    
    local_x <- x_trim[window_mask]
    local_y <- y_trim[window_mask]
    
    valley_idx <- which.min(local_y)
    valley_x   <- local_x[valley_idx]
    valley_y   <- local_y[valley_idx]
    
    peak_mean <- mean(c(h1, h2))
    dip_strength <- 1 - (valley_y / peak_mean)
    
    if (dip_strength < dip_cutoff) next
    
    return(round(valley_x, 3))
  }
  
  unimodal_threshold(exp_vector, marker_name)
}
