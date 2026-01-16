#!/usr/bin/env R

########### ~~~~ Automated threshold (flow) ~~~~ ###########  

# Load required packages 
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(LaplacesDemon))

### Automatically computes a threshold for a marker expression vector using unimodal/bimodal heuristics
# exp_vector: Numeric vector of expression values for a single marker
# Returns a single numeric threshold that divides the positive population from the negative
auto_threshold_func <- function(exp_vector) {
  
  # If distribution is unimodal
  if (is.unimodal(exp_vector, min.size = 0.05) == TRUE) {
    
    q4 <- quantile(exp_vector, 0.95) 
    names(q4) <- NULL
    
    q3 <- quantile(exp_vector, 0.75)
    names(q3) <- NULL
    
    qt_8 <- quantile(exp_vector, 0.85)
    names(qt_8) <- NULL
    
    if (q4 < 0.5 && q3 < 0.1) {
      
      val <- quantile(exp_vector, 0.98)
      names(val) <- NULL
      return(round(val, 3))
    }
    
    else if (q3 >= 0.1) {
      
      if (q4 > 0.5 && q4 < 1) {
        val <- round((q4 + q3) / 2, 3)
      }
      else if (q4 >= 1 && q4 < 1.5) {
        val <- round((q4 + q3) / 2, 3)
      }
      else {
        if (q4 >= 1.5 && qt_8 >= 1) {
          val <- round((qt_8 + q3) / 2, 3)
          names(val) <- NULL
        }
        else {
          val <- quantile(exp_vector, 0.85)
          names(val) <- NULL
        }
      }
      return(round(val, 3))
    }
    
    else {
      val <- quantile(exp_vector, 0.75) + 0.05
      names(val) <- NULL
      return(round(val, 3))
    }
  }
  
  # Otherwise, use kernel density to find peaks and valley between them
  else {
    
    dens <- stats::density(exp_vector) 
    
    med_x <- abs(median(dens$x)) 
    
    q3 <- abs(quantile(exp_vector, 0.75))
    names(q3) <- NULL
    
    ix01_y <- which.max(dens$y)
    peak_neg <- dens$x[ix01_y] 
    
    upper_mask <- dens$x > q3
    
    if (!any(upper_mask)) {
      
      val <- quantile(exp_vector, 0.75)
      names(val) <- NULL
      
      return(round(val, 3))
    }
    upper_val <- dens$y[upper_mask]
    
    y_idx_local <- which.max(upper_val)
    y_idx <- which(upper_mask)[y_idx_local]
    peak_pos <- dens$x[y_idx] 
    
    # Finding the valley between negative and positive peaks
    if (peak_pos != peak_neg) {

      window_mask <- dens$x < peak_pos & dens$x > peak_neg
      if (!any(window_mask)) {
        valley <- round(median(exp_vector), 3)
        names(valley) <- NULL
        return(round(valley, 3))
      }
      
      min_y <- min(dens$y[window_mask]) 

      y_trough <- which(dens$y == min_y)[1] 
      valley <- dens$x[y_trough] 
      
      return(round(valley, 3)) 
    }
    
    else { 
      val <- quantile(exp_vector, 0.75)
      names(val) <- NULL
      return(round(val, 3))
    }
  }
}

