
#!/usr/bin/env R

########### ~~~~ Required Libraries ~~~~ ###########

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggridges))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(LaplacesDemon))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(flowCore))
suppressPackageStartupMessages(library(FlowSOM))
suppressPackageStartupMessages(library(ConsensusClusterPlus))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(Rtsne))
suppressPackageStartupMessages(library(umap))
suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(progress))
suppressPackageStartupMessages(library(crayon))
suppressPackageStartupMessages(library(patchwork))

########### ~~~~ Plot Colors ~~~~ ###########


sample_colours <- c("#8ad6cc","#47857c","#f6b250","#e35988",
                    RColorBrewer::brewer.pal(name = "Paired",n = 12),
                    RColorBrewer::brewer.pal(name = "Dark2",n = 8))

batchColours <- c("#5C8E39","#0072B2","#FF69B4","#E1A100","#19E7E7","#5E4FA2",
                  "#35B779","#00bbf9","#FF700A","#787878","#C24141","#3c096c",
                  "#FF5349","#00FF00","#0000FF","#83d632","#482878","#FF00FF",
                  "#008080","#FFC033","#c3acf4","#8AD6CC","#006400","#D88A45",
                  "#4169E1","#FF69B4","#9932CC","#6495ED","#91B0B9","#AE3A3A",
                  "#3288BD","#B91F48","#486BAF","#88CFA4","#9E0142","#4CA5B1",
                  "#66C2A5","#FF9500","#8338ec","#ABDDA4","#003049","#669BBC",
                  "#C8E99E","#ffbe0b","#588157","#f4a261","#6A3D9A","#ff4d6d")

clustColors <- c("#ED553B","#20639B","#F6d55C","#3CAEA3","#FEA8ED",
                 "#A8FAA4","#FCFC97","#B6056A","#0D44BA","#880ED4",
                 "#DE8297","#B8E0FF","#FF1919","#CE8CF8","#69794E",
                 "#71A6B8","#B20245","#236A3A","#C576F6","#FCB254",
                 "#559050","#CA9BB0","#C78D40","#9EB994","#24BAB9",
                 "#EE4483","#F25A25","#93C47D","#C2272D","#38761D",
                 "#C898B1","#98C8C7","#C8AF98","#668BAD","#98B1C8",
                 "#1CD103","#4C03BD","#E6E50E","#C500CB","#00ADE3",
                 "#8ad6cc","#47857c","#f6b250","#e35988","#482878",
                 "#5C8E39","#0072B2","#FF69B4","#E1A100","#787878")


########### ~~~~ Data Transformation ~~~~ ########### 

### Creates a transformed dataframe (arcsinh or logicle) from FCS files
# transformation_class: Either "arcsinh" or "logicle" specifying the transformation to apply
# cf: Cofactor used for arcsinh transformation (default = 6000)
# fcs_info: A dataframe describing files (must include columns 'file', 'batch', 'control', and optionally 'sample_id')
# tag_list: A character vector of fluorphore/metal tag names present in the FCS files
# marker_list: A character vector of protein marker channel names present in the FCS files
# marker_tbl: A table mapping tags to markers (columns specified by tag_colm and marker_colm)
# tag_colm: Column name in marker_tbl containing tag names (e.g., PnN)
# marker_colm: Column name in marker_tbl containing marker/protein names (e.g., PnS)
# file_dir: Directory where the FCS files are located
# output_dir: Directory where the transformed CSV will be saved
create_transformed_dataframe_function <- function(transformation_class, 
                                                  cf = 6000, 
                                                  fcs_info, 
                                                  tag_list, 
                                                  marker_list,
                                                  marker_tbl, 
                                                  tag_colm, 
                                                  marker_colm, 
                                                  file_dir, 
                                                  output_dir) {
  
  # Use arcsinh or logicle transformer based on user choice
  if (transformation_class == "arcsinh") {
    
    ffdf <- create_arcsinh_transformed_dataframe(cf = cf, fcs_info = fcs_info, tag_list = tag_list, 
                                                 marker_list = marker_list, marker_tbl = marker_tbl, 
                                                 tag_colm = tag_colm, marker_colm = marker_colm, 
                                                 output_dir = output_dir, file_dir = file_dir)
    
  }
  
  if (transformation_class == "logicle") {
    
    ffdf <- create_logicle_transformed_dataframe(fcs_info = fcs_info, tag_list = tag_list, 
                                                 marker_list = marker_list, marker_tbl = marker_tbl, 
                                                 tag_colm = tag_colm, marker_colm = marker_colm, 
                                                 output_dir = output_dir, file_dir = file_dir)
    
  }
  
  # Return transformed dataframe
  return(ffdf)
}

### Performs arcsinh transformation of FCS files and saves the transformed dataframe to CSV
# cf: Cofactor for arcsinh transformation (e.g., 6000 for flow cytometry)
# fcs_info: A dataframe with files metadata (columns: file, batch, control, optionally sample_id)
# tag_list: Character vector of tag/fluorophore names (alternative to marker_list)
# marker_list: Character vector of protein marker channel names (alternative to tag_list)
# marker_tbl: Mapping table from tags to markers (columns specified by tag_colm and marker_colm)
# tag_colm: Column name in marker_tbl containing tags (e.g., PnN)
# marker_colm: Column name in marker_tbl containing marker names (e.g., PnS)
# file_dir: Directory path to the FCS files
# output_dir: Directory where the transformed CSV will be saved
create_arcsinh_transformed_dataframe <- function(cf, 
                                                 fcs_info, 
                                                 tag_list, 
                                                 marker_list,
                                                 marker_tbl, 
                                                 tag_colm, 
                                                 marker_colm,
                                                 file_dir, 
                                                 output_dir) {
  
  # Style for progress bar text
  colrSet <- crayon::make_style("#539DDD")
  
  # Create progress bar
  pBar <- progress::progress_bar$new(
    format = paste0(colrSet("Transforming FCS files"),
                    " [", colrSet(":bar"), "] ",
                    colrSet(":percent"), " | File: :current/:total ",
                    "| Elapsed: :elapsed | ETA: :eta"),
    total = nrow(fcs_info),
    clear = FALSE,
    width = 80,
    complete = "=",
    incomplete = "-"
  )
  
  # Normalize column names to lower case
  colnames(fcs_info) <- tolower(colnames(fcs_info))
  
  # Build full paths to FCS files
  fcs_info[, "path"] <- file.path(file_dir, fcs_info[, "file"])
  
  # Read the first file to inspect available channels
  tmp1 <- flowCore::read.FCS(fcs_info[1, "path"], transformation = FALSE, truncate_max_range = FALSE)
  flowcolnames <- flowCore::colnames(tmp1)
  
  # Detect missing tags/markers in the FCS file channels
  tag_match <- setdiff(tag_list, flowcolnames)
  mark_match <- setdiff(marker_list, flowcolnames)
  
  # Create sample_id column if one doesn't exist
  if (is.null(fcs_info[, "sample_id"])) {
    fcs_info[, "sample_id"] <- paste0(fcs_info[, "batch"], "_", fcs_info[, "control"])
  }
  
  # Initialize transformed expression dataframe
  expr_df <- data.frame()
  
  # Case 1: tags are present in the FCS files
  if (all(tag_list %in% flowcolnames) == TRUE) {
    for (i in 1:nrow(fcs_info)) {
      
      filename <- fcs_info[i, "file"]
      
      # Update progress bar
      pBar$tick(tokens = list(current = i))
      
      # Read flow set without internal transformations
      myfile <- flowCore::read.flowSet(fcs_info[i, "path"], 
                                       transformation = FALSE, 
                                       truncate_max_range = FALSE)
      # Apply arcsinh transformation on tag_list
      myfile <- flowCore::fsApply(myfile, function(x, cofactor = asinh_scale) {
        expr <- flowCore::exprs(x)
        expr <- asinh(expr[, tag_list] / cf) ### use list of fluorphores/metal tags 
        flowCore::exprs(x) <- expr
        x
      })
      
      # Extract expressions and attach filename
      tmp_mat <- flowCore::fsApply(myfile, flowCore::exprs)
      tmp_mat <- data.frame(tmp_mat)
      tmp_mat$file <- filename
      expr_df <- rbind(expr_df, tmp_mat)
    }
    
    # Name columns and map tags to marker names via marker_tbl
    colnames(expr_df) <- c(tag_list, "file")
    
    # Match metal/fluorophore names to panel file column
    match_names <- match(names(expr_df), marker_tbl[, tag_colm])
    
    # Rename to protein channel names
    names(expr_df)[!is.na(match_names)] <- marker_tbl[, marker_colm][na.omit(match_names)]
    expr_df <- data.frame(expr_df)
    
  }
  # Case 2: markers are present in the FCS files
  else if (all(marker_list %in% flowcolnames) == TRUE) {
    for (i in 1:nrow(fcs_info)) {
      
      filename <- fcs_info[i, "file"]
      
      # Update progress bar
      pBar$tick(tokens = list(current = i))
      
      # Read flow set, no internal transformations
      myfile <- flowCore::read.flowSet(fcs_info[i, "path"], 
                                       transformation = FALSE, 
                                       truncate_max_range = FALSE)
      
      # Apply arcsinh transformation on marker_list
      myfile <- flowCore::fsApply(myfile, function(x, cofactor = asinh_scale) {
        expr <- flowCore::exprs(x)
        expr <- asinh(expr[, marker_list] / cf) # Use list of protein channels 
        flowCore::exprs(x) <- expr
        x
      })
      
      # Extract expressions and attach filename
      tmp_mat <- flowCore::fsApply(myfile, flowCore::exprs)
      tmp_mat <- data.frame(tmp_mat)
      tmp_mat$file <- filename
      expr_df <- rbind(expr_df, tmp_mat)
    }
    
    # Name columns for markers + file
    colnames(expr_df) <- c(marker_list, "file")
  }
  # Case 3: neither tags nor markers fully match file channels
  else {
    stop("Names in fcs files do not match with the ones in the marker table")
    if (!is.null(tag_match)) {
      mm <- match(marker_tbl$PnN, tag_match)
      cat("marker =", marker_tbl[, "PnS"][na.omit(mm)], ", tag =", marker_tbl[, "PnN"][na.omit(mm)], "\n")
    }
    else if (!is.null(mark_match)) {
      mm <- match(marker_tbl$PnS, mark_match)
      cat("marker =", marker_tbl[, "PnS"][na.omit(mm)], ", tag =", marker_tbl[, "PnN"][na.omit(mm)], "\n")
    }
  }
  
  # Match file column to get batch, control, sample id information
  match_cols <- match(expr_df$file, fcs_info$file)
  expr_df$sample_id <- fcs_info$sample_id[match_cols]
  expr_df$batch <- fcs_info$batch[match_cols]
  expr_df$control <- fcs_info$control[match_cols]
  
  # Build output CSV path
  file_path <- file.path(output_dir, "transformed_dataframe_of_control_samples.csv")
  
  # Save transformed dataframe to CSV
  tryCatch({
    readr::write_csv(expr_df, file = file_path)
    cat("Transformed data saved as", file_path)
  }, error = function(e) {
    stop("Error saving table as CSV: ", e$message)
  })
  
  # Return the transformed dataframe
  return(expr_df)
}

### Builds a logicle transformation list for specified channels of a flowFrame
# ff: A flowCore::flowFrame object containing expression data
# channels: Character vector of channel names to transform
# w: Width of the linear region (default = 0.5)
# m: Total number of decades (default = 4.5)
# a: Additional negative decades (default = 0)
# Returns a flowCore::transformList to be applied on the given channels
build_logicle_transformer <- function(ff, channels,
                                      w = 0.5,  # width of linear region
                                      m = 4.5,  # total decades
                                      a = 0) {  # extra negative decades
  # Extract expressions from the flowFrame
  expr <- flowCore::exprs(ff)
  
  # Per-channel T = max positive value 
  t_vals <- apply(expr[, channels, drop = FALSE], 2,
                  function(z) max(z[is.finite(z)], na.rm = TRUE))
  
  # Build logicle transform per channel with channel-specific T
  tf_list <- lapply(seq_along(channels), function(i) {
    flowCore::logicleTransform(w = w, t = t_vals[i], m = m, a = a)
  })
  names(tf_list) <- channels
  
  # Return transform list for application via flowCore::transform
  flowCore::transformList(channels, tf_list)
}

### Performs logicle transformation of FCS files and saves the transformed dataframe to CSV
# cf: Cofactor parameter kept for interface compatibility (unused in logicle)
# fcs_info: A dataframe with files metadata (columns: file, batch, control, optionally sample_id)
# tag_list: Character vector of tag/fluorophore names (alternative to marker_list)
# marker_list: Character vector of protein marker channel names (alternative to tag_list)
# marker_tbl: Mapping table from tags to markers (columns specified by tag_colm and marker_colm)
# tag_colm: Column name in marker_tbl containing tags (e.g., PnN)
# marker_colm: Column name in marker_tbl containing marker names (e.g., PnS)
# file_dir: Directory path to the FCS files
# output_dir: Directory where the transformed CSV will be saved
create_logicle_transformed_dataframe <- function(cf,  # kept for interface compatibility (unused)
                                                 fcs_info,
                                                 tag_list,
                                                 marker_list,
                                                 marker_tbl,
                                                 tag_colm,
                                                 marker_colm,
                                                 file_dir,
                                                 output_dir) {
  
  # Style for progress bar text
  colrSet <- crayon::make_style("#539DDD")
  
  # Create progress bar
  pBar <- progress::progress_bar$new(
    format = paste0(colrSet("Transforming FCS files (logicle)"),
                    " [", colrSet(":bar"), "] ",
                    colrSet(":percent"), " | File: :current/:total ",
                    "| Elapsed: :elapsed | ETA: :eta"),
    total = nrow(fcs_info),
    clear = FALSE,
    width = 80,
    complete = "=",
    incomplete = "-"
  )
  
  # Normalize column names to lower case
  colnames(fcs_info) <- tolower(colnames(fcs_info))
  
  # Build full paths to FCS files
  fcs_info[, "path"] <- file.path(file_dir, fcs_info[, "file"])
  
  # Read the first file to inspect available channels
  tmp1 <- flowCore::read.FCS(fcs_info[1, "path"],
                             transformation = FALSE,
                             truncate_max_range = FALSE)
  
  flowcolnames <- flowCore::colnames(tmp1)
  
  # Detect missing tags/markers in the FCS file channels
  tag_match  <- setdiff(tag_list,    flowcolnames)
  mark_match <- setdiff(marker_list, flowcolnames)
  
  # Create sample_id column if one doesn't exist
  if (is.null(fcs_info[, "sample_id"])) {
    fcs_info[, "sample_id"] <- paste0(fcs_info[, "batch"], "_", fcs_info[, "control"])
  }
  
  # Initialize transformed expression dataframe
  expr_df <- data.frame()
  
  # Case 1: tag_list matches 
  if (all(tag_list %in% flowcolnames) == TRUE) {
    
    for (i in 1:nrow(fcs_info)) {
      filename <- fcs_info[i, "file"]
      
      # Update progress bar
      pBar$tick(tokens = list(current = i))
      
      # Read flow set without internal transformations
      myfile <- flowCore::read.flowSet(fcs_info[i, "path"],
                                       transformation = FALSE,
                                       truncate_max_range = FALSE)
      
      # Build logicle transform for this file on tag_list channels
      lgcl <- build_logicle_transformer(myfile[[1]], channels = tag_list,
                                        w = 0.5, m = 4.5, a = 0)
      
      # Apply logicle transform
      myfile <- flowCore::transform(myfile, lgcl)
      
      # Extract only transformed tag_list columns 
      tmp_mat <- flowCore::fsApply(
        myfile,
        function(x) {
          expr <- flowCore::exprs(x)
          expr[, tag_list, drop = FALSE]
        }
      )
      
      tmp_mat <- data.frame(tmp_mat)
      tmp_mat$file <- filename
      expr_df <- rbind(expr_df, tmp_mat)
    }
    
    # Name columns and map tags to marker names via marker_tbl
    colnames(expr_df) <- c(tag_list, "file")
    
    # Match metal/fluorophore names to panel file column
    match_names <- match(names(expr_df), marker_tbl[, tag_colm])
    
    # Rename to protein channel names
    names(expr_df)[!is.na(match_names)] <- marker_tbl[, marker_colm][na.omit(match_names)]
    expr_df <- data.frame(expr_df)
    
  }
  
  # Case 2: marker_list matches 
  else if (all(marker_list %in% flowcolnames) == TRUE) {
    
    for (i in 1:nrow(fcs_info)) {
      filename <- fcs_info[i, "file"]
      
      # Update progress bar
      pBar$tick(tokens = list(current = i))
      
      # Read flow set without internal transformations
      myfile <- flowCore::read.flowSet(fcs_info[i, "path"],
                                       transformation = FALSE,
                                       truncate_max_range = FALSE)
      
      # Build logicle transform for this file on marker_list channels
      lgcl <- build_logicle_transformer(myfile[[1]], channels = marker_list,
                                        w = 0.5, m = 4.5, a = 0)
      
      # Apply logicle transform
      myfile <- flowCore::transform(myfile, lgcl)
      
      # Extract only transformed marker_list columns
      tmp_mat <- flowCore::fsApply(
        myfile,
        function(x) {
          expr <- flowCore::exprs(x)
          expr[, marker_list, drop = FALSE]
        }
      )
      
      tmp_mat <- data.frame(tmp_mat)
      tmp_mat$file <- filename
      expr_df <- rbind(expr_df, tmp_mat)
    }
    
    # Name columns for markers + file
    colnames(expr_df) <- c(marker_list, "file")
  }
  
  # Case 3: mismatch
  else {
    stop("Names in fcs files do not match with the ones in marker table")
    
    if (!is.null(tag_match)) {
      mm <- match(marker_tbl$PnN, tag_match)
      cat("marker =", marker_tbl[, "PnS"][na.omit(mm)],
          ", tag =", marker_tbl[, "PnN"][na.omit(mm)], "\n")
    } else if (!is.null(mark_match)) {
      mm <- match(marker_tbl$PnS, mark_match)
      cat("marker =", marker_tbl[, "PnS"][na.omit(mm)],
          ", tag =", marker_tbl[, "PnN"][na.omit(mm)], "\n")
    }
  }
  
  # Match file column to get batch, control, sample id information
  match_cols <- match(expr_df$file, fcs_info$file)
  expr_df$sample_id <- fcs_info$sample_id[match_cols]
  expr_df$batch <- fcs_info$batch[match_cols]
  expr_df$control <- fcs_info$control[match_cols]
  
  # Build output CSV path
  file_path <- file.path(output_dir, "logicle_transformed_dataframe_of_control_samples.csv")
  
  # Save transformed dataframe to CSV
  tryCatch({
    readr::write_csv(expr_df, file = file_path)
    cat("Transformed data saved as", file_path)
  }, error = function(e) {
    stop("Error saving table as CSV: ", e$message)
  })
  
  # Return the transformed dataframe
  return(expr_df)
}

########### ~~~~ Marker Threshold Calculation ~~~~ ########### 

### Subsamples cells from a dataframe based on specified strategy
# df: A dataframe containing the transformed cytometry data
# sampling_type: Set as "overall_min", "batch_min", "fixed", or "none"
#   - "overall_min": Equalizes cell counts across all samples using the global minimum
#   - "batch_min": Equalizes cell counts within each batch using the batch-specific minimum
#   - "fixed": Samples a fixed number of cells per sample (requires n_cells)
#   - "none": Returns the dataframe unchanged (no subsampling)
# n_cells: Integer specifying number of cells per sample (only for "fixed")
# seed: Integer seed for reproducibility of random sampling
# Returns a dataframe with subsampled rows according to the chosen strategy
sample_cells <- function(df,
                         sampling_type = c("overall_min", "batch_min", "fixed", "none"),
                         n_cells = NULL,
                         seed = 250) {
  
  sampling_type <- match.arg(sampling_type)
  
  # Validate input dataframe and required columns
  if (!is.data.frame(df)) stop("`df` must be a data.frame.")
  if (!"sample_id" %in% names(df)) stop("df must contain a column named 'sample_id'.")
  if (anyNA(df$sample_id)) stop("df$sample_id contains NA; please fix before sampling.")
  
  # Additional checks for batch_min strategy
  if (sampling_type == "batch_min") {
    if (!"batch" %in% names(df)) stop("df must contain a column named 'batch' for sampling_type = 'batch_min'.")
    if (anyNA(df$batch)) stop("df$batch contains NA; please fix before sampling.")
  }
  
  # Set seed for reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # If no sampling requested, return original dataframe
  if (sampling_type == "none") { return(df) }
  
  # Split row indices by sample_id
  idx_by_sample <- split(seq_len(nrow(df)), df$sample_id)
  
  # Compute per-sample sizes
  sample_sizes <- vapply(idx_by_sample, length, integer(1))
  if (length(sample_sizes) == 0) return(df)
  
  # Determine target number of cells per sample based on sampling_type
  n_per_sample <- switch(
    sampling_type,
    "overall_min" = as.integer(min(sample_sizes)),  # global minimum across all samples
    "fixed" = {
      if (is.null(n_cells)) stop("`n_cells` must be provided when sampling_type = 'fixed'.")
      if (!is.numeric(n_cells) || length(n_cells) != 1 || is.na(n_cells) || n_cells <= 0)
        stop("`n_cells` must be a single positive integer.")
      as.integer(n_cells)
    },
    "batch_min" = NULL  # handled separately below
  )
  
  # Compute per-sample target counts for batch_min or replicate for other strategies
  if (sampling_type == "batch_min") {
    # Map each sample_id to its batch
    sample_batch <- vapply(idx_by_sample, function(idx) df$batch[idx[1]], FUN.VALUE = df$batch[1])
    
    # Compute minimum sample size within each batch
    batch_min <- tapply(sample_sizes, sample_batch, min)
    
    # Assign per-sample target based on its batch minimum
    n_per_sample_vec <- as.integer(batch_min[as.character(sample_batch)])
  } else {
    # Same target for all samples
    n_per_sample_vec <- rep.int(as.integer(n_per_sample), length(idx_by_sample))
    names(n_per_sample_vec) <- names(idx_by_sample)
  }
  
  # Perform sampling within each sample_id (without replacement)
  # If sample has fewer cells than target, take all
  sampled_idx <- unlist(
    Map(function(idx, n_target) {
      n_take <- min(as.integer(n_target), length(idx))
      sample(idx, size = n_take, replace = FALSE)
    }, idx_by_sample, n_per_sample_vec),
    use.names = FALSE
  )
  
  # Return subsampled dataframe
  df[sampled_idx, , drop = FALSE]
}

### Subsamples a fixed number of rows per batch from a dataframe
# df: A dataframe containing all observations
# batch_list: A character vector of batch IDs 
# batch_colm: The column name in df that contains batch identifiers
# num: Number of cells to downsample to
# seed: Integer seed for reproducibility of random sampling
# Returns a dataframe with an equal number of randomly sampled rows per batch
subsample_batch <- function(df,
                            batch_list,
                            batch_colm, 
                            num,
                            seed) {
  
  # Define number of rows to sample per batch
  num <- num
  
  # Set seed for reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Initialize vector to store sampled row indices
  subsamp <- c()
  
  # Loop through each batch in batch_list
  for (fi in 1:length(batch_list)) {
    # Find indices of rows belonging to the current batch
    tmp <- which(df[, batch_colm] == batch_list[fi])
    
    # Randomly sample 'num' indices from this batch
    subsamp <- c(subsamp, tmp[sample(1:length(tmp), num)])
  }
  
  # Subset the dataframe using sampled indices
  rdf <- df[subsamp, ]
  
  # Return the subsampled dataframe
  return(rdf)
}

subsample_batch <- function(df,
                            batch_list,
                            batch_colm,
                            num,
                            seed) {
  
  # Check if df is a data frame
  stopifnot(is.data.frame(df))
  
  # Check if batch column names are in the data frame
  if (!batch_colm %in% names(df)) {
    stop(sprintf("Column '%s' not found in df.", batch_colm))
  }
  
  # Keep only requested batches 
  if (!missing(batch_list) && length(batch_list) > 0) {
    df <- df[df[[batch_colm]] %in% batch_list, , drop = FALSE]
  }
  
  if (nrow(df) == 0L) {
    warning("subsample_batch(): input df is empty after filtering.")
    return(df)
  }
  
  # Set seed
  set.seed(seed)
  
  # Split by batch and sample per batch
  by_batch <- split(df, df[[batch_colm]])
  
  sampled <- lapply(by_batch, function(d) {
    n <- nrow(d)
    size <- min(num, n)  
    if (n == 0L) return(d[integer(0), , drop = FALSE])
    idx <- sample.int(n, size = size, replace = FALSE)
    d[idx, , drop = FALSE]
  })
  
  out <- do.call(rbind, sampled)
  rownames(out) <- NULL
  out
}

### Checks unimodality of marker distributions across control samples
# df: A dataframe containing marker expression columns and a 'control' column
# markers: A character vector of marker column names 
# control_list: A character vector of control sample IDs 
# Returns a list with two character vectors:
# unimodal_markers: markers that are unimodal across all controls
# bimodal_markers: markers that are not unimodal in at least one control
check_unimodality <- function(df, markers, control_list) {
  
  # Initialize collections for unimodal and bimodal markers
  unimodal_markers <- character()
  bimodal_markers  <- character()
  
  # Iterate over each marker
  for (m in markers) {
    
    # Assume unimodal until proven otherwise 
    marker_is_unimodal_everywhere <- TRUE  
    
    # Iterate over each control sample
    for (ctrl in control_list) {
      
      # Subset data for the current control
      tmp <- df[df$control == ctrl, , drop = FALSE]
      
      # Safely extract the current marker vector by name
      x <- tmp[[m]]
      
      # Drop NAs
      x <- x[!is.na(x)]
      
      # If there is no data for this control/marker, skip this control
      if (length(x) < 2) next
      
      # Check unimodality for this control 
      is_uni <- is.unimodal(x, min.size = 0.05)
      
      # If any control is not unimodal
      if (!isTRUE(is_uni)) {
        marker_is_unimodal_everywhere <- FALSE
        break  # No need to check other controls for this marker
      }
    }
    
    # Classify the marker based on checks across all controls
    if (marker_is_unimodal_everywhere) {
      unimodal_markers <- c(unimodal_markers, m)
    } else {
      bimodal_markers  <- c(bimodal_markers, m)
    }
  }
  
  # Ensure uniqueness (just in case)
  unimodal_markers <- unique(unimodal_markers)
  bimodal_markers  <- unique(bimodal_markers)
  
  # Return a named list of results
  list(
    unimodal_markers = unimodal_markers,
    bimodal_markers  = bimodal_markers
  )
}

### Computes automated marker thresholds across controls and saves them to CSV
# df: A dataframe containing marker expression columns and metadata (including 'control' and batch column)
# markers: A character vector of marker column names
# control_list: A character vector of control sample IDs 
# batch_list: A character vector of batch IDs 
# batch_colm: The column name in df that contains batch identifiers (used by subsample_batch)
# subsample: Either "yes" or "no"; if "yes", subsamples per batch before threshold computation
# seed: Integer seed for reproducibility of random sampling
# num: Number of cells to downsample to
# dataset_type: Either "flow" or "mass"; selects the appropriate threshold function
# output_dir: Directory where the CSV will be saved
automated_threshold_DF <- function(df,
                                   markers,
                                   control_list,
                                   batch_list,
                                   batch_colm,
                                   subsample = "yes",
                                   seed = 250,
                                   num = 20000,
                                   dataset_type = "flow",
                                   output_dir) {
  
  # Allowed option values
  valid_options <- c("yes", "no")
  data_options <- c("flow", "mass")
  
  # Normalize user inputs
  subsample <- tolower(subsample)
  dataset_type <- tolower(dataset_type)
  
  # Validate 'subsample' input
  if (!subsample %in% valid_options) {
    stop("Invalid input: please specify 'yes' or 'no'")
  }
  
  # Validate 'dataset_type' input
  if (!dataset_type %in% data_options) {
    stop("Invalid input: please specify if dataset is 'flow' or 'mass'")
  }
  
  # Initialize the output dataframe
  auto_cutoff <- data.frame()
  
  # Ensure provided control_list exists in df
  if (!all(control_list %in% unique(df$control))) {
    stop("Check if dataframe has all unique control samples")
  }
  
  # Progress bar styling
  colrSet <- crayon::make_style("#539DDD")
  
  # Total iterations for progress bar tracking
  total_iter <- length(markers) * length(control_list)
  current_iter <- 0
  
  # Create progress bar
  pBar <- progress::progress_bar$new(
    format = paste0(
      colrSet("Calculating automated marker threshold"),
      " [", colrSet(":bar"), "] ",
      colrSet(":percent"), " | Marker in ", length(control_list),
      " control samples: :current/:total ",
      "| Elapsed: :elapsed | ETA: :eta"
    ),
    total = total_iter,
    clear = FALSE,
    width = 80,
    complete = "=",
    incomplete = "-"
  )
  
  # Branch on dataset type to select appropriate threshold function
  if (dataset_type == "flow") {
    for (i in 1:length(control_list)) {
      
      # Filter dataframe to current control sample
      tmp <- df[(df$control == control_list[i]), ]
      
      # Optional, per-batch subsampling
      if (subsample == "yes") {
        tmp <- data.frame(subsample_batch(df = tmp, batch_list = batch_list, batch_colm = batch_colm, num = num, seed = seed))
      }
      else {
        tmp <- tmp
      }
      
      # Compute thresholds per marker for this control
      for (m in 1:length(markers)) {
        
        current_iter <- current_iter + 1
        
        # Update progress bar
        pBar$tick(tokens = list(current = current_iter))
        
        # Use threshold function
        myfile <- data.frame(marker = markers[m], cutoff = auto_threshold_func(tmp[, m]))
        myfile$control <- as.character(control_list[i])
        auto_cutoff <- rbind(auto_cutoff, myfile)
      }
      
      # Message per control
      cat("Marker thresholds for control", i, "of", length(control_list), "calculated", "\n")
    }
  }
  else {
    for (i in 1:length(control_list)) {
      
      # Filter dataframe to current control sample
      tmp <- df[(df$control == control_list[i]), ]
      
      # Optional per-batch subsampling
      if (subsample == "yes") {
        tmp <- data.frame(subsample_batch(df = tmp, batch_list = batch_list, batch_colm = batch_colm, num = num, seed = seed))
      }
      else {
        tmp <- tmp
      }
      
      # Compute thresholds per marker for this control
      for (m in 1:length(markers)) {
        
        current_iter <- current_iter + 1
        
        # Update progress bar
        pBar$tick(tokens = list(current = current_iter))
        
        # Use mass cytometry-based threshold function
        myfile <- data.frame(marker = markers[m], cutoff = auto_threshold_func_w_cytof(tmp[, m]))
        myfile$control <- as.character(control_list[i])
        auto_cutoff <- rbind(auto_cutoff, myfile)
      }
      
      # Informative message per control
      cat("marker thresholds for control", i, "of", length(control_list), "calculated", "\n")
    }
  }
  
  # Aggregate across controls: median cutoff per marker
  auto_cutoff <- data.frame(auto_cutoff %>% dplyr::group_by(marker) %>% 
                              dplyr::summarise(cutoff = median(cutoff)))
  
  # Round cutoff values and order by input marker list
  auto_cutoff[, "cutoff"] <- round(auto_cutoff[, "cutoff"], 3)
  auto_cutoff <- auto_cutoff %>% arrange(match(marker, markers))
  
  # Build output CSV path
  file_path <- file.path(output_dir, "automated_marker_thresholds.csv")
  
  # Save the cutoff table to CSV
  tryCatch({
    readr::write_csv(auto_cutoff, file = file_path)
    cat("Marker threshold table saved as", file_path)
  }, error = function(e) {
    stop("Error saving table as CSV: ", e$message)
  })
  
  # Return the dataframe of automated thresholds
  return(auto_cutoff)
}

########### ~~~~ Initial visualizations (# Cells, MDS, UMAP) ~~~~ ########### 

### Generates a PNG showing the number of cells per batch across control samples
# df: A dataframe containing at least the batch and control columns
# batch_colm: The column name in df that contains batch identifiers
# control_colm: The column name in df that contains control sample identifiers
# batch_list: A character vector of batch IDs 
# control_list: A character vector of control sample IDs 
# output_dir: Directory where the PNG will be saved
# axis_size: Numeric font size for axis labels/titles in the bar plots
# width, height: Numeric values specifying the dimensions of the output plot in inches
generate_number_of_cells_plot <- function(df,
                                          batch_colm,
                                          control_colm,
                                          batch_list,
                                          control_list,
                                          output_dir,
                                          axis_size = 20,
                                          width = 11.5,
                                          height = 8) {
  
  # File name for output PNG
  file_name <- "number_cells_per_batch_per_control.png"
  
  # Build full file path
  file_path <- file.path(output_dir, file_name)
  
  # Create sample_id if one doesn't exist
  if (is.null(df[, "sample_id"])) {
    df[, "sample_id"] <- paste0(df[, batch_colm], "_", df[, control_colm]) 
  }
  
  # Compute the number of cells per sample_id
  ncells <- data.frame(table(df[, "sample_id"]))
  colnames(ncells) <- c("sample_id", "Freq")
  
  # Match to append batch and control metadata columns
  mm <- match(ncells$sample_id, df$sample_id)
  ncells$batch <- df$batch[mm]
  ncells$control <- df$control[mm]
  ncells$Freq <- as.numeric(ncells$Freq)
  
  # Order factor levels based on provided lists
  ncells$control <- factor(ncells$control, levels = control_list)
  ncells$batch <- factor(ncells$batch, levels = batch_list)
  
  # Number of facet rows (half the number of control samples)
  num_cols <- as.numeric(length(control_list) / 2)
  
  # If you'd like to ensure an odd number is rounded appropriately, uncomment:
  # if (num_cols %% 2 != 0) {
  #   num_cols <- as.numeric((length(control_list) + 1) / 2)
  # }
  
  # Build the bar plot with facets per control sample
  nplot <- ggplot(ncells, aes(x = Freq, y = batch)) + theme_bw() + 
    geom_bar(colour = "#121518", fill = "#7A7A78", position = "dodge", stat = "identity", width = 0.8) +
    facet_wrap(~control, nrow = num_cols) +
    scale_y_discrete(limits = rev(batch_list)) +
    ylab("Batch") + xlab("Number of Cells") +
    ggtitle("Number of cells, per batch per sample") +
    scale_x_continuous(
      limits = c(0, max(ncells$Freq)),
      expand = expansion(mult = c(0, .02)),
      labels = scales::label_comma()) +
    theme(axis.text.y = element_text(size = axis_size + 5, colour = "#4B4B4B"),
          axis.text.x = element_text(size = axis_size - 6),
          axis.title = element_text(size = axis_size + 10, colour = "#000000"),
          strip.text = element_text(size = axis_size + 2, face = 'bold'),
          plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
          panel.spacing = unit(1, "lines"))
  
  # Save the plot to file
  ggplot2::ggsave(file_path, plot = nplot, width = width, height = height, dpi = 400)
  
  # Notify user of completion
  cat("Number of cells plot saved as", file_path, "\n")
}

### Generates a PNG of MDS plots per control sample using median marker expression
# df: A dataframe containing marker expression columns and metadata (batch, control, stimulation)
# batch_colm: The column name in df that contains batch identifiers
# control_colm: The column name in df that contains control sample identifiers
# marker_list: A character vector of marker column names
# sample_colours: A named vector of colours for control samples
# output_dir: Directory where the PNG will be saved
# axis_size: Numeric font size for axis labels/titles in the MDS plot
# width, height: Numeric values specifying the dimensions of the output plot in inches
generate_MDS_plot <- function(df,
                              batch_colm,
                              control_colm,
                              marker_list,
                              sample_colours,
                              output_dir,
                              axis_size = 22,
                              width = 12,
                              height = 12) {
  
  # File name for output PNG
  file_name <- "MDS.png"
  
  # Build full file path
  file_path <- file.path(output_dir, file_name)
  
  # Create sample_id if one doesn't exist
  if (is.null(df[, "sample_id"])) {
    df[, "sample_id"] <- paste0(df[, batch_colm], "_", df[, control_colm]) 
  }
  
  # Unique lists for ordering and labeling
  batch_list <- as.character(unique(df[, batch_colm])) ## unique batch names
  control_list <- as.character(unique(df[, control_colm])) ## unique control sample names
  
  # Calculate median marker expression per sample_id
  expr_median <- data.frame(df[, marker_list], sample_id = df[, "sample_id"]) #
  expr_median <- expr_median %>% dplyr::group_by(sample_id) %>% 
    dplyr::summarize_all(list(median))
  
  # Transpose median marker expression matrix for MDS (features in rows)
  exp_med_tbl <- t(expr_median[, -1])
  colnames(exp_med_tbl) <- expr_median$sample_id
  
  # Compute MDS dimensions (no plotting)
  mds <- limma::plotMDS(exp_med_tbl, plot = FALSE)
  
  # Build a dataframe for plotting MDS results
  mdsdf <- data.frame(MDS1 = mds$x, MDS2 = mds$y,
                      sample_id = colnames(exp_med_tbl))
  
  # Match to append batch/control/stimulation metadata
  mm <- match(mdsdf$sample_id, df$sample_id)
  mdsdf$batch <- df$batch[mm]
  mdsdf$control <- df$control[mm]
  mdsdf$stimulation <- df$stimulation[mm]
  
  # Order factors for consistent facet/legend display
  mdsdf$control <- factor(mdsdf$control, levels = control_list)
  mdsdf$batch <- factor(mdsdf$batch, levels = batch_list)
  
  # Determine plotting limits (square window based on min/max across axes)
  if (min(mdsdf$MDS1) < min(mdsdf$MDS2)) {
    min_lim <- as.numeric(round(min(mdsdf$MDS1), 1))
  }
  else {
    min_lim <- as.numeric(round(min(mdsdf$MDS2), 1))
  }
  
  if (max(mdsdf$MDS1) > max(mdsdf$MDS2)) {
    max_lim <- as.numeric(round(max(mdsdf$MDS1), 1))
  }
  else {
    max_lim <- as.numeric(round(max(mdsdf$MDS2), 1))
  }
  
  # Break sequence (not used directly in plot; kept for potential grid/scale control)
  brk <- seq(min_lim, max_lim, by = 0.1)
  
  # Build the MDS scatter plot
  mds_plot <- ggplot(mdsdf, aes(x = MDS1, y = MDS2, color = control)) +
    geom_point(size = 2) +
    geom_label_repel(aes(label = batch, fontface = 'bold', size = 5.2), 
                     max.overlaps = Inf, , show.legend = F) + # 
    scale_color_manual("Control Sample", values = sample_colours) +
    coord_cartesian(xlim = c(min_lim, max_lim), ylim = c(min_lim, max_lim)) +
    guides(color = guide_legend(title = "Sample")) +
    labs(x = "MDS 1", y = "MDS 2") + 
    ggtitle("MDS plot, colored by sample and labeled by batch") +
    theme_bw() +
    theme(axis.text = element_text(size = axis_size),
          axis.title = element_text(size = axis_size + 4),
          legend.text = element_text(size = axis_size + 2),
          legend.title = element_text(size = axis_size + 2),
          plot.title = element_text(size = axis_size + 2, face = "bold", hjust = 0.5))
  
  # Save the plot to file
  ggplot2::ggsave(file_path, plot = mds_plot, width = width, height = height, dpi = 500)
  
  # Notify user of completion
  cat("MDS plot saved as", file_path, "\n")
}

### Subsets a fixed number of rows per sample_id from a dataframe
# df: A dataframe containing a "sample_id" column
# samps: A character vector of sample_id values 
# num: Number of rows to sample per sample_id (default = 20000)
# Returns a dataframe with an equal number of randomly sampled rows per sample_id
subsetting <- function(df,
                       samps,
                       num = 20000) {
  
  # Initialize vector for sampled row indices
  subsamp <- c()
  
  # Loop through each requested sample_id
  for (fi in 1:length(samps)) {
    # Find indices of rows belonging to the current sample_id
    tmp <- which(df[, "sample_id"] == samps[fi])
    # Randomly sample 'num' indices from this sample_id
    subsamp <- c(subsamp, tmp[sample(1:length(tmp), num)])
  }
  
  # Subset the dataframe using sampled indices
  rdf <- df[subsamp, ]
  
  # Return the subsetted dataframe
  return(rdf)
}

### Generates a PNG UMAP of control samples using selected markers
# df: A dataframe containing marker expression columns and metadata (batch, control, and optionally sample_id)
# batch_colm: The column name in df that contains batch identifiers
# control_colm: The column name in df that contains control sample identifiers
# marker_list: A character vector of marker column names 
# batch_list: A character vector of batch IDs 
# control_list: A character vector of control sample IDs 
# sample_ids: A character vector of sample_id values 
# batch_colours: A named vector of colours for batches (names must match batch_list)
# output_dir: Directory where the PNG will be saved
# axis_size: Numeric font size for axis labels/titles in the UMAP plot
# width, height: Numeric values specifying the dimensions of the output plot in inches
generate_UMAP <- function(df,
                          batch_colm,
                          control_colm,
                          marker_list,
                          batch_list,
                          control_list,
                          sample_ids,
                          batch_colours,
                          output_dir,
                          axis_size = 22,
                          width = 12,
                          height = 12) {
  
  # Inform user that UMAP computation may be slow
  cat("This process may take a long time, please be patient\n")
  
  # Output file name and full path
  file_name <- "UMAP.png"
  
  file_path <- file.path(output_dir, file_name)
  
  # Create sample_id column if one doesn't exist
  if (is.null(df[, "sample_id"])) {
    df[, "sample_id"] <- paste0(df[, batch_colm], "_", df[, control_colm]) 
  }
  
  # Determine the number of rows to sample per sample_id (bounded by 2500)
  nums <- as.numeric(round(min(table(df$sample_id)) / length(batch_list)))
  if (nums > 2500) {
    nums <- 2500
  }
  
  # Subsample dataframe to create UMAP input (balanced across sample_ids)
  subdf <- data.frame(subsetting(df, sample_ids, nums))
  
  # Capture vector of sample_id values for alignment check
  samp_filt <- subdf[, "sample_id"]
  
  cat("Calculating UMAP dimensions...\n")
  
  # Compute UMAP dimensions based on selected markers
  out_umap <- umap::umap(subdf[, marker_list], config = umap.defaults)
  dims_umap <- out_umap$layout
  colnames(dims_umap) <- c("UMAP_1", "UMAP_2")
  
  # Sanity check: UMAP rows should match number of sampled cells
  stopifnot(nrow(dims_umap) == length(samp_filt))
  
  # Build dataframe for plotting (add sample_id and a type tag)
  dims_umap <- cbind(as.data.frame(dims_umap), sample_id = samp_filt, type = "UMAP")
  
  # Match to append batch and control metadata; set factor levels for consistency
  mm <- match(dims_umap$sample_id, subdf$sample_id)
  dims_umap$batch <- subdf$batch[mm]
  dims_umap$control <- subdf$control[mm]
  dims_umap$batch <- factor(dims_umap$batch, levels = batch_list)
  dims_umap$control <- factor(dims_umap$control, levels = control_list)
  
  # Ensure batch_colours are named by batch_list for consistent mapping
  names(batch_colours) <- batch_list
  
  # Minimal dataframe used for plotting
  umapdf <- subset(dims_umap, select = c(UMAP_1, UMAP_2, control, batch))
  
  cat("Generating UMAP plot...\n")
  
  # Build UMAP plot colored by batch
  umap_plot <- ggplot(umapdf, aes(x = UMAP_1, y = UMAP_2, colour = batch)) +
    geom_point(alpha = 0.7, size = 0.1) +
    scale_colour_manual("Batch", values = batch_colours) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    ggtitle("UMAP of cells, colored by batch") +
    theme_bw() +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    theme(axis.text = element_text(size = axis_size),
          axis.title = element_text(size = axis_size + 4),
          legend.text = element_text(size = axis_size + 2),
          legend.title = element_text(size = axis_size + 2),
          plot.title = element_text(size = axis_size + 2, face = "bold", hjust = 0.5))
  
  # Save the plot to disk
  ggplot2::ggsave(file_path, plot = umap_plot, width = width, height = height, dpi = 400)
  
  # Notify user of completion
  cat("UMAP plot saved as", file_path, "\n")
}

########### ~~~~ IQR Metrics ~~~~ ########### 

### Generates a dataframe of negative MFI per sample_id for each marker
# df: A data frame with rows for individual cells and columns including marker intensities
# markers: A character vector of marker names
# cutf: A 2-column object (e.g., data.frame or matrix) whose second column contains cutoffs per marker
# The order of cutoffs is assumed to align with 'markers'
# sample_col_name: Column name in 'df' holding the control/sample ID 
# batch_col_name: Column name in 'df' holding the batch ID 
# min_cells: If TRUE, drop sample_ids with <100 negative cells for a given marker
neg_mfi_df_all <- function(df,
                           markers,
                           cutf,
                           sample_col_name,
                           batch_col_name,
                           min_cells = TRUE) {
  
  # Rename columns
  names(df)[names(df) == sample_col_name] <- "control"
  names(df)[names(df) == batch_col_name]  <- "batch"
  
  # Ensure sample_id exists
  if (!("sample_id" %in% names(df))) {
    df$sample_id <- paste0(df$batch, "_", df$control)
  }
  
  # Convert cutoff
  cutv <- as.numeric(cutf[, 2])
  
  # Output container
  out <- vector("list", length(markers))
  drop_log <- list()
  
  # Collect dropped sample-marker combinations (for summary printing)
  dropped_samples <- data.frame(
    sample_id = character(),
    marker    = character(),
    neg_Freq  = integer(),
    stringsAsFactors = FALSE
  )
  
  for (m in seq_along(markers)) {
    
    mn     <- markers[m]
    cutoff <- cutv[m]
    
    # Extract numeric vector
    expr <- suppressWarnings(as.numeric(df[[mn]]))
    
    # Filter negatives
    keep     <- expr < cutoff
    tmp_expr <- expr[keep]
    tmp_sid  <- df$sample_id[keep]
    
    if (length(tmp_expr) == 0) next
    
    # Count negative cells per sample_id 
    freq  <- table(tmp_sid)
    small <- names(freq[freq < 100])
    
    # Compute medians 
    med <- tapply(tmp_expr, tmp_sid, median)
    
    # Remove small groups (and record them for summary)
    if (min_cells && length(small) > 0) {
      dropped_samples <- rbind(
        dropped_samples,
        data.frame(
          sample_id = small,
          marker    = mn,
          neg_Freq  = as.integer(freq[small]),
          stringsAsFactors = FALSE
        )
      )
      med <- med[!(names(med) %in% small)]
      drop_log[[mn]] <- small
    }
    
    if (length(med) == 0) next
    
    # Expand metadata
    sid <- names(med)
    idx <- match(sid, df$sample_id)
    
    out[[m]] <- data.frame(
      sample_id = sid,
      value     = as.numeric(med),
      control   = df$control[idx],
      batch     = df$batch[idx],
      marker    = mn,
      stringsAsFactors = FALSE
    )
  }
  
  ndf <- do.call(rbind, out)
  if (!is.null(ndf) && nrow(ndf) > 0) {
    ndf$type <- "-MFI"
  } else {
    ndf <- data.frame(sample_id = character(),
                      value     = numeric(),
                      control   = character(),
                      batch     = character(),
                      marker    = character(),
                      type      = character(),
                      stringsAsFactors = FALSE)
  }
  
  # Summary of dropped items (<100 negative cells)
  cat("\nSummary of dropped sample-marker combinations (<100 negative cells):\n")
  if (nrow(dropped_samples) > 0) {
    print(dropped_samples[, c("sample_id", "marker", "neg_Freq")])
  } else {
    cat("No sample-marker combinations were dropped for <100 negative cells.\n")
  }
  
  return(ndf)
}

### Generates a dataframe of positive MFI per sample_id for each marker
# df: A data frame with rows for individual cells and columns including marker intensities
# markers: A character vector of marker names
# cutf: A 2-column object (e.g., data.frame or matrix) whose second column contains cutoffs per marker
# The order of cutoffs is assumed to align with 'markers'
# sample_col_name: Column name in 'df' holding the control/sample ID 
# batch_col_name: Column name in 'df' holding the batch ID 
# min_cells: If TRUE, drop sample_ids with <100 positive cells for a given marker
pos_mfi_df_all <- function(df,
                           markers,
                           cutf,
                           sample_col_name,
                           batch_col_name,
                           min_cells = TRUE) {
  
  # Rename columns to standard names
  names(df)[names(df) == sample_col_name] <- "control"
  names(df)[names(df) == batch_col_name]  <- "batch"
  
  # Ensure sample_id exists
  if (!("sample_id" %in% names(df))) {
    df$sample_id <- paste0(df$batch, "_", df$control)
  }
  
  # Convert cutoff vector (assumes alignment with `markers`)
  cutv <- as.numeric(cutf[, 2])
  
  # Output container
  out <- vector("list", length(markers))
  drop_log <- list()
  
  # Collect dropped sample-marker combinations (for summary printing)
  dropped_samples <- data.frame(
    sample_id = character(),
    marker    = character(),
    pos_Freq  = integer(),
    stringsAsFactors = FALSE
  )
  
  for (m in seq_along(markers)) {
    
    mn     <- markers[m]
    cutoff <- cutv[m]
    
    # Extract numeric expression vector for the marker
    expr <- suppressWarnings(as.numeric(df[[mn]]))
    
    # Keep positive cells (>= cutoff)
    keep     <- expr >= cutoff
    tmp_expr <- expr[keep]
    tmp_sid  <- df$sample_id[keep]
    
    # If no positive cells for this marker, skip
    if (length(tmp_expr) == 0) next
    
    # Count positive cells per sample_id
    freq  <- table(tmp_sid)
    small <- names(freq[freq < 100])
    
    # Median of positive expression per sample_id
    med <- tapply(tmp_expr, tmp_sid, median)
    
    # Optionally drop sample_ids with <100 positive cells
    if (min_cells && length(small) > 0) {
      # Record dropped items and their positive counts for summary
      dropped_samples <- rbind(
        dropped_samples,
        data.frame(
          sample_id = small,
          marker    = mn,
          pos_Freq  = as.integer(freq[small]),
          stringsAsFactors = FALSE
        )
      )
      med <- med[!(names(med) %in% small)]
      drop_log[[mn]] <- small
    }
    
    # If nothing remains after dropping, skip
    if (length(med) == 0) next
    
    # Expand metadata from original df
    sid <- names(med)
    idx <- match(sid, df$sample_id)
    
    out[[m]] <- data.frame(
      sample_id = sid,
      value     = as.numeric(med),
      control   = df$control[idx],
      batch     = df$batch[idx],
      marker    = mn,
      stringsAsFactors = FALSE
    )
  }
  
  # Bind per-marker results
  ndf <- do.call(rbind, out)
  if (!is.null(ndf) && nrow(ndf) > 0) {
    ndf$type <- "+MFI"
  } else {
    ndf <- data.frame(sample_id = character(),
                      value     = numeric(),
                      control   = character(),
                      batch     = character(),
                      marker    = character(),
                      type      = character(),
                      stringsAsFactors = FALSE)
  }
  
  # Summary of dropped items (<100 positive cells)
  cat("\nSummary of dropped sample-marker combinations (<100 positive cells):\n")
  if (nrow(dropped_samples) > 0) {
    # Print selected columns in a tidy form
    print(dropped_samples[, c("sample_id", "marker", "pos_Freq")])
  } else {
    cat("No sample-marker combinations were dropped for <100 positive cells.\n")
  }
  
  return(ndf)
}

### Generates a dataframe of percent positive cells per sample_id for each marker
# df: A data frame with rows for individual cells and columns including marker intensities
# markers: A character vector of marker names
# cutf: A 2-column object (data frame or matrix) whose second column contains cutoffs per marker.
# The order of cutoffs is assumed to align with 'markers'
# sample_col_name: Column name in 'df' holding the control/sample ID 
# batch_col_name: Column name in 'df' holding the batch ID 
percent_pos_df_all <- function(df,
                               markers,
                               cutf,
                               sample_col_name,
                               batch_col_name) {
  
  # Standardize column names
  names(df)[names(df) == sample_col_name] <- "control"
  names(df)[names(df) == batch_col_name]  <- "batch"
  
  # Create sample_id (control_batch)
  df <- cbind(df, sample_id = paste0(df[, "control"], "_", df[, "batch"]))
  
  # Total cells per sample_id
  ncdf <- data.frame(table(df$sample_id))
  names(ncdf)[names(ncdf) == "Var1"] <- "sample_id"
  
  # Output accumulator (preserves original initialization behavior)
  perc_df <- data.frame(sample_id = NA, pos_Freq = NA, control = NA, batch = NA,
                        tot_cells = NA, value = NA, marker = NA)
  
  if (!is.null(df[, "control"])) {
    for (m in 1:length(markers)) {
      mn <- markers[m]
      
      # Filter to positive cells for this marker (>= cutoff)
      tmp <- data.frame(expression = df[, mn],
                        sample_id = df[, "sample_id"],
                        control   = df[, "control"],
                        batch     = df[, "batch"],
                        marker    = mn) %>%
        dplyr::filter(expression >= cutf[m, 2])
      
      # Count positive cells per sample_id
      pMark_freq <- data.frame(table(tmp$sample_id))
      colnames(pMark_freq) <- c("sample_id", "pos_Freq")
      
      # Recover control and batch info
      mm <- match(pMark_freq$sample_id, tmp$sample_id)
      pMark_freq$control <- tmp$control[mm]
      pMark_freq$batch   <- tmp$batch[mm]
      
      # Compute percent positive using total cells per sample_id
      mm <- match(pMark_freq$sample_id, ncdf$sample_id)
      pMark_freq$tot_cells <- ncdf$Freq[mm]
      pMark_freq <- pMark_freq %>% dplyr::mutate(value = ((pos_Freq / tot_cells) * 100))
      pMark_freq$marker <- mn
      
      # Append to output
      perc_df <- rbind(perc_df, pMark_freq)
    }
    
    # Drop the initial NA row and keep only the required columns
    perc_df <- perc_df[-1, ]
    perc_df <- data.frame(subset(perc_df, select = c(sample_id, value, control, batch, marker)))
    perc_df$type <- "%pos"
    
    return(perc_df)
  } else {
    print("Please provide control information of each batch and protein channel")
  }
}


### IQR outlier flagging function

### Detects outliers in a numeric vector using the IQR rule
# x: A numeric vector of values
# Returns a logical vector indicating TRUE for outliers and FALSE otherwise
detect_outlier <- function(x) { 
  
  iqr <- 1.5
  Q1 <- quantile(x, 0.25)
  Q3 <- quantile(x, 0.75)
  
  return(x > Q3 + (iqr * IQR(x)) | x < Q1 - (iqr * IQR(x)))
}

### Flags outlier batches per marker/control using the IQR rule
# df: A data frame containing rows of summarized values per cell/sample with columns including:
# 'type' (metric type), 'marker', 'value', 'batch', and a sample/control identifier
# markers: A character vector of marker names
# sample_col_name: Column name in 'df' holding the control/sample ID 
# batch_col_name: Column name in 'df' holding the batch ID 
# Returns a data frame identical to input rows for the selected markers, with an added column:
# 'flagged_batch': "none" if value is NA or not an outlier; otherwise the batch ID
iqr_metric_function <- function(df,
                                markers,
                                sample_col_name,
                                batch_col_name) {
  
  # Standardize column names
  names(df)[names(df) == sample_col_name] <- "control"
  names(df)[names(df) == batch_col_name]  <- "batch"
  
  # Preserve metric type (assumes a single type present; keeps original value)
  type <- as.character(unique(df$type))
  
  iqr15_metric <- list()
  
  # Per-marker outlier detection and flagging
  for (m in markers) {
    iqr15_metric[[m]] <- list()
    iqr15_metric[[m]] <- df %>%
      dplyr::filter(marker == m) %>%
      dplyr::group_by(control) %>%
      dplyr::mutate(
        flagged_batch = ifelse(
          is.na(value), "none",
          ifelse(detect_outlier(value), paste0(batch), "none")
        )
      )
  }
  
  # Unlist into a single data frame
  MFIdf <- do.call(rbind.data.frame, iqr15_metric)
  rownames(MFIdf) <- NULL
  MFIdf$type <- as.character(type)
  MFIdf <- data.frame(MFIdf)
  
  return(MFIdf)
}

### IQR metrics dataframe for every marker, control sample, and batch
# neg_df: A dataframe of per-sample negative MFI medians (type = "-MFI")
# pos_df: A dataframe of per-sample positive MFI medians (type = "+MFI")
# percent_pos_df: A dataframe of percent positive cells per sample (type = "%pos")
# markers: A character vector of marker names to include in the IQR flagging
# sample_col_name: Column name in 'df' holding the control/sample ID 
# batch_col_name: Column name in 'df' holding the batch ID 
# output_dir: Directory where the CSV of IQR flags across markers/controls/batches will be saved
iqr_dataframe <- function(neg_df,
                          pos_df,
                          percent_pos_df,
                          markers,
                          sample_col_name,
                          batch_col_name,
                          output_dir) {
  
  # Combine the three metric dataframes
  df <- data.frame(rbind(neg_df, pos_df, percent_pos_df))
  type_list <- c("-MFI", "+MFI", "%pos")
  
  iqr_df <- data.frame()
  
  # Compute IQR flags per metric type and marker, grouped by control
  for (t in 1:length(type_list)) {
    metric <- type_list[t]
    tmpdf <- df[(df$type == type_list[t]), ]
    
    # Flag outliers using IQR rule via iqr_metric_function (adds 'flagged_batch')
    iqr_tmp <- data.frame(iqr_metric_function(df = tmpdf,
                                              markers = markers,
                                              sample_col_name = sample_col_name,
                                              batch_col_name = batch_col_name))
    
    # Binary batch effect indicator: 1 if flagged, else 0 (NA values -> 0)
    iqr_tmp <- iqr_tmp %>%
      dplyr::mutate(batch_effect = ifelse(is.na(value), 0,
                                          ifelse(flagged_batch == "none", 0, 1)))
    
    # Plotting color helper: NA -> light grey; flagged -> red; not flagged -> grey
    iqr_tmp <- iqr_tmp %>%
      dplyr::mutate(cols = ifelse(is.na(value), "#F5F5F5",
                                  ifelse(batch_effect == 1, "#8E000C", "#AEAEAE")))
    
    iqr_df <- rbind(iqr_df, iqr_tmp)
  }
  
  # Save the combined flags table
  filename <- paste0("IQR_metrics_flags.csv")
  file_path <- file.path(output_dir, filename)
  
  # Reformat dataframe for CSV file
  iqr_output <- iqr_df %>%
    dplyr::select(batch, control, marker, type, value, batch_effect, flagged_batch) %>%
    dplyr::rename(
      above_threshold = batch_effect
    ) %>%
    dplyr::mutate(
      above_threshold = as.logical(as.integer(above_threshold)),
      type = dplyr::recode(type,
                           "+MFI"   = "Positive MFI",
                           "-MFI"   = "Negative MFI",
                           "%pos"   = "Percent positive"))
  
  tryCatch({
    readr::write_csv(iqr_output, file = file_path)
    cat("IQR metrics table saved as", file_path, "\n")
  }, error = function(e) {
    stop("Error saving table as CSV: ", e$message)
  })
  
  return(iqr_df)
}

### Generates an IQR-based boxplot for a marker with flexable y-axis limits
# df: A data frame containing rows for multiple markers and metric types, including columns:
# 'type' (metric type, e.g., "%pos", "MFI"), 'marker', 'value', 'batch', 'control', 'cols', and 'flagged_batch'
# type: A character string specifying the metric type to filter by ("%pos", "+MFI", "-MFI")
# marker: A single marker name to plot (must exist in df$marker)
# control_list: A character vector of control sample IDs 
# batch_list: A character vector of batch IDs
# colours: A named vector of colors for batches (names should match batch_list)
# ctrl_labs: Optional, a vector of labels for controls (used on x-axis). If NULL, control_list is shown
# seed: Integer seed for reproducibility of jittering in the plot
iqr_boxplot <- function(df,
                        type,
                        marker,
                        control_list,
                        batch_list,
                        colours,
                        ctrl_labs = NULL,
                        seed) {
  
  # Filter by type and marker
  df_mark <- df[df$type == type & df$marker == marker, ]
  
  # Exclude batches with < 100 cells
  df_mark <- df_mark[df_mark$cols != "#F5F5F5", ]
  
  # Convert "none" to NA for flagged batches
  df_mark$flagged_batch <- ifelse(df_mark$flagged_batch == "none", NA, df_mark$flagged_batch)
  
  # Set factor levels
  df_mark$control <- factor(df_mark$control, levels = control_list)
  df_mark$batch <- factor(df_mark$batch, levels = batch_list)
  
  names(colours) <- batch_list
  
  if(!is.null(ctrl_labs)) names(ctrl_labs) <- control_list
  
  has_flags <- any(!is.na(df_mark$flagged_batch))
  
  # Apply horizontal jitter to points
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  df_mark$jitter_x <- as.numeric(df_mark$control) + runif(nrow(df_mark), -0.15, 0.15)
  
  # Determine y-range for vertical nudge
  y_range <- max(df_mark$value, na.rm = TRUE) - min(df_mark$value, na.rm = TRUE)
  
  nudge_amount <- 0.05 * y_range  
  
  # Plot
  p <- ggplot(df_mark, aes(x = control, y = value)) +
    ggtitle(paste0("IQR flags for ", type, " values in ", marker, "\n")) +
    geom_boxplot(outlier.color = NA, colour = "#000000", fill = "#F2F2F2", width = 0.6) +
    geom_point(aes(x = jitter_x, fill = batch), size = 3.8, pch = 21, colour = "#4D4D4D") +
    geom_label_repel(
      data = df_mark[!is.na(df_mark$flagged_batch), ],
      aes(x = jitter_x, y = value, label = flagged_batch, color = batch),
      size = 6.8,
      show.legend = FALSE,
      force = 0.1,
      force_pull = 2,
      max.time = Inf,
      direction = "y",
      nudge_y = nudge_amount,
      segment.size = 0.4) +
    scale_fill_manual("Batch", values = colours) +
    { if (has_flags) scale_color_manual(values = colours) } +
    guides(fill = guide_legend(override.aes = list(shape = 22, size = 5, colour = "black"))) +
    xlab("Control Sample") +
    ylab(type) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 26, angle = 90, vjust = 0.5),
      axis.text.y = element_text(size = 25),
      axis.title = element_text(size = 28),
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 19, face = "bold"),
      legend.text = element_text(size = 17)
    )
  
  # Apply custom x labels if provided
  if(!is.null(ctrl_labs)) {
    p <- p + scale_x_discrete(labels = ctrl_labs)
  }
  
  return(p)
}

### Generates an IQR-based boxplot for a marker with fixed y-axis limits
# df: A data frame containing rows for multiple markers and metric types, including columns:
# 'type' (metric type, e.g., "%pos", "MFI"), 'marker', 'value', 'batch', 'control', 'cols', and 'flagged_batch'
# type: A character string specifying the metric type to filter by ("%pos", "+MFI", "-MFI")
# marker: A single marker name to plot (must exist in df$marker)
# control_list: A character vector of control sample IDs 
# batch_list: A character vector of batch IDs 
# colours: A named vector of colors for batches (names should match batch_list)
# ctrl_labs: Optional, a vector of labels for controls (used on x-axis). If NULL, control_list is shown
# y_limits: Optional, a numeric vector c(min, max) for fixed y-axis range. 
# If NULL, defaults to global min/max across all markers and batches for the selected 'type'
# seed: Integer seed for reproducibility of jittering in the plot
iqr_boxplot_fixed_coords <- function(df,
                                     type,
                                     marker,
                                     control_list,
                                     batch_list,
                                     colours,
                                     ctrl_labs = NULL,
                                     y_limits = "auto",
                                     seed) {
  
  # Filter by type only for computing global limits
  df_type <- df[df$type == type, ]
  df_type$value <- as.numeric(df_type$value)
  
  # Global min/max across all markers & controls for this type
  if (is.character(y_limits) && y_limits == "auto") {
    lim1 <- min(df_type$value, na.rm = TRUE)
    lim2 <- max(df_type$value, na.rm = TRUE)
  } else {
    lim1 <- min(y_limits)
    lim2 <- max(y_limits)
  }
  
  # Now filter by marker
  df_mark <- df_type[df_type$marker == marker, ]
  
  # Exclude batches with <100 cells
  df_mark <- df_mark[df_mark$cols != "#F5F5F5", ]
  
  # Convert "none" to NA for flagged batches
  df_mark$flagged_batch <- ifelse(df_mark$flagged_batch == "none", NA, df_mark$flagged_batch)
  
  # Set factor levels
  df_mark$control <- factor(df_mark$control, levels = control_list)
  df_mark$batch <- factor(df_mark$batch, levels = batch_list)
  names(colours) <- batch_list
  if (!is.null(ctrl_labs)) names(ctrl_labs) <- control_list
  has_flags <- any(!is.na(df_mark$flagged_batch))
  
  # Horizontal jitter
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  df_mark$jitter_x <- as.numeric(df_mark$control) + runif(nrow(df_mark), -0.15, 0.15)
  
  # Vertical nudge for labels
  y_range <- lim2 - lim1
  nudge_amount <- 0.05 * y_range  
  
  # Plot
  p <- ggplot(df_mark, aes(x = control, y = value)) +
    ggtitle(paste0("IQR flags for ", type, " values in ", marker, "\n")) +
    geom_boxplot(outlier.color = NA, colour = "#000000", fill = "#F2F2F2", width = 0.6) +
    geom_point(aes(x = jitter_x, fill = batch), size = 3.8, pch = 21, colour = "#4D4D4D") +
    geom_label_repel(
      data = df_mark[!is.na(df_mark$flagged_batch), ],
      aes(x = jitter_x, y = value, label = flagged_batch, color = batch),
      size = 6.8, show.legend = FALSE,
      force = 0.1, force_pull = 2,
      max.time = Inf, direction = "y",
      nudge_y = nudge_amount, segment.size = 0.4
    ) +
    scale_fill_manual("Batch", values = colours) +
    { if (has_flags) scale_color_manual(values = colours) } +
    guides(fill = guide_legend(override.aes = list(shape = 22, size = 5, colour = "black"))) +
    xlab("Control Sample") + ylab(type) +
    scale_y_continuous(limits = c(lim1, lim2)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 26, angle = 90, vjust = 0.5),
      axis.text.y = element_text(size = 25),
      axis.title = element_text(size = 28),
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 19, face = "bold"),
      legend.text = element_text(size = 17)
    )
  
  if (!is.null(ctrl_labs)) p <- p + scale_x_discrete(labels = ctrl_labs)
  return(p)
}

### Generates a PDF of IQR-based boxplots (one per marker) using flexable y-axis limits
# iqr.plt: A data frame containing rows for multiple markers and metric types
# Data frame should include 'type' (metric type), 'marker', 'value', 'batch', 'control', 'cols', and 'flagged_batch'
# type: A character string specifying the metric type to filter by ("%pos", "+MFI", "-MFI")
# markers: A character vector of marker names to plot
# control_list: A character vector of control sample IDs defining the x-axis order
# batch_list: A character vector of batch IDs defining legend order and color mapping
# colours: A named vector of colors for batches (names should match batch_list)
# output_dir: Directory where the PDF file will be saved
# width, height: Numeric values specifying the dimensions of the output plot in inches
# ctrl_labs: Optional, a vector of labels for controls (used on x-axis). If NULL, control_list is shown
# seed: Integer seed for reproducibility of jittering in the plot
generate_IQR_metric_boxplots <- function(iqr.plt,
                                         type,
                                         markers,
                                         control_list,
                                         batch_list,
                                         colours,
                                         output_dir,
                                         width = 11, 
                                         height = 10,
                                         ctrl_labs = NULL,
                                         seed) {
  # Define the PDF file path
  file_name <- paste0("boxplots_",type,"_free_axes.pdf")
  
  if (type == "%pos") {
    file_name <- sapply(file_name, function(f)gsub("%pos","percent_pos",f))
  }
  
  pdf_file_path <- file.path(output_dir, file_name)
  
  # Open the PDF device
  pdf(pdf_file_path, width = width, height = height, useDingbats = FALSE)
  
  # Loop through the markers and create IQR metric boxplots with flags
  for (i in 1:length(markers)) {
    
    mk <- markers[i]
    
    # Generate boxplot with free axes
    plot <- iqr_boxplot(df = iqr.plt,
                        type = type,
                        marker = mk,
                        batch_list = batch_list,
                        colours = colours, 
                        control_list = control_list,
                        ctrl_labs = ctrl_labs,
                        seed = seed)
    print(plot)
  }
  
  # Close the PDF device
  dev.off()
  
  # Notify user of completion
  cat(type, "boxplots (flexable/free axes) saved as", pdf_file_path, "\n")
  
}

### Generates a PDF of IQR-based boxplots (one per marker) using fixed y-axis limits
# iqr.plt: A data frame containing rows for multiple markers and metric types
# Data frame should include: 'type' (metric type), 'marker', 'value', 'batch', 'control', 'cols', and 'flagged_batch'
# type: A character string specifying the metric type to filter by ("%pos", "+MFI", "-MFI")
# markers: A character vector of marker names to plot; each marker will get a page/plot in the PDF
# control_list: A character vector of control sample IDs defining the x-axis order
# batch_list: A character vector of batch IDs defining legend order and color mapping
# colours: A named vector of colors for batches (names should match batch_list)
# output_dir: Directory where the PDF file will be saved
# width, height: Numeric values specifying the dimensions of the output plot in inches
# ctrl_labs: Optional, a vector of labels for controls (used on x-axis). If NULL, control_list is shown
# y_limits: Optional, a numeric vector c(min, max) for fixed y-axis range. 
# If NULL, defaults to global min/max across all markers and batches for the selected 'type'
# seed: Integer seed for reproducibility of jittering in the plot
generate_fixed_IQR_metric_boxplots <- function(iqr.plt,
                                               type,
                                               markers,
                                               control_list,
                                               batch_list,
                                               colours,
                                               output_dir,
                                               width = 11, 
                                               height = 10,
                                               ctrl_labs = NULL,
                                               y_limits = "auto",
                                               seed) {
  
  # Define the PDF file path
  file_name <- paste0("boxplots_",type,"_fixed_axes.pdf")
  
  if (type == "%pos") {
    file_name <- gsub("%pos", "percent_pos", file_name)
  }
  pdf_file_path <- file.path(output_dir, file_name)
  
  # Open the PDF device (ensure it closes on exit even if errors occur)
  pdf(pdf_file_path, width = width, height = height, useDingbats = FALSE)
  on.exit(dev.off(), add = TRUE)
  
  # Loop through the markers and create IQR metric boxplots with fixed axes
  for (i in seq_along(markers)) {
    mk <- markers[i]
    
    plot <- iqr_boxplot_fixed_coords(df = iqr.plt,
                                     type = type,
                                     marker = mk,
                                     control_list = control_list,
                                     batch_list = batch_list,
                                     colours = colours,
                                     ctrl_labs = ctrl_labs,
                                     y_limits = y_limits,
                                     seed = seed)  
    print(plot)
  }
  
  # Notify user of completion
  cat(type, "boxplots (fixed axes) saved as", pdf_file_path, "\n")
  
  invisible(pdf_file_path)
}

########### ~~~~ EMD Metric ~~~~ ########### 

### Creates a molten dataframe of marker expressions by sample_id
# df: A dataframe containing marker columns and metadata (batch and control columns)
# sample_col_name: The column name in df that contains control/sample identifiers
# batch_col_name: The column name in df that contains batch identifiers
# markers: A character vector of marker column names 
# Returns a long-format dataframe with columns 'sample_id', 'marker', 'expression', 'batch', and 'control'
reshape_df <- function(df,
                       sample_col_name,
                       batch_col_name,
                       markers) {
  
  # Construct sample_id by concatenating batch and control/sample columns
  df[, "sample_id"] <- paste0(df[, batch_col_name], "_", df[, sample_col_name])
  
  # Melt the marker expression columns into long format
  molten_df <- melt(
    data.frame(df[, markers], sample_id = df[, "sample_id"]),
    id.vars = "sample_id", variable.name = "marker",
    value.name = "expression"
  )
  
  # Match back to df to append batch and control columns
  mm <- match(molten_df$sample_id, df$sample_id)
  molten_df[, "batch"] <- df[, batch_col_name][mm]
  molten_df[, "control"] <- df[, sample_col_name][mm]
  
  # Return dataframe
  return(molten_df)
}

### Subsamples a fixed number of rows per marker across selected sample_ids
# cdf: A long-format dataframe (e.g., from reshape_df) containing 'marker' and 'sample_id' columns
# samps: A character vector of sample_id values 
# markers: A character vector of marker names
# num: Number of rows to sample per marker per sample_id (default = 20000)
# seed: Integer seed for reproducibility of random sampling
# Returns a dataframe concatenating subsampled rows across all specified markers
subsample_df <- function(cdf,
                         samps,
                         markers,
                         num = 20000,
                         seed) {
  
  # Set seed for reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Initialize output dataframe
  subdf <- data.frame()
  
  # Loop over markers and perform subsampling using helper 'subsetting'
  for (i in 1:length(markers)) {
    tmp <- cdf[(cdf$marker == markers[i]), ]
    tmpdf <- subsetting(tmp, samps, num)
    subdf <- rbind(subdf, tmpdf)
  }
  
  # Return dataframe
  return(subdf)
}

### Computes Earth Mover's Distance (EMD) matrices per control across batches
# df: A long-format dataframe with columns 'expression', 'control', and 'batch'
# binSize: Numeric bin width 
# control_list: A character vector of control sample IDs 
# batch_list: A character vector of batch IDs
# seed: Integer seed for reproducibility of random sampling
# Returns a nested list of EMD matrices 
EMD_calc <- function(df,
                     binSize,
                     control_list,
                     batch_list,
                     seed) {
  
  # Select upper and lower limits for binning 
  lim1 <- (min(df[, "expression"]) - binSize)
  min_lim <- round((min(df[, "expression"]) - binSize))
  if (lim1 < min_lim) {
    min_lim <- round((lim1 - 1))
  }
  
  lim2 <- (max(df[, "expression"]) + binSize)
  max_lim <- round((max(df[, "expression"]) + binSize))
  if (lim2 > max_lim) {
    max_lim <- round((lim2 + 1))
  }
  
  # Create list of binned distribution per control and batch
  distr <- list()
  
  for (tcs in control_list) {
    distr[[tcs]] <- list()
    for (bat in batch_list) {
      distr[[tcs]][[bat]] <- df %>%
        dplyr::filter(control == tcs, batch == bat) %>%
        dplyr::select(expression) %>%
        apply(2, function(x) {
          
          # Bin the data into densities 
          bins <- seq(min_lim, max_lim, by = binSize)
          if (length(x) == 0) {
            rep(0, times = length(bins) - 1)
          } else {
            graphics::hist(x, breaks = bins, plot = FALSE)$density
          }
        })
    }
  }
  
  # Compute pairwise EMD distances between batch distributions within each control
  distances <- list()
  
  for (tcs in control_list) {
    distances[[tcs]] <- list()
    distances[[tcs]] <- matrix(
      NA, nrow = length(batch_list), ncol = length(batch_list),
      dimnames = list(batch_list, batch_list)
    )
    for (i in seq_along(batch_list)[-length(batch_list)]) {
      batch1 <- batch_list[i]
      A <- matrix(distr[[tcs]][[batch1]])
      for (j in seq(i + 1, length(batch_list))) {
        batch2 <- batch_list[j]
        B <- matrix(distr[[tcs]][[batch2]])
        distances[[tcs]][batch1, batch2] <- emdist::emd2d(A, B)
      }
    }
  }
  
  # Return nested list of EMD matrices per control
  return(distances)
}

### Builds an EMD matrix list for every marker across control samples and batches
# df: A dataframe containing columns 'marker', 'expression', 'control', 'batch', and 'sample_id'
# samps: A character vector of sample_id values 
# markers: A character vector of marker names 
# batch_list: A character vector of batch IDs 
# control_list: A character vector of control sample IDs
# num: Number of rows to sample per marker per sample_id (default = 20000)
# seed: Integer seed for reproducibility of random sampling
# Returns EMD matrices
EMD_list <- function(df,
                     samps,
                     markers,
                     batch_list,
                     control_list,
                     num = 20000,
                     seed) {
  
  # Subsample the dataframe across markers and sample_ids for efficiency
  subdf <- data.frame(subsample_df(cdf = df, samps = samps, markers = markers, num = 20000, seed = seed))
  
  # Reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Bin size for density histograms used in EMD calculation
  binSize <- 0.1
  
  # Style for progress bar text
  colrSet <- crayon::make_style("#539DDD")
  total_iter <- length(markers)
  
  # Create progress bar
  pBar <- progress::progress_bar$new(
    format = paste0(colrSet("Calculating pairwise EMDs for each marker and control"),
                    " [", colrSet(":bar"), "] ",
                    colrSet(":percent"), " | Marker: :current/:total ",
                    "| Elapsed: :elapsed | ETA: :eta"),
    total = total_iter,
    clear = FALSE,
    width = 80,
    complete = "=",
    incomplete = "-"
  )
  
  # Output list, EMD matrices per marker and control
  emds_list <- list()
  
  # Iterate over markers and compute EMD matrices per control
  for (m in markers) {
    
    emds_list[[m]] <- list()
    tmp <- subdf[(subdf$marker == m), ]
    emds_list[[m]] <- EMD_calc(df = tmp, binSize = binSize,
                               control_list = control_list, batch_list = batch_list, seed = seed)
    
    # Update progress bar
    pBar$tick()
  }
  
  # Return nested list of EMD matrices
  return(emds_list)
}

### Computes median EMD values per batch for a given marker and control
# batch_list: A character vector of batch IDs t
# distf: An EMD matrix 
# control_list: A control sample ID (character)
# marker: The marker name 
# Returns a dataframe with columns 'batch', 'value' (median EMD), 'control', and 'marker'
median_vals <- function(batch_list,
                        distf,
                        control_list,
                        marker) {
  
  # EMD matrix to long format 
  # Remove NAs
  df <- melt(distf, na.rm = TRUE)
  
  # Compute medians across both Var1 and Var2 
  
  # Median by Var1 (row dimension)
  med_Var1 <- data.frame(df %>% arrange(match(Var1, batch_list)))
  med_Var1 <- med_Var1[, -2]
  colnames(med_Var1) <- c("batch", "value")
  
  # Median by Var2 (column dimension)
  med_Var2 <- data.frame(df %>% arrange(match(Var2, batch_list)))
  med_Var2 <- med_Var2[, -1]
  colnames(med_Var2) <- c("batch", "value")
  
  # Combine and compute median per batch
  med <- data.frame(rbind(med_Var1, med_Var2))
  med <- med %>% group_by(batch) %>% dplyr::summarise(value = median(value))
  colnames(med) <- c("batch", "value")
  
  # Annotate with control and marker
  med$control <- control_list
  med$marker  <- marker
  
  # Return median values per batch
  return(med)
}

### Produces a molten dataframe of EMD values including per-batch medians
# batch_list: A character vector of batch IDs
# distf: An EMD matrix 
# control_list: The control sample ID (character) 
# Returns a long-format dataframe with original EMD entries, plus a 'median' row per batch
dist_median_emd <- function(batch_list,
                            distf,
                            control_list) {
  
  # EMD matrix to long format 
  # Remove NAs
  df <- melt(distf, na.rm = TRUE)
  df$reference <- "no"
  
  # Compute medians across both Var1 and Var2
  med_Var1 <- data.frame(df %>% dplyr::arrange(match(Var1, batch_list)))
  med_Var1 <- med_Var1[, -2]
  colnames(med_Var1) <- c("batch", "value")
  
  med_Var2 <- data.frame(df %>% dplyr::arrange(match(Var2, batch_list)))
  med_Var2 <- med_Var2[, -1]
  colnames(med_Var2) <- c("batch", "value")
  
  med <- data.frame(rbind(med_Var1, med_Var2))
  med <- med %>% group_by(batch) %>% dplyr::summarise(value = median(value))
  colnames(med) <- c("Var2", "value")
  
  # Create a row per batch marking the median
  med <- cbind(Var1 = batch_list, med)
  med$Var2 <- "median"
  med$reference <- "yes"
  
  # Combine EMD entries with per-batch medians 
  # Annotate with control
  df.plt <- rbind(df, med)
  df.plt <- cbind(df.plt, control = control_list)
  df.plt$Var2 <- factor(df.plt$Var2, levels = c(batch_list, "median"))
  
  # Return combined dataframe
  return(df.plt)
}

### Builds a flags dataframe of median EMD per batch across markers and controls, and saves to CSV
# em_dists: A nested list of EMD matrices 
# markers: A character vector of marker names 
# samps: A character vector of sample_id values (unused here; kept for interface consistency)
# control_list: A character vector of control sample IDs
# batch_list: A character vector of batch IDs 
# threshold: Numeric threshold for flagging batches showing potential batch effect
# output_dir: Directory where the CSV will be saved
# num: Number of rows to sample (unused here; kept for interface consistency)
# Returns a dataframe of median EMD per batch/control/marker with flag and color columns
emd_med_vals <- function(em_dists,
                         markers,
                         samps,
                         control_list,
                         batch_list,
                         threshold,
                         output_dir,
                         num = 20000) {
  
  # Ensure numeric threshold
  threshold <- as.numeric(threshold)
  
  # Initialize output dataframe
  emdf_all <- data.frame(batch = NA, value = NA, control = NA, marker = NA)
  
  # Iterate over markers and controls to collect medians
  for (m in 1:length(markers)) {
    mrk <- markers[m]
    for (j in 1:length(control_list)) {
      ct <- control_list[j]
      dist_list <- em_dists[[mrk]][[ct]]
      tmp_med <- data.frame(median_vals(batch_list = batch_list, distf = dist_list, 
                                        control_list = ct, marker = markers[m]))
      emdf_all <- rbind(emdf_all, tmp_med)
    }
  }
  
  # Remove initial placeholder row
  emdf_all <- emdf_all[-1, ]
  
  # Create flag and color columns based on threshold
  emdf_all <- emdf_all %>% 
    dplyr::mutate(flagged_batch = ifelse(value >= threshold, paste0(batch), "none"))
  emdf_all <- emdf_all %>% 
    dplyr::mutate(batch_effect = ifelse(value >= threshold, as.numeric(1), as.numeric(0)))
  emdf_all <- emdf_all %>% 
    dplyr::mutate(cols = ifelse(value >= threshold, "#8E000C", "#AEAEAE"))
  
  # Order rows by control, batch, and marker lists
  emdf_all <- emdf_all %>% arrange(match(control, control_list))
  emdf_all <- emdf_all %>% arrange(match(batch, batch_list))
  emdf_all <- emdf_all %>% arrange(match(marker, markers))
  
  # Annotate type
  emdf_all$type <- "EMD"
  
  # Build output CSV path and filename
  filename <- paste0("EMD_flags.csv")
  file_path <- file.path(output_dir, filename)
  
  # Reformat dataframe for CSV file
  emdf_output <- emdf_all %>%
    dplyr::select(batch, control, marker, value, batch_effect, flagged_batch) %>%
    dplyr::rename(
      EMD = value,
      above_threshold = batch_effect
    ) %>%
    dplyr::mutate(
      above_threshold = as.logical(as.integer(above_threshold))
    )
  
  # Save flags dataframe robustly to CSV
  tryCatch({
    readr::write_csv(emdf_output, file = file_path)
    cat("EMD table saved as", file_path)
  }, error = function(e) {
    stop("Error saving table as CSV: ", e$message)
  })
  
  # Return the flags dataframe
  return(emdf_all)
}

### Generates a heatmap of pairwise Earth Mover’s Distance (EMD) values for a given marker and control
# emds_list: A nested list of EMD results indexed by [[marker]][[control]]
# marker: The marker name to display 
# control: The control sample ID to display 
# batch_list: A character vector of batch IDs to order axes and compute the plotting dataframe
# axis_size: Numeric size for axis label text (default: 18)
# threshold: Numeric threshold used to define the color scale upper band (default: 5)
EMD_hmap <- function(emds_list,
                     marker,
                     control,
                     batch_list,
                     axis_size = 18,
                     threshold = 5) {
  
  # Extract the per-control distance list for the selected marker
  dist_list <- emds_list[[marker]][[control]]
  
  # Build plotting dataframe (pairwise EMDs by batch vs batch)
  df.plt <- dist_median_emd(batch_list, dist_list, control)
  rownames(df.plt) <- NULL
  
  # Sequence of fill colors (based on min/max EMD values)
  col_seq <- seq(round(min(df.plt$value)), round(max(df.plt$value) + 1))
  
  # Define color mapping and legend breaks/labels based on threshold
  if (threshold == 5) {
    col_fun_1 <- colorRamp2(c(0, 1, 2.5, threshold), c("grey", "grey75", "bisque", "firebrick3"))
    colour_heat <- col_fun_1(col_seq)
    mid <- 2.5
    lim <- as.numeric(threshold + 1)
    breaks <- c(0, 0.5, seq(1, 6))
    labls <- c("0.0", "", "1.0", "2.0", "3.0", "4.0", "5.0", "> 5")
  } else {
    col_fun_1 <- colorRamp2(c(0, 1, 2.5, threshold), c("grey", "grey75", "bisque", "firebrick3"))
    colour_heat <- col_fun_1(col_seq)
    mid <- as.numeric((threshold) / 2)
    lim <- as.numeric(threshold + 1)
    breaks <- c(0, 0.5, seq(1, lim))
    labls <- c("0.0", "", as.character(seq(1, threshold)), paste0("> ", threshold))
  }
  
  # Choose palette/scale depending on whether max exceeds threshold
  if (max(df.plt$value) >= threshold) {
    hm01 <- ggplot(data = df.plt, aes(x = Var1, y = Var2, color = reference, fill = value)) +
      geom_tile(colour = "black") +
      scale_fill_gradientn(colours = colour_heat, guide = "colourbar",
                           breaks = breaks, labels = labls) +
      labs(x = "", y = "", fill = "Earth Mover's \nDistance") +
      geom_text(aes(x = Var1, y = Var2, label = round(value, 2)), color = "black",
                fontface = "bold", size = 5.3) +
      scale_x_discrete(position = "top") +
      ggtitle(paste0("Pairwise EMDs between batches for ", marker, "\n", " in sample ", control, "\n")) +
      theme(axis.text.x = element_text(face = "bold", size = axis_size, angle = 70, vjust = 1, hjust = 0),
            axis.text.y = element_text(face = "bold", size = axis_size),
            plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
            legend.title = element_text(face = "bold", size = 19),
            legend.text = element_text(size = axis_size - 7),
            legend.key.height = unit(1.2, 'cm'))
  } else {
    
    # Legend breaks differ when the max value is below the threshold
    if (lim < 2) {
      breaks <- seq(0, lim, by = 0.0001)
    } else {
      breaks <- seq(0, lim)
    }
    hm01 <- ggplot(data = df.plt, aes(x = Var1, y = Var2, color = reference, fill = value)) +
      geom_tile(colour = "black") +
      scale_fill_gradient2(high = "firebrick3", mid = "bisque",
                           low = "grey", midpoint = mid,
                           guide = "colourbar", breaks = breaks) +
      labs(x = "", y = "", fill = "Earth Mover's \nDistance") +
      geom_text(aes(x = Var1, y = Var2, label = round(value, 2)), color = "black",
                fontface = "bold", size = 5.3) +
      scale_x_discrete(position = "top") +
      ggtitle(paste0("Pairwise EMDs between batches for ", marker, "\n", " in sample ", control, "\n")) +      
      theme(axis.text.x = element_text(face = "bold", size = axis_size, 
                                       angle = 70, vjust = 1, hjust = 0),
            axis.text.y = element_text(face = "bold", size = axis_size),
            plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
            legend.title = element_text(face = "bold", size = 19),
            legend.text = element_text(size = axis_size - 7),
            legend.key.height = unit(1.2, 'cm'))
  }
  return(hm01)
}

### Generates a PDF of EMD heatmaps across markers and controls
# emds_list: A nested list of EMD results indexed by [[marker]][[control]]
# batch_list: A character vector of batch IDs to order axes and compute the plotting dataframe
# markers: A character vector of markers 
# controls: A character vector of control sample IDs
# threshold: Numeric threshold used for the EMD color scale 
# output_dir: Directory where the PDF will be saved
# axis_size: Numeric font size for axis labels in the heatmaps 
# width, height: Numeric values specifying the dimensions of the output plot in inches
generate_EMD_heatmaps <- function(emds_list,
                                  batch_list,
                                  markers,
                                  controls,
                                  threshold,
                                  output_dir,
                                  axis_size = 18,
                                  width = 11.5,
                                  height = 8) {
  # Define the PDF file path
  file_name <- "EMD_heatmaps.pdf"
  pdf_file_path <- file.path(output_dir, file_name)
  
  # Progress bar setup
  colrSet <- crayon::make_style("#539DDD")
  total_iter <- length(markers)
  
  # Create progress bar
  pBar <- progress::progress_bar$new(
    format = paste0(colrSet("Plotting marker expression densities"),
                    " [", colrSet(":bar"), "] ",
                    colrSet(":percent"), " | Marker: :current/:total ",
                    "| Elapsed: :elapsed | ETA: :eta\n"),
    total = total_iter,
    clear = FALSE,
    width = 80,
    complete = "=",
    incomplete = "-"
  )
  
  # Open the PDF device
  pdf(pdf_file_path, width = width, height = height, useDingbats = FALSE)
  
  # Loop through markers and controls; generate and print each heatmap
  for (i in 1:length(markers)) {
    mk <- as.character(markers[i])
    for (c in 1:length(controls)) {
      ct <- controls[c]
      
      # Console status
      cat("\nProcessing marker:", mk, "control:", ct, "\n")
      
      # Build the heatmap for the given marker-control pair
      hmap <- EMD_hmap(emds_list = emds_list, marker = mk, control = ct,
                       batch_list = batch_list, axis_size = axis_size,
                       threshold = threshold)
      print(hmap)
    }
    
    # Update progress bar
    pBar$tick(tokens = list(current = i))
  }
  
  # Close the PDF device
  dev.off()
  
  # Notify user of completion
  cat("Pairwise EMD heatmaps saved as", pdf_file_path, "\n")
}

### Generates an ordered pairwise EMD boxplot for a given marker and control
# distdf: A nested list of EMD matrices indexed by [[marker]][[control]]
# batch_colours: A named vector of colours for batches (names must match batch_list)
# control: A character string specifying the control sample ID
# marker: A character string specifying the marker name
# batch_list: A character vector of batch IDs used for ordering and mapping shapes/colours
# threshold: Numeric threshold for flagging batches with high median EMDs (default = 5)
# axis_size: Numeric font size for axis labels and titles in the boxplot
EMD_boxplot <- function(distdf,
                        batch_colours,
                        control,
                        marker,
                        batch_list,
                        threshold = 5) {
  
  # Print progress message for current marker/control
  cat("\nProcessing marker:", marker, "and control:", control, "\n")
  
  # Extract EMD matrix for the given marker and control
  EMD_matrix <- distdf[[marker]][[control]]
  lower_triangle <- t(EMD_matrix)
  
  # Build a symmetric matrix by combining upper and lower triangles
  emd_square_matrix <- matrix(NA, nrow = nrow(EMD_matrix), ncol = ncol(EMD_matrix))
  emd_square_matrix[upper.tri(emd_square_matrix)] <- EMD_matrix[upper.tri(EMD_matrix)]
  emd_square_matrix[lower.tri(emd_square_matrix)] <- lower_triangle[lower.tri(lower_triangle)]
  rownames(emd_square_matrix) <- rownames(EMD_matrix)
  colnames(emd_square_matrix) <- colnames(EMD_matrix)
  
  # Convert to long format for plotting
  emd_square_matrix <- melt(emd_square_matrix)
  df1 <- emd_square_matrix %>% drop_na()
  colnames(df1) <- c("pair", "ref", "value")
  
  # Compute median EMD per batch and order batches by decreasing median
  df1_med <- data.frame(value = df1$value, batch = df1$ref) %>%
    group_by(batch) %>% dplyr::summarise_all(list(median))
  df1_med <- df1_med[order(df1_med$value, decreasing = TRUE), ]
  df1_med$batch <- as.character(df1_med$batch)
  
  # Reorder df1 to match median order
  df1 <- df1 %>% arrange(match(ref, df1_med$batch))
  df1$ref <- factor(df1$ref, levels = unique(df1_med$batch))
  df1$pair <- factor(as.character(df1$pair), levels = unique(df1$pair))
  
  # Legend order follows user inputted batch_list
  df1$batch <- factor(df1$ref, levels = batch_list)           
  
  # X-axis order follows median order 
  df1$ref_ordered <- factor(df1$ref, levels = unique(df1_med$batch))
  
  # Assign colors based on threshold
  df1_med <- data.frame(df1_med) %>%
    dplyr::mutate(cols = ifelse(value > threshold, "#FF2316", "#999999"))
  mm <- match(df1$ref, df1_med$batch)
  df1$cols <- df1_med$cols[mm]
  
  # Identify flagged batches (median > threshold)
  out_batch <- unique(as.character(df1[which(df1$cols == "#FF2316"), "ref"]))
  out_batch <- out_batch[order(match(out_batch, batch_list))]
  
  # Define shape vector for flagged batches
  shape_vect <- c(rep(c(22,23,24,25,18,8,9,11,4,13,3,14,12,10,1,5,2,7,0,6,16,15,19,17), 4))
  
  # Default shape assignment
  df1 <- cbind(df1, shapes = 21)
  
  # Assign distinct shapes to flagged batches
  if (length(out_batch) != 0) {
    for (b in 1:length(out_batch)) {
      for (i in 1:nrow(df1)) {
        if (df1[i, "pair"] == out_batch[b]) {
          df1[i, "shapes"] <- shape_vect[b]
        }
      }
    }
  }
  
  # Ensure batch_colours are named
  names(batch_colours) <- batch_list
  
  # Dynamic y-axis limit if max value >= 20
  if (max(df1$value) >= 20) {
    max_lim <- round(max(df1$value))
    if (max(df1$value) > max_lim) {
      max_lim <- max_lim + 1
    }
    
    ggplot(df1, aes(x = ref_ordered, y = value)) +
      geom_boxplot(aes(color = cols), alpha = 1.0, fill = "#F2F2F2", outlier.color = NA) +
      ggplot2::scale_color_identity(guide = "none") +   # suppress threshold color legend
      geom_point(aes(fill = batch, shape = shapes), colour = "#000000", size = 1.8,
                 alpha = 1.0, position = position_jitter(width = 0.4, height = 0)) +
      ggplot2::scale_fill_manual(
        name   = "Batch",
        values = batch_colours,
        breaks = batch_list,  
        drop   = FALSE        
      ) +
      ggplot2::scale_shape_identity(guide = "none") +
      ggplot2::geom_hline(yintercept = threshold, colour = "#FF2316", linewidth = 0.9) +
      ggplot2::labs(x = "", y = "Earth Mover's Distance",
                    title = paste0("Ordered pairwise EMDs per batch for ", marker,
                                   "\n", " in sample ", control, "\n")) +
      ggplot2::theme_bw() +
      ggplot2::scale_y_continuous(limits = c(0, max_lim)) +
      ggplot2::theme(
        axis.text.x   = element_text(size = 26, angle = 70, vjust = 0.9, hjust = 0.9),
        axis.text.y   = element_text(size = 25, vjust = 0.75),
        axis.title.y  = element_text(size = 28),
        plot.title    = element_text(size = 24, face = "bold", hjust = 0.5),
        legend.title  = element_text(size = 19, face = "bold"),
        legend.text   = element_text(size = 17),
        legend.key.size = unit(3, "mm")
      ) +
      guides(fill = guide_legend(override.aes = list(shape = 22, size = 5, colour = "black")))
    
  } else {
    
    ggplot(df1, aes(x = ref_ordered, y = value)) +
      geom_boxplot(aes(color = cols), alpha = 1.0, fill = "#F2F2F2", outlier.color = NA) +
      ggplot2::scale_color_identity(guide = "none") +   # suppress threshold color legend
      geom_point(aes(fill = batch, shape = shapes), colour = "#000000", size = 1.8,
                 alpha = 1.0, position = position_jitter(width = 0.4, height = 0)) +
      ggplot2::scale_fill_manual(
        name   = "Batch",
        values = batch_colours,
        breaks = batch_list,  
        drop   = FALSE         
      ) +
      ggplot2::scale_shape_identity(guide = "none") +
      ggplot2::geom_hline(yintercept = threshold, colour = "#FF2316", linewidth = 0.9) +
      ggplot2::labs(x = "", y = "Earth Mover's Distance",
                    title = paste0("Ordered pairwise EMDs per batch for ", marker,
                                   "\n", " in sample ", control, "\n")) +
      ggplot2::theme_bw() +
      ggplot2::scale_y_continuous(limits = c(0, 20)) +
      ggplot2::theme(
        axis.text.x   = element_text(size = 26, angle = 70, vjust = 0.9, hjust = 0.9),
        axis.text.y   = element_text(size = 25, vjust = 0.75),
        axis.title.y  = element_text(size = 28),
        plot.title    = element_text(size = 24, face = "bold", hjust = 0.5),
        legend.title  = element_text(size = 19, face = "bold"),
        legend.text   = element_text(size = 17),
        legend.key.size = unit(3, "mm")
      ) +
      guides(fill = guide_legend(override.aes = list(shape = 22, size = 5, colour = "black")))
  }
}

### Generates a PDF of ordered pairwise EMD boxplots across markers and controls
# distdf: A nested list of pairwise EMD matrices indexed by [[marker]][[control]]
# batch_list: A character vector of batch IDs used for ordering and legend mapping
# markers: A character vector of markers to iterate over
# control_list: A character vector of control sample IDs to iterate over
# batch_colours: A named vector of fill colours for batches (names must match batch_list)
# threshold: Numeric threshold used to flag batches with high median EMDs (passed to EMD_boxplot)
# output_dir: Directory where the PDF will be saved
# axis_size: Numeric font size for axis labels in the boxplots
# width, height: Numeric values specifying the dimensions of the output plot in inches
generate_EMD_boxplots <- function(distdf,
                                  batch_list,
                                  markers,
                                  control_list,
                                  batch_colours,
                                  threshold,
                                  output_dir,
                                  width = 11.5,
                                  height = 8) {
  
  # Define the PDF file path and name
  file_name <- "EMD_boxplots.pdf"
  pdf_file_path <- file.path(output_dir, file_name)
  
  # Create a color style for progress bar text
  colrSet <- crayon::make_style("#539DDD")
  total_iter <- length(markers)
  
  # Create progress bar 
  pBar <- progress::progress_bar$new(
    format = paste0(
      colrSet("Plotting marker expression densities"),
      " [", colrSet(":bar"), "] ",
      colrSet(":percent"), " | Marker: :current/:total ",
      "| Elapsed: :elapsed | ETA: :eta"
    ),
    total = total_iter,
    clear = FALSE,
    width = 80,
    complete = "=",
    incomplete = "-"
  )
  
  # Ensure batch_colours are named consistently by batch_list
  names(batch_colours) <- batch_list
  
  # Open the PDF device (useDingbats = FALSE to avoid font embedding issues)
  pdf(pdf_file_path, width = width, height = height, useDingbats = FALSE)
  
  # Loop through the markers and controls, generate boxplots, and print to the PDF device
  for (i in 1:length(markers)) {
    mk <- markers[i]
    for (c in 1:length(control_list)) {
      ct <- control_list[c]
      box_plot <- EMD_boxplot(
        distdf = distdf,
        batch_colours = batch_colours, 
        marker = mk,
        control = ct,
        batch_list = batch_list, 
        threshold = threshold
      )
      print(box_plot)
    }
    
    # Update progress bar after finishing all controls for this marker
    pBar$tick(tokens = list(current = i))
  }
  
  # Close the PDF device (writes the file to disk)
  dev.off()
  
  # Notify user of completion and location of the saved PDF
  cat("EMD boxplots saved as", pdf_file_path, "\n")
}

########### ~~~~ Summary of Flags (IQR and EMD) ~~~~ ########### 

### Automatically determines axis text size based on the number of markers/controls
# axis: Either 'x' or 'y' indicating which axis to size
# num_markers: Numeric count of markers (used when axis == 'x')
# num_controls: Numeric count of controls (used when axis == 'y')
# Returns a numeric axis text size appropriate for the plot layout
auto_axis_size <- function(axis,
                           num_markers,
                           num_controls) {
  
  # Validate axis input
  if (!axis %in% c('x','y')) {
    stop("axis must be specified - either 'x' or 'y'")
  }
  
  # Validate numeric inputs
  if (!is.numeric(num_markers) || !is.numeric(num_controls)) {
    stop("num_markers (number of markers) and num_controls (number of controls) must be numeric")
  }
  
  # Default base text size for axes
  base_size <- 18
  
  # Compute size for x-axis based on number of markers
  if (axis == 'x') {
    set <- 4 * (num_markers)
    if (set < 80) {
      axis_size <- as.numeric(base_size) * 2
    }
    else if (set < 100 && set >= 80) {
      axis_size <- as.numeric(base_size)
    }
    else {
      size_factor <- round(min(1, 100 / set), 1)
      axis_size <- round(base_size * size_factor)
    }
  }
  
  # Compute size for y-axis based on number of controls
  else {
    min_size <- as.numeric(num_controls)
    if (min_size <= 4 && min_size > 2) {
      axis_size <- as.numeric(base_size) + 1
    }
    else if (min_size <= 2) {
      axis_size <- as.numeric(base_size + (base_size / 4))
    }
    
    else {
      # Decrease size linearly as the number of controls increases beyond 4
      size_factor <- as.numeric(subtract((num_controls), 4))
      if (size_factor < base_size) { 
        axis_size <- as.numeric(subtract(base_size, size_factor))
      }
      else if (size_factor == base_size) {
        axis_size <- as.numeric(base_size / 8)
      }
      else if (size_factor > base_size && size_factor <= base_size + 2) {
        axis_size <- as.numeric(base_size / 8)
      }
      else {
        axis_size <- as.numeric(1)
      }
    }
  }
  
  # Ensure numeric output
  axis_size <- as.numeric(axis_size)
  return(axis_size)
}

### Generates a ranked heatmap of flagged batch effects across markers and controls
# emd_df: A dataframe of EMD-based flags (columns must include 'batch', 'value', 'control', 'marker', 'flagged_batch', 'batch_effect', 'cols', and 'type')
# iqr_df: A dataframe of IQR-based flags (same required columns as emd_df)
# batch_list: A character vector of batch IDs 
# controls: A character vector of control sample IDs 
# markers: A character vector of marker names 
# batch_cols: A named vector of colors for batches (names must match batch_list)
# output_dir: Optional directory where a CSV summary will be saved (default = NULL)
# min_batch_effect: Minimum number of flags per batch required to keep the batch in the heatmap (default = 0)
# min_marker_effect: Minimum number of flags per marker required to keep the marker in the heatmap (default = 0)
ranked_flagged_hmap <- function(emd_df,
                                iqr_df,
                                batch_list,
                                controls,
                                markers,
                                batch_cols,
                                output_dir = NULL,
                                min_batch_effect = 0,
                                min_marker_effect = 0) {
  
  # Conditions used to order panels/columns in the heatmap
  conditions <- c("-MFI", "+MFI", "%pos", "EMD")
  
  # Select only the necessary columns from IQR and EMD flag dataframes
  iqr_df <- subset(iqr_df, select = c(batch, value, control, marker, flagged_batch,
                                      batch_effect, cols, type))
  emd_df <- subset(emd_df, select = c(batch, value, control, marker, flagged_batch,
                                      batch_effect, cols, type))
  
  # Create a combined variable for (type_marker) to rank conditions per marker
  iqr_df$var2 <- paste0(iqr_df$type, "_", iqr_df$marker)
  emd_df$var2 <- paste0(emd_df$type, "_", emd_df$marker)
  
  # Combine into a single dataframe for the ranked heatmap
  hmdf <- data.frame(rbind(iqr_df, emd_df))
  
  # Re-arrange by list of conditions, batches, and markers
  # This will order the dataframe by each condition for each marker
  hmdf <- hmdf %>% arrange(match(type, conditions))
  hmdf <- hmdf %>% arrange(match(batch, batch_list))
  hmdf <- hmdf %>% arrange(match(marker, markers))
  
  # Sort by decreasing batch_effect 
  hmdf <- hmdf[order(hmdf$batch_effect, decreasing = TRUE), ]
  hmdf <- hmdf %>% arrange(match(batch, batch_list))
  hmdf <- hmdf %>% arrange(match(marker, markers))
  
  # Force order of controls for the plot 
  hmdf$control <- factor(hmdf$control, levels = rev(controls))
  
  ## Get the frequency of flags for all marker-condition combinations
  freq_var2 <- data.frame(table(hmdf$var2, hmdf$batch_effect))
  colnames(freq_var2) <- c("var2", "flag", "Freq")
  
  # Order by decreasing order of flags per combo
  freq_var2 <- freq_var2[order(freq_var2$flag, freq_var2$Freq, decreasing = TRUE), ]
  freq_var2$var2 <- as.character(freq_var2$var2)
  
  # Reorder hmap dataframe to match order of marker-condition flags
  hmdf <- hmdf %>% arrange(match(var2, freq_var2$var2))
  
  # Get the frequency of flags for all batches
  batch_freq <- data.frame(table(hmdf$batch, hmdf$batch_effect))
  colnames(batch_freq) <- c("batch", "flag", "Freq")
  
  # Order by frequency of flags per batch
  batch_freq <- batch_freq[order(batch_freq$flag, batch_freq$Freq, decreasing = TRUE), ]
  batch_freq$batch <- as.character(batch_freq$batch)
  
  # Reorder hmap dataframe to match order of flags of batches
  hmdf <- hmdf %>% arrange(match(batch, batch_freq$batch))
  
  # Get the frequency of flags for all panel markers
  marker_freq <- data.frame(table(hmdf$marker, hmdf$batch_effect))
  colnames(marker_freq) <- c("marker", "flag", "Freq")
  
  # Order by frequency of flags per marker
  marker_freq <- marker_freq[order(marker_freq$flag, marker_freq$Freq, decreasing = TRUE), ]
  marker_freq$marker <- as.character(marker_freq$marker)
  
  # Re-order hmap dataframe by conditions, flagged batch frequency, and flagged marker frequency
  hmdf <- hmdf %>% arrange(match(batch, batch_freq$batch))
  hmdf <- hmdf %>% arrange(match(type, conditions))
  hmdf <- hmdf %>% arrange(match(marker, marker_freq$marker))
  hmdf <- hmdf %>% arrange(match(batch, batch_freq$batch))
  
  # Re-order hmap dataframe by decreasing order of flags again for final stability
  hmdf <- hmdf[order(hmdf$batch_effect, decreasing = TRUE), ]
  hmdf <- hmdf %>% arrange(match(var2, freq_var2$var2))
  hmdf <- hmdf %>% arrange(match(batch, batch_freq$batch))
  hmdf <- hmdf %>% arrange(match(marker, marker_freq$marker))
  hmdf <- hmdf[order(hmdf$batch_effect, decreasing = TRUE), ]
  
  # Get a vector of markers with most to least flags
  mord <- unique(marker_freq$marker)
  
  # Make factors
  hmdf$marker <- factor(hmdf$marker, levels = mord)
  hmdf$batch <- factor(hmdf$batch, levels = unique(batch_freq$batch))
  
  # Force hmap plot order of conditions flagged for each marker
  var_match <- data.frame(var2 = unique(hmdf$var2))
  var_match$type <- sapply(var_match$var2, function(c) gsub("\\_\\D+\\w+", "", c))
  var_match$marker <- sapply(var_match$var2, function(c) gsub("\\W\\w+\\_", "", c))
  var_match$marker <- sapply(var_match$marker, function(c) gsub("\\w+\\_", "", c))
  
  # Order by list of conditions and order of flagged markers
  var_match <- var_match %>% arrange(match(type, conditions))
  var_match <- var_match %>% arrange(match(marker, mord))
  
  # Re-order hmap dataframe by order of marker-condition combos and decreasing order of flagged batches
  hmdf <- hmdf %>% arrange(match(var2, unique(var_match$var2)))
  hmdf <- hmdf %>% arrange(match(batch, unique(batch_freq$batch)))
  
  # Make factors
  hmdf$var2 <- factor(hmdf$var2, levels = unique(var_match$var2))
  hmdf$type <- factor(hmdf$type, levels = conditions)
  
  # Reorder named batch colors to order of frequency of batches flagged
  # This is for the strip-l function used for the heatmap
  names(batch_cols) <- batch_list
  ord_bcols <- batch_cols
  ord_bcols <- ord_bcols[order(factor(names(ord_bcols), 
                                      levels = unique(batch_freq$batch)))]
  
  # Reformat dataframe for CSV file
  hmdf_output <- hmdf %>%
    dplyr::select(batch, control, marker, type, value, batch_effect, flagged_batch) %>%
    dplyr::rename(
      above_threshold = batch_effect
    ) %>%
    dplyr::mutate(
      above_threshold = as.logical(as.integer(above_threshold)),
      type = dplyr::recode(type,
                           "+MFI"   = "Positive MFI",
                           "-MFI"   = "Negative MFI",
                           "%pos"   = "Percent positive")
    ) %>%
    dplyr::arrange(as.character(marker), as.character(batch), as.character(control), type)
  
  
  # Optionally write summary CSV to output_dir
  if (!is.null(output_dir)) {
    csv_file <- "summary_IQR_EMD_flags.csv"
    csv_file_path <- file.path(output_dir, csv_file)
    readr::write_csv(hmdf_output, file = csv_file_path)
    cat("Summary of IQR and EMD flags table saved as", csv_file_path, "\n")
  }
  
  # Get the size of x and y axis using auto_axis_size helper
  xlab_size <- auto_axis_size(axis = 'x', num_markers = length(markers),
                              num_controls = length(controls))
  
  ylab_size <- auto_axis_size(axis = 'y', num_markers = length(markers),
                              num_controls = length(controls))
  
  # Count number of flags per batch to optionally filter batches shown
  be_count <- hmdf %>%
    dplyr::group_by(batch) %>%
    dplyr::summarise(n_flag = sum(batch_effect == 1, na.rm = TRUE),
                     .groups = "drop")
  
  # min_batch_effect: Minimum number of flags per batch-control combo required to keep the batch-control combo in the heatmap (default = 0)
  if (min_batch_effect > 0) {
    
    bc_flag_count <- hmdf %>%
      dplyr::group_by(batch, control) %>%
      dplyr::summarise(
        n_flags = sum(batch_effect == TRUE, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Keep only batch–control combos meeting the threshold
    hmdf <- hmdf %>%
      dplyr::semi_join(
        dplyr::filter(bc_flag_count, n_flags >= min_batch_effect),
        by = c("batch", "control")
      )
    
    # Drop unused factor levels for batch & control
    hmdf$batch   <- factor(hmdf$batch)
    hmdf$control <- factor(hmdf$control)
  }
  
  # min_marker_effect: Minimum number of flags per marker required to keep the marker combo in the heatmap (default = 0)
  if (min_marker_effect > 0) {
    
    bc_flag_count <- hmdf %>%
      dplyr::group_by(marker) %>%
      dplyr::summarise(
        n_flags = sum(batch_effect == TRUE, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Keep only markers meeting the threshold
    hmdf <- hmdf %>%
      dplyr::semi_join(
        dplyr::filter(bc_flag_count, n_flags >= min_marker_effect),
        by = c("marker")
      )
    
    # Drop unused factor levels for markers
    hmdf$marker   <- factor(hmdf$marker)
  }
  
  # Build the heatmap plot of flags
  hmap_plot <- ggplot(hmdf, aes(x = type, y = control, fill = cols)) +
    geom_tile(colour = "black", lwd = 0.6) +
    scale_fill_identity() +
    scale_y_discrete(position = "right", expand = c(0, 0)) +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    facet_grid(batch ~ marker, scales = "free", space = "free_x", switch = "y") +
    xlab("") +
    ylab("") +
    ggtitle("Flags for IQR and EMR metrics, per batch, per sample") +
    theme(plot.title = element_text(size = xlab_size + 10, face = "bold", hjust = 0.5),
          strip.text.y = element_text(size = 21, face = "bold", colour = "black"),
          strip.text.x = element_text(size = xlab_size - 2, face = "bold"),
          strip.background.x = element_blank(),
          strip.placement = "outside",
          axis.text.x = element_text(size = xlab_size + 2.5, angle = 90,
                                     hjust = -0.3), # , hjust = -0.3, vjust = 3
          axis.text.y = element_text(size = ylab_size + 2, hjust = 0.7, vjust = 0.4),
          panel.spacing.x = unit(0.4, 'lines'),
          panel.spacing.y = unit(0.5, 'lines'),
          panel.grid = element_blank())
  
  # Assigning strip colors based on batch number
  grid.plt <- ggplot_gtable(ggplot_build(hmap_plot))
  stripr <- which(grepl('strip-l', grid.plt$layout$name))
  
  fills <- batch_cols[as.character(unique(hmdf$batch))]
  
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', grid.plt$grobs[[i]]$grobs[[1]]$childrenOrder))
    grid.plt$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k + 1
  }
  
  # Plot ranked heatmap of flags 
  grid::grid.draw(grid.plt)
}

### Computes PNG image width (in pixels) based on number of markers and conditions
# num_markers: Numeric count of markers in the panel
# Returns numeric width in pixels
image_width <- function(num_markers) {
  
  # Number of conditions plotted per marker ( -MFI, +MFI, %pos, EMD )
  num_conditions <- 4
  
  # Width grows linearly 
  size <- as.numeric((num_markers * num_conditions * 100) + 7500)
  
  # Return width
  return(size)
}

### Computes PNG image height (in pixels) based on number of batches and controls
# num_batch: Numeric count of batches 
# num_controls: Numeric count of controls 
# Returns numeric height in pixels
image_height <- function(num_batch,
                         num_controls) {
  
  # Height grows linearly
  size <- as.numeric((num_controls * num_batch * 100) + 3000)
  
  # Return height
  return(size)
}

### Generates PNG ranked flagged heatmaps
# emd_df: A dataframe of EMD-based flags (columns must include 'batch', 'value', 'control', 'marker', 'flagged_batch', 'batch_effect', 'cols', and 'type')
# iqr_df: A dataframe of IQR-based flags (same required columns as emd_df)
# batch_list: A character vector of batch IDs 
# controls: A character vector of control sample IDs 
# markers: A character vector of marker names 
# batch_cols: A named vector of colours for batches (names must match batch_list)
# output_dir: Directory where PNG images will be saved
# markers_per_plot: Number of markers to include per image panel (default = all markers)
# batches_per_plot: Number of batches to include per image panel (default = all batches)
# min_batch_effect: Minimum number of flags per batch-control combination required to keep the batch in the heatmap (default = 0)
# min_marker_effect: Minimum number of flags per marker required to keep the marker in the heatmap (default = 0)
generate_ranked_flagged_hmap <- function(emd_df,
                                         iqr_df,
                                         batch_list,
                                         controls,
                                         markers,
                                         batch_cols,
                                         output_dir,
                                         markers_per_plot = NULL,
                                         batches_per_plot = NULL,
                                         min_batch_effect = 0,
                                         min_marker_effect = 0) {
  
  # Defaults: use everything in one image (backwards compatible)
  if (is.null(markers_per_plot)) markers_per_plot  <- length(markers)
  if (is.null(batches_per_plot)) batches_per_plot  <- length(batch_list)
  
  # Cache sizes
  nm  <- length(markers)
  nb  <- length(batch_list)
  nct <- length(controls)
  
  # Basic safety checks
  if (markers_per_plot <= 0L) stop("markers_per_plot must be >= 1")
  if (batches_per_plot <= 0L) stop("batches_per_plot must be >= 1")
  
  # Keep only needed columns once 
  iqr_df <- subset(
    iqr_df,
    select = c(batch, value, control, marker, flagged_batch,
               batch_effect, cols, type)
  )
  
  # Initialize image counter
  img_counter <- 0L
  
  # Loop over marker chunks
  for (m_start in seq(1L, nm, by = markers_per_plot)) {
    m_end          <- min(m_start + markers_per_plot - 1L, nm)
    markers_chunk  <- markers[m_start:m_end]
    
    # Loop over batch chunks
    for (b_start in seq(1L, nb, by = batches_per_plot)) {
      b_end        <- min(b_start + batches_per_plot - 1L, nb)
      batch_chunk  <- batch_list[b_start:b_end]
      
      img_counter <- img_counter + 1L
      
      # Size depends on how many markers/batches in this panel
      width  <- image_width(length(markers_chunk))
      height <- image_height(length(batch_chunk), nct)
      
      # Optional: subset iqr_df (and emd_df) to speed things up for each panel
      iqr_df_chunk <- subset(
        iqr_df,
        marker %in% markers_chunk & batch %in% batch_chunk
      )
      
      emd_df_chunk <- subset(
        emd_df,
        marker %in% markers_chunk & batch %in% batch_chunk
      )
      
      # Decide whether everything fits in a single plot after any filtering
      all_markers_fit <- markers_per_plot >= length(markers)      # markers = post-filter list
      all_batches_fit <- batches_per_plot >= length(batch_list)   # batch_list = post-filter list
      
      has_marker_filter <- min_marker_effect > 0
      has_batch_filter  <- min_batch_effect  > 0
      
      # Prefix to reflect filter state in the filename
      filter_prefix <- if (has_marker_filter && has_batch_filter) {
        "filtered_markers_filtered_batches_"
      } else if (has_marker_filter) {
        "filtered_markers_"
      } else if (has_batch_filter) {
        "filtered_batches_"
      } else {
        ""
      }
      
      if (all_markers_fit && all_batches_fit) {
        # No split
        if (filter_prefix == "") {
          filename <- "summary_heatmap_IQR_EMD_flags_all_markers_all_batches.png"
        } else {
          # Filtered but still a single plot: avoid "all_markers_all_batches" wording
          filename <- paste0("summary_heatmap_IQR_EMD_flags_", filter_prefix, "single_plot.png")
        }
      } else {
        # Split
        if (all_markers_fit && !all_batches_fit) {
          filename <- sprintf(
            "summary_heatmap_IQR_EMD_flags_%sall_markers_batches%02d-%02d.png",
            filter_prefix, b_start, b_end
          )
        } else if (!all_markers_fit && all_batches_fit) {
          filename <- sprintf(
            "summary_heatmap_IQR_EMD_flags_%smarkers%02d-%02d_all_batches.png",
            filter_prefix, m_start, m_end
          )
        } else {
          filename <- sprintf(
            "summary_heatmap_IQR_EMD_flags_%smarkers%02d-%02d_batches%02d-%02d.png",
            filter_prefix, m_start, m_end, b_start, b_end
          )
        }
      }
      
      png_file_path <- file.path(output_dir, filename)
      
      # Open PNG device 
      png(filename = png_file_path, width = width, height = height, res = 400)
      
      tryCatch({
        
        # Draw ranked flagged heatmap for the current chunk
        ranked_flagged_hmap(
          emd_df          = emd_df_chunk,
          iqr_df          = iqr_df_chunk,
          batch_list      = batch_chunk,
          controls        = controls,
          markers         = markers_chunk,
          batch_cols      = batch_cols,
          output_dir      = output_dir,
          min_batch_effect = min_batch_effect,
          min_marker_effect = min_marker_effect
        )
      }, error = function(e) {
        message("Error in generating summary heatmap (image ", img_counter, "): ",
                e$message)
      }, finally = {
        
        # Close
        dev.off()
      })
      
      # Console message
      if (nm > markers_per_plot || nb > batches_per_plot) {
        
        # Multiple heatmaps (splitting occurred)
        cat("Summary of IQR and EMD flags heatmap, part", img_counter,
            "saved as", png_file_path, "\n")
      } else {
        
        # Only one heatmap (no split)
        cat("Summary of IQR and EMD flags heatmap saved as", png_file_path, "\n")
      }
    }
  }
}

########### ~~~~ Clustering Metric ~~~~ ########### 

### Uses FlowSOM unsupervised clustering to get put each cell into a cluster based on marker expression values
# data_f: A flowFrame or flowSet object containing the cytometry data
# markers: A character vector of marker names to use for clustering
# seed: Numeric value used for reproducibility in clustering steps
# meta_cluster_num: Total number of meta-clusters to generate
# samps: A vector of sample IDs to include in the analysis
# num: Numeric value specifying the number of cells to sample for clustering
# output_dir: Directory where clustering results and plots will be saved
getFlowSOM_clusters <- function(data_f, 
                                markers, 
                                seed, 
                                meta_cluster_num, 
                                samps, 
                                num, 
                                output_dir) {
  
  # Create a directory to store all the consensus plots for metaclusters
  plot_outdir <- file.path(paste0(output_dir,"/cluster_plots_k",
                                  meta_cluster_num,"_seed",seed))
  
  dir.create(plot_outdir, recursive = TRUE, showWarnings = FALSE)
  
  # Create directory and check if successful
  if (!dir.create(plot_outdir, recursive = TRUE, showWarnings = FALSE)) {
    if (!dir.exists(plot_outdir)) {
      stop("Cannot create directory: ", plot_outdir)
    }
  }
  
  # Ensure we have write permissions
  if (file.access(plot_outdir, mode = 2) != 0) {
    stop("No write permission in directory: ", plot_outdir)
  }
  
  old_dir <- getwd()
  
  tryCatch({
    
    # Subsample dataframe for clustering
    if (!is.null(seed)) {
      set.seed(seed)
    }
    
    df <- data.frame(subsetting(df = data_f, samps = samps, num = num))
    
    # Convert transformed dataframe into a flowFrame object
    ctrl_ff <- flowCore::flowFrame(data.matrix(df[,markers]))
    
    # Create FlowSOM input 
    fsomIN <- FlowSOM::ReadInput(ctrl_ff, transform = FALSE, scale = FALSE)
    
    # Build SOM clusters 
    som <- FlowSOM::BuildSOM(fsomIN)
    
    # SOM clusters information
    cell_clustering_som<-som$map$mapping[,1]
    codes <- som$map$codes
    
    # Remove the fluorochrome/metal tag names from the prettyColnames
    som[["prettyColnames"]] <- sapply(som[["prettyColnames"]],
                                      function(x)gsub("\\s\\W\\w+\\W","",x))
    
    # Call the ConsensusClusterPlus function to get meta-clusters
    mc <- ConsensusClusterPlus(t(codes), maxK = meta_cluster_num, reps = 100,
                               pItem = 0.9, pFeature = 1, title = plot_outdir, 
                               plot = "png", clusterAlg = "hc", innerLinkage = "average", 
                               finalLinkage = "average",
                               distance = "euclidean", seed = seed)
    
    # Get cluster ids and labels from the metacluster object
    code_clustering_1 <- mc[[meta_cluster_num]]$consensusClass
    cell_clustering_1 <- code_clustering_1[cell_clustering_som]
    
    df <- cbind(df, cluster = as.character(cell_clustering_1))
    
    filename <- paste0("clusters", meta_cluster_num,"_seed",seed,".csv")
    file_path <- file.path(output_dir, filename)
    
    tryCatch({
      readr::write_csv(df, file = file_path)
      cat("ConsensusClusterPlus plots saved in", plot_outdir)
      cat("\nCluster data saved as", file_path)
    }, error = function(e) {
      stop("Error saving table as CSV: ", e$message)
    })
    return(df)
  }, error = function(e) {
    stop("Error in cluster processing: ", e$message)
  }, finally = {
    setwd(old_dir)
  })
}

### Generates a bar plot of cluster sizes, ordered from largest to smallest cluster by percentage of cells
# cluster_df: A data frame containing at a cluster column (called 'cluster') for each cell and marker expression columns
# seed: Numeric value used for reproducibility in upstream steps (displayed in the plot title, included in the output file name)
# meta_cluster_num: Total number of clusters (also used in the plot title)
# output_dir: Directory where the PNG file will be saved
# width, height: Numeric values specifying the dimensions of the output plot in inches
clusterSizes_plot <- function(cluster_df, 
                              seed, 
                              meta_cluster_num, 
                              output_dir, 
                              width = 11.5, 
                              height = 8) {
  
  minclust <- 0.005 
  csizes  <- data.frame(table(cluster_df[,"cluster"]))
  
  # Percent of cells in each cluster
  csizes <- cbind(csizes, percent = csizes[,"Freq"]/nrow(cluster_df))
  csizes$xlab <- paste0(csizes[,"Var1"], " (", csizes[,"Freq"], ")")
  csizes[,"xlab"] <- factor(csizes[,"xlab"],
                            levels = as.character(csizes[rev(order(csizes[,"percent"])),"xlab"]))
  csizes[,"Var1"] <- factor(csizes[,"Var1"],
                            levels = as.character(csizes[rev(order(csizes[,"percent"])),"Var1"]))
  
  # Plot
  p <- ggplot(data = csizes, aes(x = xlab, y = percent)) +
    geom_bar(stat = "identity", aes(color = xlab), show.legend = FALSE) +
    scale_y_continuous(labels = scales::percent) +
    ggtitle(paste0("Size for ", meta_cluster_num,
                   " metaclusters (seed ", seed, ")")) +
    geom_hline(yintercept = c(minclust), colour = "#F25A25") +
    xlab("Cluster (# of cells)") + ylab("% of all analyzed cells") +
    theme(
      axis.text.x = element_text(size = 20, angle = 45, vjust = 1.0, hjust = 1.0),
      axis.text.y = element_text(size = 20),
      axis.title = element_text(size = 30, colour = "#000000"),
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
    )
  
  # File name and path
  filename <- paste0("number_cells_clusters", meta_cluster_num, "_seed", seed, ".png")
  
  png_file_path <- file.path(output_dir, filename)
  
  # Save as PNG
  tryCatch({
    ggsave(filename = png_file_path, plot = p, width = width, height = height, dpi = 300, device = "png")
    cat("Cluster size bar plot saved as", png_file_path, "\n")
  }, error = function(e) {
    stop("Error saving plot:", e$message)
  })
}

### Generates a heatmap of scaled median marker expression for each cluster, ordered by cluster number
# cluster_df: A data frame containing a 'cluster' column for each cell, a 'batch' column, and marker expression columns
# markers: A character vector of marker names to include in the heatmap
# cluster_colrs: A named vector of colors for clusters
# cell_cluster_vector: A vector assigning each cell to a cluster
# meta_cluster_num: Total number of clusters (used in the plot title)
# seed: Numeric value used for reproducibility in upstream steps (displayed in the plot title, included in the output file name)
# axis_size: Numeric value specifying font size for axis labels
# num_size: Numeric value specifying font size for numbers displayed in the heatmap cells
# output_dir: Directory where the PNG file will be saved
# width, height: Numeric values specifying the dimensions of the output plot in inches
clustering_heatmap <- function(cluster_df, 
                               markers, 
                               cluster_colrs, 
                               cell_cluster_vector, 
                               meta_cluster_num, 
                               seed, 
                               axis_size = 15, 
                               num_size = 8.4, 
                               output_dir, 
                               width = 11.5, 
                               height = 8) {
  
  # Initial input checks
  if (!is.data.frame(cluster_df)) {
    stop("cluster_df must be a data frame")
  }
  if (!all(c("cluster", "batch") %in% colnames(cluster_df))) {
    stop("cluster_df must contain 'cluster' and 'batch' columns")
  }
  if (!dir.exists(output_dir)) {
    stop("Output directory does not exist")
  }
  
  # Expression values normalized to 0-1 scale
  expr_t<-as.matrix(cluster_df[,markers])
  rng1 <- colQuantiles(expr_t, probs = c(0.01, 0.99), useNames = TRUE)
  expr01t <- t((t(expr_t) - rng1[, 1]) / (rng1[, 2] - rng1[, 1]))
  expr01t[expr01t < 0] <- 0
  expr01t[expr01t > 1] <- 1
  sample_ids<-cluster_df[,"sample_id"]
  controls<-cluster_df[,"control"]
  batches<-cluster_df[,"batch"]
  
  # Calculate the median expression and scaled median expressions
  expr_median <- data.frame(expr_t, cell_clustering = cell_cluster_vector) %>%
    group_by(cell_clustering) %>% summarize_all(list(median))
  expr01_median <- data.frame(expr01t, cell_clustering = cell_cluster_vector) %>%
    group_by(cell_clustering) %>% summarize_all(list(median))
  
  # Calculate cluster frequencies and proportions
  clustering_table <- as.numeric(table(cell_cluster_vector))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  
  # Sort the clusters with hierarchical clustering
  d <- dist(expr_median[, colnames(expr_t)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_median[, colnames(expr01t)])
  rownames(expr_heat) <- expr01_median$cell_clustering
  
  # Colors for the heatmap
  color_heat <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)
  
  # Color breaks for the legend
  legend_breaks = seq(from = 0, to = 1.2, by = 0.2)
  
  # Convert cluster names to numeric for sorting
  expr01_median$cell_clustering <- as.numeric(as.character(expr01_median$cell_clustering))
  
  # Sort rows by numeric cluster order
  expr01_median <- expr01_median[order(expr01_median$cell_clustering), ]
  expr_heat <- expr_heat[order(as.numeric(rownames(expr_heat))), ]
  
  # Update labels after sorting
  labels_row <- paste0("Cluster ", expr01_median$cell_clustering, " (", clustering_prop[order(expr01_median$cell_clustering)], "%)")
  
  # Annotation for the original clusters
  annotation_row <- data.frame(Cluster = factor(expr01_median$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  clustColors1 <- cluster_colrs[1:nlevels(annotation_row$Cluster)]
  names(clustColors1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = clustColors1)
  
  # Plot
  plt <- ComplexHeatmap::pheatmap(expr_heat, name = "Scaled Median\nExpression",
                                  color = color_heat, cluster_cols = FALSE,
                                  cluster_rows = FALSE, labels_row = labels_row,
                                  display_numbers = TRUE, number_color = "black",
                                  fontsize = axis_size, fontsize_number = num_size, 
                                  column_names_side = c("top"),
                                  legend_breaks = legend_breaks, 
                                  legend_labels = c("0","0.2","0.4","0.6","0.8","1.0",""),
                                  annotation_legend = FALSE,
                                  main = paste("Scaled median marker expressions for", meta_cluster_num, "metaclusters"))
  
  # File name and path
  filename <- paste0("heatmap_clusters", meta_cluster_num, "_seed", seed, ".png")
  
  png_file_path <- file.path(output_dir, filename)
  
  # Save as PNG
  png(png_file_path, width = width, height = height, units = "in", res = 300)
  
  print(plt)
  dev.off()
  
  # Notify user of completion
  cat("Cluster heatmap saved as", png_file_path, "\n")
}

### Generates a bar plot of batch proportions for each cluster, outlining in bold batches that exceed a specified threshold
# cluster_df: A data frame containing a 'cluster' column for each cell, a 'batch' column, and marker expression columns
# cell_cluster_vector: A vector assigning each cell to a cluster
# set_threshold1: The threshold, specifically a numeric multiplier for the expected proportion threshold (e.g., 2 for 2x expected)
# batch_list: A character vector of batch names to include in the plot
# batch_colours: A named vector of colors corresponding to batches
# seed: Numeric value used for reproducibility in upstream steps (displayed in the plot title, included in the output file name)
# output_dir: Directory where the PNG file will be saved
# width, height: Numeric values specifying the dimensions of the output plot in inches
# axis_size: Numeric value specifying font size for axis labels
proportion_of_batches_per_cluster <- function(cluster_df, 
                                              cell_cluster_vector,
                                              set_threshold1 = 2,
                                              batch_list, 
                                              batch_colours, 
                                              seed, 
                                              output_dir, 
                                              width=11.5, 
                                              height=8, 
                                              axis_size=14) {
  
  # Inital input checks
  if (!is.data.frame(cluster_df)) {
    stop("cluster_df must be a data frame")
  }
  if (!all(c("cluster", "batch") %in% colnames(cluster_df))) {
    stop("cluster_df must contain 'cluster' and 'batch'")
  }
  if (!dir.exists(output_dir)) {
    stop("Output directory does not exist")
  }
  if (!all(batch_list %in% unique(cluster_df$batch))) { 
    stop("Some batches in batch_list are not present")
  }
  
  # Number of clusters
  meta_cluster_num <- length(unique(cluster_df$cluster))
  
  # Cluster size
  csizes <- data.frame(table(cluster_df$cluster))
  csizes <- csizes %>% mutate(percent = Freq / sum(Freq))
  csizes$xlab <- paste0(csizes$Var1, " (", csizes$Freq, ")")
  csizes$xlab <- factor(csizes$xlab, levels = rev(csizes$xlab))
  
  # Compute batch proportions
  counts_table <- table(cell_cluster_vector, cluster_df$batch)
  prop_df <- data.frame(counts_table)
  colnames(prop_df) <- c("cluster", "batch", "Freq")
  
  mm <- match(prop_df$cluster, csizes$Var1)
  prop_df$clust_freq <- csizes$Freq[mm]
  prop_df$clust_props <- csizes$percent[mm]
  prop_df <- prop_df %>% mutate(proportion = (Freq / clust_freq) * 100)
  
  # Sort clusters numerically
  prop_df$cluster <- as.numeric(as.character(prop_df$cluster))
  prop_df <- prop_df[rev(order(prop_df$cluster)), ]
  
  # Labels for axis
  prop_df$xlab_new <- round(prop_df$clust_props * 100, 2)
  prop_df$xlab_new <- paste0("Cluster ", prop_df$cluster, " (", prop_df$xlab_new, "%)")
  prop_df$xlab_new <- factor(prop_df$xlab_new, levels = unique(prop_df$xlab_new))
  
  # Get threshold
  prop_threshold <- set_threshold1 * round(100 / length(batch_list))
  
  # Highlight column for bars above threshold 
  prop_df$highlight <- prop_df$proportion > prop_threshold
  
  # Assign colors 
  names(batch_colours) <- batch_list
  
  # Plot
  barplot <- ggplot(prop_df, aes(x = proportion, y = xlab_new, fill = batch)) +
    geom_bar(
      stat = "identity",
      width = 0.78,                      
      color = ifelse(prop_df$highlight, "black", NA),  
      size = 1.5                           
    ) +
    scale_fill_manual("Batch", values = batch_colours) +
    xlab("Batch Proportions") +
    ylab("") +
    ggtitle("Batch proportions per cluster") +
    scale_x_continuous(position = "bottom", expand = expansion(mult = c(0.003, 0.02))) +
    scale_y_discrete(
      position = "right",
      expand = expansion(add = 0.5) 
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = axis_size + 2, hjust = 0.2),
      axis.text.y = element_text(size = axis_size + 3),
      axis.title = element_text(size = axis_size + 3),
      legend.title = element_text(size = axis_size + 1),
      legend.text = element_text(size = axis_size),
      plot.title = element_text(size = axis_size + 3, face = "bold", hjust = 0.5)
    )
  
  # Print summary
  cat("Threshold used (anything over is flagged): ", prop_threshold, "%\n")
  
  if (sum(prop_df$highlight) == 0) {
    cat(paste0("No batches with ", set_threshold1, "x higher than the expected proportion found.\n"))
    cat("---------------------------------------------------\n")
  } else {
    cat("Batches exceeding threshold:\n")
    batch_props_filt_1 <- prop_df[prop_df$highlight, ]
    print(batch_props_filt_1[, c("cluster", "batch", "proportion")])
    cat("---------------------------------------------------\n")
  }
  
  # File name and path
  filename <- paste0("proportions_barplot_threshold", set_threshold1, "_clusters",
                     meta_cluster_num, "_seed", seed, ".png")
  
  png_file_path <- file.path(output_dir, filename)
  
  # Save as PNG
  tryCatch({
    ggsave(filename = png_file_path, plot = barplot, width = width, height = height, units = "in", dpi = 350)
    # Notify user of completion
    cat("Bar plot saved as", png_file_path, "\n")
  }, error = function(e) {
    stop("Error generating plot: ", e$message)
  })
  
  return(invisible(NULL))
}

### Generates a table of batch proportions for each cluster, with a column indicating batches that exceed a specified threshold
# cluster_df: A data frame containing a 'cluster' column for each cell, a 'batch' column, and marker expression columns
# cell_cluster_vector: A vector assigning each cell to a cluster
# set_threshold1: The threshold, specifically a numeric multiplier for the expected proportion threshold (e.g., 2 for 2x expected)
# batch_list: A character vector of batch names to include in the analysis
# seed: Numeric value used for reproducibility in upstream steps (included in the output file name)
# output_dir: Directory where the CSV file will be saved
proportion_of_batches_per_cluster_table <- function(cluster_df, 
                                                    cell_cluster_vector,
                                                    set_threshold1 = 2,
                                                    batch_list, 
                                                    seed, 
                                                    output_dir) {
  # Initial input checks
  if (!is.data.frame(cluster_df)) {
    stop("cluster_df must be a data frame")
  }
  if (!all(c("cluster", "batch") %in% colnames(cluster_df))) {
    stop("cluster_df must contain 'cluster' and 'batch' columns")
  }
  if (!dir.exists(output_dir)) {
    stop("Output directory does not exist")
  }
  if (!all(batch_list %in% unique(cluster_df$batch))) {
    stop("Some batches in batch_list are not present in the data")
  }
  
  # Number of clusters
  meta_cluster_num <- length(unique(cluster_df$cluster))
  
  # Cluster sizes
  csizes <- data.frame(table(cluster_df[,"cluster"]))
  csizes <- cbind(csizes, percent = csizes[,"Freq"]/nrow(cluster_df))
  
  # Counts and proportions
  counts_table <- table(cell_cluster_vector, cluster_df[,"batch"])
  prop_df <- data.frame(counts_table)
  colnames(prop_df) <- c("cluster","batch","Freq")
  
  # Order and add cluster info
  prop_df <- prop_df %>% arrange(match(cluster, levels(csizes$Var1)))
  mm <- match(prop_df$cluster, csizes$Var1)
  prop_df$clust_freq <- csizes$Freq[mm]
  prop_df$clust_props <- csizes$percent[mm]
  prop_df <- prop_df %>% mutate(proportion = (Freq/clust_freq)*100)
  
  # Compute threshold and add flag column
  prop_threshold_1 <- set_threshold1 * round((100 / length(batch_list)))
  prop_df <- prop_df %>%
    mutate(above_threshold = proportion > prop_threshold_1)
  
  # Order by cluster number
  prop_df$cluster <- as.numeric(as.character(prop_df$cluster)) 
  prop_df <- prop_df[order(prop_df[, 1]), ]
  
  # Print summary
  cat("Threshold used (anything over is flagged): ", prop_threshold_1, "%\n")
  
  if (sum(prop_df$above_threshold) == 0) {
    cat(paste0("No batches with ", set_threshold1, "x higher than the expected proportion found.\n"))
    at("---------------------------------------------------\n")
  } else {
    cat("Batches exceeding threshold:\n")
    # Filter rows above threshold
    batch_props_filt_1 <- prop_df[prop_df$above_threshold, ]
    print(batch_props_filt_1[, c("cluster", "batch", "proportion")])
    cat("---------------------------------------------------\n")
  }
  
  # File name and path
  filename <- paste0("proportions_table_threshold", set_threshold1, "_clusters",
                     meta_cluster_num, "_seed", seed, ".csv")
  
  file_path <- file.path(output_dir, filename)
  
  # Save as CSV
  tryCatch({
    readr::write_csv(prop_df[, c("cluster","batch","proportion","above_threshold")], 
                     file = file_path)
    cat("All proportions saved as", file_path, "\n")
  }, error = function(e) {
    stop("Error saving table as CSV: ", e$message)
  })
  
  return(prop_df[, c("cluster","batch","proportion","above_threshold")])
}

### Generates boxplots of batch proportions across control samples for clusters impacted by batch effects, highlighting batches that exceed a specified threshold
# cluster_df: A data frame containing 'cluster', 'sample_id', 'control', and 'batch' columns for each cell
# cell_cluster_vector: A vector assigning each cell to a cluster
# batch_props_df: A data frame from the previous step indicating clusters flagged for batch effects
# set_threshold2: The threshold, specifically a numeric multiplier for the expected proportion threshold (e.g., 1.5 for 1.5x expected)
# batch_list: A character vector of batch names to include in the plot
# control_list: A character vector of control sample names to include in the plot
# batch_colours: A named vector of colors corresponding to batches
# output_dir: Directory where the PNG file will be saved
# width, height: Numeric values specifying the dimensions of the output plot in inches
# axis_size: Numeric value specifying font size for axis labels
# seed: Numeric value used for reproducibility in upstream steps (included in the output file name)
proportion_of_batches_across_control_samples <- function(cluster_df, 
                                                         cell_cluster_vector,
                                                         batch_props_df, 
                                                         set_threshold2 = 1.5,
                                                         batch_list, 
                                                         control_list, 
                                                         batch_colours, 
                                                         output_dir, 
                                                         width = 11.5, 
                                                         height = 8, 
                                                         axis_size = 12, 
                                                         seed) {
  
  # Inital input checks
  if (!is.data.frame(cluster_df)) {
    stop("cluster_df must be a data frame")
  }
  if (!all(c("cluster", "sample_id") %in% colnames(cluster_df))) {
    stop("cluster_df must contain 'cluster' and 'sample_id' columns")
  }
  if (!dir.exists(output_dir)) {
    stop("Output directory does not exist")
  }
  
  # Get number of clusters
  meta_cluster_num <- length(unique(cluster_df$cluster))
  
  # File name and path
  filename <- paste0(
    "proportions_samples_boxplot_threshold", set_threshold2,
    "_clusters", meta_cluster_num,
    "_seed", seed, ".png"
  )
  
  png_file_path <- file.path(output_dir, filename)
  
  csizes <- data.frame(table(cluster_df[, "cluster"]))
  csizes <- cbind(csizes, percent = csizes[, "Freq"] / nrow(cluster_df))
  
  csizes$Var1 <- factor(
    csizes$Var1,
    levels = as.character(sort(as.numeric(as.character(csizes$Var1))))
  )
  
  csizes$xlab <- paste0("Cluster ", csizes$Var1, " (", csizes$Freq, ")")
  csizes$xlab <- factor(
    csizes$xlab,
    levels = paste0("Cluster ", levels(csizes$Var1), " (", csizes$Freq, ")")
  )
  
  sample_ids <- cluster_df[, "sample_id"]
  
  counts_table <- table(cell_cluster_vector, sample_ids)
  prop_df <- data.frame(counts_table)
  colnames(prop_df) <- c("cluster", "sample_id", "Freq")
  
  prop_df <- prop_df %>% arrange(match(cluster, unique(csizes$Var1)))
  mm <- match(prop_df$cluster, csizes$Var1)
  prop_df$clust_freq <- csizes$Freq[mm]
  prop_df$clust_prop <- csizes$percent[mm]
  prop_df <- prop_df %>% mutate(proportion = (Freq / clust_freq) * 100)
  
  mm <- match(prop_df$sample_id, cluster_df$sample_id)
  prop_df$controls <- cluster_df$control[mm]
  prop_df$batch <- cluster_df$batch[mm]
  
  prop_df$cluster <- factor(prop_df$cluster, levels = unique(prop_df$cluster))
  
  prop_df$strip_lab <- paste0("Cluster ", prop_df$cluster)
  prop_df$strip_lab <- factor(
    prop_df$strip_lab,
    levels = paste0(
      "Cluster ",
      sort(as.numeric(as.character(unique(prop_df$cluster))))
    )
  )
  
  prop_df$batch <- factor(prop_df$batch, levels = batch_list)
  prop_df$control <- factor(prop_df$control, levels = control_list)
  
  out_clust <- unique(batch_props_df$cluster[batch_props_df$above_threshold])
  prop_df <- prop_df[prop_df$cluster %in% out_clust, ]
  
  combodf <- data.frame(prop_df)
  
  prop_thr <- set_threshold2 * round((100 / length(batch_list)))
  
  flagged <- combodf[combodf$proportion > prop_thr, ]
  
  cat("Threshold used (anything over is flagged): ", prop_thr, "%\n")
  
  if (nrow(flagged) == 0) {
    cat(paste0(
      "No controls with ", set_threshold2,
      "x higher than the expected proportion found.\n"
    ))
    cat("---------------------------------------------------\n")
  } else {
    cat("Controls exceeding threshold:\n")
    print(flagged[, c("cluster", "controls", "batch", "proportion")])
    cat("---------------------------------------------------\n")
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  combodf$jitter_x <- as.numeric(combodf$control) +
    runif(nrow(combodf), -0.15, 0.15)
  
  flagged$jitter_x <- combodf$jitter_x[
    match(
      paste(flagged$cluster, flagged$controls, flagged$batch),
      paste(combodf$cluster, combodf$controls, combodf$batch)
    )
  ]
  
  y_range <- max(combodf$proportion, na.rm = TRUE) -
    min(combodf$proportion, na.rm = TRUE)
  nudge_amount <- 0.05 * y_range
  
  ctrl_labels <- setNames(as.character(seq_along(control_list)), control_list)
  
  # Factor
  combodf$batch <- factor(combodf$batch, levels = names(batch_colours) %||% batch_list)
  flagged$batch <- factor(flagged$batch, levels = levels(combodf$batch))

  # Plot
  plot <- ggplot(combodf, aes(control, proportion)) +
    geom_boxplot(
      outlier.color = NA,
      colour = "#000000",
      fill = "#F2F2F2",
      width = 0.6
    ) +
    geom_point(
      aes(x = jitter_x, fill = batch, colour = batch),  # map BOTH
      size = 3.8,
      pch = 21,
      stroke = 0.4
    ) +
    geom_label_repel(
      data = flagged,
      aes(
        x = jitter_x,
        y = proportion,
        label = batch,
        colour = batch
      ),
      fill = "white",
      size = 6.8,
      show.legend = FALSE,
      force = 0.1,
      force_pull = 2,
      max.time = Inf,
      direction = "y",
      nudge_y = nudge_amount,
      segment.size = 0.4
    ) +
    # Single scale applied to BOTH colour and fill:
    scale_colour_manual(
      name = "Batch",
      values = batch_colours,
      breaks = levels(combodf$batch),
      aesthetics = c("colour", "fill")
    ) +
    facet_wrap(~ strip_lab, nrow = round(length(levels(prop_df$cluster)) / 2)) +
    scale_x_discrete(labels = ctrl_labels) +
    xlab("Controls") +
    ylab("Proportion (%)") +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 23),
      axis.text.y = element_text(size = 25),
      axis.title = element_text(size = 25),
      strip.text.x = element_text(size = 21, face = "bold"),
      legend.title = element_text(size = 15, face = "bold"),
      legend.text = element_text(size = 12)
    )
  
  tryCatch({
    ggsave(
      filename = png_file_path,
      plot = plot,
      width = width,
      height = height,
      units = "in",
      dpi = 400
    )
    cat("Boxplots saved as", png_file_path, "\n")
  }, error = function(e) {
    stop("Error generating plot: ", e$message)
  })
}

### Generates a table of batch proportions across control samples for clusters that were over threshold 1, with a column indicating batches that exceed a specified threshold 2
# cluster_df: A data frame containing 'cluster', 'sample_id', 'control', and 'batch' columns for each cell
# cell_cluster_vector: A vector assigning each cell to a cluster
# batch_props_df: A data frame from the previous step indicating clusters flagged for batch effects
# set_threshold2: he threshold, specifically a numeric multiplier for the expected proportion threshold (e.g., 1.5 for 1.5x expected)
# batch_list: A character vector of batch names included in the analysis
# control_list: A character vector of control sample names included in the analysis
# output_dir: Directory where the CSV file will be saved
# seed: Numeric value used for reproducibility in upstream steps (included in the output file name)
proportion_of_batches_across_control_samples_table <- function(cluster_df,
                                                               cell_cluster_vector,
                                                               batch_props_df, 
                                                               set_threshold2 = 1.5,
                                                               batch_list, 
                                                               control_list, 
                                                               output_dir, 
                                                               seed) {
  # Initial input checks
  if (!is.data.frame(cluster_df)) {
    stop("cluster_df must be a data frame")
  }
  if (!all(c("cluster", "sample_id") %in% colnames(cluster_df))) {
    stop("cluster_df must contain 'cluster' and 'sample_id' columns")
  }
  if (!dir.exists(output_dir)) {
    stop("Output directory does not exist")
  }
  
  # Number of clusters
  meta_cluster_num <- length(unique(cluster_df$cluster))
  
  # Cluster sizes
  csizes <- data.frame(table(cluster_df$cluster))
  csizes <- cbind(csizes, percent = csizes$Freq / nrow(cluster_df))
  
  # Counts and proportions
  counts_table <- table(cell_cluster_vector, cluster_df$sample_id)
  prop_df <- data.frame(counts_table)
  colnames(prop_df) <- c("cluster", "sample_id", "Freq")
  
  # Add cluster info
  mm <- match(prop_df$cluster, csizes$Var1)
  prop_df$clust_freq <- csizes$Freq[mm]
  prop_df$clust_prop <- csizes$percent[mm]
  prop_df <- prop_df %>% mutate(proportion = (Freq / clust_freq) * 100)
  
  # Add control and batch info
  mm2 <- match(prop_df$sample_id, cluster_df$sample_id)
  prop_df$control <- cluster_df$control[mm2]
  prop_df$batch <- cluster_df$batch[mm2]
  
  # Only keep clusters flagged in batch_props_df
  out_clust <- unique(batch_props_df$cluster[batch_props_df$above_threshold])
  prop_df <- prop_df[prop_df$cluster %in% out_clust, ]
  
  # Compute threshold and flag
  prop_thr <- set_threshold2 * round(100 / length(batch_list))
  prop_df <- prop_df %>% mutate(above_threshold = proportion > prop_thr)
  
  # Print summary
  cat("Threshold used (anything over is flagged): ", prop_thr, "%\n")
  
  if (sum(prop_df$above_threshold) == 0) {
    cat(paste0("No controls with ", set_threshold2, "x higher than expected proportion found.\n"))
    cat("---------------------------------------------------\n")
  } else {
    cat("Controls exceeding threshold:\n")
    flagged <- prop_df[prop_df$above_threshold, ]
    print(flagged[, c("cluster", "control", "batch", "proportion")])
    cat("---------------------------------------------------\n")
  }
  
  # Order by cluster, then batch, then control 
  prop_df$cluster <- as.numeric(as.character(prop_df$cluster))
  prop_df <- prop_df %>% dplyr::arrange(cluster, batch, control)
  
  # File name and path
  filename <- paste0("proportions_samples_table_threshold", set_threshold2, "_clusters",
                     meta_cluster_num, "_seed", seed, ".csv")
  
  file_path <- file.path(output_dir, filename)
  
  tryCatch({
    readr::write_csv(prop_df[, c("cluster", "control", "batch", "proportion", "above_threshold")],
                     file = file_path)
    cat("All proportions saved as", file_path, "\n")
  }, error = function(e) {
    stop("Error saving table as CSV: ", e$message)
  })
  
  return(prop_df[, c("cluster", "control", "batch", "proportion", "above_threshold")])
}

### Generates a UMAP plot of subsampled cells colored by batch
# df: A data frame containing per-cell data with columns 'cluster', 'sample_id', 'batch', 'control', and marker columns
# marker_list: A character vector of marker column names 
# batch_list: A character vector of batch names
# batch_colours: A named character vector of colors for each batch in batch_list (names must match batch_list)
# seed: Numeric value for reproducibility; controls subsampling and UMAP initialization
# output_dir: Directory path where the PNG file will be saved
# axis_size: Numeric size for axis label text and titles (default: 15)
# width, height: Numeric values specifying the dimensions of the output plot in inches
generate_cluster_UMAP <- function(df,
                                  marker_list,
                                  batch_list,
                                  batch_colours,
                                  seed,
                                  num = 2000,
                                  output_dir,
                                  axis_size = 22,
                                  width = 12,
                                  height = 12) {
  
  # Inform the user that the process may be slow
  cat("This process may take a long time, please be patient", "\n")
  
  # Validate that the output directory exists
  if (!dir.exists(output_dir)) {
    stop("Output directory does not exist")
  }
  
  # Determine the number of meta-clusters present
  meta_cluster_num <- length(unique(df$cluster))
  
  # Build the output file name and full path
  file_name <- paste0("UMAP_clusters", meta_cluster_num, "_seed", seed, ".png")
  file_path <- file.path(output_dir, file_name)
  
  # Compute number of cells per cluster
  ncells <- data.frame(table(df$cluster))
  colnames(ncells) <- c("cluster", "Freq")
  
  # Attach cluster size to each row in df
  mm <- match(df$cluster, ncells$cluster)
  df$csize <- ncells$Freq[mm]
  
  # Set the seed for reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  subsamp <- c()
  
  # List of unique clusters to iterate over
  samps <- unique(df$cluster)
  
  # Subsample up to N cells per cluster (or fewer if the cluster has <N cells)
  for (fi in 1:length(samps)) {
    tmp <- which(df[, "cluster"] == samps[fi])
    nums <- as.numeric(unique(df[which(df[, "cluster"] == samps[fi]), "csize"]))
    subsamp <- c(subsamp, tmp[sample(1:length(tmp), min(nums, num))])
  }
  
  # Create the subsampled data frame
  subdf <- df[subsamp, ]
  subdf$samp_clust <- paste0(subdf$sample_id, "_", subdf$cluster)
  
  # Cache the sample-cluster identifier for later joins
  samp_clust <- subdf[, "samp_clust"]
  
  # Announce the UMAP computation
  cat("Calculating UMAP dimensions...", "\n")
  
  # Compute UMAP
  out_umap <- umap::umap(subdf[, marker_list], config = umap.defaults)
  
  # Extract and label the UMAP coordinates
  dims_umap <- out_umap$layout
  colnames(dims_umap) <- c("UMAP_1", "UMAP_2")
  
  stopifnot(nrow(dims_umap) == length(samp_clust))
  
  # Build a UMAP results dataframe with identifiers and type
  dims_umap <- cbind(as.data.frame(dims_umap), sample_id = samp_clust, type = "UMAP")
  
  # Map back cluster, batch, and control metadata to the UMAP points
  mm <- match(dims_umap$sample_id, subdf$samp_clust)
  dims_umap$cluster <- subdf$cluster[mm]
  dims_umap$batch <- subdf$batch[mm]
  dims_umap$control <- subdf$control[mm]
  
  # Set batch as a factor
  dims_umap$batch <- factor(dims_umap$batch, levels = batch_list)
  
  # Name the batch_colours vector to match batch levels
  names(batch_colours) <- batch_list
  
  # Keep only the columns needed for plotting
  umapdf <- subset(dims_umap, select = c(UMAP_1, UMAP_2, batch))
  
  # Plot generation
  cat("Generating UMAP plot...", "\n")
  
  # Plot
  umap_plot <- ggplot(umapdf, aes(x = UMAP_1, y = UMAP_2, colour = batch)) +
    geom_point(alpha = 0.7, size = 0.1) +
    scale_colour_manual("Batch", values = batch_colours) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    ggtitle("UMAP of cells used in FlowSOM, \ncolored by batch") +
    theme_bw() +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    theme(axis.text = element_text(size = axis_size),
          axis.title = element_text(size = axis_size + 4),
          legend.text = element_text(size = axis_size + 2),
          legend.title = element_text(size = axis_size + 2),
          plot.title = element_text(size = axis_size + 2, face = "bold", hjust = 0.5))
  
  # Save plot
  ggplot2::ggsave(file_path, plot = umap_plot, width = width, height = height, dpi = 400)
  
  # Notify user of completion
  cat("UMAP saved as", file_path, "\n")
}

########### ~~~~ Summary of Flags (IQR, EMD, and cluster) ~~~~ ########### 

### Generates a ranked heatmap of flagged batch effects across markers and controls
# emd_df: A dataframe of EMD-based flags (columns must include 'batch', 'value', 'control', 'marker', 'flagged_batch', 'batch_effect', 'cols', and 'type')
# iqr_df: A dataframe of IQR-based flags (same required columns as emd_df)
# cluster_df: A data frame containing 'cluster', 'sample_id', 'control', and 'batch' columns for each cell
# cell_cluster_vector: A vector assigning each cell to a cluster
# set_threshold1: Threshold 1 for flag 1, specifically a numeric multiplier for the expected proportion threshold (e.g., 1.5 for 1.5x expected)
# set_threshold2: Threshold 2 for flag 2, specifically a numeric multiplier for the expected proportion threshold (e.g., 1.5 for 1.5x expected)
# batch_list: A character vector of batch IDs 
# controls: A character vector of control sample IDs 
# markers: A character vector of marker names 
# batch_cols: A named vector of colors for batches (names must match batch_list)
# output_dir: Optional directory where a CSV summary will be saved (default = NULL)
# min_batch_effect: Minimum number of flags per batch-control combination required to keep the batch in the heatmap (default = 0)
# min_marker_effect: Minimum number of flags per marker required to keep the marker in the heatmap (default = 0)
ranked_flagged_hmap_all <- function(emd_df,
                                    iqr_df,
                                    cluster_df,
                                    cell_cluster_vector,
                                    set_threshold1 = 2,
                                    set_threshold2 = 1.5,
                                    batch_list,
                                    controls,
                                    markers,
                                    batch_cols,
                                    output_dir = NULL,
                                    min_batch_effect = 0,
                                    min_marker_effect = 0) {
  
  # Compute cluster proportions 
  
  # Compute overall cluster sizes for normalization
  csizes <- data.frame(table(cluster_df$cluster))
  csizes$percent <- csizes$Freq / nrow(cluster_df)
  
  # Counts
  counts_df <- as.data.frame(
    table(
      cluster = cell_cluster_vector,
      batch   = cluster_df$batch,
      control = cluster_df$control
    )
  )
  colnames(counts_df) <- c("cluster", "batch", "control", "Freq")
  
  csizes2 <- csizes
  colnames(csizes2) <- c("cluster", "clust_freq")
  mm <- match(counts_df$cluster, csizes2$cluster)
  counts_df$clust_freq <- csizes2$clust_freq[mm]
  
  # Ensure all combinations of batch × control × cluster exist
  all_batch_ctrl <- expand.grid(
    batch   = unique(cluster_df$batch),
    control = unique(cluster_df$control),
    cluster = unique(cluster_df$cluster)
  )
  
  # Merge counts and compute proportions 
  prop_df <- all_batch_ctrl %>%
    dplyr::left_join(counts_df, by = c("batch", "control", "cluster")) %>%
    dplyr::mutate(
      Freq = ifelse(is.na(Freq), 0L, Freq),
      control = as.character(control),
      batch   = as.character(batch),
      cluster = as.character(cluster)
    ) %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(
      cluster_total = sum(Freq, na.rm = TRUE),
      proportion = ifelse(cluster_total > 0,
                          (Freq / cluster_total) * 100,
                          0)
    ) %>%
    dplyr::ungroup()
  
  # Threshold 1 
  # Set threshold 1
  prop_thr1 <- set_threshold1 * (100 / length(batch_list))
  
  thr1_batch <- prop_df %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(
      cluster_total = sum(Freq, na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      batch_prop = ifelse(cluster_total > 0,
                          (Freq / cluster_total) * 100,
                          0)
    ) %>%
    dplyr::group_by(cluster, batch) %>%
    dplyr::summarise(
      batch_prop = sum(batch_prop),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      above_thr1 = batch_prop > prop_thr1
    ) %>%
    dplyr::group_by(batch) %>%
    dplyr::summarise(
      thr1_flag = any(above_thr1),
      .groups = "drop"
    )
  
  # Threshold 2 
  # Compute threshold 2 similarly to threshold 1
  prop_thr2 <- set_threshold2 * round(100 / length(batch_list))
  
  # Flag rows above Threshold 2
  prop_df <- prop_df %>%
    dplyr::mutate(above_thr2 = proportion > prop_thr2)
  
  # TRUE if any cluster exceeds threshold 2
  thr2_ctrl <- prop_df %>%
    dplyr::group_by(batch, control) %>%
    dplyr::summarise(
      thr2_flag_raw = any(above_thr2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # Only consider batches that passed Flag1
    dplyr::left_join(thr1_batch, by = "batch") %>%
    dplyr::mutate(
      thr2_flag = dplyr::case_when(
        !thr1_flag            ~ NA,     # batch didn't pass Flag1
        thr2_flag_raw == TRUE ~ TRUE,   # control passes Flag2
        TRUE                  ~ FALSE   # control fails Flag2
      )
    ) %>%
    dplyr::select(batch, control, thr2_flag)
  
  # Subset only required columns for IQR and EMD
  iqr_df <- subset(iqr_df,
                   select = c(batch, value, control, marker,
                              flagged_batch, batch_effect, cols, type))
  emd_df <- subset(emd_df,
                   select = c(batch, value, control, marker,
                              flagged_batch, batch_effect, cols, type))
  
  # Combine IQR and EMD into one dataframe for plotting
  hmdf <- dplyr::bind_rows(iqr_df, emd_df)
  
  # Add Threshold 1 marker 
  thr1_marker <- "Batch > threshold"
  
  thr1_hmdf <- hmdf %>%
    dplyr::distinct(batch, control) %>%
    dplyr::mutate(batch = as.character(batch),
                  control = as.character(control)) %>%
    dplyr::left_join(thr1_batch %>% mutate(batch = as.character(batch)), by = "batch") %>%
    dplyr::mutate(
      marker = thr1_marker,
      type = "Threshold",
      value = thr1_flag,
      batch_effect = thr1_flag,
      flagged_batch = as.character(thr1_flag),
      cols = case_when(
        thr1_flag ~ "#8E000C",
        TRUE      ~ "#AEAEAE"
      )
    ) %>%
    dplyr::select(names(hmdf))
  
  # Add Threshold 2 marker
  thr2_marker <- "Batch > threshold2"
  
  thr2_hmdf <- hmdf %>%
    dplyr::distinct(batch, control) %>%
    dplyr::mutate(batch = as.character(batch),
                  control = as.character(control)) %>%
    dplyr::left_join(thr2_ctrl, by = c("batch", "control")) %>%
    dplyr::mutate(
      marker = thr2_marker,
      type = "Threshold2",
      value = thr2_flag,
      batch_effect = thr2_flag,
      flagged_batch = as.character(thr2_flag),
      cols = dplyr::case_when(
        is.na(thr2_flag) ~ NA_character_,
        thr2_flag        ~ "#8E000C",
        TRUE             ~ "#AEAEAE"
      )
    ) %>%
    dplyr::select(names(hmdf))
  
  # Combine IQR/EMD with Threshold markers
  hmdf <- dplyr::bind_rows(hmdf, thr1_hmdf, thr2_hmdf)
  markers <- unique(c(markers, thr1_marker, thr2_marker))
  
  # Ordering & ranking
  conditions <- c("-MFI", "+MFI", "%pos", "EMD", "Threshold", "Threshold2")
  
  # Sort by type, batch, and marker
  hmdf <- hmdf %>%
    dplyr::arrange(match(type, conditions)) %>%
    dplyr::arrange(match(batch, batch_list)) %>%
    dplyr::arrange(match(marker, markers))
  
  # Prioritize flagged batches at top
  hmdf <- hmdf %>%
    dplyr::arrange(
      dplyr::desc(batch_effect),
      batch,
      control
    )
  
  # Convert controls to factor with reversed order for plotting
  hmdf$control <- factor(hmdf$control, levels = rev(controls))
  
  # Compute batch frequency of flags for reordering
  batch_freq <- data.frame(table(hmdf$batch, hmdf$batch_effect))
  colnames(batch_freq) <- c("batch", "flag", "Freq")
  batch_freq <- batch_freq[order(batch_freq$flag, batch_freq$Freq,
                                 decreasing = TRUE), ]
  
  hmdf$batch <- factor(
    hmdf$batch,
    levels = unique(as.character(batch_freq$batch))
  )
  
  # Ensure marker is character before factoring
  hmdf$marker <- as.character(hmdf$marker)
  
  # Factor markers 
  markers <- unique(c(markers, "Clusters"))
  hmdf$marker <- factor(hmdf$marker, levels = markers)
  
  # min_batch_effect: Minimum number of flags per batch-control combo required to keep the batch-control combo in the heatmap (default = 0)
  if (min_batch_effect > 0) {
    
    bc_flag_count <- hmdf %>%
      dplyr::filter(marker != "Clusters") %>%
      dplyr::group_by(batch, control) %>%
      dplyr::summarise(
        n_flags = sum(batch_effect == TRUE, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Keep only batch–control combos meeting the threshold
    hmdf <- hmdf %>%
      dplyr::semi_join(
        dplyr::filter(bc_flag_count, n_flags >= min_batch_effect),
        by = c("batch", "control")
      )
    
    # Drop unused factor levels for batch & control
    hmdf$batch   <- factor(hmdf$batch)
    hmdf$control <- factor(hmdf$control)
  }
  
  # Factor
  if (is.factor(hmdf$marker) && !"Clusters" %in% levels(hmdf$marker)) {
    levels(hmdf$marker) <- c(levels(hmdf$marker), "Clusters")
  }
  
  hmdf$marker[hmdf$type %in% c("Threshold", "Threshold2")] <- "Clusters"
  
  # min_marker_effect: Minimum number of flags per marker required
  # to keep the marker combo in the heatmap (default = 0).
  if (min_marker_effect > 0) {
    
    # Count flags per marker
    bc_flag_count <- hmdf %>%
      dplyr::filter(marker != "Clusters", !type %in% c("Threshold","Threshold2","Flag1","Flag2")) %>%
      dplyr::group_by(marker) %>%
      dplyr::summarise(n_flags = sum(batch_effect == TRUE, na.rm = TRUE), .groups = "drop")
    
    # Keep only markers meeting the threshold
    keep_markers <- bc_flag_count %>%
      dplyr::filter(n_flags >= min_marker_effect) %>%
      dplyr::pull(marker) %>%
      unique()
    
    # Filter hmdf 
    hmdf <- hmdf %>%
      dplyr::filter(marker == "Clusters" | marker %in% keep_markers)
    
    # Drop unused levels 
    hmdf$marker <- droplevels(hmdf$marker)  
  }
  
  # Rename Threshold types to Flag1/Flag2
  hmdf$type <- dplyr::recode(
    hmdf$type,
    "Threshold"  = "Flag1",
    "Threshold2" = "Flag2"
  )
  
  hmdf$type <- factor(
    hmdf$type,
    levels = c("-MFI", "+MFI", "%pos", "EMD", "Flag1", "Flag2")
  )
  
  # Reorder markers by flag frequency 
  hmdf$marker <- as.character(hmdf$marker)
  
  marker_order <- hmdf %>%
    dplyr::filter(marker != "Clusters") %>%
    dplyr::group_by(marker) %>%
    dplyr::summarise(
      n_flag = sum(batch_effect == TRUE, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(n_flag)) %>%
    dplyr::pull(marker)
  
  final_marker_level <- c(marker_order, "Clusters")
  hmdf$marker <- factor(hmdf$marker, levels = final_marker_level)
  
  # Axis sizes via helper 
  xlab_size <- auto_axis_size(axis = 'x', num_markers = length(markers),
                              num_controls = length(controls))
  
  ylab_size <- auto_axis_size(axis = 'y', num_markers = length(markers),
                              num_controls = length(controls))
  
  # Plot heatmap 
  p <- ggplot(hmdf, aes(x = type, y = control, fill = cols)) +
    geom_tile(color = "black", lwd = 0.6) +
    scale_fill_identity(na.value = "white") +
    scale_y_discrete(position = "right", expand = c(0, 0)) +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    facet_grid(batch ~ marker, scales = "free_x", space = "free_x", switch = "y") +
    theme_bw() +
    theme(
      plot.title = element_text(size = xlab_size + 10, face = "bold", hjust = 0.5),
      strip.text.y = element_text(size = 21, face = "bold", colour = "black"),
      strip.text.x = element_text(size = xlab_size - 2, face = "bold"),
      strip.background.x = element_blank(),
      strip.placement = "outside",
      axis.text.x = element_text(size = xlab_size + 3.5, angle = 90, hjust = -0.3),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = ylab_size + 3, hjust = 0.7, vjust = 0.4),
      panel.spacing.x = unit(0.4, 'lines'),
      panel.spacing.y = unit(0.5, 'lines'),
      panel.grid = element_blank()
    ) +
    ggtitle("Flags for IQR, EMD, and clustering metrics with batch thresholds")
  
  grid.plt <- ggplot_gtable(ggplot_build(p))
  
  # Assign colors to batch strips
  names(batch_cols) <- batch_list
  ord_bcols <- batch_cols[levels(hmdf$batch)]
  
  stripr <- which(grepl("strip-l", grid.plt$layout$name))
  k <- 1
  for (i in stripr) {
    j <- which(grepl("rect", grid.plt$grobs[[i]]$grobs[[1]]$childrenOrder))
    grid.plt$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- ord_bcols[k]
    k <- k + 1
  }
  
  # Draw final heatmap
  grid::grid.draw(grid.plt)
}

### Generates PNG ranked flagged heatmaps
# emd_df: A dataframe of EMD-based flags (columns must include 'batch', 'value', 'control', 'marker', 'flagged_batch', 'batch_effect', 'cols', and 'type')
# iqr_df: A dataframe of IQR-based flags (same required columns as emd_df)
# cluster_df: A data frame containing 'cluster', 'sample_id', 'control', and 'batch' columns for each cell
# cell_cluster_vector: A vector assigning each cell to a cluster
# set_threshold1: Threshold 1 for flag 1, specifically a numeric multiplier for the expected proportion threshold (e.g., 1.5 for 1.5x expected)
# set_threshold2: Threshold 2 for flag 2, specifically a numeric multiplier for the expected proportion threshold (e.g., 1.5 for 1.5x expected)
# batch_list: A character vector of batch IDs 
# controls: A character vector of control sample IDs 
# markers: A character vector of marker names 
# batch_cols: A named vector of colours for batches (names must match batch_list)
# output_dir: Directory where PNG images will be saved
# markers_per_plot: Number of markers to include per image panel (default = all markers)
# batches_per_plot: Number of batches to include per image panel (default = all batches)
# min_batch_effect: Minimum number of flags per batch required to keep the batch in the heatmap
# min_marker_effect: Minimum number of flags per marker required to keep the marker in the heatmap (default = 0)
generate_ranked_flagged_hmap_all <- function(emd_df,
                                             iqr_df,
                                             cluster_df,
                                             cell_cluster_vector,
                                             set_threshold1 = 2,
                                             set_threshold2 = 1.5,
                                             batch_list,
                                             controls,
                                             markers,
                                             batch_cols,
                                             output_dir,
                                             markers_per_plot = NULL,
                                             batches_per_plot = NULL,
                                             min_batch_effect = 0,
                                             min_marker_effect = 0) {
  
  # Defaults: use everything in one image (backwards compatible)
  if (is.null(markers_per_plot)) markers_per_plot  <- length(markers)
  if (is.null(batches_per_plot)) batches_per_plot  <- length(batch_list)
  
  # Cache sizes
  nm  <- length(markers)
  nb  <- length(batch_list)
  nct <- length(controls)
  
  # Basic safety checks
  if (markers_per_plot <= 0L) stop("markers_per_plot must be >= 1")
  if (batches_per_plot <= 0L) stop("batches_per_plot must be >= 1")
  
  # Keep only needed columns once 
  iqr_df <- subset(
    iqr_df,
    select = c(batch, value, control, marker, flagged_batch,
               batch_effect, cols, type)
  )
  
  # Initialize image counter
  img_counter <- 0L
  
  # Loop over marker chunks
  for (m_start in seq(1L, nm, by = markers_per_plot)) {
    m_end          <- min(m_start + markers_per_plot - 1L, nm)
    markers_chunk  <- markers[m_start:m_end]
    
    # Loop over batch chunks
    for (b_start in seq(1L, nb, by = batches_per_plot)) {
      b_end        <- min(b_start + batches_per_plot - 1L, nb)
      batch_chunk  <- batch_list[b_start:b_end]
      
      img_counter <- img_counter + 1L
      
      # Size depends on how many markers/batches in this panel
      width  <- image_width(length(markers_chunk))
      height <- image_height(length(batch_chunk), nct)
      
      # Optional: subset iqr_df (and emd_df) to speed things up for each panel
      iqr_df_chunk <- subset(
        iqr_df,
        marker %in% markers_chunk & batch %in% batch_chunk
      )
      
      emd_df_chunk <- subset(
        emd_df,
        marker %in% markers_chunk & batch %in% batch_chunk
      )
      
      # # Construct file name
      # if (markers_per_plot >= length(markers) && batches_per_plot >= length(batch_list)) {
      #   # No split: all markers and all batches in one plot
      #   filename <- "summary_heatmap_IQR_EMD_cluster_flags_all_markers_all_batches.png"
      # } else {
      #   # Split: use ranges in the name
      #   filename <- sprintf(
      #     "summary_heatmap_IQR_EMD_cluster_flags_markers%02d-%02d_batches%02d-%02d.png",
      #     m_start, m_end, b_start, b_end
      #   )
      # }
      
      
      # Decide whether everything fits in a single plot after any filtering
      all_markers_fit <- markers_per_plot >= length(markers)      # markers = post-filter list
      all_batches_fit <- batches_per_plot >= length(batch_list)   # batch_list = post-filter list
      
      has_marker_filter <- min_marker_effect > 0
      has_batch_filter  <- min_batch_effect  > 0
      
      # Prefix to reflect filter state in the filename
      filter_prefix <- if (has_marker_filter && has_batch_filter) {
        "filtered_markers_filtered_batches_"
      } else if (has_marker_filter) {
        "filtered_markers_"
      } else if (has_batch_filter) {
        "filtered_batches_"
      } else {
        ""
      }
      
      if (all_markers_fit && all_batches_fit) {
        # No split
        if (filter_prefix == "") {
          filename <- "summary_heatmap_IQR_EMD_clusters_flags_all_markers_all_batches.png"
        } else {
          # Filtered but still a single plot: avoid "all_markers_all_batches" wording
          filename <- paste0("summary_heatmap_IQR_EMD_clusters_flags_", filter_prefix, "single_plot.png")
        }
      } else {
        # Split
        if (all_markers_fit && !all_batches_fit) {
          filename <- sprintf(
            "summary_heatmap_IQR_EMD_clusters_flags_%sall_markers_batches%02d-%02d.png",
            filter_prefix, b_start, b_end
          )
        } else if (!all_markers_fit && all_batches_fit) {
          filename <- sprintf(
            "summary_heatmap_IQR_EMD_clusters_flags_%smarkers%02d-%02d_all_batches.png",
            filter_prefix, m_start, m_end
          )
        } else {
          filename <- sprintf(
            "summary_heatmap_IQR_EMD_clusters_flags_%smarkers%02d-%02d_batches%02d-%02d.png",
            filter_prefix, m_start, m_end, b_start, b_end
          )
        }
      }
      
      png_file_path <- file.path(output_dir, filename)
      
      # Open PNG device 
      png(filename = png_file_path, width = width, height = height, res = 400)
      
      tryCatch({
        
        # Draw ranked flagged heatmap for the current chunk
        ranked_flagged_hmap_all(
          emd_df          = emd_df_chunk,
          iqr_df          = iqr_df_chunk,
          cluster_df      = cluster_df,
          cell_cluster_vector = cell_cluster_vector,
          set_threshold1 = set_threshold1,
          set_threshold2 = set_threshold2,
          batch_list      = batch_chunk,
          controls        = controls,
          markers         = markers_chunk,
          batch_cols      = batch_cols,
          output_dir      = output_dir,
          min_batch_effect = min_batch_effect,
          min_marker_effect = min_marker_effect
        )
      }, error = function(e) {
        message("Error in generating summary heatmap (image ", img_counter, "): ",
                e$message)
      }, finally = {
        
        # Close
        dev.off()
      })
      
      # Console message
      if (nm > markers_per_plot || nb > batches_per_plot) {
        
        # Multiple heatmaps (splitting occurred)
        cat("Summary of IQR, EMD, and cluster flags heatmap, part", img_counter,
            "saved as", png_file_path, "\n")
      } else {
        
        # Only one heatmap (no split)
        cat("Summary of IQR, EMD, and cluster flags heatmap saved as", png_file_path, "\n")
      }
    }
  }
}


