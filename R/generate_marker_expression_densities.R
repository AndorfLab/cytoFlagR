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
suppressPackageStartupMessages(library(progress))
suppressPackageStartupMessages(library(crayon))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggpmisc))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(grid))

########### ~~~~ Marker Expression Density/Histogram Plot ~~~~ ###########

### Subsamples rows for reference density plots by sampling a fixed or minimum number of cells per sample
# df: A dataframe containing at least columns 'sample_id' and 'batch'
# samps: A character vector of sample_id values to include in the subsampling
# ref_subsampling_size: Integer used when sampling_type == "fixed" (number of cells per sample)
# seed: Integer seed for reproducibility of random sampling
# Returns a dataframe (subset of df) containing the concatenated subsamples across all 'samps'
subsample_for_density_plots_refs <- function(df, samps, ref_subsampling_size, seed) {
  
  # Get number of cells per sample and set sample size as min. cells per batch
  ncells <- data.frame(table(df$sample_id))
  colnames(ncells) <- c("sample_id", "Freq")
  
  # Map each sample_id to its batch in df
  mm <- match(ncells$sample_id, df$sample_id)
  ncells$batch <- df$batch[mm]
  
  # Compute per-batch minimum number of cells among its samples
  ncells <- ncells %>%
    dplyr::group_by(batch) %>%
    dplyr::mutate(min_per_batch = min(Freq))
  
  # Annotate df with min_per_batch by sample_id
  mm <- match(df$sample_id, ncells$sample_id)
  df$min_cells <- ncells$min_per_batch[mm]
  
  # Set seed for reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Container for sampled indices
  subsamp <- c()
  
  # Takes the number of cells from the file that has the smallest number of cells
  for (fi in 1:length(samps)) {
    
    # Indices for current sample_id
    tmp <- which(df[, "sample_id"] == samps[fi])
    
    # Determine number of cells to take:
    # default: minimum across all files (global min of df$min_cells)
    # if sampling_type == "sample_min": use per-sample min_cells (here implemented as global min)
    # if sampling_type == "fixed": use ref_subsampling_size
    mc <- min(unique(df$min_cells))  # a fixed number, the minimum in any file, applied to all files
   
    if (sampling_type == "sample_min") { mc <- min(unique(df$min_cells)) }
    if (sampling_type == "fixed")      { mc <- ref_subsampling_size }
    
    # Sample 'mc' indices from current sample_id
    subsamp <- c(subsamp, tmp[sample(1:length(tmp), mc)])
  }
  
  # Subset df to sampled rows
  rdf <- df[subsamp, ]
  
  # Return the subsampled dataframe
  return(rdf)
}

### Generates density plots of marker expression with an automated cutoff line
# df: A dataframe containing columns 'expression', 'marker', 'batch', 'control', and 'sample_id' if ref subsampling is used.
# marker: A character string indicating the markers
# batch_list: A character vector listing batches
# control_list: A character vector listing control groups
# cutf: Numeric value for the vertical cutoff line
# axis_size: Numeric base size for axis text and titles (default = NULL)
# ref_subsampling: Either "no" or another string
# if "no", duplicates df as reference. Otherwise, creates a reference by subsampling via subsample_for_density_plots_refs().
# samps: A character vector of sample_id values used only when ref_subsampling != "no"
# ref_subsampling_size: Integer passed to subsample_for_density_plots_refs() when using fixed sampling
# seed: Integer seed for reproducibility of random sampling
# Returns a ggplot object with densities, faceted by control, and colored by reference status
density_plots_w_thresholds <- function(df,
                                       marker,
                                       batch_list,
                                       control_list,
                                       cutf,
                                       axis_size = NULL,
                                       ref_subsampling,
                                       samps = samps,
                                       ref_subsampling_size,
                                       seed) {
  
  if (ref_subsampling == "no") {
    
    # Main plotting df
    mdf <- data.frame(df)
    mdf$reference <- "no"
    
    # Create reference by duplicating df and tagging as 'ref' batch
    ggd_bg <- mdf
    ggd_bg$batch <- "ref"
    ggd_bg$reference <- "yes"
    
    # Combine layers and set factor levels for control/batch
    ggd_plot <- rbind(mdf, ggd_bg)
    ggd_plot$control <- factor(ggd_plot$control, levels = control_list)
    ggd_plot$batch   <- factor(ggd_plot$batch,   levels = c(rev(batch_list), "ref"))
    
  } else {
    
    # Use minimal amount of cells for sampling 
    ggd_bg <- subsample_for_density_plots_refs(df, samps = samps, ref_subsampling_size = ref_subsampling_size, seed = seed)
    ggd_bg$min_cells <- NULL  
    
    # Main plotting df
    mdf <- data.frame(df)
    mdf$reference <- "no"
    
    # Tag reference layer
    ggd_bg$batch     <- "ref"
    ggd_bg$reference <- "yes"
    
    # Combine layers and set factor levels for control/batch
    ggd_plot <- rbind(mdf, ggd_bg)
    ggd_plot$control <- factor(ggd_plot$control, levels = control_list)
    ggd_plot$batch   <- factor(ggd_plot$batch,   levels = c(rev(batch_list), "ref"))
  }
  
  # Density plots with cutoff line 
  # Faceting by control
  ggplot() +
    ggridges::geom_density_ridges(
      data = ggd_plot,
      aes(x = expression, y = batch, color = reference, fill = reference),
      alpha = 0.3,
      show.legend = FALSE
    ) +
    ggtitle(paste0("Expression densities for ", marker, " with automated cutoff")) +
    geom_vline(xintercept = cutf, linewidth = 1.0, colour = "#121518") +
    facet_wrap(~ control, nrow = 1) +
    xlab("Expression") +
    ylab("Batch") +
    ggridges::theme_ridges(center_axis_labels = TRUE) +
    theme(
      axis.text.x  = element_text(size = axis_size),
      axis.text.y  = element_text(size = axis_size + 7),
      axis.title   = element_text(size = axis_size + 9),
      plot.title   = element_text(size = axis_size + 9, face = "bold", hjust = 0.5),
      strip.text   = element_text(size = axis_size + 5, face = "bold"),
      legend.position = "none"
    )
}
  
### Creates a long dataframe of marker expression values that can be used for density plots
# df: A dataframe containing marker columns and metadata (must include 'sample_id')
# sample_col_name: The column name in df that contains control/sample identifiers
# batch_col_name: The column name in df that contains batch identifiers
# markers: A character vector of marker column names
# Returns a long-format dataframe with columns 'sample_id', 'marker', 'expression', 'batch', and 'control'
reshape_df <- function(df, 
                       sample_col_name, 
                       batch_col_name, 
                       markers) {
  
  # Melt marker columns into long format with expression values
  molten_df <- melt(data.frame(df[,markers],sample_id = df[,"sample_id"]),
                    id.vars = "sample_id",variable.name = "marker",
                    value.name = "expression")
  
  # Match sample_id to append batch and control labels
  mm <- match(molten_df$sample_id,df$sample_id)
  molten_df[,"batch"] <- df[,batch_col_name][mm]
  molten_df[,"control"] <- df[,sample_col_name][mm]
  
  # Return the molten dataframe
  return(molten_df)
}         
         
### Plots and saves marker expression density plots to a multi-page PDF
# marker_list: A character vector of marker names to plot
# df: A wide-form dataframe containing marker columns and metadata. Must include 'batch' and 'control' columns
# batch_list: A character vector listing batches 
# control_list: A character vector listing control groups 
# auto_cutoffs: A dataframe with columns 'marker' and 'cutoff' providing numeric cutoff per marker
# batch_colm: The column name in df representing batches 
# control_colm: The column name in df representing control/sample identifiers
# output_dir: Directory where the PDF file will be written
# wd, ht: Numeric width and height of the PDF in inches
# axis_size: Numeric base size for plot axes and titles (default = 14)
# file_name: Base name for the output PDF file (default = "marker_expression_densities")
# ref_subsampling: Either "no" or "yes"
# samps: A character vector of sample_id values used when ref_subsampling != "no" 
# ref_subsampling_size: Integer specifying number of cells per sample for reference subsampling (default = 1000)
# seed: Integer seed for reproducibility of random sampling
# Returns a PDF file containing density plots for all markers in marker_list
write_density_plot_per_marker <- function(marker_list, 
                                          df, 
                                          batch_list, 
                                          control_list, 
                                          auto_cutoffs, 
                                          batch_colm, 
                                          control_colm, 
                                          output_dir, 
                                          wd = 12, 
                                          ht = 11, 
                                          axis_size = 14, 
                                          file_name = "marker_expression_densities",
                                          ref_subsampling, 
                                          samps, 
                                          ref_subsampling_size = 1000,
                                          seed) {
  
  # Inform user that the process may take a while
  cat("This process may take a long time, please be patient\n")
  
  # Reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Define a color style for progress bar text
  colrSet <- crayon::make_style("#539DDD")
  
  # Determine total number of markers to plot
  total_iter <- length(marker_list)
  
  # Create progress bar
  pBar <- progress::progress_bar$new(
    format = paste0(colrSet("Plotting marker expression densities"),
                    " [", colrSet(":bar"), "] ",
                    colrSet(":percent"), " | Marker: :current/:total ",
                    "| Elapsed: :elapsed | ETA: :eta"),
    total = total_iter,
    clear = FALSE,
    width = 80,
    complete = "=",
    incomplete = "-"
  )
  
  # Build PDF file name and path
  file_name <- paste0(file_name, ".pdf")
  pdf_file_path <- file.path(output_dir, file_name)
  
  # Ensure sample_id column exists; create if missing
  if (is.null(df[, "sample_id"])) {
    df[, "sample_id"] <- paste0(df[, batch_colm], "_", df[, control_colm]) 
  }
  
  # Reshape dataframe for plotting (long format)
  cdf <- data.frame(reshape_df(df, control_colm, batch_colm, marker_list))
  
  # Open the PDF device for multi-page output
  pdf(pdf_file_path, width = wd, height = ht, useDingbats = FALSE)
  
  # Loop through each marker and generate plots
  for (i in 1:length(marker_list)) {
    
    # Current marker
    mk <- marker_list[i]
    
    # Update progress bar
    pBar$tick(tokens = list(current = i))
    
    # Subset dataframe for the current marker
    mdf <- cdf[(cdf$marker == marker_list[i]), ]
    mdf$marker <- as.character(mdf$marker)
    
    # Retrieve automated cutoff for current marker
    cutf <- as.numeric(auto_cutoffs[which(auto_cutoffs$marker == mk), "cutoff"])
    
    # Build the density plot with thresholds and optional reference subsampling
    plot <- density_plots_w_thresholds(df = mdf,
                                       marker = mk,
                                       batch_list = batch_list,
                                       control_list = control_list,
                                       cutf = cutf,
                                       axis_size = axis_size,
                                       ref_subsampling = ref_subsampling,
                                       samps = samps,
                                       ref_subsampling_size = ref_subsampling_size,
                                       seed = seed)
    
    # Print the plot to the current PDF page
    print(plot)
    
    # Notify progress in console
    cat("\n", i, "of", length(marker_list), "markers plotted", "\n")
  }
  
  # Close the PDF device
  dev.off()
  
  # Notify user of completion
  cat("Density plot saved as", pdf_file_path, "\n")
}

########### ~~~~ Marker Expression Biaxial Dotplots ~~~~ ###########

### Computes 2D kernel density values for paired vectors
# x: Numeric vector for the x-axis values
# y: Numeric vector for the y-axis values
# ...: Additional arguments passed to MASS::kde2d (e.g., n, h, lims)
# Returns a numeric vector of density values aligned with input (x, y) pairs
get_density <- function(x, 
                        y, 
                        ...) {
  
  # Compute 2D kernel density over a grid
  dens <- MASS::kde2d(x, y, ...)
  
  # Find x grid indices for each x
  ix <- findInterval(x, dens$x)
  
  # Find y grid indices for each y
  iy <- findInterval(y, dens$y)
                     
  # Bind index pairs to sample from the density matrix
  ii <- cbind(ix, iy)
  
  # Return density values 
  return(dens$z[ii])
}              
                     
### Generates biaxial dotplots faceted by batch and control
# df: A dataframe with columns 'marker', 'ref', 'density', 'batch', and 'control'
# ref_marker: A character string for the y-axis label representing the reference marker
# select_marker: A character string for the x-axis label representing the selected marker
# colrs: A vector of colors 
# cutf: Numeric cutoff value for a vertical line on the x-axis
# axis_size: Numeric base size for axis text and titles (default = NULL)
# Returns a biaxial dotplot ggplot
biax_grid_plot <- function(df,
                           ref_marker, 
                           select_marker, 
                           colrs,
                           cutf, 
                           axis_size=NULL) {
  
  # Build the dot plot with density coloring and facets
  ggplot() +
    geom_point(data = df, aes(x = marker, y = ref, color = density), size=0.1) +
    facet_grid(batch~control) +
    geom_vline(xintercept = cutf, colour="#222528", linewidth=1.5) + # Marker cutoff
    xlab(select_marker)+
    ylab(ref_marker)+
    scale_color_gradientn(colours = colrs) +
    ggtitle(paste0("Biaxial density dot plots of ",select_marker," against ",ref_marker)) +
    theme_bw() +
    theme(axis.text.y = element_text(size = axis_size+2),
          axis.text.x = element_text(size = axis_size),
          axis.title = element_text(size = axis_size+6,face = "bold"),
          strip.text.y = element_text(size = axis_size+6,face = "bold"),
          strip.text.x = element_text(size = axis_size+7, face = "bold"),
          legend.title = element_text(size = (axis_size-3)),
          legend.text = element_text(size = (axis_size-5)),
          plot.title = element_text(size = 24, face = "bold", hjust = 0.5))
}
  
### Subsamples rows per sample for biaxial dot plot creation
# df: A dataframe containing at least columns 'sample_id' and 'batch'
# samps: A character vector of sample_id values 
# sampling_type: Either "batch_min" or "overall_min", each specifying subsampling strategy
# seed: Integer seed for reproducibility of random sampling
# Returns a dataframe containing the subsampled rows across all 'samps'
subsample_for_biax_plot <- function(df, 
                                    samps, 
                                    sampling_type, 
                                    seed){
  
  # Get number of cells per sample and set sample size as min. cells per batch
  ncells <- data.frame(table(df$sample_id)) 
  colnames(ncells) <- c("sample_id","Freq")
  
  # Map sample_id to batch in df
  mm <- match(ncells$sample_id,df$sample_id)
  ncells$batch <- df$batch[mm]
  ncells <- ncells %>% group_by(batch) %>% mutate(min_per_batch = min(Freq))
  
  # Annotate df with min_per_batch by sample_id
  mm <- match(df$sample_id,ncells$sample_id)
  df$min_cells <- ncells$min_per_batch[mm]
  
  # Set seed for reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Initialize vector of sampled indices
  subsamp <- c()
  
  # Takes the number of cells from the file that has the smallest number of cells
  for(fi in 1:length(samps)){
    tmp <- which(df[,"sample_id"]==samps[fi])
    
    # Choose sample size based on sampling_type
    if (sampling_type == "batch_min") {mc <- unique(df[which(df[,"sample_id"]==samps[fi]),"min_cells"])}
    if (sampling_type == "overall_min") {mc <- min(unique(df$min_cells))} # A fixed number, the minimum in any file, applued to all files}
    if (sampling_type != "overall_min" | sampling_type == "batch_min") {stop("Invalid sampling_type selection")}
   
     # Sample indices for this sample_id
    subsamp <- c(subsamp, tmp[sample(1:length(tmp), mc)])
  }
  
  # Subset df to sampled rows
  rdf <- df[subsamp,]
  
  # Return the subsampled dataframe
  return(rdf)
}

### Generates and saves biaxial density dot plots of markers against a reference marker
# df: A dataframe containing marker columns and metadata (requires 'batch' and 'control')
# ref_marker: A character string specifying the reference marker 
# marker_list: A character vector of markers to plot 
# control_list: A character vector specifying control order
# batch_list: A character vector specifying batch order
# samps: A character vector of sample_id values
# auto_cutoffs: A dataframe with columns 'marker' and 'cutoff'
# batch_col_name: The column name in df representing batches
# control_col_name: The column name in df representing control/sample identifiers
# wd, ht: Numeric values specifying the dimensions of the output plot in inches
# axis_size: Numeric base size for axis text and titles (default = 20)
# output_dir: Directory where PNG files will be saved 
# Returns: One PNG per marker (excluding the reference marker), saved to `output_dir`
draw_biaxial_dotplots <- function(df, 
                                  ref_marker, 
                                  marker_list, 
                                  control_list, 
                                  batch_list, 
                                  samps, 
                                  auto_cutoffs, 
                                  batch_col_name, 
                                  control_col_name, 
                                  wd = 15, 
                                  ht = 20, 
                                  axis_size = 20, 
                                  output_dir){
  
  # Check
  if(missing(output_dir)) {
    stop("output_dir parameter is required")
  }
  
  if (!dir.exists(output_dir)) {
    stop("Output directory does not exist: ", output_dir)
  }
  
  # Harmonize control column name to 'control'
  if(control_col_name!="control") {
    names(df)[names(df)==control_col_name] <- "control"
  }
  
  # Harmonize batch column name to 'batch'
  if(batch_col_name!="batch") {
    names(df)[names(df)==batch_col_name] <- "batch"
  }
  
  # Get sample id column if one doesn't exist
  if(is.null(df[,"sample_id"])) {
    df[,"sample_id"] <- paste0(df[,"batch"], "_", df[,"control"])
  }
  
  # Inform user that the process may take a while
  cat("This process may take a long time, please be patient\n")
  
  # Progress bar styling and setup
  colrSet <- crayon::make_style("#539DDD")
  total_iter <- length(marker_list) - 1
  pBar <- progress::progress_bar$new(
    format = paste0(colrSet("Plotting density dotplots"),
                    " [", colrSet(":bar"), "] ",
                    colrSet(":percent"), " | Marker: :current/:total ",
                    "| Elapsed: :elapsed | ETA: :eta"),
    total = total_iter,
    clear = FALSE,
    width = 80,
    complete = "=",
    incomplete = "-"
  )
  
  # Working dataframe
  cdf <- df
  
  # Color vector to represent dense regions 
  densColrs <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
  
  # Exclude the reference marker from the plotting list
  marker_list <- marker_list[marker_list != ref_marker]
  
  # Loop through each marker to create plots
  for (i in 1:length(marker_list)) {
    mk <- marker_list[i]
    
    # Build output file path (PNG)
    pdf_file_path <- file.path(
      output_dir, 
      paste0("density_dotplot_", mk, "_vs_", ref_marker, "_reference_marker.png")
    )
    
    # Columns required for plotting
    colms.to_use <- c(ref_marker, mk, "batch", "control")
    
    # Update progress bar
    pBar$tick(tokens = list(current = i))
    
    # Retrieve automated cutoff for current marker
    cutf <- as.numeric(auto_cutoffs[which(auto_cutoffs$marker == mk), "cutoff"])
    
    # Prepare plotting dataframe and compute densities
    mdf <- cdf[, names(cdf) %in% colms.to_use]
    mdf$density <- get_density(mdf[, ref_marker], mdf[, mk], n = 100) ### get densities
    colnames(mdf) <- c("ref", "marker", "batch", "control", "density")
    mdf$control <- factor(mdf$control, levels = control_list)
    mdf$batch   <- factor(mdf$batch,   levels = batch_list)
    
    # Generate dotplot
    plotp <- biax_grid_plot(df = mdf,
                            ref_marker = ref_marker,
                            select_marker = mk,
                            colrs = densColrs, 
                            cutf = cutf,
                            axis_size = axis_size)
    
    # Save plot as PNG
    tryCatch({
      png(pdf_file_path, width = wd, height = ht, units = "in", res = 300)
      print(plotp)
      dev.off()
      
      cat("\n", i, "of", length(marker_list), "markers plotted", "\n")
      
    }, error = function(r) {
      message("Error in plotting densities: ", r$message)
    })
  }
  
  # Notify user of completion
  cat("Biaxial plots saved in", output_dir, "\n")
}
