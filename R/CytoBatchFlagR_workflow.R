# ======================================================================
# Load libraries
# ======================================================================

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
suppressPackageStartupMessages(library(emdist))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(LaplacesDemon))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(flowCore))
suppressPackageStartupMessages(library(FlowSOM))
suppressPackageStartupMessages(library(ConsensusClusterPlus))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(umap))
suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(progress))
suppressPackageStartupMessages(library(crayon))
suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ggpmisc))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(shinyjs))
suppressPackageStartupMessages(library(shinyWidgets))
suppressPackageStartupMessages(library(shinycssloaders))
suppressPackageStartupMessages(library(bslib))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(rstudioapi))
suppressPackageStartupMessages(library(grid))

source("CytoBatchFlagR_functions.R")
source("auto_threshold_func.R")
source("editMarkerTable.R")
source("generate_marker_expression_densities.R")
source("selectMarkerList.R")
source("referenceMarkerSelection.R")

# ======================================================================
# Step 1: Pre-processing and Visual Assessment of the Data
# ======================================================================

# Select the folder that contains the FCS files, FCS metadata file, and marker information file. Make sure this folder is unzipped.

# For more information about the format these files should be in, please go here: https://github.com/AndorfLab/CytoBatchFlagR/wiki/Input-Files

# Folder where CytoBatchFlagR and related files are stored
work.dir <- getwd() 

# Folder where all FCS files, marker information file, and FCS metadata file are stored
data.dir <- rstudioapi::selectDirectory(caption = "Select FCS file folder", 
                                      label = "Select", 
                                      path = getActiveProject()) 

# Load in the FCS metadata file and marker information file.

# Read FCS metadata/information file
fcs.info <- data.frame(readr::read_delim(file.path(data.dir,"fcs_file_info.csv"))) 

# Read the marker and fluorophore/metal tag information file
marker_tbl <- data.frame(readr::read_delim(file.path(data.dir, "panel_markers.csv" ))) 
# Replace "-" with "."
marker_tbl$PnS <- gsub("-", ".", marker_tbl$PnS)

# Store the tag and marker lists, as well as their column names
tag_list <- marker_tbl[,"PnN"]
markers <- marker_tbl[,"PnS"]

# Replace "-" with "."
markers <- gsub("-", ".", markers)

# Next, set the output directory. All outputted files and plots will go here. If using the code below, this directory will be made within the same directory this Rmd file is currently located.

# Name of output directory - it will include the date/time this line is run
output.dir <- file.path(paste0("output_", format(Sys.time(), "%Y.%m.%d_%H%M%S")))

# This is the folder where your outputs are stored
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE) 

# Here, you should transform the data if you haven't already done so. This tool contains a function for the data to be transformed with an arcsinh function.

# For the arcsinh transformation, you must set a cofactor. Cofactors commonly used in the field are:
#   - Mass Cytometry (CyTOF) datasets: 5 
# - Conventional Flow Cytometry datasets: 150
# - Spectral Flow Cytometry datasets: 6000

# Create a transformed dataframe of all FCS files
transformation_class <- "arcsinh" 

# Set a cofactor for arcsinh transformation
cf <- 6000

# Transform
ffdf <- create_transformed_dataframe_function(transformation_class = transformation_class, 
                                              cf = cf, 
                                              fcs_info = fcs.info, 
                                              tag_list = tag_list, 
                                              marker_list = markers, 
                                              marker_tbl = marker_tbl,
                                              tag_colm = "PnN", 
                                              marker_colm = "PnS",
                                              output_dir = output.dir, 
                                              file_dir = data.dir)

# Set the names of the controls, batches, and sample IDs. If this is in the uploaded FCS information file, you can read that in with the code below. 

# Control names
controls <- as.character(unique(fcs.info$control)) 

# Batch names
batches <- as.character(unique(fcs.info$batch)) 

# Sample ID names
sample_ids <- as.character(unique(fcs.info$sample_id))

# Now, plots can be created as an initial visual assessment of the data. 

# The code below will create 3 plots:
# - Plot 1 shows the number of cells per sample
# - Plot 2 is an MDS plot
# - Plot 3 is an UMAP plot

# Plot 1: number of cells per sample
# Axis size, plot width, and plot height are customizable
generate_number_of_cells_plot(df = ffdf, 
                              batch_colm = "batch", 
                              control_colm = "control", 
                              batch_list = batches, 
                              control_list = controls, 
                              output_dir = output.dir, 
                              axis_size = 20, 
                              width = 11.5, 
                              height = 8)

# Plot 2: MDS plot across samples 
# Axis size, plot width, and plot height are customizable
generate_MDS_plot(df = ffdf, 
                  batch_colm = "batch", 
                  control_colm = "control",
                  marker_list = markers, 
                  sample_colours = sample_colours, 
                  output_dir = output.dir, 
                  axis_size = 24, 
                  width = 12, 
                  height = 12)

# Plot 3: UMAP of controls annotated by batch colors
# Axis size, plot width, and plot height are customizable
generate_UMAP(df = ffdf, 
              batch_colm = "batch", 
              control_colm = "control", 
              marker_list = markers, 
              batch_list = batches, 
              control_list = controls, 
              sample_ids = sample_ids, 
              batch_colours = batchColours, 
              output_dir = output.dir, 
              axis_size = 24, 
              width = 12, 
              height = 12)

# Finally, select the markers you want to use for the IQR-based and EMD-based assessments. Once the markers are selected, click 'Done'.

# Use this interface to select markers you want to use for the IQR-based assessment
markers_1 <- selectMarkerList(marker_tbl, 
                              column_name = "PnS")

# ======================================================================
# Step 2: Interquartile Range (IQR)-Based Assessment
# ======================================================================

# This section handles the IQR-based assessment. 

# First, the tool generates automated marker thresholds for each marker in the dataset. These thresholds split the data into positive and negative groups.

# Make sure you set: 
# - dataset_type = "flow" for spectral or conventional flow cytometry datasets
# - dataset_type = "mass" for mass cytometry (CyTOF) datasets

# Depending on how large your data is, this process could take a long time to run. You may want to subsample cells by setting subsample = "yes" and specifying the number of cells by setting num. 

# Subsample: change to "yes" if you want CytoBatchFlagR to subsample your data to calculate the marker thresholds
# If you set subsample = "yes", the seed and number of cells per batch to downsample to (num) should also be set
auto_cutoffs <- data.frame(automated_threshold_DF(df = ffdf, 
                                                  markers = markers_1,
                                                  control_list = controls, 
                                                  batch_list = batches,
                                                  batch_colm = "batch", 
                                                  subsample = "no", 
                                                  seed = 450,
                                                  num = 20000,
                                                  output_dir = output.dir))

# Next, density plots and biaxial plots can be created to visualize the data and help determine if the automated marker thresholds are sufficient or if they should be overriden. 

# Prior to making these plots, you may want to downsample your data because it can be time-consuming to use all cells when your dataset is large. 

# There are 4 sampling_type options:
#   - "none": No subsampling
# - "overall_min": Take the same number of cells equal to the smallest sample size across all samples
# - "fixed": Take exactly n_cells cells per sample. If a sample has fewer than n_cells, it takes all available cells for that sample
# - "batch_min": For each batch, compute the minimum sample size within that batch. Then for every sample in that batch, take that batch’s minimum count. This equalizes cell counts within each batch.

# When sampling_type = "fixed", you must have n_cells set. n_cells should be an integer specifying how many cells to take per sample.

# The subsampling done here will be used when making the density and biaxial plots  

# Set sampling type
# Set as "fixed", "overall_min", "batch_min", or "none" 
sampling_type <-  "overall_min"  

# This will only be used if "fixed" was chosen for sampling_type
n_cells <- 15000

# Show original number of rows and columns of the cytometry data
print(dim(ffdf))

# Subsample data
ffdf_sampled <- sample_cells(df = ffdf, 
                             sampling_type = sampling_type, 
                             n_cells = n_cells)

# Show post-subsampled number of rows and columns of the cytometry data
print(dim(ffdf_sampled))

# Now, you can plot the marker expression densities across controls and batches with the automated cutoff shown. This process can be very slow and may take a while to generate density plots, depending on data size.

# To further speed up the process, ref_subsampling can be set to "yes", which means the reference will be subsampled. Concurrently, the ref_subsampling_size can be set to a number, specifying number of cells per sample for reference subsampling. The seed can also be set to ensure reproducability. If ref_subsampling = "no", the reference will contain all cells.

# If your data has a large amount batches, you may want to split the data into multiple plots to improve readability. The batches_per_plot parameter can be used in this cases. 
# batches_per_plot indicates the maximum number of batches to include in each plot. For example, if you have 12 batches and batches_per_plot = 5, then for each marker you will have 3 plots (plot 1: batches 1-5, plot 2: batches 6-10, plot 3: batches 11-12)
# batches_per_plot = NULL is used if you want to include all batches.

# Once these plots are outputted, examine the densities to help determine if the automated cutoffs sufficiently split the marker expression values into a negative and positive population. 

# Plot and save densities to pdf
# Axis size, plot width, plot height, output file name, and seed are customizable
write_density_plot_per_marker(marker_list = markers_1, 
                              df = ffdf_sampled, 
                              batch_list = batches,
                              control_list = controls, 
                              auto_cutoffs = auto_cutoffs,
                              batch_colm = "batch", 
                              control_colm = "control",
                              output_dir = output.dir, 
                              wd = 11, 
                              ht = 10, 
                              axis_size = 14,
                              file_name = "marker_expression_densities_original_cutoffs",
                              ref_subsampling = "yes", 
                              samps = sample_ids, 
                              ref_subsampling_size = 20000,
                              seed = 450,
                              batches_per_plot = NULL)

# Based on your assessment of the data, you may want to edit some marker threshold values.

# You can run the function below to optionally edit and update marker thresholds. Make sure you close the app after downloading a file with the new, edited threshold values.

# Launch interface to update your marker thresholds
editMarkerTable(auto_cutoffs)

# After downloading the file, the new cutoffs are saved and assigned as updated_auto_cutoffs.

# If you want to further examine the new automated cutoffs, you can optionally choose to re-plot marker expression densities across controls and batches with the new automated cutoffs.

# Axis size, plot width, and plot height are customizable
write_density_plot_per_marker(marker_list = markers_1, 
                              df = ffdf_sampled, 
                              batch_list = batches, 
                              control_list = controls, 
                              auto_cutoffs = updated_auto_cutoffs, 
                              batch_colm = "batch", 
                              control_colm = "control", 
                              output_dir = output.dir, 
                              wd = 11, 
                              ht = 10, 
                              axis_size = 14, 
                              file_name = "marker_expression_densities_updated_cutoffs", 
                              ref_subsampling = "yes", 
                              samps = sample_ids,
                              seed = 450)

# Next, you can generate biaxial dot plots against a reference marker.

# This process can be very slow and may take a while to generate plots, depending on data size.

# Choose a reference marker using the pop-up interface
ref_marker <- referenceMarkerSelection(marker_tbl, 
                                       column_name = "PnS")

# Pick a sampling_type for the biaxial plots: "overall_min" or "batch_min"
# "overall_min": pick the minimum number of cells of any file as sampling size (quicker plot)
# "batch_min": the minimum number of cells of a file in a batch as sample size (very slow) 
# Axis size, plot width, and plot height are customizable
draw_biaxial_dotplots(df = ffdf_sampled, 
                      ref_marker = ref_marker, 
                      marker_list = markers_1, 
                      control_list = controls, 
                      batch_list = batches, 
                      batch_col_name = "batch", 
                      control_col_name = "control", 
                      samps = sample_ids, 
                      auto_cutoffs = updated_auto_cutoffs, 
                      wd = 12, 
                      ht = 11, 
                      axis_size = 18, 
                      output_dir = output.dir)

# Using the thresholds, get the negative population (-MFI), positive population (+MFI), and percent positive (%pos) cells for all markers and controls.

# -MFI
nMFIdf <- data.frame(neg_mfi_df_all(df = ffdf, 
                                    markers = markers_1, 
                                    cutf = updated_auto_cutoffs, 
                                    sample_col_name = "control", 
                                    batch_col_name = "batch", 
                                    min_cells = TRUE))

# +MFI
pMFIdf <- data.frame(pos_mfi_df_all(df = ffdf, 
                                    markers = markers_1, 
                                    cutf = updated_auto_cutoffs, 
                                    sample_col_name = "control", 
                                    batch_col_name = "batch", 
                                    min_cells = TRUE))

# %pos
posdf <- data.frame(percent_pos_df_all(df = ffdf, 
                                       markers = markers_1, 
                                       cutf = updated_auto_cutoffs, 
                                       sample_col_name = "control", 
                                       batch_col_name = "batch"))

# Now, you can generate a dataframe that contains flags generated from all 3 IQR metrics above (-MFI, +MFI, and %pos).

iqrdf <- iqr_dataframe(neg_df = nMFIdf, 
                       pos_df = pMFIdf, 
                       percent_pos_df = posdf, 
                       markers = markers_1, 
                       sample_col_name = "control", 
                       batch_col_name = "batch", 
                       output_dir = output.dir)

# Optional: you can change the names of controls to 1, 2, 3, 4, etc.
control_labels <- as.character(seq(1,length(controls)))

# Assign control sample name to labels
names(control_labels) <- controls

# The resulting data can be visualized. The flags generated by the IQR-based metrics are shown in a boxplot. These flags are marked with the name of that batch next to the flagged point. In contrast, a point without a label indicates it was not flagged.

# The boxplots below have free axes, meaning the axes will change with each boxplot depending on the underlying data. 

# Seed can be set to ensure the jitter pattern (offset location of the points) in the plot is the same each time. If you want the jitter to be random each time you run the function, set seed = NULL.

# If using different labels for controls, use ctrl_labs = control_labels, 
# Otherwise, use ctrl_labs = controls while generating the IQR boxplots below

# Plot width, plot height, and seed are customizable

# -MFI
# Boxplots with flags 
generate_IQR_metric_boxplots(iqr.plt = iqrdf, 
                             type = "-MFI", 
                             markers = markers_1,
                             ctrl_labs = controls, 
                             colours = batchColours, 
                             batch_list = batches, 
                             control_list = controls, 
                             output_dir = output.dir, 
                             width = 11, 
                             height = 10,
                             seed = 450)

# +MFI
# Boxplots with flags
generate_IQR_metric_boxplots(iqr.plt = iqrdf, 
                             type = "+MFI", 
                             markers = markers_1, 
                             ctrl_labs = controls, 
                             colours = batchColours,
                             batch_list = batches, 
                             control_list = controls, 
                             output_dir = output.dir, 
                             width = 11, 
                             height = 10,
                             seed = 450)

# %pos
# Boxplots with flags
generate_IQR_metric_boxplots(iqr.plt = iqrdf, 
                             type = "%pos", 
                             markers = markers_1,
                             ctrl_labs = controls, 
                             colours = batchColours, 
                             batch_list = batches, 
                             control_list = controls, 
                             output_dir = output.dir, 
                             width = 11, 
                             height = 10,
                             seed = 450)

# The boxplots below have fixed axes, meaning the axes are set to be the same across all boxplots per metric. This is done to have a standard comparison with the same axes across all the markers.

# Seed can be set to ensure the jitter pattern (offset location of the points) in the plot is the same each time. If you want the jitter to be random each time you run the function, set seed = NULL.

# If using different labels for controls, use ctrl_labs = control_labels, 
# Otherwise, use ctrl_labs = controls while generating the IQR boxplots below

# Plot width, plot height, and seed are customizable
# y_limits can be set at "auto", which sets the minimum and maximum respectively as the minimum and maximum values across controls and markers.
# Alternatively, these values can be set with c(min, max). For example, c(0, 10) would set the minimum as 0 and the maximum as 10 across all markers for that metric. 

# -MFI
# Boxplots with flags 
generate_fixed_IQR_metric_boxplots(iqr.plt = iqrdf, 
                                   type = "-MFI", 
                                   markers = markers_1, 
                                   colours = batchColours, 
                                   control_list = controls, 
                                   batch_list = batches, 
                                   output_dir = output.dir, 
                                   width = 11, 
                                   height = 10,
                                   ctrl_labs = controls,
                                   y_limits = "auto",
                                   seed = 450)

# +MFI
# Boxplots with flags 
generate_fixed_IQR_metric_boxplots(iqr.plt = iqrdf, 
                                   type = "+MFI", 
                                   markers = markers_1, 
                                   colours = batchColours, 
                                   control_list = controls, 
                                   batch_list = batches, 
                                   output_dir = output.dir, 
                                   width = 11, 
                                   height = 10,
                                   ctrl_labs = controls,
                                   y_limits = "auto",
                                   seed = 450)

# %pos
# Boxplots with flags 
generate_fixed_IQR_metric_boxplots(iqr.plt = iqrdf, 
                                   type = "%pos", 
                                   markers = markers_1, 
                                   colours = batchColours, 
                                   control_list = controls, 
                                   batch_list = batches, 
                                   output_dir = output.dir, 
                                   width = 11, 
                                   height = 10, 
                                   ctrl_labs = controls,
                                   y_limits = "auto",
                                   seed = 450)

# Lastly, you may want to remove some objects that are loaded into R to save space before moving on to the next step.

rm(unimodal_and_bimodal_markers)
rm(auto_cutoffs)
rm(updated_auto_cutoffs)
rm(ffdf_sampled)
rm(nMFIdf)
rm(pMFIdf)
rm(posdf)

# ======================================================================
# Step 3: Earth Mover’s Distance (EMD)-Based Assessment
# ======================================================================

# This section handles the EMD-based assessment. 

# First, the data needs to be prepared prior to calculating the EMD.

cdf <- data.frame(reshape_df(df = ffdf, 
                             sample_col_name = "control", 
                             batch_col_name = "batch", 
                             markers = markers_1))

# Next, get the pairwise EMD between markers for each sample.

# To speed up this process, subsampling can be done prior to computing EMD. num is the number of cells to subsample per marker per sample, while seed is set to make all subsampling reproducible. If you want the subsampling to be random each time you run the function, set seed = NULL.

# Get EMDs
# num and seed can be set to a number
pairwise_EMDs_per_marker <- EMD_list(df = cdf, 
                                     samps = sample_ids, 
                                     markers = markers_1, 
                                     batch_list = batches, 
                                     control_list = controls, 
                                     num = 20000,
                                     seed = 450)

# Optionally, save EMD matrix list
# saveRDS(pairwise_EMDs_per_marker, 
#        file = file.path(output.dir,"pairwise_EMDs_across_all_markers_and_controls.rds"))

# If you later want to read that file back in, use:
# readRDS(file = file.path(output.dir,"pairwise_EMDs_across_all_markers_and_controls.rds"))

# The code below uses heatmaps to visualize EMDs of pairs of batches for every marker and control.

# The default threshold is 5, but this can be changed based on your data. 

# Define EMD threshold here
emd_threshold <- 5

# Axis size, plot width, and plot height are customizable
generate_EMD_heatmaps(emds_list = pairwise_EMDs_per_marker, 
                      batch_list = batches, 
                      markers = markers_1, 
                      controls = controls, 
                      threshold = emd_threshold, 
                      output_dir = output.dir, 
                      axis_size = 18, 
                      width = 11.5, 
                      height = 8)

# Next, create boxplots of pairwise EMDs ordered by median EMDs per batch. 

# In the produced figures, the red horizontal line represents the user-defined threshold. The box itself is outlined in red when the batch's median EMD is greater then the threshold.

# Plot width and plot height are customizable
generate_EMD_boxplots(distdf = pairwise_EMDs_per_marker, 
                      batch_list = batches, 
                      markers = markers_1, 
                      control_list = controls, 
                      batch_colours = batchColours, 
                      threshold = emd_threshold, 
                      output_dir = output.dir,
                      width = 11.5, 
                      height = 8)

# Finally, create a results dataframe of median EMD values and the corresponding generated flags. 

emdf <- data.frame(emd_med_vals(em_dists = pairwise_EMDs_per_marker, 
                              markers = markers_1, 
                              samps = sample_ids, 
                              control_list = controls, 
                              batch_list = batches, 
                              threshold = emd_threshold, 
                              output_dir = output.dir, 
                              num = 20000))

# Again, you may want to remove some objects that are loaded into R to save space before moving on to the next step.

rm(cdf)
rm(pairwise_EMDs_per_marker)

# ======================================================================
# Step 4: Summary of Results (for Steps 2 and 3)
# ======================================================================

# Once Steps 2 and 3 have been completed, a summary of flags from the IQR-based and EMD-based metrics can be created. The following function creates both a table and a heatmap depicting the flags.

# This plot is ranked by most to least flagged batches and markers across metrics. Flagged boxes appear in red, non-flagged boxes in grey, and NA values are shown in white.

# If your data has a large amount of markers and/or batches, you may want to split the data into multiple plots to improve readability. The markers_per_plot and batches_per_plot parameters can be used in these cases. 

# markers_per_plot and batches_per_plot indicate the maximum number of markers and batches to include in each plot, respectively. For example if you set markers_per_plot = 20 and the panel includes 49 markers, then 3 images will be created. Plot 1 and 2 will have 20 markers each and Plot 3 will have 9 markers. If you set them at NULL, then all markers/batches will be included in the outputted plot. 

# Additionally, min_batch_effect can be used to filter out batch-control combinations that contain few flags. A value of 0 means no filtering. Every flagged batch, even with very small effect sizes, will be shown. A value of 3 means all batch-controls that have less than 3 flags are removed. This can be done to focus the heatmap on batch-controls that have increased risk of batch effects. Likewise, min_marker_effect can be used to filter out markers in the same manner. A value of 5 would mean markers with less than 5 flags across the IQR-based and EMD-based metrics are filtered out. 

# Large heatmap with all the flags
generate_ranked_flagged_hmap(emd_df = emdf,
                             iqr_df = iqrdf,
                             batch_list = batches, 
                             controls = controls, 
                             markers = markers_1,
                             batch_cols = batchColours, 
                             output_dir = output.dir,
                             markers_per_plot = NULL,
                             batches_per_plot = NULL,
                             min_batch_effect = 0,
                             min_marker_effect = 0) 

# ======================================================================
# Step 5: Unsupervised Clustering-Based Assessment
# ======================================================================

# The final metric is a clustering-based assessment of the dataset. This can be used to check if problematic batches impact downstream cell clusters.

# First, choose the markers you want to include for clustering. That can either be done by:
# - Selecting the markers in a pop-up window, or
# - If a 'status' column was included in the marker file, then filtering by values in that column

# Select the markers used for clustering in this pop-up window
# Click 'Done' once markers are chosen
pheno_markers <- selectMarkerList(marker_tbl, column_name = "PnS")

# Alternatively, use this if you want to select all 'phenotype/lineage' markers
# pheno_markers <- as.character(marker_tbl[which(marker_tbl$status=="phenotype"),"PnS"])

# Once the markers have been selected, set some parameters. This should include the number of cells to randomly downsample to, the number of metaclusters to make in the unsupervised clustering algorithm, and thresholds 1 and 2. 

# The thresholds should be a numeric multiplier for the expected proportion threshold (e.g., 2 for 2x expected). The default for threshold 1 is 2 and the default for threshold 2 is 1.5.

# To make the cells that are subsampled reproducible, set the seed to a number. If you want the subsampling to be random each time you run the function, set seed <- NULL.

# Number of cells to downsample
num <- 50000

# Number of clusters to make in FlowSOM
meta_cluster_num <- 9

# Seed for reproducability 
seed <- 150

# Set threshold 1
# If the contribution of cells from a batch in a given cluster is
# greater than [threshold1] times the expected proportion based on the total number of batches, these batches are flagged
set_threshold1 <- 2.0

# Set threshold 2
# If a cluster had 1 or more batches that were flagged in part 1 (where threshold 1 was used), then part 2 is done for that cluster(s)
# Here, a batch is flagged if its contribution is greater than [threshold2] times the expected proportion in that cluster for the given control 
set_threshold2 <- 1.5

# Next, run the FlowSOM unsupervised clustering algorithm, which will put cells into clusters based on their marker expression values. 

# FlowSOM is a popular clustering algorithm made for cytometry data. Read more about it here: https://doi.org/10.1002/cyto.a.22625 

# This function also uses ConsensusClusterPlus output plots to aid in determining the optimal number of clusters (k). Read more about it here: https://pubmed.ncbi.nlm.nih.gov/20427518/ 

# The ConsensusClusterPlus plots are saved in the folder that begins with cluster_plots. 

# Uses FlowSOM unsupervised clustering to get put each cell into a cluster based on marker expression values
filter_df <- getFlowSOM_clusters(data_f = ffdf, 
                                 markers = pheno_markers, 
                                 seed = seed,
                                 meta_cluster_num = meta_cluster_num, 
                                 samps = sample_ids, 
                                 num = num, 
                                 output_dir = output.dir)

# Gets the cluster for each cell from the FlowSOM clustering above
cell_cluster_vector <- filter_df[,"cluster"]

# Now a bar plot showing cluster sizes (number and percent of cells) can be created.

# The horizontal red line represents 0.5% of all analyzed cell and can be used to indicate the very small clusters that are below this line.

# Generate cluster sizes plot
# Axis size, plot width, and plot height are customizable
clusterSizes_plot(cluster_df = filter_df, 
                  seed = seed, 
                  meta_cluster_num = meta_cluster_num,
                  output_dir = output.dir, 
                  width = 11.5, 
                  height = 8)

# The scaled median marker expression values for the clusters can be shown in a heatmap.

# Generates a heatmap of scaled median marker expression for each cluster, ordered by cluster number
# axis_size: Numeric value specifying font size for axis labels
# num_size: Numeric value specifying font size for numbers displayed in the heatmap cells
# width, height: Numeric values specifying the dimensions of the output plot in inches
clustering_heatmap(cluster_df = filter_df, 
                   markers = pheno_markers, 
                   cluster_colrs = clustColors, 
                   cell_cluster_vector = cell_cluster_vector, 
                   seed = seed, 
                   meta_cluster_num = meta_cluster_num, 
                   axis_size = 15, 
                   num_size = 8.4, 
                   output_dir = output.dir, 
                   width = 11.5, 
                   height = 8)

# Next, a bar plot and text file is generated that depict batch proportions per cluster and indicate batches that exceed the threshold. This uses threshold 1, which was set above.

# In the bar plot, the batches that are above the threshold are designated with a bold, black outline. 

# Additionally, the batches that are above the threshold are designated in the console and in the CSV file (TRUE in the above_threshold column).

# Bar plot of batch proportions per cluster; list of batches > 2x expected proportion saved

# Batch colors
colr_batch <- batchColours[1:length(batches)]

# Generates a bar plot of batch proportions for each cluster, outlining in bold batches that exceed a specified threshold
# width, height: Numeric values specifying the dimensions of the output plot in inches
# axis_size: Numeric value specifying font size for axis labels
proportion_of_batches_per_cluster(cluster_df = filter_df, 
                                  cell_cluster_vector = cell_cluster_vector, 
                                  set_threshold1 = set_threshold1, 
                                  batch_list = batches, 
                                  batch_colours = batchColours, 
                                  seed = seed, 
                                  output_dir = output.dir, 
                                  width = 11.5, 
                                  height = 8, 
                                  axis_size = 18)

# Generates a table of batch proportions for each cluster, with a column indicating batches that exceed a specified threshold
batch_props <- proportion_of_batches_per_cluster_table(cluster_df = filter_df,
                                                        cell_cluster_vector = cell_cluster_vector,
                                                        set_threshold1 = set_threshold1,
                                                        batch_list = batches, 
                                                        seed = seed, 
                                                        output_dir = output.dir)

# Now, the second assessment that examines batch proportions per cluster across controls is done. This uses threshold 2, designated above. Boxplots and a text file of results are generated with the code below. 

# In the plot, flags that signify a batch above the threshold are marked with the name of that batch next to the flagged point. In contrast, a point without a label indicates it was not flagged. Only clusters that have at least one batch above the first threshold, threshold 1, are included in the outputted plot. 

# Additionally, the batches that are above the threshold are designated in the console and in the CSV file (TRUE in the above_threshold column).

# If the code above (that used threshold 1) revealed no batches impacting clustering, skip this step and go directly to plotting the UMAP of clusters (next section).

# A seed can be set to ensure the jitter pattern (offset location of the points) in the plot is the same each time. If you want the jitter to be random each time you run the function, set seed = NULL.

# Boxplots of batch proportions per cluster across controls

# Generates boxplots of batch proportions across control samples for clusters impacted by batch effects, highlighting batches that exceed a specified threshold
# axis_size: Numeric value specifying font size for axis labels
# width, height: Numeric values specifying the dimensions of the output plot in inches
proportion_of_batches_across_control_samples(cluster_df = filter_df,
                                             cell_cluster_vector = cell_cluster_vector,
                                             batch_props_df = batch_props,
                                             set_threshold2 = set_threshold2, 
                                             batch_list = batches, 
                                             control_list = controls, 
                                             batch_colours = batchColours, 
                                             output_dir = output.dir, 
                                             seed = seed, 
                                             axis_size = 12, 
                                             width = 11.5, 
                                             height = 8)

# Generates a table of batch proportions across control samples for clusters that were over threshold 1, with a column indicating batches that exceed a specified threshold 2
bath_controls_props <- proportion_of_batches_across_control_samples_table(cluster_df = filter_df, 
                                                                          cell_cluster_vector = cell_cluster_vector, 
                                                                          set_threshold2 = set_threshold2, 
                                                                          batch_props_df = batch_props, 
                                                                          batch_list = batches, control_list = controls, 
                                                                          output_dir = output.dir, 
                                                                          seed = seed)

# Finally, a UMAP that colors cells by batch is made as a final visualization of the clustering assessment.

# Set num to the number of cells to downsample to per clusters. If a cluster has < num cells, all cells will be included. To make the cells that are subsampled reproducible, set the seed to a number. If you want the subsampling to be random each time you run the function, set seed = NULL.

# UMAP of clusters annotated by batch colours

# Axis size, plot width, plot height, seed, and number of cells (num) are customizable
generate_cluster_UMAP(df = filter_df, 
                      marker_list = pheno_markers, 
                      batch_list = batches, 
                      batch_colours = batchColours, 
                      seed = seed,
                      num = 2000,
                      output_dir = output.dir, 
                      axis_size = 24, 
                      width = 12, 
                      height = 12)

# ======================================================================
# Step 6: Summary of Results (for Steps 2, 3, and 5)
# ======================================================================

# Once Steps 2, 3, and 5 are completed, all IQR-based, EMD-based, and cluster-based flags have been generated. As a result, an overall summary of flags from all assessments can be created.

# The following function creates a heatmap depicting all flags. This conveys the same information from Step 4, with the addition of the clustering results.

# This plot is ranked by most to least flagged batches and markers across metrics. Flagged boxes appear in red, non-flagged boxes in grey, and NA values are shown in white.

# If your data has a large amount of markers and/or batches, you may want to split the data into multiple plots to improve readability. The markers_per_plot and batches_per_plot parameters can be used in these cases. 

# markers_per_plot and batches_per_plot indicate the maximum number of markers and batches to include in each plot, respectively. For example if you set markers_per_plot = 20 and the panel includes 49 markers, then 3 images will be created. Plot 1 and 2 will have 20 markers each and Plot 3 will have 9 markers. If you set them at NULL, then all markers/batches will be included in the outputted plot. 

# Additionally, min_batch_effect can be used to filter out batch-control combinations that contain few flags. A value of 0 means no filtering. Every flagged batch, even with very small effect sizes, will be shown. A value of 3 means all batch-controls that have less than 3 flags are removed. This can be done to focus the heatmap on batch-controls that have increased risk of batch effects. Likewise, min_marker_effect can be used to filter out markers in the same manner. A value of 5 would mean markers with less than 5 flags across the IQR-based and EMD-based metrics are filtered out. 

# Ranked, summarized heatmap of all assessment metrics

# Large heatmap with all the flags
generate_ranked_flagged_hmap_all(emd_df = emdf,
                                 iqr_df = iqrdf,
                                 cluster_df = filter_df, 
                                 cell_cluster_vector = cell_cluster_vector,
                                 set_threshold1 = set_threshold1,
                                 set_threshold2 = set_threshold2,
                                 batch_list = batches, 
                                 controls = controls, 
                                 markers = markers_1,
                                 batch_cols = batchColours, 
                                 output_dir = output.dir,
                                 markers_per_plot = NULL,
                                 batches_per_plot = NULL,
                                 min_batch_effect = 0,
                                 min_marker_effect = 0) 
