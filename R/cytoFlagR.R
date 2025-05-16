
###############################################################################################
####################################~~~~ Run cytoFlagR ~~~~####################################
###############################################################################################
##### Installation of dependencies #####
### load this function to check for missing packages and install them
source("package_installer.R")
### from CRAN project
("https://cran.r-project.org/")

### packages needed
requiredPackages<-c("dplyr","scales","tidyr","reshape2","readr","matrixStats","readxl",
                    "ggplot2","ggpubr","ggridges","MASS","RColorBrewer","cowplot",
                    "ggrepel","emdist","circlize","gridExtra","stats",
                    "LaplacesDemon","pheatmap","umap","progress","crayon","patchwork",
                    "ggpmisc","viridis","tidyverse","shiny","shinyjs","DT","bslib","shinyjs",
                    "cluster","rstudioapi","grid","shinycssloaders","shinyWidgets")
### install
package_installer(requiredPackages)

### check if the packages can be loaded
lapply(requiredPackages, require, character.only = TRUE)

### from Bioconductor, install Bioconductor Manager (BiocManager) and then install Bioconductor packages
("https://www.bioconductor.org/install/")

checkBiocManager_install()

### packages needed
required_BioconductorPackages<-c("flowCore","FlowSOM","ComplexHeatmap","limma","ConsensusClusterPlus")

### install
BioC_package_installer(required_BioconductorPackages)

### check if the packages can be loaded
lapply(required_BioconductorPackages, require, character.only = TRUE)

##### load all packages #####

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
suppressPackageStartupMessages(library(shinycssloaders))
suppressPackageStartupMessages(library(bslib))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(rstudioapi))
suppressPackageStartupMessages(library(grid))

##### source all required functions #####
source("cytoFlagR_R_functions.R")
source("auto_threshold_func.R")
source("auto_threshold_func_cytof.R")
source("editMarkerTable.R")
source("generate_marker_expression_densities.R")
source("selectMarkerList.R")
source("referenceMarkerSelection.R")

##### Set directory to read files #####
### folder where cytoFlagR and related files are stored
work.dir<-getwd()
### folder where all FCS files, marker information file, and fcs meta-data file exist
data.dir<-rstudioapi::selectDirectory(caption = "Select FCS file folder", label = "Select", 
                                      path = getActiveProject()) 

##### Load fcs meta-data file, marker information file #####
### fcs information file
fcs.info<-data.frame(readr::read_delim(file.path(data.dir,"fcs_file_info.csv")))

### marker and fluorophore/metal tag information file
marker_tbl<-data.frame(readr::read_delim(file.path(data.dir,"panel_markers.tsv")))

### store the tag and marker lists as well as their column names
tag_list<-marker_tbl[,"PnN"]
markers<-marker_tbl[,"PnS"]

############################################
#                                          #
#                  STEP 1                  #
#                                          #
############################################

### set your own cofactor depending on your dataset type
cf<-6000
### this is the folder where your outputs are stored
output.dir<-file.path(paste0("output_",format(Sys.time(), "%Y.%m.%d_%H%M%S")))
dir.create(output.dir,recursive = TRUE, showWarnings = FALSE) 

### create a transformed dataframe of all FCS files

ffdf<-create_transformed_dataframe(cf = cf, fcs_info = fcs.info, tag_list = tag_list, 
                                   marker_list = markers, marker_tbl = marker_tbl, 
                                   tag_colm = "PnN", marker_colm = "PnS", 
                                   output_dir = output.dir, file_dir = data.dir)

### you may input your own character list if you want the names to be ordered. For example: controls<-as.character("control_1","control_2","control_3")
controls<-as.character(unique(fcs.info$control)) 

### you may input your own character list if you want the names to be ordered. For example: batches<-as.character("batch_1","batch_2",...,"batch_10")
batches<-as.character(unique(fcs.info$batch)) 

sample_ids<-as.character(unique(fcs.info$sample_id))

#### number of cells per sample ####
#### axis size, plot width, and plot height are customizable
generate_number_of_cells_plot(df = ffdf, batch_colm = "batch", control_colm = "control", 
                              batch_list = batches, control_list = controls, output_dir = output.dir, 
                              axis_size = 20, width = 22.5, height = 14.6)

#### MDS plot across samples #####
#### axis size, plot width, and plot height are customizable
generate_MDS_plot(df = ffdf, batch_colm = "batch", control_colm = "control", 
                  marker_list = markers,sample_colours = sample_colours, output_dir = output.dir, 
                  axis_size = 15, width = 6.2, height = 6)

#### UMAP of controls annotated by batch colours ####
#### axis size, plot width, and plot height are customizable
generate_UMAP(df = ffdf, batch_colm = "batch", control_colm = "control", 
              marker_list = markers, batch_list = batches, control_list = controls, 
              sample_ids = sample_ids, batch_colours = batchColours, output_dir = output.dir, 
              axis_size = 11, width = 5, height = 5)

############################################
#                                          #
#                  STEP 2                  #
#                                          #
############################################

### Use this to interface to deselect markers you don't want to evaluate and 
### save the rest to calculate marker cutoffs:
markers_1<-selectMarkerList(marker_tbl, column_name = "PnS")

#### dataset type: change to "flow" if using spectral or conventional flow data
#### subsample: change to "yes" if you want cytoFlagR to sample your data 
### to calculate the marker thresholds
#### calculate the automated marker thresholds ####
auto_cutoffs<-data.frame(automated_threshold_DF(df = ffdf, markers = markers_1, 
                                                control_list = controls, batch_list = batches, 
                                                batch_colm = "batch", subsample = "no", 
                                                dataset_type = "mass",output_dir = output.dir))

### plot and save densities to pdf
#### axis size, plot width, and plot height are customizable
#### generate marker expression densities with automated thresholds ####
write_density_plot_per_marker(marker_list = markers_1, df = ffdf, batch_list = batches, 
                              control_list = controls, auto_cutoffs = auto_cutoffs, 
                              batch_colm = "batch", control_colm = "control", 
                              output_dir = output.dir, wd = 12, ht = 11, axis_size = 14,
                              file_name = "panel_marker_expression_densities")

#### launch interface to update your marker thresholds and save your automated thresholds ####

editMarkerTable(auto_cutoffs)

### new cutoffs saved as a file and assigned as: updated_auto_cutoffs

#### axis size, plot width, and plot height are customizable
#### generate marker expression densities with updated automated thresholds ####
write_density_plot_per_marker(marker_list = markers_1, df = ffdf, batch_list = batches, 
                              control_list = controls, auto_cutoffs = updated_auto_cutoffs, 
                              batch_colm = "batch", control_colm = "control", 
                              output_dir = output.dir, wd = 12, ht = 11, axis_size = 14, 
                              file_name = "marker_expression_densities_with_updated_cutoffs")

### choose a reference marker using interface:
ref_marker<-referenceMarkerSelection(marker_tbl, column_name = "PnS")

#### generate marker expression density biaxial dotplots with updated automated thresholds ####
#### axis size, plot width, and plot height are customizable
draw_biaxial_dotplots(df = ffdf, ref_marker = ref_marker, marker_list = markers_1, 
                      control_list = controls, batch_list = batches, batch_col_name = "batch", 
                      control_col_name = "control", samps = sample_ids, 
                      auto_cutoffs = updated_auto_cutoffs, wd = 5000, 
                      ht = 6000, axis_size = 18, output_dir = output.dir)

#### Get negative populations, positive populations, and 
### percent positive cells of all markers and all controls ####
nMFIdf<-data.frame(neg_mfi_df_all(df = ffdf, markers = markers_1, cutf = updated_auto_cutoffs, 
                                  sample_col_name = "control", batch_col_name = "batch", 
                                  min_cells = TRUE))
pMFIdf<-data.frame(pos_mfi_df_all(df = ffdf, markers = markers_1, cutf = updated_auto_cutoffs, 
                                  sample_col_name = "control", batch_col_name = "batch", 
                                  min_cells = TRUE))
posdf<-data.frame(percent_pos_df_all(df = ffdf, markers = markers_1, cutf = updated_auto_cutoffs, 
                                     ctrls = controls, sample_col_name = "control", 
                                     batch_col_name = "batch"))

#### Generate flags of IQR assessment metrics to find problem batches and markers ####
iqrdf<-iqr_dataframe(neg_df = nMFIdf,pos_df = pMFIdf, percent_pos_df = posdf, 
                     markers = markers_1, sample_col_name = "control", 
                     batch_col_name = "batch", output_dir = output.dir)

### optional: you can change the names of controls to 1, 2, 3, 4, etc.
control_labels<-as.character(seq(1,length(controls)))
names(control_labels)<-controls ### assign control sample name to labels

### if using different labels for controls, use ctrl_labs = control_labels, 
### else use ctrl_labs = controls while generating IQR boxplots below:

#### Generate boxplot of flagged batches for each marker (variable axes) ####
#### axis size, plot width, and plot height are customizable
#### boxplots with flags for -MFI metric ####

generate_IQR_metric_boxplots(iqr.plt = iqrdf, type = "-MFI", markers = markers_1, 
                             ctrl_labs = controls, colours = batchColours, 
                             batch_list = batches, control_list = controls, output_dir = output.dir, 
                             width = 11, height = 10)
#### boxplots with flags for +MFI metric ####
generate_IQR_metric_boxplots(iqr.plt = iqrdf, type = "+MFI", markers = markers_1, 
                             ctrl_labs = controls, colours = batchColours, 
                             batch_list = batches, control_list = controls, output_dir = output.dir, 
                             width = 11, height = 10)
#### boxplots with flags for %pos metric ####
generate_IQR_metric_boxplots(iqr.plt = iqrdf, type = "%pos", markers = markers_1, 
                             ctrl_labs = controls, colours = batchColours, 
                             batch_list = batches, control_list = controls, output_dir = output.dir, 
                             width = 11, height = 10)

#### Generate boxplot of flagged batches for each marker (fixed axes) ####
#### axis size, plot width, and plot height are customizable
#### boxplots with flags for -MFI metric ####
generate_fixed_IQR_metric_boxplots(iqr.plt = iqrdf, type = "-MFI", markers = markers_1, 
                                   ctrl_labs = controls,colours = batchColours, 
                                   batch_list = batches, control_list = controls, 
                                   output_dir = output.dir, width = 11, height = 10)
#### boxplots with flags for +MFI metric ####
generate_fixed_IQR_metric_boxplots(iqr.plt = iqrdf, type = "+MFI", markers = markers_1, 
                                   ctrl_labs = controls,colours = batchColours, 
                                   batch_list = batches, control_list = controls, 
                                   output_dir = output.dir, width = 11, height = 10)
#### boxplots with flags for %pos metric ####
generate_fixed_IQR_metric_boxplots(iqr.plt = iqrdf, type = "%pos", markers = markers_1, 
                                   ctrl_labs = controls,colours = batchColours, 
                                   batch_list = batches, control_list = controls, 
                                   output_dir = output.dir, width = 11, height = 10)


############################################
#                                          #
#                  STEP 3                  #
#                                          #
############################################
#### EMD-based assessment ####
#### prepare data for EMD calculations
cdf<-data.frame(reshape_df(df = ffdf, sample_col_name = "control", 
                           batch_col_name = "batch", markers = markers_1))

#### calculate pairwise EMD per marker for each control sample across batches
pairwise_EMDs_per_marker<-EMD_list(df = cdf, samps = sample_ids, markers = markers_1, 
                                   batch_list = batches, control_list = controls, num = 20000)

#### optional: save EMD matrix list
saveRDS(pairwise_EMDs_per_marker, 
        file = file.path(output.dir,"pairwise_EMDs_across_all_markers_and_controls.rds"))

#### you would read it in with: 
### readRDS(file = file.path(output.dir,"pairwise_EMDs_across_all_markers_and_controls.rds"))

### define EMD threshold here:
emd_threshold<-3

#### Pairwise EMD heatmaps ####
#### axis size, plot width, and plot height are customizable
generate_EMD_heatmaps(emds_list = pairwise_EMDs_per_marker, batch_list = batches, 
                      markers = markers_1, controls = controls, threshold = emd_threshold, 
                      output_dir = output.dir, axis_size = 18, width = 13, height = 12)

#### ordered pairwise EMD boxplots ####
#### axis size, plot width, and plot height are customizable
generate_EMD_boxplots(distdf = pairwise_EMDs_per_marker, batch_list = batches, 
                      markers = markers_1, control_list = controls, batch_colours = batchColours, 
                      threshold = emd_threshold, output_dir = output.dir, axis_size = 18, 
                      width = 11.5, height = 10.5)

#### Dataframe of median EMD values per marker, per control, per batch and flags are generated ####
emdf<-data.frame(emd_med_vals(em_dists = pairwise_EMDs_per_marker, markers = markers_1, 
                              samps = sample_ids, control_list = controls, 
                              batch_list = batches, threshold = emd_threshold, 
                              output_dir = output.dir, num = 20000))

#### Remove re-shaped dataframe after EMD assessment is complete
rm(cdf)

############################################
#                                          #
#                  STEP 4                  #
#                                          #
############################################

##### Ranked, summarized heatmap of flags across all assessment metrics #####
generate_ranked_flagged_hmap(emd_df = emdf, iqr_df = iqrdf, batch_list = batches, 
                             controls = controls, markers = markers_1, 
                             batch_cols = batchColours, output_dir = output.dir)

############################################
#                                          #
#                  STEP 5                  #
#                                          #
############################################

##### Clustering-based assessment of dataset #####
#### Launch the interface to select markers to use for clustering, close the app to save list
### get the phenotype/lineage markers used for clustering
pheno_markers<-selectMarkerList(marker_tbl, column_name = "PnS")

##### Generate unique clusters using FlowSOM #####
### default num = 20000, you may change this value along with meta_cluster_num and seed
meta_cluster_num<-20
seed<-250
filter_df<-getFlowSOM_clusters(data_f = ffdf, markers = pheno_markers, seed = seed, 
                               meta_cluster_num = meta_cluster_num, samps = sample_ids, num = 20000, 
                               output_dir = output.dir)

### get cluster vector
cell_cluster_vector<-filter_df[,"cluster"]

##### Generate barplot of cluster sizes #####
#### axis size, plot width, and plot height are customizable
clusterSizes_plot(cluster_df = filter_df, seed = seed, meta_cluster_num = meta_cluster_num, 
                  output_dir = output.dir, width = 5, height = 4)

##### Generate heatmap of marker expressions across clusters #####
#### axis size, plot width, and plot height are customizable
clustering_heatmap(cluster_df = filter_df, markers = pheno_markers, cluster_colrs = clustColors, 
                   cell_cluster_vector = cell_cluster_vector, seed = seed, 
                   meta_cluster_num = meta_cluster_num, axis_size = 15, num_size = 8.4, 
                   output_dir = output.dir,width = 11.5, height = 8)

##### Barplot of batch proportions per cluster to determine impact of batches on clustering #####

#### axis size, plot width, and plot height are customizable
proportion_of_batches_per_cluster(cluster_df = filter_df, cell_cluster_vector = cell_cluster_vector, 
                                  batch_list = batches, batch_colours = batchColours, seed = seed, 
                                  output_dir = output.dir, width = 3000, height = 1900, 
                                  axis_size = 12)

#### save the list of batches that impact clustering
batch_props<-proportion_of_batches_per_cluster_table(cluster_df = filter_df, 
                                                     cell_cluster_vector = cell_cluster_vector, 
                                                     batch_list = batches, seed = seed, 
                                                     output_dir = output.dir)

##### Boxplots of batch proportions per cluster across controls #####
#### axis size, plot width, and plot height are customizable
proportion_of_batches_across_control_samples(cluster_df = filter_df, 
                                             cell_cluster_vector = cell_cluster_vector, 
                                             batch_props_df = batch_props, batch_list = batches, 
                                             control_list = controls, batch_colours = batchColours, 
                                             output_dir = output.dir, seed = seed,
                                             axis_size = 12, width = 3800, height = 5800)

##### UMAP of clusters annotated by batch colours #####
#### axis size, plot width, and plot height are customizable
generate_cluster_UMAP(df = filter_df, marker_list = pheno_markers, batch_list = batches, 
                      batch_colours = batchColours, seed = seed, output_dir = output.dir, 
                      axis_size = 15, width = 5, height = 5)

###############################################################################################

