#!/usr/bin/env R

###############################################################################################
#############################~~~~ cytoFlagR functions ~~~~#############################
###############################################################################################

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

###############################################################################################

sample_colours<-c("#8ad6cc","#47857c","#f6b250","#e35988",
                  RColorBrewer::brewer.pal(name = "Paired",n = 12),
                  RColorBrewer::brewer.pal(name = "Dark2",n = 8))

batchColours<-c("#5C8E39","#0072B2","#FF69B4","#E1A100","#19E7E7","#5E4FA2",
                "#35B779","#6A3D9A","#FF700A","#787878","#C24141","#3c096c",
                "#FF5349","#00FF00","#0000FF","#83d632","#482878","#FF00FF",
                "#008080","#FFC033","#c3acf4","#8AD6CC","#006400","#D88A45",
                "#4169E1","#FF69B4","#9932CC","#6495ED","#91B0B9","#AE3A3A",
                "#3288BD","#B91F48","#486BAF","#88CFA4","#9E0142","#4CA5B1",
                "#66C2A5","#FF9500","#8338ec","#ABDDA4","#003049","#669BBC",
                "#C8E99E","#ffbe0b","#588157","#f4a261","#00bbf9","#ff4d6d")

clustColors<-c("#ED553B","#20639B","#F6d55C","#3CAEA3","#FEA8ED",
               "#A8FAA4","#FCFC97","#B6056A","#0D44BA","#880ED4",
               "#DE8297","#B8E0FF","#FF1919","#CE8CF8","#69794E",
               "#71A6B8","#B20245","#236A3A","#C576F6","#FCB254",
               "#559050","#CA9BB0","#C78D40","#9EB994","#24BAB9",
               "#EE4483","#F25A25","#93C47D","#C2272D","#38761D",
               "#C898B1","#98C8C7","#C8AF98","#668BAD","#98B1C8",
               "#1CD103","#4C03BD","#E6E50E","#C500CB","#00ADE3",
               "#8ad6cc","#47857c","#f6b250","#e35988","#482878",
               "#5C8E39","#0072B2","#FF69B4","#E1A100","#787878")


###############################################################################################
############################~~~~ cytoFlagR data transformation ~~~~############################
##########################~~~~ and marker threshold calculation ~~~~###########################
###############################################################################################

#### perform arcsinh transformation of fcs files
create_transformed_dataframe<-function(cf,fcs_info,tag_list,marker_list,
                                       marker_tbl,tag_colm,marker_colm,file_dir,output_dir){
  
  colrSet<-crayon::make_style("#539DDD")
  
  ### create progress bar
  pBar<-progress::progress_bar$new(
    
    format = paste0(colrSet("Transforming FCS files"),
                    " [",colrSet(":bar"),"] ",
                    colrSet(":percent")," | File: :current/:total ",
                    "| Elapsed: :elapsed | ETA: :eta"),
    total = nrow(fcs_info),
    clear = FALSE,
    width = 80,
    complete = "=",
    incomplete = "-"
  )
  
  colnames(fcs_info)<-tolower(colnames(fcs_info))
  
  fcs_info[,"path"]<-file.path(file_dir, fcs_info[,"file"])
  
  tmp1<-flowCore::read.FCS(fcs_info[1,"path"], transformation = FALSE, truncate_max_range = FALSE)
  flowcolnames<-flowCore::colnames(tmp1)
  
  tag_match<-setdiff(tag_list,flowcolnames)
  mark_match<-setdiff(marker_list,flowcolnames)
  
  ### get sample id column if one doesn't exist
  if(is.null(fcs_info[,"sample_id"])) {
    fcs_info[,"sample_id"]<-paste0(fcs_info[,"batch"],"_",fcs_info[,"control"])
  }
  
  expr_df<-data.frame()
  
  if(all(tag_list %in% flowcolnames) == TRUE) {
    for (i in 1:nrow(fcs_info)) {
      
      filename<-fcs_info[i,"file"]
      
      ### update progress bar
      pBar$tick(tokens = list(current = i))
      
      myfile<-flowCore::read.flowSet(fcs_info[i,"path"], 
                                     transformation = FALSE, 
                                     truncate_max_range = FALSE)
      myfile<-flowCore::fsApply(myfile, function(x, cofactor=asinh_scale){
        expr<-flowCore::exprs(x)
        expr<-asinh(expr[,tag_list]/cf) ### use list of fluorphores/metal tags 
        flowCore::exprs(x)<-expr
        x
      })
      tmp_mat<-flowCore::fsApply(myfile, flowCore::exprs)
      tmp_mat<-data.frame(tmp_mat)
      tmp_mat$file<-filename
      expr_df<-rbind(expr_df,tmp_mat)
    }
    colnames(expr_df)<-c(tag_list,"file")
    #### match metal/fluorophore names to panel file column
    match_names<-match(names(expr_df),marker_tbl[,tag_colm])
    
    #### rename to protein channel names
    names(expr_df)[!is.na(match_names)]<-marker_tbl[,marker_colm][na.omit(match_names)]
    expr_df<-data.frame(expr_df)
    
  }
  else if(all(marker_list %in% flowcolnames) == TRUE){
    for (i in 1:nrow(fcs_info)) {
      
      filename<-fcs_info[i,"file"]
      
      ### update progress bar
      pBar$tick(tokens = list(current = i))
      
      myfile<-flowCore::read.flowSet(fcs_info[i,"path"], 
                                     transformation = FALSE, 
                                     truncate_max_range = FALSE)
      myfile<-flowCore::fsApply(myfile, function(x, cofactor=asinh_scale){
        expr<-flowCore::exprs(x)
        expr<-asinh(expr[,marker_list]/cf) ### use list of protein channels 
        flowCore::exprs(x)<-expr
        x
      })
      tmp_mat<-flowCore::fsApply(myfile, flowCore::exprs)
      tmp_mat<-data.frame(tmp_mat)
      tmp_mat$file<-filename
      expr_df<-rbind(expr_df,tmp_mat)
    }
    colnames(expr_df)<-c(marker_list,"file")
  }
  else {
    stop("names in fcs files do not match with the ones in marker table")
    if(!is.null(tag_match)) {
      mm<-match(marker_tbl$PnN,tag_match)
      cat("marker =", marker_tbl[,"PnS"][na.omit(mm)], ", tag =", marker_tbl[,"PnN"][na.omit(mm)],"\n")
    }
    else if(!is.null(mark_match)){
      mm<-match(marker_tbl$PnS,mark_match)
      cat("marker =", marker_tbl[,"PnS"][na.omit(mm)], ", tag =", marker_tbl[,"PnN"][na.omit(mm)],"\n")
    }
  }
  
  #### match file column to get batch, control, sample id information
  match_cols<-match(expr_df$file,fcs_info$file)
  expr_df$sample_id<-fcs_info$sample_id[match_cols]
  expr_df$batch<-fcs_info$batch[match_cols]
  expr_df$control<-fcs_info$control[match_cols]
  
  file_path <- file.path(output_dir, "transformed_dataframe_of_control_samples.csv")
  
  tryCatch({
    readr::write_delim(expr_df, file = file_path)
    message("Table saved as: ", file_path)
  }, error = function(e) {
    stop("Error saving table as CSV: ", e$message)
  })
  
  return(expr_df)
}


subsample_batch<-function(df,batch_list,batch_colm) {
  
  num<-20000
  set.seed(3000)
  
  subsamp<-c()
  
  for(fi in 1:length(batch_list)){
    tmp<-which(df[,batch_colm]==batch_list[fi])
    subsamp<-c(subsamp, tmp[sample(1:length(tmp), num)])
  }
  rdf<-df[subsamp,]
  return(rdf)
}

##### get the automated marker thresholds for all panel markers
automated_threshold_DF<-function(df,markers,control_list,batch_list,
                                 batch_colm,subsample="yes",dataset_type="flow",output_dir) {
  
  valid_options <- c("yes", "no")
  data_options <- c("flow", "mass")
  
  subsample<-tolower(subsample)
  dataset_type<-tolower(dataset_type)
  
  if(!subsample %in% valid_options) {
    stop("Invalid input: please specify 'yes' or 'no'")
  }
  
  if(!dataset_type %in% data_options) {
    stop("Invalid input: please specify if dataset is 'flow' or 'mass'")
  }
  
  auto_cutoff<-data.frame()
  
  if(!all(control_list %in% unique(df$control))) {
    stop("check if dataframe has all unique control samples!")
  }
  
  colrSet<-crayon::make_style("#539DDD")
  
  total_iter<-length(markers)*length(control_list)
  current_iter<-0
  
  ### create progress bar
  pBar<-progress::progress_bar$new(
    
    format = paste0(colrSet("Calculating automated marker threshold"),
                    " [",colrSet(":bar"),"] ",
                    colrSet(":percent")," | Marker in ",length(control_list),
                    " control samples: :current/:total ",
                    "| Elapsed: :elapsed | ETA: :eta"),
    total = total_iter,
    clear = FALSE,
    width = 80,
    complete = "=",
    incomplete = "-"
  )
  
  if (dataset_type == "flow") {
    for (i in 1:length(control_list)) {
      
      tmp<-df[(df$control==control_list[i]),]
      if(subsample == "yes") {
        tmp<-data.frame(subsample_batch(df = tmp,batch_list = batch_list,batch_colm = batch_colm))
      }
      else {
        tmp<-tmp
      }
      
      for (m in 1:length(markers)) {
        
        current_iter<-current_iter+1
        ### update progress bar
        pBar$tick(tokens = list(current = current_iter))
        
        myfile<-data.frame(marker = markers[m], cutoff = auto_threshold_func(tmp[,m]))
        myfile$control<-as.character(control_list[i])
        auto_cutoff<-rbind(auto_cutoff,myfile)
      }
      cat("marker thresholds for control",i,"of",length(control_list),"calculated","\n")
    }
  }
  else {
    for (i in 1:length(control_list)) {
      
      tmp<-df[(df$control==control_list[i]),]
      if(subsample == "yes") {
        tmp<-data.frame(subsample_batch(df = tmp,batch_list = batch_list,batch_colm = batch_colm))
      }
      else {
        tmp<-tmp
      }
      
      for (m in 1:length(markers)) {
        
        current_iter<-current_iter+1
        ### update progress bar
        pBar$tick(tokens = list(current = current_iter))
        
        myfile<-data.frame(marker = markers[m], cutoff = auto_threshold_func_w_cytof(tmp[,m]))
        myfile$control<-as.character(control_list[i])
        auto_cutoff<-rbind(auto_cutoff,myfile)
      }
      cat("marker thresholds for control",i,"of",length(control_list),"calculated","\n")
    }
  }
  auto_cutoff<-data.frame(auto_cutoff %>% dplyr::group_by(marker) %>% 
                            dplyr::summarise(cutoff = median(cutoff)))
  auto_cutoff[,"cutoff"]<-round(auto_cutoff[,"cutoff"],3)
  auto_cutoff<-auto_cutoff %>% arrange(match(marker,markers))
  
  file_path <- file.path(output_dir, "automated_marker_threshold_values.csv")
  
  tryCatch({
    readr::write_delim(auto_cutoff, file = file_path)
    message("Table saved as: ", file_path)
  }, error = function(e) {
    stop("Error saving table as CSV: ", e$message)
  })
  
  return(auto_cutoff)
}

###############################################################################################
######################~~~~ cytoFlagR initial visualization functions ~~~~######################
###############################################################################################

generate_number_of_cells_plot<-function(df,batch_colm,control_colm,batch_list,
                                        control_list,output_dir,axis_size=20,
                                        width=22.5,height=14.6){
  
  file_name<-"number_of_cells_per_batch_across_control_samples.png"
  
  file_path<-file.path(output_dir, file_name)
  
  ### get sample id colm if one doesn't exist
  if(is.null(df[,"sample_id"])) {
    df[,"sample_id"]<-paste0(df[,batch_colm],"_",df[,control_colm]) 
  }
  
  ### get the number of cells per sample per batch
  ncells<-data.frame(table(df[,"sample_id"]))
  colnames(ncells)<-c("sample_id","Freq")
  
  ### match columns to append batch and control information
  mm<-match(ncells$sample_id,df$sample_id)
  ncells$batch<-df$batch[mm]
  ncells$control<-df$control[mm]
  ncells$Freq<-as.numeric(ncells$Freq)
  
  ncells$control<-factor(ncells$control,levels = control_list)
  ncells$batch<-factor(ncells$batch,levels = batch_list)
  
  num_cols<-as.numeric(length(control_list)/2)
  
  # if(num_cols %% 2 != 0) {
  #   num_cols<-as.numeric((length(control_list)+1)/2)
  # }
  
  nplot<-ggplot(ncells, aes(x = Freq, y = batch)) + theme_bw() + 
    geom_bar(colour="#121518", fill= "#7A7A78", position="dodge", stat="identity", width=0.8)+
    facet_wrap(~control, nrow = num_cols)+
    ylab("Batch") + xlab("Number of Cells")+
    scale_x_continuous(limits = c(0, as.numeric(max(ncells$Freq))), 
                       expand = expansion(mult = c(0, .02)))+
    theme(axis.text.y = element_text(size = axis_size+5,colour = "#4B4B4B"),
          axis.text.x = element_text(size = axis_size),
          axis.title = element_text(size = axis_size+10,colour = "#000000"),
          strip.text = element_text(size = axis_size+12,face='bold'))
  
  ### save the plot
  ggplot2::ggsave(file_path,plot = nplot,width = width,height = height,dpi = 400)
  
  # Notify user of completion
  cat("Number of cells plot saved as", file_path, "\n")
}

generate_MDS_plot<-function(df,batch_colm,control_colm,marker_list,
                            sample_colours,output_dir,axis_size=15,
                            width=8,height=7){
  
  file_name<-"MDS_per_control_sample.png"
  
  file_path<-file.path(output_dir, file_name)
  
  ### get sample id colm if one doesn't exist
  if(is.null(df[,"sample_id"])) {
    df[,"sample_id"]<-paste0(df[,batch_colm],"_",df[,control_colm]) 
  }
  
  batch_list<-as.character(unique(df[,batch_colm])) ## unique batch names
  control_list<-as.character(unique(df[,control_colm])) ## unique control sample names
  
  ### calculate median marker expression matrix
  expr_median<-data.frame(df[,marker_list],sample_id = df[,"sample_id"])
  expr_median<-expr_median %>%  dplyr::group_by(sample_id) %>% 
    dplyr::summarize_all(list(median))
  ### transpose median marker expression matrix 
  exp_med_tbl<-t(expr_median[,-1])
  colnames(exp_med_tbl) <- expr_median$sample_id
  
  ### calculate MDS dimensions
  mds<-limma::plotMDS(exp_med_tbl, plot = FALSE)
  
  mdsdf<-data.frame(MDS1 = mds$x, MDS2 = mds$y,
                    sample_id = colnames(exp_med_tbl))
  ### match columns to append batch and control information
  mm<-match(mdsdf$sample_id,df$sample_id)
  mdsdf$batch<-df$batch[mm]
  mdsdf$control<-df$control[mm]
  mdsdf$stimulation<-df$stimulation[mm]
  
  mdsdf$control<-factor(mdsdf$control,levels = control_list)
  mdsdf$batch<-factor(mdsdf$batch,levels = batch_list)
  
  if(min(mdsdf$MDS1) < min(mdsdf$MDS2)) {
    min_lim<-as.numeric(round(min(mdsdf$MDS1),1))
  }
  else {
    min_lim<-as.numeric(round(min(mdsdf$MDS2),1))
  }
  
  if(max(mdsdf$MDS1) > max(mdsdf$MDS2)) {
    max_lim<-as.numeric(round(max(mdsdf$MDS1),1))
  }
  else {
    max_lim<-as.numeric(round(max(mdsdf$MDS2),1))
  }
  
  brk<-seq(min_lim,max_lim,by=0.1)
  
  mds_plot<-ggplot(mdsdf, aes(x = MDS1, y = MDS2, color = control)) +
    geom_point(size = 2) +
    geom_label_repel(aes(label = batch,fontface = 'bold',size=5.2), 
                     max.overlaps = Inf, , show.legend = F) + # 
    scale_color_manual("Control Sample",values = sample_colours) +
    coord_cartesian(xlim = c(min_lim,max_lim),ylim = c(min_lim,max_lim)) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size=4.3))) +
    theme_bw() +
    theme(axis.text = element_text(size = axis_size),
          axis.title = element_text(size = axis_size+2),
          legend.title = element_text(size = axis_size+2),
          legend.text = element_text(size = axis_size+2))
  
  ### save the plot
  ggplot2::ggsave(file_path,plot = mds_plot,width = width,height = height,dpi = 500)
  
  # Notify user of completion
  cat("MDS saved as", file_path, "\n")
}

#### subsampling function ####
subsetting<-function(df,samps,num=20000) {
  
  subsamp<-c()
  
  for(fi in 1:length(samps)){
    tmp<-which(df[,"sample_id"]==samps[fi])
    subsamp<-c(subsamp, tmp[sample(1:length(tmp), num)])
  }
  
  rdf<-df[subsamp,]
  return(rdf)
}

generate_UMAP<-function(df,batch_colm,control_colm,marker_list,batch_list,
                        control_list,sample_ids,batch_colours,output_dir,
                        axis_size=11,width=5,height=5){
  
  cat("This process will take a long time, please be patient\n")
  file_name<-"UMAP_of_control_samples.png"
  
  file_path<-file.path(output_dir, file_name)
  
  ### get sample id colm if one doesn't exist
  if(is.null(df[,"sample_id"])) {
    df[,"sample_id"]<-paste0(df[,batch_colm],"_",df[,control_colm]) 
  }
  
  nums<-as.numeric(round(min(table(df$sample_id))/length(batch_list)))
  
  if(nums > 2500) {
    nums<-2500
  }
  
  ### subsample dataframe to create umap
  subdf<-data.frame(subsetting(df,sample_ids,nums))
  
  #### get character vector of all cells denoted by sample_id
  samp_filt<-subdf[,"sample_id"]
  cat("calculating UMAP dimensions...\n")
  ### get UMAP dimensions
  out_umap<-umap::umap(subdf[,marker_list], config = umap.defaults)
  
  dims_umap<-out_umap$layout
  colnames(dims_umap) <- c("UMAP_1", "UMAP_2")
  
  stopifnot(nrow(dims_umap)==length(samp_filt))
  
  dims_umap <- cbind(as.data.frame(dims_umap), sample_id = samp_filt, type = "UMAP")
  
  ### match columns to append batch and control information
  mm<-match(dims_umap$sample_id,subdf$sample_id)
  dims_umap$batch<-subdf$batch[mm]
  dims_umap$control<-subdf$control[mm]
  dims_umap$batch<-factor(dims_umap$batch,levels = batch_list)
  dims_umap$control<-factor(dims_umap$control,levels = control_list)
  
  names(batch_colours)<-batch_list
  
  umapdf<-subset(dims_umap,select = c(UMAP_1,UMAP_2,control,batch))
  
  cat("Generating UMAP plot...\n")
  umap_plot<-ggplot(umapdf, aes(x = UMAP_1, y = UMAP_2, colour = batch)) + 
    geom_point(alpha = 0.3,size=0.1) + 
    scale_colour_manual("Batch",values = batch_colours) + 
    labs(x = "UMAP_1", y = "UMAP_2") + 
    theme_bw() + 
    guides(color = guide_legend(override.aes = list(size=4))) +
    theme(aspect.ratio = 1, axis.text.x = element_text(size = axis_size),
          axis.text.y = element_text(size = axis_size,vjust = 0.75),
          axis.title = element_text(size = axis_size+3,colour = "#4d4d4d"),
          legend.key.size = unit(axis_size-3, "mm"))
  
  ### save the plot
  ggplot2::ggsave(file_path,plot = umap_plot,width = width,height = height,dpi = 400)
  
  # Notify user of completion
  cat("UMAP saved as", file_path, "\n")
}

###############################################################################################
#############################~~~~ cytoFlagR metrics functions ~~~~#############################
###############################################################################################

#### negative MFI dataframe function ####
neg_mfi_df_all<-function(df,markers,cutf,sample_col_name,batch_col_name,min_cells = TRUE) {
  
  names(df)[names(df)==sample_col_name]<-"control"
  names(df)[names(df)==batch_col_name]<-"batch"
  
  ncdf<-data.frame(table(df$sample_id))
  names(ncdf)[names(ncdf)=="Var1"]<-"sample_id"
  
  ### get sample id column if one doesn't exist
  if(is.null(df[,"sample_id"])) {
    df[,"sample_id"]<-paste0(df[,"batch"],"_",df[,"control"])
  }
  
  ndf<-data.frame(sample_id=NA,value=NA,control=NA,batch=NA,marker=NA)
  
  if(!(is.null(df[,'control']))) {
    for (m in 1:length(markers)) {
      mn<-markers[m]
      tmp<-data.frame(expression = df[,mn], sample_id = df[,"sample_id"],
                      marker = mn) %>% dplyr::filter(expression < cutf[m,2])
      
      ### get number of positive cells
      tmp_freq<-data.frame(table(tmp$sample_id))
      colnames(tmp_freq)<-c("sample_id","neg_Freq")
      
      ### get the batch and control info
      mm<-match(tmp_freq$sample_id,tmp$sample_id)
      tmp_freq$control<-tmp$control[mm]
      tmp_freq$batch<-tmp$batch[mm]
      tmp_freq$marker<-mn

      ### calculate MFI of positive population
      tmp_median<-data.frame(tmp %>% group_by(sample_id) %>% 
                               dplyr::summarise(value = median(expression)))
      
      mm<-match(tmp_median$sample_id,df$sample_id)
      tmp_median$control<-df$control[mm]
      tmp_median$batch<-df$batch[mm]
      tmp_median$marker<-markers[m]
      
      ### find batches with number cells < 100 
      cells_below_threshold<-tmp_freq[tmp_freq$pos_Freq < 100, ]
      
      if(min_cells == TRUE) {
        if(nrow(cells_below_threshold) > 0) {
          for(i in 1:nrow(cells_below_threshold)) {
            idx<-which(tmp_median$sample_id == cells_below_threshold$sample_id[i] & 
                         tmp_median$marker == cells_below_threshold$marker[i] &
                         tmp_median$control == cells_below_threshold$control[i] & 
                         tmp_median$batch == cells_below_threshold$batch[i])
            if(length(idx) > 0) {
              tmp_median$sample_id[idx]<-NA
              tmp_median$value[idx]<-NA
              tmp_median$control[idx]<-NA
              tmp_median$batch[idx]<-NA
              tmp_median$marker[idx]<-NA
            }
          }
        }
      }
      ndf<-rbind(ndf,tmp_median)
    }
    ndf<-ndf[-1,]
    ndf$type<-"-MFI"
    return(ndf)
  }
  else {
    stop("Please provide control information of each batch and protein channel")
  }
}

#### positive MFI dataframe function ####
pos_mfi_df_all<-function(df,markers,cutf,sample_col_name,batch_col_name,min_cells = TRUE) {
  
  names(df)[names(df)==sample_col_name]<-"control"
  names(df)[names(df)==batch_col_name]<-"batch"
  
  ### get sample id column if one doesn't exist
  if(is.null(df[,"sample_id"])) {
    df[,"sample_id"]<-paste0(df[,"batch"],"_",df[,"control"]) 
  }
  
  ncdf<-data.frame(table(df$sample_id))
  names(ncdf)[names(ncdf)=="Var1"]<-"sample_id"
  
  pdf<-data.frame(sample_id=NA,value=NA,control=NA,batch=NA,marker=NA)
  
  if(!(is.null(df[,'control']))) {
    for (m in 1:length(markers)) {
      mn<-markers[m]
      tmp<-data.frame(expression = df[,mn], sample_id = df[,"sample_id"], 
                      marker = mn) %>% dplyr::filter(expression >= cutf[m,2])
      
      ### get number of positive cells
      tmp_freq<-data.frame(table(tmp$sample_id))
      colnames(tmp_freq)<-c("sample_id","pos_Freq")
      
      ### get the batch and control info
      mm<-match(tmp_freq$sample_id,tmp$sample_id)
      tmp_freq$control<-tmp$control[mm]
      tmp_freq$batch<-tmp$batch[mm]
      tmp_freq$marker<-mn
      
      ### calculate MFI of positive population
      tmp_median<-data.frame(tmp %>% group_by(sample_id) %>% 
                        dplyr::summarise(value = median(expression)))
      
      mm<-match(tmp_median$sample_id,df$sample_id)
      tmp_median$control<-df$control[mm]
      tmp_median$batch<-df$batch[mm]
      tmp_median$marker<-markers[m]
      
      ### find batches with number cells < 100 
      cells_below_threshold<-tmp_freq[tmp_freq$pos_Freq < 100, ]
      
      if(min_cells == TRUE) {
        if(nrow(cells_below_threshold) > 0) {
          for(i in 1:nrow(cells_below_threshold)) {
            idx<-which(tmp_median$sample_id == cells_below_threshold$sample_id[i] & 
                         tmp_median$marker == cells_below_threshold$marker[i] &
                         tmp_median$control == cells_below_threshold$control[i] & 
                         tmp_median$batch == cells_below_threshold$batch[i])
            if(length(idx) > 0) {
              tmp_median$sample_id[idx]<-NA
              tmp_median$value[idx]<-NA
              tmp_median$control[idx]<-NA
              tmp_median$batch[idx]<-NA
              tmp_median$marker[idx]<-NA
            }
          }
        }
      }
      pdf<-rbind(pdf,tmp_median)
    }
    pdf<-pdf[-1,]
    pdf$type<-"+MFI"
    return(pdf)
  }
  
  else {
    stop("Please provide control information of each batch and protein channel")
  }
}

#### percent positive cells dataframe function ####
percent_pos_df_all<-function(df,markers,cutf,ctrls,sample_col_name,batch_col_name) {
  
  names(df)[names(df)==sample_col_name]<-"control"
  names(df)[names(df)==batch_col_name]<-"batch"
  
  df<-cbind(df, sample_id = paste0(df[,'control'],"_",df[,'batch']))
  
  ncdf<-data.frame(table(df$sample_id))
  names(ncdf)[names(ncdf)=="Var1"]<-"sample_id"
  
  perc_df<-data.frame(sample_id=NA,pos_Freq=NA,control=NA,batch=NA,tot_cells=NA,
                      value=NA,marker=NA)
  
  if(!(is.null(df[,'control']))) {
    for (m in 1:length(markers)) {
      mn<-markers[m]
      tmp<-data.frame(expression = df[,mn], sample_id = df[,"sample_id"], 
                      control = df[,'control'], batch = df[,'batch'],
                      marker = mn) %>% dplyr::filter(expression >= cutf[m,2])
      
      pMark_freq<-data.frame(table(tmp$sample_id))
      colnames(pMark_freq)<-c("sample_id","pos_Freq")
      
      ### get the batch and control info
      mm<-match(pMark_freq$sample_id,tmp$sample_id)
      pMark_freq$control<-tmp$control[mm]
      pMark_freq$batch<-tmp$batch[mm]
      
      ### get the percent of positive cells
      mm<-match(pMark_freq$sample_id,ncdf$sample_id)
      pMark_freq$tot_cells<-ncdf$Freq[mm]
      pMark_freq<-pMark_freq %>% dplyr::mutate(value = ((pos_Freq/tot_cells)*100))
      pMark_freq$marker<-mn
      
      perc_df<-rbind(perc_df,pMark_freq)
    }
    perc_df<-perc_df[-1,]
    perc_df<-data.frame(subset(perc_df,select=c(sample_id,value,control,batch,marker)))
    perc_df$type<-"%pos"
    return(perc_df)
  }
  
  else {
    print("Please provide control information of each batch and protein channel")
  }
}

#### IQR outlier flagging function ####

detect_outlier <- function(x) {
  iqr<-1.5
  Q1<-quantile(x, 0.25)
  Q3<-quantile(x, 0.75)
  return(x > Q3 + (iqr*IQR(x)) | x < Q1 - (iqr*IQR(x)))
}

iqr_metric_function<-function(df,markers,sample_col_name,batch_col_name){
  
  names(df)[names(df)==sample_col_name]<-"control"
  names(df)[names(df)==batch_col_name]<-"batch"
  
  type<-as.character(unique(df$type))

  iqr15_metric<-list()
  
  for (m in markers) {
    iqr15_metric[[m]]<-list()
    iqr15_metric[[m]]<-df %>% 
      dplyr::filter(marker == m) %>% 
      group_by(control) %>%
      dplyr::mutate(flagged_batch = ifelse(is.na(value),"none",
                                           ifelse(detect_outlier(value),
                                                  paste0(batch),"none")))
  }
  
  #### unlist into dataFrames ####
  MFIdf<-do.call(rbind.data.frame, iqr15_metric)
  rownames(MFIdf)<-NULL
  MFIdf$type<-as.character(type)
  MFIdf<-data.frame(MFIdf)
  
  return(MFIdf)
  
}

#### IQR metrics dataframe for every marker, control sample and batch ####
iqr_dataframe<-function(neg_df,pos_df,percent_pos_df,markers,sample_col_name,batch_col_name,output_dir) {
  
  df<-data.frame(rbind(neg_df,pos_df,percent_pos_df))
  type_list<-c("-MFI","+MFI","%pos")
  
  iqr_df<-data.frame()
  
  for (t in 1:length(type_list)) {
    metric<-type_list[t]
    tmpdf<-df[(df$type==type_list[t]),]
    iqr_tmp<-data.frame(iqr_metric_function(df=tmpdf,markers=markers,
                                            sample_col_name=sample_col_name,
                                            batch_col_name=batch_col_name))
    iqr_tmp<-iqr_tmp %>% 
      dplyr::mutate(batch_effect = ifelse(is.na(value),0,
                                          ifelse(flagged_batch=="none",
                                                 0,1)))
    iqr_tmp<-iqr_tmp %>% 
      dplyr::mutate(cols = ifelse(is.na(value),"#F5F5F5",
                                  ifelse(batch_effect==1,"#8E000C","#AEAEAE")))
    
    iqr_df<-rbind(iqr_df,iqr_tmp)
  }
  filename<-paste0("IQR_metrics_flags_across_markers_controls_batches.csv")
  file_path <- file.path(output_dir, filename)
  
  tryCatch({
    readr::write_delim(iqr_df, file = file_path)
    message("Table saved as: ", file_path)
  }, error = function(e) {
    stop("Error saving table as CSV: ", e$message)
  })
  
  return(iqr_df)
}

iqr_boxplot<-function(df,type,marker,colours,batch_list,control_list,ctrl_labs=NULL) {
  
  df<-df[(df$type==type),]
  
  ### select marker
  df_mark<-df[(df$marker==marker),]
  ### exclude batches that are assigned as having < 100 cells
  df_mark<-df_mark[!(df_mark$cols=="#F5F5F5"),]
  ### assign non-outliers as NA
  df_mark[,"flagged_batch"]<-sapply(df_mark[,"flagged_batch"], 
                                    function (x) gsub("none",NA,x))
  df_mark$control<-factor(df_mark$control,levels = control_list)
  df_mark$batch<-factor(df_mark$batch,levels = batch_list)
  
  names(colours)<-batch_list
  names(ctrl_labs)<-control_list
  
  if(!is.null(ctrl_labs)) {
    #### plot the IQR boxplot ###
    
    ggplot(df_mark, aes(control, value)) +
      ggtitle(paste0("IQR flags for ",type," values in ",marker,"\n")) +
      geom_boxplot(outlier.color = NA,colour = "#000000",fill="#F2F2F2",width=0.6) +
      geom_point(aes(fill = batch), size = 3.8, pch = 21, colour = "#4D4D4D", 
                 position = position_jitter(width = 0.3), key_glyph = "rect") +
      scale_fill_manual("Batch",values = colours) +
      geom_label_repel(data = df_mark[complete.cases(df_mark$flagged_batch),],
                       aes(color = flagged_batch, fontface = 'bold', label = flagged_batch),
                       show.legend = F, max.overlaps = Inf, size=6.8) +
      scale_color_manual(values = colours) +
      scale_x_discrete(label = ctrl_labs) +
      xlab("Control Sample")+
      ylab(paste0(type))+
      theme_bw() +
      theme(axis.text.x = element_text(size = 26,angle = 90, vjust = 0.5), # , hjust = 0.9
            axis.text.y = element_text(size = 25),
            axis.title = element_text(size = 28),
            plot.title = element_text(size = 22,face = "bold"),
            legend.title = element_text(size = 19,face = "bold"), 
            legend.text = element_text(size = 17))
  }
  else {
    #### plot the IQR boxplot ###
    
    ggplot(df_mark, aes(control, value)) +
      ggtitle(paste0("IQR flags for ",type," values in ",marker,"\n")) +
      geom_boxplot(outlier.color = NA,colour = "#000000",fill="#F2F2F2",width=0.6) +
      geom_point(aes(fill = batch), size = 3.8, pch = 21, colour = "#4D4D4D", 
                 position = position_jitter(width = 0.3), key_glyph = "rect") +
      scale_fill_manual("Batch",values = colours) +
      geom_label_repel(data = df_mark[complete.cases(df_mark$flagged_batch),],
                       aes(color = flagged_batch, fontface = 'bold', label = flagged_batch),
                       show.legend = F, max.overlaps = Inf, size=6.8) +
      scale_color_manual(values = colours) +
      scale_x_discrete(label = control_list) +
      xlab("Control Sample")+
      ylab(paste0(type))+
      theme_bw() +
      theme(axis.text.x = element_text(size = 26,angle = 90, vjust = 0.5), # , hjust = 0.9
            axis.text.y = element_text(size = 25),
            axis.title = element_text(size = 28),
            plot.title = element_text(size = 22,face = "bold"),
            legend.title = element_text(size = 19,face = "bold"), 
            legend.text = element_text(size = 17))
  }
}


iqr_boxplot_fixed_coords<-function(df,type,marker,control_list,batch_list,colours,ctrl_labs=NULL) {
  
  df<-df[(df$type==type),]
  
  lim1<-round(min(df$value),3)
  lim2<-round(max(df$value),3)

  ### select marker
  df_mark<-df[(df$marker==marker),]
  ### exclude batches that are assigned as having < 100 cells
  df_mark<-df_mark[!(df_mark$cols=="#F5F5F5"),]
  ### assign non-outliers as NA
  df_mark[,"flagged_batch"]<-sapply(df_mark[,"flagged_batch"], 
                                    function (x) gsub("none",NA,x))
  
  df_mark$control<-factor(df_mark$control,levels = control_list)
  df_mark$batch<-factor(df_mark$batch,levels = batch_list)
  
  names(colours)<-batch_list
  names(ctrl_labs)<-control_list

  if(!is.null(ctrl_labs)) {
    #### plot the IQR boxplot with fixed axes ###
    
    ggplot(df_mark, aes(control, value)) +
      ggtitle(paste0("IQR flags for ",type," values in ",marker,"\n")) +
      geom_boxplot(outlier.color = NA,colour = "#000000",fill="#F2F2F2",width=0.6) +
      geom_point(aes(fill = batch), size = 3.8, pch = 21, colour = "#4D4D4D", 
                 position = position_jitter(width = 0.3), key_glyph = "rect") +
      scale_fill_manual("Batch",values = colours) +
      geom_label_repel(data = df_mark[complete.cases(df_mark$flagged_batch),],
                       aes(color = flagged_batch, fontface = 'bold', label = flagged_batch),
                       show.legend = F, max.overlaps = Inf, size=6.8) +
      scale_color_manual(values = colours) +
      xlab("Control Sample")+
      ylab(paste0(type))+
      scale_x_discrete(label = ctrl_labs) +
      scale_y_continuous(limits = c(lim1,lim2)) +
      theme_bw() +
      theme(axis.text.x = element_text(size = 26,angle = 90, vjust = 0.5),
            axis.text.y = element_text(size = 25),
            axis.title = element_text(size = 28),
            plot.title = element_text(size = 22,face = "bold"),
            legend.title = element_text(size = 19,face = "bold"), 
            legend.text = element_text(size = 17))
  }
  else {
    #### plot the IQR boxplot with fixed axes ###
    
    ggplot(df_mark, aes(control, value)) +
      ggtitle(paste0("IQR flags for ",type," values in ",marker,"\n")) +
      geom_boxplot(outlier.color = NA,colour = "#000000",fill="#F2F2F2",width=0.6) +
      geom_point(aes(fill = batch), size = 3.8, pch = 21, colour = "#4D4D4D", 
                 position = position_jitter(width = 0.3), key_glyph = "rect") +
      scale_fill_manual("Batch",values = colours) +
      geom_label_repel(data = df_mark[complete.cases(df_mark$flagged_batch),],
                       aes(color = flagged_batch, fontface = 'bold', label = flagged_batch),
                       show.legend = F, max.overlaps = Inf, size=6.8) +
      scale_color_manual(values = colours) +
      xlab("Control Sample")+
      ylab(paste0(type))+
      scale_x_discrete(label = control_list) +
      scale_y_continuous(limits = c(lim1,lim2)) +
      theme_bw() +
      theme(axis.text.x = element_text(size = 26,angle = 90, vjust = 0.5),
            axis.text.y = element_text(size = 25),
            axis.title = element_text(size = 28),
            plot.title = element_text(size = 22,face = "bold"),
            legend.title = element_text(size = 19,face = "bold"), 
            legend.text = element_text(size = 17))
  }
}



generate_IQR_metric_boxplots<-function(iqr.plt,type,markers,colours,
                                       batch_list,control_list,output_dir,
                                       width = 11, height = 10,ctrl_labs) {
  # Define the PDF file path
  file_name<-paste0("boxplots_of_all_markers_for_",type,"_metric.pdf")
  if (type == "%pos") {
    file_name<-sapply(file_name, function(f)gsub("%pos","percent_pos",f))
  }
  pdf_file_path<-file.path(output_dir, file_name)
  
  # Open the PDF device
  pdf(pdf_file_path, width = width, height = height, useDingbats = FALSE)
  
  # Loop through the markers and create IQR metric boxplots with flags
  for (i in 1:length(markers)) {
    mk<-markers[i]
    # generate boxplot with free axes
    plot<-iqr_boxplot(df = iqr.plt,type = type,marker = mk,batch_list = batch_list,
                      colours = colours, control_list = control_list,ctrl_labs=ctrl_labs)
    print(plot)
  }
  # Close the PDF device
  dev.off()
  
  # Notify user of completion
  cat("PDF saved as", pdf_file_path, "\n")
}

generate_fixed_IQR_metric_boxplots<-function(iqr.plt,type,markers,
                                             control_list,colours,batch_list,output_dir,
                                             width = 11, height = 10,ctrl_labs) {
  
  # Define the PDF file path
  file_name<-paste0("fixed_axes_boxplots_of_all_markers_for_",type,"_metric.pdf")
  if (type == "%pos") {
    file_name<-sapply(file_name, function(f)gsub("%pos","percent_pos",f))
  }
  
  pdf_file_path<-file.path(output_dir, file_name)
  
  # Open the PDF device
  pdf(pdf_file_path, width = width, height = height, useDingbats = FALSE)
  
  # Loop through the markers and create IQR metric boxplots with flags
  for (i in 1:length(markers)) {
    mk<-markers[i]
    # generate boxplot with fixes axes
    plot<-iqr_boxplot_fixed_coords(df = iqr.plt,type = type,marker = mk,
                                   colours = colours,control_list = control_list,
                                   batch_list = batch_list,ctrl_labs=ctrl_labs)
    print(plot)
  }
  # Close the PDF device
  dev.off()
  
  # Notify user of completion
  cat("PDF saved as", pdf_file_path, "\n")
}

#### create molten transformed dataframe ####
reshape_df<-function(df,sample_col_name,batch_col_name,markers) {
  
  df[,"sample_id"]<-paste0(df[,batch_col_name],"_",df[,sample_col_name])
  
  molten_df<-melt(data.frame(df[,markers],sample_id = df[,"sample_id"]),
                  id.vars = "sample_id",variable.name = "marker",
                  value.name = "expression")
  mm<-match(molten_df$sample_id,df$sample_id)
  molten_df[,"batch"]<-df[,batch_col_name][mm]
  molten_df[,"control"]<-df[,sample_col_name][mm]
  
  return(molten_df)
}

#### subsample entire dataframe ####
subsample_df<-function(cdf,samps,markers,num = 20000) {
  set.seed(3000)
  subdf<-data.frame()
  
  for (i in 1:length(markers)) {
    tmp<-cdf[(cdf$marker==markers[i]),]
    tmpdf<-subsetting(tmp,samps,num)
    subdf<-rbind(subdf,tmpdf)
  }
  return(subdf)
}

EMD_calc<-function(df,binSize,control_list,batch_list){
  
  ### select a upper and lower limits for binning (required by EMD)
  lim1<-(min(df[,"expression"]) - binSize)
  min_lim<-round((min(df[,"expression"]) - binSize))
  if(lim1 < min_lim) {
    min_lim<-round((lim1-1))
  }
  
  lim2<-(max(df[,"expression"]) + binSize)
  max_lim<-round((max(df[,"expression"]) + binSize))
  if(lim2>max_lim) {
    max_lim<-round((lim2+1))
  }
  
  # Create list of binned distribution matrices
  distr<-list()
  
  for (tcs in control_list) {
    distr[[tcs]]<-list()
    for (bat in batch_list) {
      distr[[tcs]][[bat]]<-df %>% 
        dplyr::filter(control == tcs, batch == bat) %>%
        dplyr::select(expression) %>%
        apply(2, function(x){
          ### bin the data
          bins<-seq(min_lim, max_lim, by = binSize)
          if(length(x)==0) {
            rep(0, times = length(bins) - 1)
          }
          else {
            graphics::hist(x, breaks = bins, plot = FALSE)$density
          }
        })
    }
  }
  # compute EMD
  distances<-list()
  
  for (tcs in control_list) {
    distances[[tcs]]<-list()
    distances[[tcs]]<-matrix(NA, nrow = length(batch_list), 
                             ncol = length(batch_list), dimnames = list(batch_list,batch_list))
    for (i in seq_along(batch_list)[-length(batch_list)]) {
      batch1<-batch_list[i]
      A<-matrix(distr[[tcs]][[batch1]])
      for (j in seq(i + 1, length(batch_list))) {
        batch2 <- batch_list[j]
        B<-matrix(distr[[tcs]][[batch2]])
        distances[[tcs]][batch1, batch2] <- emdist::emd2d(A,B)
      }
    }
  }
  return(distances)
}

#### EMD matrix list for every marker, control sample and batch ####
EMD_list<-function(df,samps,markers,batch_list,control_list,num = 20000) {
  
  subdf<-data.frame(subsample_df(cdf = df, samps = samps, markers = markers, num = 20000))
  
  set.seed(3000)
  binSize<-0.1
  
  colrSet<-crayon::make_style("#539DDD")
  total_iter<-length(markers)
  ### create progress bar
  pBar<-progress::progress_bar$new(
    
    format = paste0(colrSet("Calculating pairwise EMDs for each marker and control"),
                    " [",colrSet(":bar"),"] ",
                    colrSet(":percent")," | Marker: :current/:total ",
                    "| Elapsed: :elapsed | ETA: :eta"),
    total = total_iter,
    clear = FALSE,
    width = 80,
    complete = "=",
    incomplete = "-"
  )
  
  emds_list<-list()
  
  for (m in markers) {
    
    emds_list[[m]]<-list()
    tmp<-subdf[(subdf$marker==m),]
    emds_list[[m]]<-EMD_calc(df = tmp,binSize = binSize,
                             control_list = control_list,batch_list = batch_list)
    pBar$tick()
  }
  return(emds_list)
}

median_vals<-function(batch_list,distf,control_list,marker) {
  
  df<-melt(distf,na.rm = TRUE)
  ### medians
  
  med_Var1<-data.frame(df %>% arrange(match(Var1,batch_list)))
  med_Var1<-med_Var1[,-2]
  colnames(med_Var1)<-c("batch","value")
  med_Var2<-data.frame(df %>% arrange(match(Var2,batch_list)))
  med_Var2<-med_Var2[,-1]
  colnames(med_Var2)<-c("batch","value")
  
  med<-data.frame(rbind(med_Var1,med_Var2))
  med<-med %>% group_by(batch) %>% dplyr::summarise(value = median(value))
  colnames(med)<-c("batch","value")
  med$control<-control_list
  med$marker<-marker
  return(med)
}

#### molten median function ####
dist_median_emd<-function(batch_list,distf,control_list) {
  
  df<-melt(distf,na.rm = TRUE)
  df$reference<-"no"
  
  ### medians
  
  med_Var1<-data.frame(df %>% dplyr::arrange(match(Var1,batch_list)))
  med_Var1<-med_Var1[,-2]
  colnames(med_Var1)<-c("batch","value")
  med_Var2<-data.frame(df %>% dplyr::arrange(match(Var2,batch_list)))
  med_Var2<-med_Var2[,-1]
  colnames(med_Var2)<-c("batch","value")
  
  med<-data.frame(rbind(med_Var1,med_Var2))
  med<-med %>% group_by(batch) %>% dplyr::summarise(value = median(value))
  colnames(med)<-c("Var2","value")
  med<-cbind(Var1 = batch_list,med)
  
  med$Var2<-"median"
  med$reference<-"yes"
  
  df.plt<-rbind(df,med)
  df.plt<-cbind(df.plt,control = control_list)
  df.plt$Var2<-factor(df.plt$Var2,
                      levels = c(batch_list,"median"))
  return(df.plt)
}

#### get median EMD metric flags dataframe ####
emd_med_vals<-function(em_dists,markers,samps,control_list,batch_list,threshold,
                       output_dir,num = 20000) {
  
  threshold<-as.numeric(threshold)
  binSize<-0.1
  
  emdf_all<-data.frame(batch=NA,value=NA,control=NA,marker=NA)
  
  
  for (m in 1:length(markers)) {
    mrk<-markers[m]
    for(j in 1:length(control_list)) {
      ct<-control_list[j]
      dist_list<-em_dists[[mrk]][[ct]]
      tmp_med<-data.frame(median_vals(batch_list = batch_list, distf = dist_list, 
                                      control_list = ct,marker = markers[m]))
      emdf_all<-rbind(emdf_all,tmp_med)
    }
  }
  emdf_all<-emdf_all[-1,]
  
  emdf_all<-emdf_all %>% 
    dplyr::mutate(flagged_batch = ifelse(value>=threshold,paste0(batch),"none"))
  emdf_all<-emdf_all %>% 
    dplyr::mutate(batch_effect = ifelse(value>=threshold,as.numeric(1),as.numeric(0)))
  emdf_all<-emdf_all %>% 
    dplyr::mutate(cols = ifelse(value>=threshold,"#8E000C","#AEAEAE"))
  
  emdf_all<-emdf_all %>% arrange(match(control,control_list))
  emdf_all<-emdf_all %>% arrange(match(batch,batch_list))
  emdf_all<-emdf_all %>% arrange(match(marker,markers))
  
  emdf_all$type<-"EMD"
  
  filename<-paste0("median_EMDs_per_batch_across_markers_and_controls_flags.csv")
  file_path <- file.path(output_dir, filename)
  
  tryCatch({
    readr::write_delim(emdf_all, file = file_path)
    message("Table saved as: ", file_path)
  }, error = function(e) {
    stop("Error saving table as CSV: ", e$message)
  })
  return(emdf_all)
}

EMD_hmap<-function(emds_list,marker,control,batch_list,axis_size=18,threshold=5) {
  
  dist_list<-emds_list[[marker]][[control]]
  
  df.plt<-dist_median_emd(batch_list,dist_list,control)
  rownames(df.plt)<-NULL
  
  col_seq<-seq(round(min(df.plt$value)),round(max(df.plt$value)+1))
  
  if(threshold==5) {
    col_fun_1 = colorRamp2(c(0,1,2.5,threshold), c("grey","grey75","bisque","firebrick3"))
    colour_heat<-col_fun_1(col_seq)
    mid<-2.5
    lim<-as.numeric(threshold+1)
    breaks<-c(0,0.5,seq(1,6))
    labls<-c("0.0","","1.0","2.0","3.0","4.0","5.0","> 5")
  }
  else {
    col_fun_1 = colorRamp2(c(0,1,2.5,threshold), c("grey","grey75","bisque","firebrick3"))
    colour_heat<-col_fun_1(col_seq)
    mid<-as.numeric((threshold)/2)
    lim<-as.numeric(threshold+1)
    breaks<-c(0,0.5,seq(1,lim))
    labls<-c("0.0","",as.character(seq(1,threshold)),paste0("> ",threshold))
  }

  if(max(df.plt$value) >= threshold) {
    hm01<-ggplot(data = df.plt, aes(x = Var1, y = Var2, color = reference, fill = value)) +
      geom_tile(colour="black") +
      scale_fill_gradientn(colours = colour_heat,guide = "colourbar",
                           breaks=breaks,labels=labls) +
      labs(x = "", y = "", fill = "Earth Mover's \n Distance") +
      geom_text(aes(x = Var1, y = Var2, label = round(value, 3)), color = "black", 
                fontface = "bold", size = 5.3) +
      scale_x_discrete(position = "top") +
      ggtitle(paste0("Pairwise EMDs for ",marker," in control sample: ",control,"\n")) +
      theme(axis.text.x = element_text(face="bold", size = axis_size, 
                                       angle = 70, vjust = 1, hjust = 0),
            axis.text.y = element_text(face="bold", size = axis_size),
            plot.title = element_text(size = 22,face = "bold"),
            legend.title = element_text(face="bold", size = axis_size-5),
            legend.text = element_text(size = axis_size-7),legend.key.height = unit(1.2,'cm'))
  }
  else {
    if(lim < 2) {
      breaks<-seq(0,lim,by=0.0001)
    }
    else {
      breaks<-seq(0,lim)
    }
    hm01<-ggplot(data = df.plt, aes(x = Var1, y = Var2, color = reference, fill = value)) +
      geom_tile(colour="black") +
      scale_fill_gradient2(high = "firebrick3",mid = "bisque",
                           low = "grey",midpoint = mid,
                           guide = "colourbar",breaks=breaks) +
      labs(x = "", y = "", fill = "Earth Mover's \n Distance") +
      geom_text(aes(x = Var1, y = Var2, label = round(value, 4)), color = "black", 
                fontface = "bold", size = 5.3) +
      scale_x_discrete(position = "top") +
      ggtitle(paste0("Pairwise EMDs for ",marker," in control sample: ",control,"\n")) +
      theme(axis.text.x = element_text(face="bold", size = axis_size, 
                                       angle = 70, vjust = 1, hjust = 0),
            axis.text.y = element_text(face="bold", size = axis_size),
            plot.title = element_text(size = 22,face = "bold"),
            legend.title = element_text(face="bold", size = axis_size-5),
            legend.text = element_text(size = axis_size-7),legend.key.height = unit(1.2,'cm'))
  }
  return(hm01)
}


generate_EMD_heatmaps<-function(emds_list,batch_list,markers,controls,
                                threshold,output_dir,axis_size=18,width=13,height=12) {
  # Define the PDF file path
  file_name<-"EMD_heatmap_for_all_markers.pdf"
  pdf_file_path<-file.path(output_dir, file_name)
  
  colrSet<-crayon::make_style("#539DDD")
  total_iter<-length(markers)
  ### create progress bar
  pBar<-progress::progress_bar$new(
    
    format = paste0(colrSet("Plotting marker expression densities"),
                    " [",colrSet(":bar"),"] ",
                    colrSet(":percent")," | Marker: :current/:total ",
                    "| Elapsed: :elapsed | ETA: :eta"),
    total = total_iter,
    clear = FALSE,
    width = 80,
    complete = "=",
    incomplete = "-"
  )
  
  # Open the PDF device
  pdf(pdf_file_path, width = width, height = height, useDingbats = FALSE)
  
  # Loop through the markers and create IQR metric boxplots with flags
  for (i in 1:length(markers)) {
    mk<-as.character(markers[i])
    for(c in 1:length(controls)) {
      ct<-controls[c]
      
      cat("\nProcessing marker:", mk, "control:", ct, "\n")
      
      hmap<-EMD_hmap(emds_list=emds_list,marker = mk,control = ct,
                     batch_list = batch_list,axis_size = axis_size,
                     threshold = threshold)
      print(hmap)
    }
    ### update progress bar
    pBar$tick(tokens = list(current = i))
  }
  # Close the PDF device
  dev.off()
  
  # Notify user of completion
  cat("PDF saved as", pdf_file_path, "\n")
}


EMD_boxplot<-function(distdf,batch_colours,control,marker,
                      batch_list,threshold=5,axis_size=18) {
  
  cat("Processing marker:", marker, "and control:", control, "\n")
  
  EMD_matrix<-distdf[[marker]][[control]]
  lower_triangle<-t(EMD_matrix)
  
  emd_square_matrix<-matrix(NA, nrow = nrow(EMD_matrix), ncol = ncol(EMD_matrix))
  emd_square_matrix[upper.tri(emd_square_matrix)] <- EMD_matrix[upper.tri(EMD_matrix)]
  emd_square_matrix[lower.tri(emd_square_matrix)] <- lower_triangle[lower.tri(lower_triangle)]
  rownames(emd_square_matrix)<-rownames(EMD_matrix)
  colnames(emd_square_matrix)<-colnames(EMD_matrix)
  emd_square_matrix<-melt(emd_square_matrix)
  
  df1<-emd_square_matrix %>% drop_na()
  colnames(df1)<-c("pair","ref","value")
  
  df1_med<-data.frame(value = df1$value, batch = df1$ref) %>% 
    group_by(batch) %>% dplyr::summarise_all(list(median))
  df1_med<-df1_med[order(df1_med$value, decreasing = TRUE),]
  df1_med$batch<-as.character(df1_med$batch)
  
  #### reorder dataframe to match the order of median controls
  df1<-df1 %>% arrange(match(ref,df1_med$batch))
  df1$ref<-factor(df1$ref, levels = unique(df1_med$batch))
  df1$pair<-as.character(df1$pair)
  df1$pair<-factor(df1$pair, levels = unique(df1$pair))
  
  ### designate colour to batches with median EMDs > 5 
  df1_med<-data.frame(df1_med) %>% dplyr::mutate(cols = ifelse(value > threshold,
                                                               "#FF2316","#999999"))
  mm<-match(df1$ref,df1_med$batch)
  df1$cols<-df1_med$cols[mm]
  
  #### get a list of flagged batches
  out_batch<-unique(as.character(df1[which(df1$cols=="#FF2316"),"ref"]))
  out_batch<-out_batch[order(match(out_batch,batch_list))]
  
  ### point shape vector
  shape_vect<-c(22,23,24,25,18,8,9,11,4,13,3,14,12,10,1,5,2,7,0,6,16,15,19,17,
                22,23,24,25,18,8,9,11,4,13,3,14,12,10,1,5,2,7,0,6,16,15,19,17,
                22,23,24,25,18,8,9,11,4,13,3,14,12,10,1,5,2,7,0,6,16,15,19,17,
                22,23,24,25,18,8,9,11,4,13,3,14,12,10,1,5,2,7,0,6,16,15,19,17)
  
  df1<-cbind(df1,shapes = 21)
  
  if(length(out_batch)!=0) {
    for(b in 1:length(out_batch)) {
      for (i in 1:nrow(df1)) {
        if(df1[i,"pair"]==out_batch[b]) {
          df1[i,"shapes"]<-shape_vect[b]
        }
      }
    }
  }
  
  names(batch_colours)<-batch_list
  
  if(max(df1$value) >= 20) {
    max_lim<-round(max(df1$value))
    if(max(df1$value) > max_lim) {
      max_lim<-max_lim+1
    }
    
    ggplot(df1, aes(ref, value)) +
      geom_boxplot(aes(color=cols),alpha = 1.0, fill="#F2F2F2", outlier.color = NA) +
      scale_color_identity()+
      geom_point(aes(fill=pair,shape = shapes), colour = "#000000",size = 1.8, 
                 alpha = 1.0, position = position_jitter(width = 0.4),
                 key_glyph = "rect") + # key_glyph = "rect"
      scale_fill_manual("Batch",values = batch_colours) +
      scale_shape_identity()+
      geom_hline(yintercept = threshold,colour="#FF2316",linewidth=0.9) +
      xlab("")+ylab("Earth Mover's Distance") + 
      theme_bw() +
      scale_y_continuous(limits = c(0,max_lim)) +
      ggtitle(paste0("Ordered pairwise EMDs per batch for ",marker,
                     "\n"," in control sample: ",control,"\n")) +
      theme(axis.text.x = element_text(size = axis_size,angle = 70, vjust = 0.9, hjust = 0.9),
            axis.text.y = element_text(size = axis_size,vjust = 0.75),
            axis.title.y = element_text(size = axis_size+4),
            plot.title = element_text(size = 22),
            legend.title = element_text(size = axis_size/2,face = "bold"), 
            legend.text = element_text(size = (axis_size/2)-2),
            legend.key.size = unit(3,'mm'))
  }
  else {
    ggplot(df1, aes(ref, value)) +
      geom_boxplot(aes(color=cols),alpha = 1.0, fill="#F2F2F2", outlier.color = NA) +
      scale_color_identity()+
      geom_point(aes(fill=pair,shape = shapes), colour = "#000000",size = 1.8, 
                 alpha = 1.0, position = position_jitter(width = 0.4),
                 key_glyph = "rect") + # key_glyph = "rect"
      scale_fill_manual("Batch",values = batch_colours) +
      scale_shape_identity()+
      geom_hline(yintercept = threshold,colour="#FF2316",linewidth=0.9) +
      xlab("")+ylab("Earth Mover's Distance") + 
      theme_bw() +
      scale_y_continuous(limits = c(0,20)) +
      ggtitle(paste0("Ordered pairwise EMDs per batch for ",marker,
                     "\n"," in control sample: ",control,"\n")) +
      theme(axis.text.x = element_text(size = axis_size,angle = 70, vjust = 0.9, hjust = 0.9),
            axis.text.y = element_text(size = axis_size,vjust = 0.75),
            axis.title.y = element_text(size = axis_size+4),
            plot.title = element_text(size = 22),
            legend.title = element_text(size = axis_size/2,face = "bold"), 
            legend.text = element_text(size = (axis_size/2)-2),
            legend.key.size = unit(3,'mm'))
  }
}

generate_EMD_boxplots<-function(distdf,batch_list,markers,control_list,batch_colours,
                                threshold,output_dir,axis_size=18,width=11.5,height=10.5) {
  # Define the PDF file path
  file_name<-"EMD_ordered_median_boxplot_for_all_markers.pdf"
  pdf_file_path<-file.path(output_dir, file_name)
  
  colrSet<-crayon::make_style("#539DDD")
  total_iter<-length(markers)
  ### create progress bar
  pBar<-progress::progress_bar$new(
    
    format = paste0(colrSet("Plotting marker expression densities"),
                    " [",colrSet(":bar"),"] ",
                    colrSet(":percent")," | Marker: :current/:total ",
                    "| Elapsed: :elapsed | ETA: :eta"),
    total = total_iter,
    clear = FALSE,
    width = 80,
    complete = "=",
    incomplete = "-"
  )
  
  names(batch_colours)<-batch_list
  
  # Open the PDF device
  pdf(pdf_file_path, width = height, height = height, useDingbats = FALSE)

  # Loop through the markers and create IQR metric boxplots with flags
  for (i in 1:length(markers)) {
    mk<-markers[i]
    for(c in 1:length(control_list)) {
      ct<-control_list[c]
      box_plot<-EMD_boxplot(distdf=distdf, batch_colours = batch_colours, 
                            marker = mk, control = ct, batch_list = batch_list, 
                            threshold = threshold, axis_size = axis_size)
      print(box_plot)
    }
    ### update progress bar
    pBar$tick(tokens = list(current = i))
  }
  # Close the PDF device
  dev.off()
    
  # Notify user of completion
  cat("PDF saved as", pdf_file_path, "\n")
}

##### set axis size function #####
auto_axis_size<-function(axis,num_markers,num_controls) {
  
  if(!axis %in% c('x','y')) {
    stop("axis must be specified - either 'x' or 'y'")
  }
  if(!is.numeric(num_markers) || !is.numeric(num_controls)) {
    stop("num_markers (number of markers) and num_controls (number of controls) must be numeric")
  }
  
  base_size<-18 ## default axis size
  
  if(axis =='x') {
    set<-4*(num_markers)
    if(set < 80) {
      axis_size<-as.numeric(base_size)*2
    }
    else if (set < 100 && set >= 80) {
      axis_size<-as.numeric(base_size)
    }
    else {
      size_factor<-round(min(1,100/set),1)
      axis_size<-round(base_size*size_factor)
    }
  }
  else {
    min_size<-as.numeric(num_controls)
    if (min_size <= 4 && min_size > 2) {
      axis_size<-as.numeric(base_size)+1
    }
    else if (min_size <= 2) {
      axis_size<-as.numeric(base_size+(base_size/4))
    }
    else {
      size_factor<-as.numeric(subtract((num_controls),4))
      if(size_factor < base_size) { 
        axis_size<-as.numeric(subtract(base_size,size_factor))
      }
      else if (size_factor == base_size) {
        axis_size<-as.numeric(base_size/8)
      }
      else if (size_factor > base_size && size_factor <= base_size+2) {
        axis_size<-as.numeric(base_size/8)
      }
      else {
        axis_size<-as.numeric(1)
      }
    }
  }
  axis_size<-as.numeric(axis_size)
  return(axis_size)
}

ranked_flagged_hmap<-function(emd_df,iqr_df,batch_list,controls,markers,batch_cols,output_dir=NULL) {
  
  conditions<-c("-MFI","+MFI","%pos","EMD")
  
  ### select columns to use
  iqr_df<-subset(iqr_df,select = c(batch,value,control,marker,flagged_batch,
                                   batch_effect,cols,type))
  emd_df<-subset(emd_df,select = c(batch,value,control,marker,flagged_batch,
                                   batch_effect,cols,type))
  
  #### column of unique conditions per marker
  iqr_df$var2<-paste0(iqr_df$type,"_",iqr_df$marker)
  emd_df$var2<-paste0(emd_df$type,"_",emd_df$marker)
  
  #### combine them for ranked hmap
  hmdf<-data.frame(rbind(iqr_df,emd_df))
  
  #### re-arrange by list of conditions, batches, and markers
  #### this will order the dataframe by each condition for each marker 
  
  hmdf<-hmdf %>% arrange(match(type,conditions))
  hmdf<-hmdf %>% arrange(match(batch,batch_list))
  hmdf<-hmdf %>% arrange(match(marker,markers))
  
  #### order in decreasing order of flags, while maintaining marker & batch order
  hmdf<-hmdf[order(hmdf$batch_effect, decreasing = TRUE),]
  hmdf<-hmdf %>% arrange(match(batch,batch_list))
  hmdf<-hmdf %>% arrange(match(marker,markers))
  
  #### force order of controls for the plot
  hmdf$control<-factor(hmdf$control, levels = rev(controls))
  
  ## get the frequency of flags for all marker-condition combinations
  freq_var2<-data.frame(table(hmdf$var2,hmdf$batch_effect))
  colnames(freq_var2)<-c("var2","flag","Freq")
  #### order by decreasing order of flags per combo
  freq_var2<-freq_var2[order(freq_var2$flag, freq_var2$Freq, decreasing = TRUE),]
  freq_var2$var2<-as.character(freq_var2$var2)
  
  #### reorder hmap dataframe to match order of marker-condition flags
  hmdf<-hmdf %>% arrange(match(var2,freq_var2$var2))
  
  ## get the frequency of flags for all batches
  batch_freq<-data.frame(table(hmdf$batch,hmdf$batch_effect))
  colnames(batch_freq)<-c("batch","flag","Freq")
  #### order by frequency of flags per batch
  batch_freq<-batch_freq[order(batch_freq$flag, batch_freq$Freq, decreasing = TRUE),]
  batch_freq$batch<-as.character(batch_freq$batch)
  
  #### reorder hmap dataframe to match order of flags of batches
  hmdf<-hmdf %>% arrange(match(batch,batch_freq$batch))
  
  ## get the frequency of flags for all panel markers
  marker_freq<-data.frame(table(hmdf$marker,hmdf$batch_effect))
  colnames(marker_freq)<-c("marker","flag","Freq")
  #### order by frequency of flags per marker
  marker_freq<-marker_freq[order(marker_freq$flag, marker_freq$Freq, decreasing = TRUE),]
  marker_freq$marker<-as.character(marker_freq$marker)
  
  #### re-order hmap dataframe by conditions, flagged batch frequency, and 
  #### flagged marker frequency
  hmdf<-hmdf %>% arrange(match(batch,batch_freq$batch))
  hmdf<-hmdf %>% arrange(match(type,conditions))
  hmdf<-hmdf %>% arrange(match(marker,marker_freq$marker))
  hmdf<-hmdf %>% arrange(match(batch,batch_freq$batch))
  
  #### re-order hmap dataframe by decreasing order of flags
  hmdf<-hmdf[order(hmdf$batch_effect, decreasing = TRUE),]
  hmdf<-hmdf %>% arrange(match(var2,freq_var2$var2))
  hmdf<-hmdf %>% arrange(match(batch,batch_freq$batch))
  hmdf<-hmdf %>% arrange(match(marker,marker_freq$marker))
  hmdf<-hmdf[order(hmdf$batch_effect, decreasing = TRUE),]
  
  ### get a vector of markers with most to least flags
  mord<-unique(marker_freq$marker)
  
  #### factorize order of markers in the hmap dataframe for the plot
  hmdf$marker<-factor(hmdf$marker, levels = mord)
  
  #### factorize order of batches in the hmap dataframe for the plot
  hmdf$batch<-factor(hmdf$batch, levels = unique(batch_freq$batch))
  
  #### force hmap plot order of conditions flagged for each marker
  var_match<-data.frame(var2=unique(hmdf$var2))
  var_match$type<-sapply(var_match$var2,function(c)gsub("\\_\\D+\\w+","",c))
  var_match$marker<-sapply(var_match$var2,function(c)gsub("\\W\\w+\\_","",c))
  var_match$marker<-sapply(var_match$marker,function(c)gsub("\\w+\\_","",c))
  
  #### order by list of conditions and order of flagged markers
  var_match<-var_match %>% arrange(match(type,conditions))
  var_match<-var_match %>% arrange(match(marker,mord))
  
  #### re-order hmap dataframe by order of marker-condition combos 
  #### decreasing order of flagged batches
  hmdf<-hmdf %>% arrange(match(var2,unique(var_match$var2)))
  hmdf<-hmdf %>% arrange(match(batch,unique(batch_freq$batch)))
  #### factorize order of markers and conditions 
  #### in the hmap dataframe for the plot
  hmdf$var2<-factor(hmdf$var2,levels = unique(var_match$var2))
  # hmdf$type<-factor(hmdf$type,levels = unique(var_match$type))
  hmdf$type<-factor(hmdf$type,levels = conditions)
  
  #### reorder named batch colors to order of freq. of batches flagged
  #### this is for the strip-l function used for the heatmap
  names(batch_cols)<-batch_list
  ord_bcols<-batch_cols
  ord_bcols<-ord_bcols[order(factor(names(ord_bcols), 
                                    levels = unique(batch_freq$batch)))]
  
  if (!is.null(output_dir)) {
    csv_file<-"summary_of_flags_all_metrics.csv"
    csv_file_path<-file.path(output_dir, csv_file)
    readr::write_delim(hmdf,file = csv_file_path)
    cat("Summary of flags saved as:",csv_file_path,"\n")
  }
  
  ### get the size of x and y axis
  xlab_size<-auto_axis_size(axis = 'x',num_markers = length(markers),
                            num_controls = length(controls))
  ylab_size<-auto_axis_size(axis = 'y',num_markers = length(markers),
                            num_controls = length(controls))
  
  hmap_plot<-ggplot(hmdf, aes(x=type,y=control,fill=cols))+
    # ggtitle(paste0("Ranked summary of IQR 1.5 and median EMD flags ",
    #                "across batches, markers, and controls"))+
    geom_tile(colour="black", lwd = 0.6)+
    scale_fill_identity()+
    scale_y_discrete(position = "right", expand = c(0,0))+
    scale_x_discrete(position = "top", expand = c(0,0))+
    facet_grid(batch~marker, scales = "free",space = "free_x", switch = "y")+
    xlab("")+
    ylab("")+
    theme(plot.title = element_text(size = 25, face = "bold"),
          strip.text.y = element_text(size = 21, face = "bold",colour = "black"),
          strip.text.x = element_text(size = xlab_size+4, face = "bold"),
          strip.background.x = element_blank(),
          strip.placement = "outside",
          axis.text.x = element_text(size = xlab_size+2.5, angle = 90,
                                     hjust = -0.3), # , hjust = -0.3,vjust = 3
          axis.text.y = element_text(size = ylab_size+2,hjust = 0.7,vjust = 0.4),
          panel.spacing.x = unit(0.4,'lines'),
          panel.spacing.y = unit(0.5,'lines'),
          panel.grid = element_blank())
  
  ### assigning strip colours based on batch number
  grid.plt <- ggplot_gtable(ggplot_build(hmap_plot))
  stripr <- which(grepl('strip-l', grid.plt$layout$name))
  fills <- ord_bcols
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', grid.plt$grobs[[i]]$grobs[[1]]$childrenOrder))
    grid.plt$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  
  #### plot ranked heatmap of flags with coloured facet strip 
  grid::grid.draw(grid.plt)
}
##### png image dimensions functions #####
image_width<-function(num_markers) {
  num_conditions<-4 ### +/-MFI, %pos, EMD assessments
  size<-as.numeric((num_markers*num_conditions*100)+7500)
}

image_height<-function(num_batch,num_controls) {
  size<-as.numeric((num_controls*num_batch*100)+3000)
}


generate_ranked_flagged_hmap<-function(emd_df,iqr_df,batch_list,controls,markers,batch_cols,output_dir){
  
  nm<-as.numeric(length(markers))
  nb<-as.numeric(length(batch_list))
  nct<-as.numeric(length(controls))
  
  # width<-(image_width(nm))/400
  # height<-(image_height(nb, nct))/400
  
  width<-image_width(nm)
  height<-image_height(nb, nct)
  
  iqr_df<-subset(iqr_df,select = c(batch,value,control,marker,flagged_batch,
                                 batch_effect,cols,type))
  
  filename<-"Ranked_summary_heatmap_of_flags_across_metrics.png"
  png_file_path<-file.path(output_dir, filename)
  
  ### open the PNG device
  png(filename = png_file_path, width = width, height = height, res = 400)
  tryCatch({
    ranked_flagged_hmap(emd_df = emd_df, iqr_df = iqr_df, batch_list = batch_list,
  controls = controls, markers = markers, batch_cols = batch_cols,output_dir = output_dir)
  }, error = function(e) {
    message("Error in generating the summary heatmap: ", e$message)
  }, finally = {
    dev.off()
  })
  
  # Notify user of completion
  cat("Plot saved as", png_file_path, "\n")
}

getFlowSOM_clusters<-function(data_f, markers, seed, meta_cluster_num, samps, 
                              num, output_dir) {
  
  ## this would create a directory to store all the consensus plots for metaclusters
  plot_outdir<-file.path(paste0(output_dir,"/clustering/new_consensus_plots_k",
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
  
  old_dir<-getwd()
  # setwd(plot_outdir)
  # if (!setwd(plot_outdir)) {
  #   stop("Cannot set working directory to: ", plot_outdir)
  # }
  
  tryCatch({
    
    ### subsample dataframe for clustering
    set.seed(1234)
    df<-data.frame(subsetting(df = data_f, samps = samps, num = num))
    
    set.seed(seed)
    ### convert transformed dataframe into a flowFrame object
    ctrl_ff<-flowCore::flowFrame(data.matrix(df[,markers]))
    
    #### create FlowSOM input ####
    fsomIN<-FlowSOM::ReadInput(ctrl_ff, transform = FALSE, scale = FALSE)
    
    #### Build SOM clusters ####
    som<-FlowSOM::BuildSOM(fsomIN)
    
    ### SOM clusters information
    cell_clustering_som<-som$map$mapping[,1]
    codes<-som$map$codes
    
    ## remove the fluorochrome/metal tag names from the prettyColnames
    som[["prettyColnames"]]<-sapply(som[["prettyColnames"]],
                                    function(x)gsub("\\s\\W\\w+\\W","",x))
    
    ### call the ConsensusClusterPlus function to get meta-clusters
    mc<-ConsensusClusterPlus(t(codes), maxK = meta_cluster_num, reps = 100,
                             pItem = 0.9, pFeature = 1, title = file.path(plot_outdir,"consensus"), 
                             plot = "png", clusterAlg = "hc", innerLinkage = "average", 
                             finalLinkage = "average",
                             distance = "euclidean", seed = seed)
    
    ## Get cluster ids and labels from the metacluster object
    code_clustering_1<-mc[[meta_cluster_num]]$consensusClass
    cell_clustering_1<-code_clustering_1[cell_clustering_som]
    
    df<-cbind(df, cluster = as.character(cell_clustering_1))
    
    filename<-paste0("clustered_dataframe_k",meta_cluster_num,"_clusters_seed",seed,".csv")
    file_path <- file.path(output_dir, filename)
    
    tryCatch({
      readr::write_delim(df, file = file_path)
      message("cluster data saved as: ", file_path)
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

clusterSizes_plot<-function(cluster_df, seed, meta_cluster_num, output_dir, width = 5, height = 4) {
  
  minclust<-0.005 ## 0.5%
  csizes <- data.frame(table(cluster_df[,"cluster"]))
  ### percent of cells in each cluster
  csizes<-cbind(csizes, percent = csizes[,"Freq"]/nrow(cluster_df))
  csizes$xlab<-paste0(csizes[,"Var1"]," (",csizes[,"Freq"],")")
  csizes[,"xlab"] <- factor(csizes[,"xlab"],
                            levels=as.character(csizes[rev(order(csizes[,"percent"])),"xlab"]))
  csizes[,"Var1"] <- factor(csizes[,"Var1"],
                            levels=as.character(csizes[rev(order(csizes[,"percent"])),"Var1"]))
  
  p<-ggplot(data=csizes, aes(x=xlab, y=percent)) +
    geom_bar(stat="identity",aes(color=xlab),show.legend = F) +
    scale_y_continuous(labels = scales::percent)+
    ggtitle(paste0("cluster sizes for: ",meta_cluster_num,
                   " metaclusters of controls","\n (seed ",seed,")")) +
    geom_hline(yintercept=c(minclust), colour = "#F25A25")+
    xlab("") + ylab("% of all analyzed cells")+
    theme(axis.text.x = element_text(angle = 45, vjust=1.0, hjust=1.0))
          
  filename<-paste0("clusterSizes_ordered_k",meta_cluster_num,"_seed",seed,".pdf")
  pdf_file_path<-file.path(output_dir, filename)
  
  tryCatch({
    ggsave(filename = pdf_file_path, plot = p, width = width, height = height, device = "pdf")
    # Notify user of completion
    cat("PDF saved as", pdf_file_path, "\n")
  }, error = function(e) {
    stop("Error saving plot:",e$message)
  })
}

clustering_heatmap<-function(cluster_df, markers, cluster_colrs, cell_cluster_vector, 
                                  meta_cluster_num, seed, axis_size = 15, num_size = 8.4, 
                                  output_dir, width = 11.5, height = 8) {
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
  ### color breaks for the legend
  legend_breaks = seq(from = 0, to = 1.2, by = 0.2)
  ### labels for the row with cluster proportions
  labels_row <- paste0(expr01_median$cell_clustering, " (", clustering_prop, "%)")
  
  # Annotation for the original clusters
  annotation_row <- data.frame(Cluster = factor(expr01_median$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  clustColors1 <- cluster_colrs[1:nlevels(annotation_row$Cluster)]
  names(clustColors1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = clustColors1)
  
  plt<-ComplexHeatmap::pheatmap(expr_heat, name = "Scaled Median\nExpression",
                           color = color_heat, cluster_cols = FALSE,
                           cluster_rows = FALSE, labels_row = labels_row,
                           display_numbers = TRUE, number_color = "black",
                           fontsize = axis_size, fontsize_number = num_size, 
                           column_names_side = c("top"),
                           legend_breaks = legend_breaks, #annotation_row = annotation_row,
                           legend_labels = c("0","0.2","0.4","0.6","0.8","1.0",""),
                           annotation_legend = FALSE)
  
  filename<-paste0("clustering_heatmap_of_k",meta_cluster_num,"clusters_seed",seed,".pdf")
  pdf_file_path<-file.path(output_dir, filename)
  
  pdf(pdf_file_path, width = width, height = height, useDingbats=F)
  print(plt)
  dev.off()
  
  # Notify user of completion
  cat("Cluster Heatmap saved as", pdf_file_path, "\n")
}

proportion_of_batches_per_cluster<-function(cluster_df, cell_cluster_vector, 
                                            batch_list, batch_colours, seed, 
                                            output_dir, width=3000, height=1900, 
                                            axis_size=12) {
  
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
  
  ### get number of clusters
  meta_cluster_num<-length(unique(cluster_df$cluster))
  
  
  filename<-paste0("batch_proportions_per_cluster_barplot_",
                   meta_cluster_num,"_clusters_seed",seed,".png")
  png_file_path<-file.path(output_dir, filename)
  
  csizes <- data.frame(table(cluster_df[,"cluster"]))
  ### percent of cells in each cluster
  csizes<-cbind(csizes, percent = csizes[,"Freq"]/nrow(cluster_df))
  csizes$xlab<-paste0(csizes[,"Var1"]," (",csizes[,"Freq"],")")
  csizes[,"xlab"] <- factor(csizes[,"xlab"],
                            levels=as.character(csizes[rev(order(csizes[,"percent"])),"xlab"]))
  csizes[,"Var1"] <- factor(csizes[,"Var1"],
                            levels=as.character(csizes[rev(order(csizes[,"percent"])),"Var1"]))
  
  batch_vector<-cluster_df[,"batch"]
  
  counts_table<-table(cell_cluster_vector,batch_vector)
  prop_df<-data.frame(counts_table)
  colnames(prop_df)<-c("cluster","batch","Freq")
  prop_df<-prop_df %>% arrange(match(cluster,levels(csizes$Var1)))
  mm<-match(prop_df$cluster,csizes$Var1)
  prop_df$clust_freq<-csizes$Freq[mm]
  prop_df$clust_props<-csizes$percent[mm]
  prop_df<-prop_df %>% mutate(proportion = (Freq/clust_freq)*100)
  
  prop_df$cluster<-factor(prop_df$cluster, levels = unique(prop_df$cluster))
  
  names(batch_colours)<-batch_list
  
  prop_df$xlab_new<-(prop_df$clust_props)*100
  prop_df$xlab_new<-round(prop_df$xlab_new,2)
  prop_df$xlab_new<-paste0(prop_df$cluster," (",prop_df$xlab_new,"%)")
  prop_df$xlab_new<-factor(prop_df$xlab_new,levels = rev(unique(prop_df$xlab_new)))
  
  barplot<-ggplot(prop_df, aes(x = proportion, y = xlab_new, fill = batch)) +
    geom_bar(stat = "identity",width = 0.92, key_glyph='rect') +
    scale_fill_manual("Batch",values = batch_colours) +
    xlab("Batch Proportions") +
    ylab("Cluster") +
    scale_x_continuous(position = 'bottom',expand = expansion(mult = c(0, .02)))+
    scale_y_discrete(position = 'right')+
    theme_bw() +
    theme(axis.text.x = element_text(size = axis_size,hjust = 0.2),
          axis.text.y = element_text(size = axis_size),
          axis.title = element_text(size = axis_size+8,colour = "#000000"),
          legend.title = element_text(size = axis_size+2,face = "bold"), 
          legend.text = element_text(size = axis_size-2))
  
  prop_threshold_1<-2*round((100/length(batch_list)))
  batch_props_filt_1<-prop_df[which(prop_df$proportion > prop_threshold_1),]
  
  if(nrow(batch_props_filt_1) == 0) {
    warning("No data points above threshold. Cannot create table!")
    png(png_file_path, width = width, height = height, res = 350)
    tryCatch({
      print(barplot)
      dev.off()
      cat("Plot saved as:",png_file_path,"\n")
    }, error = function(r){
      stop("Error in generating the plot: ",r$message)
    })
  }
  
  batch_props_filt_1<-batch_props_filt_1 %>%
    arrange(match(cluster,levels(prop_df$cluster)))
  batch_props_filt_1$batch<-as.character(batch_props_filt_1$batch)
  ### open the PNG device
  png(filename = png_file_path, width = width, height = height, res = 350)
  tryCatch({
    print(barplot)
    dev.off()
    # Notify user of completion
    cat("Plot saved as:", png_file_path, "\n")
  }, error = function(e) {
    dev.off()
    stop("Error in generating the plot: ", e$message, 
            "\nAttempted path: ", png_file_path)
  })
}

proportion_of_batches_per_cluster_table<-function(cluster_df, cell_cluster_vector, 
                                            batch_list, seed, output_dir) {
  
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
  
  ### get number of clusters
  meta_cluster_num<-length(unique(cluster_df$cluster))
  
  filename<-paste0("batches_in_higher_than_expected_proportions_across_",
                   meta_cluster_num,"_clusters_seed",seed,".csv")
  file_path <- file.path(output_dir, filename)
  
  csizes <- data.frame(table(cluster_df[,"cluster"]))
  ### percent of cells in each cluster
  csizes<-cbind(csizes, percent = csizes[,"Freq"]/nrow(cluster_df))
  csizes$xlab<-paste0(csizes[,"Var1"]," (",csizes[,"Freq"],")")
  csizes[,"xlab"] <- factor(csizes[,"xlab"],
                            levels=as.character(csizes[rev(order(csizes[,"percent"])),"xlab"]))
  csizes[,"Var1"] <- factor(csizes[,"Var1"],
                            levels=as.character(csizes[rev(order(csizes[,"percent"])),"Var1"]))
  
  batches<-cluster_df[,"batch"]
  
  counts_table<-table(cell_cluster_vector,batches)
  prop_df<-data.frame(counts_table)
  colnames(prop_df)<-c("cluster","batch","Freq")
  prop_df<-prop_df %>% arrange(match(cluster,levels(csizes$Var1)))
  mm<-match(prop_df$cluster,csizes$Var1)
  prop_df$clust_freq<-csizes$Freq[mm]
  prop_df$clust_props<-csizes$percent[mm]
  prop_df<-prop_df %>% mutate(proportion = (Freq/clust_freq)*100)
  prop_df$cluster<-factor(prop_df$cluster, levels = unique(prop_df$cluster))
  
  prop_df$xlab_new<-(prop_df$clust_props)*100
  prop_df$xlab_new<-round(prop_df$xlab_new,2)
  prop_df$xlab_new<-paste0(prop_df$cluster," (",prop_df$xlab_new,"%)")
  prop_df$xlab_new<-factor(prop_df$xlab_new,levels = rev(unique(prop_df$xlab_new)))

  prop_threshold_1<-2*round((100/length(batch_list)))
  batch_props_filt_1<-prop_df[which(prop_df$proportion > prop_threshold_1),]
  
  if(nrow(batch_props_filt_1) == 0) {
    print("batches with 2x higher than expected proportion not found!")
    return(invisible(NULL))
  }
  else {
    batch_props_filt_1<-batch_props_filt_1 %>% 
      arrange(match(cluster,levels(prop_df$cluster)))
    batch_props_filt_1$batch<-as.character(batch_props_filt_1$batch)
    batch_props_filt_1$bat_clust<-paste0(batch_props_filt_1$batch,
                                         "_",batch_props_filt_1$cluster)
    
    tryCatch({
      readr::write_delim(batch_props_filt_1, file = file_path)
      message("Table saved as: ", file_path)
    }, error = function(e) {
      stop("Error saving table as CSV: ", e$message)
    })
    return(batch_props_filt_1)
  }
}

proportion_of_batches_across_control_samples<-function(cluster_df, cell_cluster_vector,
                                                       batch_props_df, batch_list, control_list, 
                                                       batch_colours, output_dir, width=3800, height=5800, 
                                                       axis_size=12, seed) {
  
  if (!is.data.frame(cluster_df)) {
    stop("cluster_df must be a data frame")
  }
  if (!all(c("cluster", "sample_id") %in% colnames(cluster_df))) {
    stop("cluster_df must contain 'cluster' and 'sample_id' columns")
  }
  if (!dir.exists(output_dir)) {
    stop("Output directory does not exist")
  }

  ### get number of clusters
  meta_cluster_num<-length(unique(cluster_df$cluster))
  
  filename<-paste0("proportions_of_cells_per_control_across_",
                   meta_cluster_num,"_clusters_seed",seed,".png")
  png_file_path<-file.path(output_dir, filename)
  
  csizes <- data.frame(table(cluster_df[,"cluster"]))
  ### percent of cells in each cluster
  csizes<-cbind(csizes, percent = csizes[,"Freq"]/nrow(cluster_df))
  csizes$xlab<-paste0(csizes[,"Var1"]," (",csizes[,"Freq"],")")
  csizes[,"xlab"] <- factor(csizes[,"xlab"],
                            levels=as.character(csizes[rev(order(csizes[,"percent"])),"xlab"]))
  csizes[,"Var1"] <- factor(csizes[,"Var1"],
                            levels=as.character(csizes[rev(order(csizes[,"percent"])),"Var1"]))
  
  sample_ids<-cluster_df[,"sample_id"]
  
  counts_table<-table(cell_cluster_vector,sample_ids)
  prop_df<-data.frame(counts_table)
  colnames(prop_df)<-c("cluster","sample_id","Freq")
  
  prop_df<-prop_df %>% arrange(match(cluster,unique(csizes$Var1)))
  mm<-match(prop_df$cluster,csizes$Var1)
  prop_df$clust_freq<-csizes$Freq[mm]
  prop_df$clust_prop<-csizes$percent[mm]
  prop_df<-prop_df %>% mutate(proportion = (Freq/clust_freq)*100)

  mm<-match(prop_df$sample_id,cluster_df$sample_id)
  prop_df$controls<-cluster_df$control[mm]
  prop_df$batch<-cluster_df$batch[mm]
  
  prop_df$cluster<-factor(prop_df$cluster, levels = unique(prop_df$cluster))
  
  prop_df$strip_lab<-paste0("Cluster ",prop_df$cluster)
  prop_df$batch<-factor(prop_df$batch,levels = batch_list)
  prop_df$control<-factor(prop_df$control,levels = control_list)
  
  prop_df$strip_lab<-factor(prop_df$strip_lab,
                              levels = unique(prop_df$strip_lab))
  
  out_clust<-as.character(unique(batch_props_df$cluster))
  
  #### select the clusters which are impacted by batches
  combodf<-data.frame(prop_df[(prop_df$cluster %in% out_clust),])
  
  prop_thr<-1.5*round((100/length(batch_list)))
  ctrl_props_filt_all<-combodf
  ctrl_props_filt_all$bat_clust<-paste0(ctrl_props_filt_all$batch,"_",
                                        ctrl_props_filt_all$cluster)
  ctrl_props_filt_all<-ctrl_props_filt_all[(ctrl_props_filt_all$bat_clust %in% 
                                              unique(batch_props_df$bat_clust)),]
  ctrl_props_filt_all<-ctrl_props_filt_all %>% group_by(controls) %>% 
    dplyr::filter(proportion > prop_thr)
  
  num_rows<-round(length(levels(prop_df$cluster))/2)
  
  tc<-as.character(seq(1,length(control_list)))
  names(tc)<-control_list
  
  ### open the PNG device
  png(filename = png_file_path, width = width, height = height, res = 400)
  tryCatch({
    print(ggplot(prop_df, aes(control, proportion)) +
            geom_boxplot(outlier.color = NA,width = 0.5) +
            # stat_compare_means() +
            geom_point(aes(fill = batch), alpha = 0.8, pch=21, size =2.6, 
                       position = position_jitter(width = 0.3), key_glyph='rect') +
            scale_fill_manual("Batch",values = batch_colours) +
            geom_label_repel(data = ctrl_props_filt_all,
                             aes(color = batch, fontface = 'bold', label = batch),
                             show.legend = F, max.overlaps = Inf, size=6.4) +
            scale_color_manual(values = batch_colours)+
            facet_wrap(~ strip_lab, nrow = num_rows) +
            scale_x_discrete(label = tc) +
            xlab("Controls")+
            ylab("Proportions") +
            theme_bw() +
            theme(axis.text.x = element_text(size = 23),
                  axis.text.y = element_text(size = 25),
                  axis.title = element_text(size = 25),
                  strip.text.x = element_text(size = 21,face = "bold"),
                  legend.title = element_text(size = 15,face = "bold"), 
                  legend.text = element_text(size = 12)))
  }, error = function(e) {
    message("Error in generating the plot: ", e$message)
  }, finally = {
    dev.off()
  })
  # Notify user of completion
  cat("Plot saved as", png_file_path, "\n")
  
}

generate_cluster_UMAP<-function(df,marker_list,batch_list,
                                batch_colours,seed,output_dir,
                                axis_size=15,width=5,height=5){
  
  cat("This process will take a long time, please be patient","\n")
  
  if (!dir.exists(output_dir)){
    stop("Output directory does not exist")
  }
  ### get number of clusters
  meta_cluster_num<-length(unique(df$cluster))
  
  ### create file name and path
  file_name<-paste0("UMAP_of_k",meta_cluster_num,"_clusters_seed",seed,".png")
  
  file_path<-file.path(output_dir, file_name)
  
  # nums<-as.numeric(min(table(df$cluster)))
  
  ### number of cells per cluster
  ncells<-data.frame(table(df$cluster))
  colnames(ncells)<-c("cluster","Freq")
  
  mm<-match(df$cluster,ncells$cluster)
  df$csize<-ncells$Freq[mm]
  
  ### subsample dataframe to create umap
  set.seed(seed)
  
  subsamp <- c()
  
  samps<-unique(df$cluster)
  

  for(fi in 1:length(samps)){
    tmp <- which(df[,"cluster"]==samps[fi])
    nums<-as.numeric(unique(df[which(df[,"cluster"]==samps[fi]),"csize"]))
    subsamp <- c(subsamp, tmp[sample(1:length(tmp), min(nums, 2000))])
  }
  
  subdf<-df[subsamp,]
  subdf$samp_clust<-paste0(subdf$sample_id,"_",subdf$cluster)
  
  samp_clust<-subdf[,"samp_clust"]
  
  cat("calculating UMAP dimensions...\n")
  
  ### get UMAP dimensions
  out_umap<-umap::umap(subdf[,marker_list], config = umap.defaults)
  
  dims_umap<-out_umap$layout
  colnames(dims_umap) <- c("UMAP_1", "UMAP_2")
  
  stopifnot(nrow(dims_umap)==length(samp_clust))
  
  dims_umap <- cbind(as.data.frame(dims_umap), sample_id = samp_clust, type = "UMAP")
  
  ### match columns to append batch and control information
  mm<-match(dims_umap$sample_id,subdf$samp_clust)
  dims_umap$cluster<-subdf$cluster[mm]
  dims_umap$batch<-subdf$batch[mm]
  dims_umap$control<-subdf$control[mm]
  dims_umap$batch<-factor(dims_umap$batch,levels = batch_list)
  
  names(batch_colours)<-batch_list
  
  umapdf<-subset(dims_umap,select = c(UMAP_1,UMAP_2,batch))
  
  cat("Generating UMAP plot...\n")
  
  umap_plot<-ggplot(umapdf, aes(x = UMAP_1, y = UMAP_2, colour = batch)) + 
    geom_point(alpha = 0.7,size=0.1) + 
    scale_colour_manual("Batch",values = batch_colours) + 
    labs(x = "UMAP_1", y = "UMAP_2") + 
    theme_bw() + 
    guides(color = guide_legend(override.aes = list(size=4))) +
    theme(aspect.ratio = 1, axis.text.x = element_text(size = axis_size+7),
          axis.text.y = element_text(size = axis_size+7,vjust = 0.75),
          axis.title = element_text(size = axis_size+10,colour = "#222528"),
          legend.title = element_text(size = (axis_size-3)),
          legend.text = element_text(size = (axis_size-5)))
  
  ### save the plot
  ggplot2::ggsave(file_path,plot = umap_plot,width = width,height = height,dpi = 400)
  
  # Notify user of completion
  cat("UMAP saved as", file_path, "\n")
}

###############################################################################################

