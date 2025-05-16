#!/usr/bin/env R
###############################################################################################
#         function to determine the trheshold that separates the negative and positive        #
#         density distributions of markers in a cytometry dataset                             #
#         in order to flag any batch effects detected                                         #
###############################################################################################
### load required packages
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(LaplacesDemon))

###############################################################################################

auto_threshold_func_w_cytof<-function(exp_vector) {
  if(is.unimodal(exp_vector,min.size = 0.05)==TRUE) {
    q4<-quantile(exp_vector,0.95) ## if distribution is unimodal, assume 95%ile as cut-off
    names(q4)<-NULL
    med_x<-quantile(exp_vector,0.5)
    names(med_x)<-NULL
    q3<-quantile(exp_vector,0.75)
    names(q3)<-NULL
    qt_85<-quantile(exp_vector,0.85)
    names(qt_85)<-NULL
    qt_80<-quantile(exp_vector,0.8)
    names(qt_80)<-NULL
    if(q3==0 && q3 < 0.05) {
      val<-q4
      return(round(val,3))
    }
    else if(q4 < 0.5 && q3 < 0.1) {
      val<-quantile(exp_vector,0.98)
      names(val)<-NULL
      return(round(val,3))
    }
    else if(q3 >= 1.0) {
      val<-quantile(exp_vector,0.55)
      names(val)<-NULL
      return(round(val,1))
    }
    else if(med_x <= 0.001) {
      val<-quantile(exp_vector,0.75)
      names(val)<-NULL
      return(round(val,3))
    }
    else if(q3 >= 0.1 && q3 < 1.0) {
      if (q4 > 0.5 && q4 < 1) {
        val<-round((q4+q3)/2,3)
        return(round(val,3))
      }
      else if (q4 >= 1 && q4 < 1.5) {
        val<-round((q4+q3)/2,3)
        return(round(val,3))
      }
      else if (q4 >= 1.5 && qt_85 >= 1) {
        val<-round((qt_85+q3)/2,3)
        return(round(val,3))
      }
      else {
        val<-quantile(exp_vector,0.85)
        names(val)<-NULL
        return(round(val,3))
      }
      # return(round(val,3))
    }
    else {
      val<-quantile(exp_vector,0.75)+0.05
      names(val)<-NULL
      return(round(val,3))
    }
    # return(round(val,3))
  }
  else {
    dens<-stats::density(exp_vector) ### get kernel density of marker using default params
    med_x<-abs(median(dens$x)) ### find the median of the distribution
    q3<-abs(quantile(exp_vector,0.75))
    names(q3)<-NULL
    ix01_y<-which.max(dens$y)
    peak_neg<-dens$x[ix01_y] ### first peak value
    upper_val<-dens$y[dens$x > q3] ### find exp. distribution values greater than 75%ile
    y_idx<-which(dens$y == max(upper_val)) ### find index of max y-value to find 2nd peak
    peak_pos<-dens$x[y_idx] ### second peak value
    
    #### finding the valley
    
    if(peak_pos != peak_neg) {
      min_y<-min(dens$y[dens$x < peak_pos & dens$x > peak_neg]) ### find the min. y-value
      if(min_y==Inf) {
        valley<-round(median(exp_vector),3)
        names(valley)<-NULL
      }
      else {
        y_trough<-which(dens$y == min_y) ### find index of min. y-value for valley
        valley<-dens$x[y_trough] ### pos-neg expr. distribution cut-off for marker
      }
      return(round(valley,3)) 
    }
    else { ### if positive and negative peaks are the same, assume 75%ile as cut-off
      val<-quantile(exp_vector,0.75)
      names(val)<-NULL
      return(round(val,3))
    }
  }
  # else { ### if it is neither unimodal nor bimodal, return 50%ile as cutoff
  #   val<-quantile(exp_vector,0.5)
  #   names(val)<-NULL
  #   return(round(val,3))
  # }
}

###############################################################################################
