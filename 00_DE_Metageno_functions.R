###########################################################################
# Project           : MRM Health - Statistical support for MH002-UC-201 - PLX-23-35515 (1545)
# Program name      : 99_DE_Metageno_functions.R
# Developed in      : R 4.2.2
# Purpose           : Set of functions used in the analysis of metagenomic data
# Inputs            : -
#   
# Outputs           : -
#
# Revision History  :
#
#   Version   Date        Author       Revision 
#   -------   ---------   ----------   ---------------------------------
#   1.0       29092023    Pablo Villegas  Creation
#
# Reviewed by       : 
# Reference number  : NA                                      
# Validation Level  : NA
#
# Declaration of Confidentiality:
#   - I declare that this program may contain or contain confidential
#     information and that it cannot be shared as it.
###########################################################################

#######################################################
##                   Check libraries                 ##
#######################################################
# Package names
packages <- c("ggplot2", "tidyverse", "data.table", "DESeq2", "gridExtra", 
              "biomaRt", "BiocParallel", "vsn", "GOfuncR", "heatmap.2x", "ggrepel",
              "clusterProfiler","org.Hs.eg.db", "grid", "gridExtra", "plotly", "ggdendro", "DOSE")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if(any(installed_packages == FALSE)){
  install.packages(packages[!installed_packages])
}else{
  message("All packages are installed")
}

#######################################################
##                   Load Libraries                  ##
#######################################################
# note: this assumes you are using R version 3.6
library(tidyverse)
library(data.table)
library(DESeq2)
library(gridExtra)
library(biomaRt)
library(BiocParallel)
library(vsn)
library(GOfuncR)
library(ggplot2)
library(heatmap.2x)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(grid)
library(gridExtra)
library(plotly)
library(ggdendro)
library(DOSE)
library(ashr)
set.seed(100)

#######################################################
##              TFL conversion function              ##
#######################################################

plot.to.tfl <- function(p.ggplot,
                        scale.format.x = 1,
                        scale.mar.up = 1,
                        scale.mar.right = 1,
                        scale.mar.bottom = 1,
                        scale.mar.left = 1,
                        tfl.header = "", 
                        tfl.foot.1, tfl.foot.2, tfl.foot.3, tfl.foot.4,
                        outDir = "", 
                        prefix = ""){
  
    title <- tfl.header
    foot1 <- tfl.foot.1
    foot2 <- tfl.foot.2
    foot3 <- tfl.foot.3
    foot4 <- tfl.foot.4
    
    # output as a pdf
    output <- file.path(outDir, paste0(prefix, ".pdf"))
    pdf(output, height=8.5, width=11)
    
    print(p.ggplot +
            theme(axis.text.y = element_text(size=10), 
                  axis.text.x = element_text(size=11*scale.format.x), 
                  axis.title.x = element_text(size=12), 
                  axis.title.y = element_text(size=11),
                  plot.title = element_text(hjust = 0.5),
                  plot.margin=unit(c(50*scale.mar.up,
                                     50*scale.mar.right,
                                     50*scale.mar.bottom,
                                     50*scale.mar.left),"mm")))
    
    grid.text(title, x = .08, y = .93, just = c('left', 'top'), gp = gpar(fontsize = 12, col = 'black'))
    grid.text(paste0(foot1, '\n', foot2, '\n', foot3, '\n', foot4), x = .08, y = .06, just = c('left', 'bottom'), gp = gpar(fontsize = 10, col = 'black'))
    dev.off()
}

#######################################################
##            Generate Table RTF format              ##
#######################################################

table.rtf <- function(tab, 
                      rtf.header = title, 
                      rtf.foot.1 = foot1, 
                      rtf.foot.2 = foot2, 
                      rtf.foot.3 = foot3, 
                      rtf.foot.4 = foot4,
                      colW = 1.5,
                      outDir = out.dir, 
                      prefix = out.prefix){
  
  output <- file.path(outDir, paste0(prefix, "_table.rtf"))
  
  rtf <<- RTF(output, width=11, height=8.5, font.size=12, omi=rep(0.75,ncol(tab)))
  addText(rtf, "\\deff {\\fonttbl {\\f5 Times;}}\\n")
  
  startParagraph(rtf)
  
  addText(rtf,bold=F,paste("\\f5\\fs22", rtf.header, "{\\pindtabqr}","", "\n", sep=" "))
  
  endParagraph(rtf)
  
  addTable(rtf, tab, row.names = F, font.size = 11, 
           col.widths=colW,
           col.justify = rep("L",ncol(tab)))
  
  addText(rtf, "\\deff {\\fonttbl {\\f5 Times;}}\\n")
  addText(rtf, "\\n {\\pindtabql}")
  
  startParagraph(rtf)
  addText(rtf, paste("\\f5\\fs20", rtf.foot.1, "\n", sep=" "))
  addText(rtf, paste("\\f5\\fs20", rtf.foot.2, "\n", sep=" "))
  addText(rtf, paste("\\f5\\fs20", rtf.foot.3, "\n", sep=" "))
  addText(rtf, paste("\\f5\\fs20", rtf.foot.4, "\n", sep=" "))
  endParagraph(rtf)
  done(rtf)
}

#######################################################
##        Hypothesis testing Metageno pipe           ##
#######################################################

hypothesisTesting <- function(df=clin.data.biomarker_ra,analysis="RNR",bmk="CRG", res_phenotype, target_taxa, scoreName, labelScore){
  
  scoreName <- c(scoreName)
  labelScore <- c(labelScore)
  
  RA_df_analysis_corr <- lapply(seq(length(scoreName)),function(x){
    df %>% 
      dplyr::select(SUBJID,sample_name,Treatment,response=all_of(res_phenotype), score = any_of(scoreName[x])) %>% 
      mutate(measure = labelScore[x])
  }) %>% rbindlist()
  
  if(grepl("RNR",analysis)){
    RA_df_analysis_corr_Ref <- RA_df_analysis_corr %>%
      dplyr::filter(response == "NR") 
    
    RA_df_analysis_corr_Treat <- RA_df_analysis_corr %>%
      dplyr::filter(response == "R")
  }else if(grepl("PT",analysis)){
    RA_df_analysis_corr_Ref <- RA_df_analysis_corr %>%
      dplyr::filter(Treatment == "PBO") 
    
    RA_df_analysis_corr_Treat <- RA_df_analysis_corr %>%
      dplyr::filter(Treatment == "MH002")
  }
  
  if(nrow(RA_df_analysis_corr_Ref) == 0 | nrow(RA_df_analysis_corr_Treat) == 0){
    
    naval <- rep(NA,length(scoreName))
    napval <- naval
    names(naval) <- paste0("stat.",scoreName)
    names(napval) <- paste0("pvalue.",scoreName)
    
    RA_df_analysis.test_res <- cbind(data.frame(Analysis = analysis, Phenotype = res_phenotype, Taxa = target_taxa, Biomarker = bmk),
          t(data.frame(naval)),
          t(data.frame(napval)))
    
  }else{
    
    wilcox_all <- lapply(seq(length(scoreName)),function(x){
      wilcox.test(RA_df_analysis_corr_Ref %>% dplyr::filter(measure == labelScore[x]) %>% pull(score),
                  RA_df_analysis_corr_Treat %>% dplyr::filter(measure == labelScore[x]) %>% pull(score), exact = FALSE)
    })
    
    names(wilcox_all) <- labelScore
    
    RA_df_analysis.test_res <- data.frame()
    for(i in labelScore){
      resW <- wilcox_all[[i]]
      resW <- c(resW$statistic, resW$p.value)
      names(resW) <- paste0(c("stat","pvalue"))
      
      i <- case_when(i == "BMK" ~ bmk,
                     i != "BMK" ~ i)
      
      dftmp <- data.frame(Analysis =analysis,Phenotype = res_phenotype, Taxa = target_taxa, Score = i)
      dftmp <- cbind(dftmp,t(resW))
      
      RA_df_analysis.test_res <- rbind(RA_df_analysis.test_res, dftmp)
    }

  }
  return(list(RA_df_analysis_corr,RA_df_analysis.test_res))
}

#######################################################
##        Correlation testing Metageno pipe          ##
#######################################################

corTesting <- function(df,analysis,score.1,score.2, comparison,res_phenotype,target_taxa,scoreName,labelScore){
  
  scoreName <- c(scoreName)
  labelScore <- c(labelScore)
  
  RA_df_analysis_corr <- lapply(seq(length(scoreName)),function(x){
    df %>% 
      dplyr::select(SUBJID,Treatment,response=all_of(res_phenotype), score = any_of(scoreName[x])) %>% 
      mutate(measure = labelScore[x])
  }) %>% rbindlist()
  
  if(grepl("RNR",analysis)){
    poplabRef <- "NR"
    RA_df_analysis_corr_Ref <- RA_df_analysis_corr %>%
      dplyr::filter(response == poplabRef) 
    
    poplabTar <- "R"
    RA_df_analysis_corr_Treat <- RA_df_analysis_corr %>%
      dplyr::filter(response == poplabTar)
  }else if(grepl("PT",analysis)){
    poplabRef <- "PBO"
    RA_df_analysis_corr_Ref <- RA_df_analysis_corr %>%
      dplyr::filter(Treatment == poplabRef) 
    
    poplabTar <- "MH002"
    RA_df_analysis_corr_Treat <- RA_df_analysis_corr %>%
      dplyr::filter(Treatment == poplabTar)
  }
  
  ## Reference population
  test.res <- try(cor.test(RA_df_analysis_corr_Ref %>% dplyr::filter(measure == comparison) %>% pull(score),
                           RA_df_analysis_corr_Ref %>% dplyr::filter(measure == score.2) %>% pull(score), exact = TRUE),silent = T)
  
  if(is.character(test.res)){
    test.res <- "Error: Corr. not computed"
    dftmp <- data.frame(Analysis =analysis,Phenotype = res_phenotype, Taxa = target_taxa, Score.1 = score.1, Score.2 = score.2, Pop = poplabRef)
    dftmp_ref <- cbind(dftmp,data.frame(stat = test.res,pval = test.res))
  }else{
    dftmp <- data.frame(Analysis =analysis,Phenotype = res_phenotype, Taxa = target_taxa, Score.1 = score.1, Score.2 = score.2, Pop = poplabRef)
    dftmp_ref <- cbind(dftmp,data.frame(stat = test.res$statistic,pval = test.res$p.value)) %>%
      mutate(stat = round(stat,digits = 2),
             pval = round(pval,digits = 2))
  }

  ## Target population
  test.res <- try(cor.test(RA_df_analysis_corr_Treat %>% dplyr::filter(measure == comparison) %>% pull(score),
                           RA_df_analysis_corr_Treat %>% dplyr::filter(measure == score.2) %>% pull(score), exact = TRUE),silent = T)
  
  if(is.character(test.res)){
    test.res <- "sn"
    dftmp <- data.frame(Analysis =analysis,Phenotype = res_phenotype, Taxa = target_taxa, Score.1 = score.1, Score.2 = score.2, Pop = poplabTar)
    dftmp_tar <- cbind(dftmp,data.frame(stat = test.res,pval = test.res))
  }else{
    dftmp <- data.frame(Analysis =analysis,Phenotype = res_phenotype, Taxa = target_taxa, Score.1 = score.1, Score.2 = score.2, Pop = poplabTar)
    dftmp_tar <- cbind(dftmp,data.frame(stat = test.res$statistic,pval = test.res$p.value)) %>%
      mutate(stat = round(stat,digits = 2),
             pval = round(pval,digits = 2))
  }
  
  cor.res <- rbind(dftmp_ref, dftmp_tar) 
  
  if(grepl("PT",analysis)){
    cor.res$Phenotype <- NULL
  }

  return(cor.res)

}

#######################################################
##               Reformat RA dataframe               ##
#######################################################

reformatRAdf <- function(df, taxa, meta, add.col){
  
  RA_df_analysis <- df %>% as.data.frame() %>%
    dplyr::filter(rownames(.) %in% taxa) %>%
    t() %>% as.data.frame() %>%
    dplyr::filter(rownames(.) %in% meta$Pairing) %>%
    rownames_to_column("SUBJID") %>%
    pivot_longer(cols = colnames(.)[-1]) %>% 
    group_by(SUBJID) %>%
    dplyr::summarise(mean_ra = mean(value)) %>%
    left_join(meta %>% 
                dplyr::select(Pairing,Chao1,Shannon,Pielou,Beta_Div_value,GMHI,all_of(add.col)), by = c("SUBJID" = "Pairing"))
  
  trans_scores <- bestNormalize::yeojohnson(as.numeric(RA_df_analysis$mean_ra))
  trans_Chao1 <- bestNormalize::yeojohnson(as.numeric(RA_df_analysis$Chao1))
  trans_Shannon <- bestNormalize::yeojohnson(as.numeric(RA_df_analysis$Shannon))
  trans_Pielou <- bestNormalize::yeojohnson(as.numeric(RA_df_analysis$Pielou))
  trans_Beta <- bestNormalize::yeojohnson(as.numeric(RA_df_analysis$Beta_Div_value))
  trans_gmhi <- bestNormalize::yeojohnson(as.numeric(RA_df_analysis$GMHI))
  
  RA_df_analysis <- RA_df_analysis %>% 
    mutate(mean_ra_trans = trans_scores$x.t) %>%
    mutate(div_Chao1 = trans_Chao1$x.t) %>%
    mutate(div_Shannon = trans_Shannon$x.t) %>%
    mutate(div_Pielou = trans_Pielou$x.t) %>%
    mutate(div_Beta = trans_Beta$x.t) %>%
    mutate(div_gmhi = trans_gmhi$x.t)
  
  return(RA_df_analysis)
}

#######################################################
##         Reformat RA dataframe (genus/family)      ##
#######################################################

reformatRAdf_genusFam <- function(df, taxa, meta, add.col){
  RA_df_analysis <- df %>%
    dplyr::filter(rownames(.) %in% taxa) %>% 
    t() %>% as.data.frame() %>% 
    dplyr::filter(rownames(.) %in% meta$Pairing) %>%
    mutate_all(.funs = function(x){bestNormalize::yeojohnson(x)$x.t}) %>%
    rownames_to_column("SUBJID") %>%
    left_join(meta %>% dplyr::select(Pairing,Chao1,Shannon,Pielou,Beta_Div_value,GMHI,all_of(add.col)), by = c("SUBJID" = "Pairing"))
  
  trans_Chao1 <- bestNormalize::yeojohnson(as.numeric(RA_df_analysis$Chao1))
  trans_Shannon <- bestNormalize::yeojohnson(as.numeric(RA_df_analysis$Shannon))
  trans_Pielou <- bestNormalize::yeojohnson(as.numeric(RA_df_analysis$Pielou))
  trans_Beta <- bestNormalize::yeojohnson(as.numeric(RA_df_analysis$Beta_Div_value))
  trans_gmhi <- bestNormalize::yeojohnson(as.numeric(RA_df_analysis$GMHI))
  
  RA_df_analysis <- RA_df_analysis %>% 
    mutate(div_Chao1 = trans_Chao1$x.t) %>%
    mutate(div_Shannon = trans_Shannon$x.t) %>%
    mutate(div_Pielou = trans_Pielou$x.t) %>%
    mutate(div_Beta = trans_Beta$x.t) %>%
    mutate(div_gmhi = trans_gmhi$x.t)
  
  return(RA_df_analysis)
}

#######################################################
##      Merge RA with biomarker clinical data        ##
#######################################################

formatBmkandRA <- function(clin.df, ra.df, analysis, bmk){
  ## Select biomarker 
  clin.data.biomarker2 <- clin.df %>% 
    dplyr::filter(PARAMCD == bmk) 
  
  trans_scores <- bestNormalize::yeojohnson(as.numeric(clin.data.biomarker2$AVAL))
  
  clin.data.biomarker2 <- clin.data.biomarker2 %>% 
    mutate(AVAL = trans_scores$x.t)
  
  clin.data.biomarker2 <- clin.data.biomarker2 %>%
    left_join(resp.annot, by = c("SUBJID","AVISIT"="timepoint")) %>%
    dplyr::filter(!is.na(Subject.ID)) %>%
    dplyr::rename(timepoint = AVISIT)
  
  ## Merge with relative abundance data
  if(analysis == "PT.Rand"){
    clin.data.biomarker_ra <- clin.data.biomarker2 %>%
      dplyr::filter(timepoint == "Week_8") %>%
      dplyr::filter(SUBJID %in% ra.df$SUBJID) %>%
      left_join(ra.df, by = c("SUBJID"))
  }else if(analysis == "PT.Ext"){
    clin.data.biomarker_ra <- clin.data.biomarker2 %>%
      dplyr::filter(timepoint == "Week_8") %>%
      dplyr::filter(SUBJID %in% ra.df$SUBJID) %>%
      left_join(ra.df, by = c("SUBJID"))
  }else if(analysis == "RNR.Rand"){
    clin.data.biomarker_ra <- clin.data.biomarker2 %>%
      dplyr::filter(timepoint == "Week_8") %>%
      dplyr::filter(Treatment == "MH002") %>%
      dplyr::filter(SUBJID %in% ra.df$SUBJID) %>%
      left_join(ra.df, by = c("SUBJID"))
  }else if(analysis == "RNR.Ext"){
    clin.data.biomarker_ra <- clin.data.biomarker2 %>%
      dplyr::filter(timepoint == "Week_8") %>%
      dplyr::filter(Treatment == "MH002") %>%
      dplyr::filter(SUBJID %in% ra.df$SUBJID) %>%
      left_join(ra.df, by = c("SUBJID"))
  }
  
  return(clin.data.biomarker_ra)
}

#######################################################
##              Summarize data by groups             ##
#######################################################
## This function was taken from 
## http://www.cookbook-r.com/Manipulating_data/Summarizing_data/

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- plyr::ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- dplyr::rename(datac, c(score = "mean"))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  datac[is.na(datac)] <- 0 ## If only one observation, make se and sd 0
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  # ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  # datac$ci <- datac$se * ciMult
  
  return(datac)
}




