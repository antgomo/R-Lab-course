###########################################################################
# Project           : MRM Health - Statistical support for MH002-UC-201 - PLX-23-35515 (1545)
# Program name      : 01_Transcriptomics_Correlation_biomarkers.R
# Developed in      : R 4.2.2
# Purpose           : 
#
# Inputs            : DE analysis results
#                   : ADaM scores of interest 
#   
# Outputs           : RESULTS/01_transcriptomics/INT/eResponse/NR_vs_R_050923_[TRT/noTRT]_Baseline_heatmap2.pdf
#                   : RESULTS/01_transcriptomics/INT/eResponse/NR_vs_R_050923_[TRT/noTRT]_Week_8_heatmap2.pdf
#                   : RESULTS/01_transcriptomics/INT/eResponse/NR_vs_R_050923_[TRT/noTRT]_Biomark_Dist.pdf
#                   : RESULTS/01_transcriptomics/INT/eResponse/NR_vs_R_050923_[TRT/noTRT]_eResponse_corr.summary_table.rtf
#                   : RESULTS/01_transcriptomics/INT/eResponse/NR_vs_R_050923_[TRT/noTRT]_eResponse_summary_table.rtf
#                   : RESULTS/01_transcriptomics/INT/eResponse/NR_vs_R_050923_[TRT/noTRT]_eResponse_wilcox.summary_table.rtf
#                   : RESULTS/01_transcriptomics/INT/eResponse/NR_vs_R_050923_[TRT/noTRT]_eResponse_wilcox_RvsNR.summary_table.rtf
#                   : RESULTS/01_transcriptomics/INT/eResponse/NR_vs_R_050923_[TRT/noTRT]_eResponse_[ANXA2/CLDN1/MMP2/MUC2/NANOG/POU5F1/top20down/top20up]_Corr.pdf

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
##            Clean environment variables            ##
#######################################################
rm(list = ls())

#######################################################
##                   Check libraries                 ##
#######################################################
# Package names
packages <- c("ggplot2", "tidyverse", "data.table", "DESeq2", "gridExtra", 
              "biomaRt", "BiocParallel", "vsn", "GOfuncR", "GSVA", "tximport")

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
library(GSVA)
library(tximport)
library(rtf)
library(RColorBrewer)
# library(xlsx)
library(haven)
set.seed(100)

#######################################################
##                 Paths and variables               ##
#######################################################
qc.status <- "DEV"
analysis.desig <- "INT" # INT, Baseline, Week_8
timepointFiltINT <- "Week_8" # Baseline, Week_8
res_phenotype <- "eResponse" # cResponse, eResponse, cRemission
trt_covar <- FALSE # TRUE, FALSE

base.dir <- "/datadrive/clients/MRM/MH002-UC-201"
prg.dir <- file.path(base.dir, qc.status, "PROGRAMS", "transcriptomics")
prg.file <- tail(unlist(strsplit(split = "/", rstudioapi::getSourceEditorContext()$path)), 1)
in.dir <- file.path(base.dir, "SRC_DATA")
out.dir <- file.path(base.dir, qc.status, "RESULTS", "01_transcriptomics")
resp.annot <- fread(file.path(base.dir, qc.status, "DERIVED_DATA", "Responders_endoscopic_analyses_table.csv")) %>%
  dplyr::rename(Response = IRT) %>%
  dplyr::select(Subject.ID,Response) %>%
  unique()
resp.annot.transl <- fread(file.path(base.dir, qc.status, "DERIVED_DATA", "MRM_responders.csv"))
resp.annot <- resp.annot %>% 
  left_join(resp.annot.transl %>% 
              dplyr::select(sample_name=SAMPLE_ID,Subject.ID=SUBJ_ID,cResponse,cRemission,eResponse), by = c("Subject.ID")) %>%
  dplyr::filter(!is.na(sample_name)) %>%
  mutate(Response = gsub("-","_",Response))
  
in.counts.file <- file.path(in.dir, "Batch.2023-04-24_MH002-UC-201_30066", "rsem_gcounts_matrix.txt")
sample.annot <- file.path(in.dir, "Batch.2023-04-24_MH002-UC-201_30066", "sample_annotation.txt")

## Genes of interest
# Set 1. Intestinal barrier function
set1 <- c("CLDN1","CLDN2","CLDN3",
          "CLDN4", "TJP1", "OCLN",
          "MLCK", "MUC2","MUC3",
          "MUC5AC", "TFF3")

# Set 2. Epithelial restitution and Wound healing
set2 <- c("ANXA2", "CDH1", "CDX2",
          "GJA1", "MKI-67", "MMP2",
          "MMP29", "NANOG", "PCNA",
          "POU5F1", "RND3", "TLR2",
          "TLR4", "TLR5", "TLR9",
          "VIL1", "WISP1") 

# Set 3. Immune modulation
set3 <- c("L12p70", "IFNG", "IL4", 
          "IL1B", "IL17A", "IL21",
          "IL22", "IL9")

#######################################################
##                   Load functions                  ##
#######################################################
source(file.path(prg.dir, "00_DE_functions.R"))

#######################################################
##                    Reading data                   ##
#######################################################

count.matrix <- fread(in.counts.file, header = TRUE) %>% 
  column_to_rownames("V1") %>%
  mutate_all(as.integer) %>%
  as.matrix() 

colnames(count.matrix) <- substr(colnames(count.matrix),24,31)
rownames(count.matrix) <- gsub("\\..*","",rownames(count.matrix))

clin.data <- fread(sample.annot, header = TRUE)
colnames(clin.data) <- gsub(" ","_",colnames(clin.data))

clin.data <- clin.data %>%
  rename(SampleID = "Sample_ID_&_Suffix",
         phenotype = Treatment) %>%
  mutate(Timepoint = gsub(" ","_",Timepoint)) %>%
  mutate(treatTime = paste0(phenotype, "_", Timepoint)) %>%
  mutate(sample_name = substr(SampleID,11,18)) %>%
  left_join(resp.annot %>%
              dplyr::select(sample_name,SUBJ_ID=Subject.ID,Response,cResponse,cRemission,eResponse)) %>%
  column_to_rownames("sample_name")

colnames(clin.data) <- gsub(" ","_",colnames(clin.data))

## Matrix columns and clin.data rows must be in the same order.
## Filter sample IDs in clin.data and count matrix
## Sort colnames in matrix and rownames in clin.data
clin.data <- clin.data %>%
  dplyr::filter(rownames(.) %in% colnames(count.matrix)) %>%
  arrange(match(rownames(.), colnames(count.matrix)))

#######################################################
##                  Data retrieving                  ##
#######################################################
# load biomart data to get information about gene symbols and protein-coding status later on;
# make sure correct species is being queried (below query is for human)
listEnsembl()
ensembl <- useEnsembl(biomart = "ensembl", 
                      dataset = "hsapiens_gene_ensembl", 
                      mirror = "useast")

# Data Input
genemap <- fread(file.path(in.dir, "mock_data", "100_genesmapped2symbols.csv"), header = TRUE, stringsAsFactors = FALSE) %>%
  dplyr::select(-V1) %>%
  mutate(gene.id = gsub("\\..*","", gene_id))

# Get gene ids to add to results later on
gene.ids <- getBM(attributes = c("ensembl_gene_id",'gene_biotype',"start_position","end_position"),
                  filters = 'ensembl_gene_id',
                  values = rownames(count.matrix),
                  mart = ensembl)

gene.ids <- gene.ids %>%
  rename(gene.id = ensembl_gene_id) %>%
  left_join(genemap, by = c("gene.id")) %>%
  mutate(gene_length = abs(start_position - end_position))

#######################################################
##                  Analysis setting                 ##
#######################################################
## Analysis variables
problematic_samples <- c("UDI_0059","UDI_0048")
output.design <- analysis.desig 

## Filter clin.data by post.baseline 
## and problematic samples
if(analysis.desig == "INT"){
  clin.data.pBL <- clin.data %>%
    dplyr::filter(!rownames(.) %in% problematic_samples)
}

## Filter count.matrix samples
count.matrix.pBL <- count.matrix[,rownames(clin.data.pBL)]

## Remove duplicated genes due to name conversion
count.matrix.pBL <- count.matrix.pBL[!duplicated(rownames(count.matrix.pBL)),]
gene.ids <- gene.ids[!duplicated(gene.ids$gene.id),]

## Create sampleTable object to annotate deseqdataset
sampleTable <- data.frame(phenotype = factor(clin.data.pBL$Response),
                          cResponse = factor(ifelse(clin.data.pBL$cResponse==1,"R","NR")),
                          cRemission = factor(ifelse(clin.data.pBL$cRemission==1,"R","NR")),
                          eResponse = factor(ifelse(clin.data.pBL$eResponse==1,"R","NR")),
                          Treatment = factor(clin.data.pBL$phenotype),
                          timepoint = factor(clin.data.pBL$Timepoint),
                          SUBJID = factor(clin.data.pBL$SUBJ_ID),
                          sample_name = colnames(count.matrix.pBL))
rownames(sampleTable) <- sampleTable$sample_name

## Relevel factors as reference
sampleTable$timepoint<-relevel(sampleTable$time,ref="Baseline")
sampleTable$cResponse<-relevel(sampleTable$cResponse,ref="NR")
sampleTable$cRemission<-relevel(sampleTable$cRemission,ref="NR")
sampleTable$eResponse<-relevel(sampleTable$eResponse,ref="NR")
sampleTable$phenotype <- relevel(sampleTable$phenotype,ref="Responder")

if(analysis.desig == "INT"){
  if(trt_covar){
    
    ## Create design matrix
    design_f <- model.matrix(formula(paste0(paste("~", res_phenotype),paste0("+Treatment"),
                                     paste("+", res_phenotype),paste0(":SUBJID"),
                                     paste("+", res_phenotype),paste0(":timepoint"))), sampleTable)
    
    ## Avoid levels without samples (fix collinearity)
    all.zero <- apply(design_f, 2, function(x) all(x==0))

    ## Remove columns with zero
    idx <- which(all.zero)
    design_f <- design_f[,-idx]
    design_f <- design_f[,!(colnames(design_f) %in% nonEstimable(design_f))]
    
    # design_f <- paste0("~ ",res_phenotype," + timepoint + Treatment + ",res_phenotype,":timepoint")
    reslab <- "TRT"
  }else{
    
    sampleTable <- sampleTable %>% dplyr::filter(Treatment=="MH002")
    count.matrix.pBL <- count.matrix.pBL[,colnames(count.matrix.pBL) %in% sampleTable$sample_name ]

    ## Create design matrix
    design_f <- model.matrix(formula(paste0(paste("~", res_phenotype),
                                     paste("+", res_phenotype),paste0(":SUBJID"),
                                     paste("+", res_phenotype),paste0(":timepoint"))), sampleTable)
    
    ## Avoid levels without samples (fix collinearity)
    all.zero <- apply(design_f, 2, function(x) all(x==0))
    
    ##Remove columns with zero
    idx <- which(all.zero)
    design_f <- design_f[,-idx]
    design_f <- design_f[,!(colnames(design_f) %in% nonEstimable(design_f))]
    
    # design_f <- paste0("~ ",res_phenotype," + timepoint + ",res_phenotype,":timepoint")
    reslab <- "noTRT"
  }
}

out.dir <- file.path(out.dir, output.design)
if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}

#######################################################
##                  Data formatting                  ##
#######################################################
## Prepare count matrix
count.matrix.pBL <- round(count.matrix.pBL)
mode(count.matrix.pBL) <- "integer"

## Check genes unexpressed
zero_unexpressed <- (apply(count.matrix.pBL, 1, max) == 0)

## Remove them from the matrix
count.matrix.pBL <- count.matrix.pBL[!zero_unexpressed,]

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(count.matrix.pBL, 
                              sampleTable, 
                              design_f)

# Estimate size fators and dispersions
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

# Pre-filtering step to remove genes with low read counts. An additional
# more restrictive filtering process is performed in the result() function
min.counts <- 5
keep <- rowSums(counts(dds, normalized=TRUE) >= min.counts) >= min(table(sampleTable[,res_phenotype]))  
dds <- dds[keep,]

#######################################################
##                     Perform DE                    ##
#######################################################
res0 <- performDE(data = "dds",
                  phenotype = res_phenotype,
                  design = design_f,
                  test = "Wald", 
                  contrast = c("R.timepointWeek_8","NR.timepointWeek_8"),
                  annotation = gene.ids, 
                  padj.thr = 0.05,
                  ncores = 4)

#######################################################
##               Export results setting              ##
#######################################################
if(analysis.desig == "INT"){
  timepoint <- "Baseline/Week_8"
  timepointFilt <- timepointFiltINT 
  out.prefix <- paste0("NR_vs_R_050923", "_", reslab)
  sampleIDsFilt <- sampleTable %>% 
    dplyr::filter(timepoint == timepointFilt) %>%
    pull(sample_name)
}

######################################################
##               Extract biomarker scores           ##
######################################################
clin.data.biomarker <- read_sas(file.path(in.dir,"ADaM_dataset","adlb.sas7bdat")) %>%
  dplyr::filter(PARAMCD %in% c("CRP","TNFALPH","IFNGAMM","IL6","IL8","IL10"),
                RANDFL == "Y",
                ANL01FL == "Y",
                grepl("Baseline|Week 8", AVISIT),
                BASETYPE == "PERIOD 1") %>% 
  dplyr::select(USUBJID,SUBJID,AVAL,AVISIT,PARAMCD) %>%
  mutate(AVISIT = if_else(grepl("Baseline",AVISIT),"Baseline","Week_8"))

clin.data.biomarker.feccalp <- read_sas(file.path(in.dir,"ADaM_dataset","admb.sas7bdat")) %>%
  dplyr::filter(PARAMCD %in% c("FECCALP"),
                RANDFL == "Y",
                ANL01FL == "Y",
                grepl("Baseline|Week 8", AVISIT),
                BASETYPE == "PERIOD 1") %>% 
  dplyr::select(USUBJID,SUBJID,AVAL,AVISIT,PARAMCD) %>%
  mutate(AVISIT = if_else(grepl("Baseline",AVISIT),"Baseline","Week_8"))

clin.data.biomarker.hist <- read_sas(file.path(in.dir,"ADaM_dataset","adhs.sas7bdat")) %>%
  dplyr::filter(PARAMCD %in% c("NANCTSV","GEBOES","ROBARTS"),
                RANDFL == "Y",
                ANL01FL == "Y",
                grepl("Baseline|Week 8", AVISIT)) %>% 
  dplyr::select(USUBJID,SUBJID,AVAL,AVISIT,PARAMCD) %>%
  mutate(AVISIT = if_else(grepl("Baseline",AVISIT),"Baseline","Week_8"))

clin.data.biomarker.mes <- read_sas(file.path(in.dir,"ADaM_dataset","admes.sas7bdat")) %>%
  dplyr::filter(PARAMCD %in% c("MES"),
                RANDFL == "Y",
                ANL01FL == "Y",
                grepl("Baseline|Week 8", AVISIT)) %>% 
  dplyr::select(USUBJID,SUBJID,AVAL,AVISIT,PARAMCD) %>%
  mutate(AVISIT = if_else(grepl("Baseline",AVISIT),"Baseline","Week_8"))

clin.data.biomarker <- rbind(clin.data.biomarker, clin.data.biomarker.feccalp, 
                             clin.data.biomarker.hist, clin.data.biomarker.mes)

## Plot biomarker values
biomark_hist <- ggplot() + geom_histogram(clin.data.biomarker, mapping = aes(x = AVAL, fill = AVISIT), colour = "black") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) +
  facet_wrap(vars(PARAMCD), scales = "free")

title <- paste0("5.3.1 Primary Objective: Analysis of Biomarkers associated with the mechanistic effects of MH002 at ",timepointFilt,". Histogram.",
                "\nITT population (N = ",length(dds$sample_name),")", 
                "\nBiomarker Statistical Analysis Plan for MH002 in subjects with mild to moderate Ulcerative Colitis")
foot1 <- paste0("\nNote: This plot shows the distribution of the different biomarker scores used in the analysis")
foot2 <- paste0("\nProgram Path: ", file.path(prg.dir, prg.file))
foot3 <- paste("Status: Under review", sep=" ")
foot4 <- paste("Date: ",Sys.Date(), sep=" ")

out.resp <- file.path(out.dir,res_phenotype)
if(!dir.exists(out.resp)){dir.create(out.resp)}

plot.to.tfl(p.ggplot = biomark_hist,
            scale.format.x = 1,
            scale.mar.up = 0.8,
            scale.mar.right = 0.2,
            scale.mar.bottom = 0.8,
            scale.mar.left = 0.2,
            tfl.header = title, 
            tfl.foot.1 = foot1, 
            tfl.foot.2 = foot2, 
            tfl.foot.3 = foot3, 
            tfl.foot.4 = foot4,
            outDir = out.resp, 
            prefix = paste0(out.prefix,"_Biomark_Dist"))

######################################################
##              Save list of results                ##
######################################################

save(res0, file = file.path(out.resp, paste0("RNR_", res_phenotype, "_", trt_covar, "_DEG_res.RData")))

######################################################
##                    Correlations                  ##
######################################################
## Filter results
## clean of pseudogenes
library(EnsDb.Hsapiens.v86)
ts <- transcripts(EnsDb.Hsapiens.v86)
ts.pseudo <- grep("pseudo|antisense|lincRNA", unique(ts$tx_biotype), val=TRUE)
ts <- as.data.frame(ts)

## gene genes id pseudo
genes.pseudo <- ts %>%
  dplyr::filter(tx_biotype %in% ts.pseudo) %>%
  pull(gene_id) %>% unique()

## Filter out pseudo genes
my.res <- res0$DE_res %>%
  dplyr::filter(!(gsub("\\..*","",gene.id) %in% genes.pseudo))

## Filter out non-interesting genes
my.res <- my.res %>%
  dplyr::filter(!grepl("^RP|^CTD|^CTB|^AC|^ENSG|^CTC|^KB|^XX|^y|^AL|^RN|^CTA|^Y_",gene_name)) %>%
  mutate_if(is.numeric,round,digits=5)

##UP
sorted.up <- my.res %>% 
  arrange(desc(log2FoldChange)) %>% 
  dplyr::filter(log2FoldChange > 0)

##DOWN
sorted.down <- my.res %>% 
  arrange(log2FoldChange) %>% 
  dplyr::filter(log2FoldChange < 0)

## Select interesting genes
intGenes <- c("MMP2","POU5F1","CLDN1","MUC2","NANOG","ANXA2")

tableResults <- data.frame()
cor.df.all <- data.frame()
wilcox.res.all <- data.frame()
wilcox.res.all_RvsNR <- data.frame()
for(top_genes in c("top20up","top20down",intGenes)){
  if(top_genes == "top20up"){
    filtGenesToPlot <- sorted.up %>% 
      dplyr::filter(log2FoldChange > 2) %>% 
      pull(gene.id)
    
    biomarkerLab <- "Top20"
    
    transDf <- res0$DE_norm_counts[filtGenesToPlot,] %>% t() 
    
    for(x in colnames(transDf)){
      trans_v <- bestNormalize::yeojohnson(transDf[,x], standardize = T)
      transDf[,x] <- trans_v$x.t
    }
    
    norm.expr.corr <- transDf %>% t() %>%
      as.data.frame() %>%
      rownames_to_column("gene.id") %>%
      pivot_longer(cols = colnames(.)[-1], names_to = "sample_name", values_to = "NormCounts")
    
  }else if(top_genes == "top20down"){
    filtGenesToPlot <- sorted.down %>% 
      dplyr::filter(log2FoldChange < 2) %>%
      pull(gene.id)
    
    biomarkerLab <- "Top20"
    
    transDf <- res0$DE_norm_counts[filtGenesToPlot,] %>% t() 
    
    for(x in colnames(transDf)){
      trans_v <- bestNormalize::yeojohnson(transDf[,x], standardize = T)
      transDf[,x] <- trans_v$x.t
    }
    
    norm.expr.corr <- transDf %>% t() %>%
      as.data.frame() %>%
      rownames_to_column("gene.id") %>%
      pivot_longer(cols = colnames(.)[-1], names_to = "sample_name", values_to = "NormCounts")
  }else{
    filtGenesToPlot <- gene.ids %>%
      dplyr::filter(gene_name %in% intGenes) %>%
      dplyr::select(gene.id,gene_name) %>%
      mutate(gene.id = gsub("\\..*","",gene.id))

    biomarkerLab <- "Gene"
    
    geneExtract <- filtGenesToPlot %>% dplyr::filter(gene_name == top_genes) %>% pull(gene.id)
    
    transDf <- res0$DE_norm_counts %>% t() 
    
    for(x in colnames(transDf)){
      trans_v <- bestNormalize::yeojohnson(transDf[,x], standardize = T)
      transDf[,x] <- trans_v$x.t
    }
    
    norm.expr.corr <- transDf %>% t() %>%
      as.data.frame() %>%
      rownames_to_column("gene.id") %>%
      pivot_longer(cols = colnames(.)[-1], names_to = "sample_name", values_to = "NormCounts") %>%
      dplyr::filter(gene.id == geneExtract)
    
    if(nrow(norm.expr.corr) == 0){next}
  }
  
  ## Incorporate subject ID
  norm.expr.corr <- norm.expr.corr %>%
    left_join(resp.annot %>% dplyr::select(SUBJID="Subject.ID","sample_name")) %>%
    left_join(sampleTable %>% dplyr::select(sample_name,timepoint,all_of(res_phenotype)))
  
  listplots <- list()
  for(i in unique(clin.data.biomarker$PARAMCD)){
    
    ## Select biomarker 
    clin.data.biomarker2 <- clin.data.biomarker %>% 
      dplyr::filter(PARAMCD == i) 
    
    ## Calculate z scores on biomarker scores
    #z_scores <- (clin.data.biomarker2$AVAL-mean(clin.data.biomarker2$AVAL))/sd(clin.data.biomarker2$AVAL)
    
    trans_scores <- bestNormalize::yeojohnson(as.numeric(clin.data.biomarker2$AVAL))
    
    clin.data.biomarker2 <- clin.data.biomarker2 %>% 
      mutate(AVAL = trans_scores$x.t)
    
    ## Merge with expression data and select "Responders"
    clin.data.biomarker2_R_NR <- clin.data.biomarker2 %>% 
      left_join(norm.expr.corr, by = c("SUBJID", "AVISIT"="timepoint")) %>%
      dplyr::filter(!is.na(sample_name)) %>%
      as.data.frame() %>%
      dplyr::rename(timepoint = AVISIT)
    
    clin.data.biomarker2 <- clin.data.biomarker2 %>% 
      left_join(norm.expr.corr, by = c("SUBJID", "AVISIT"="timepoint")) %>%
      as.data.frame() %>%
      dplyr::filter(get(res_phenotype) == "R") %>%
      dplyr::rename(timepoint = AVISIT)
    
    ## Export table with correlation results
    if(!grepl("top",top_genes)){
      expResTable <- clin.data.biomarker2 %>%
        dplyr::mutate(gene.id = top_genes) %>%
        dplyr::rename(bmk_trans = AVAL,
                      normExp_trans = NormCounts,
                      Biomarker = PARAMCD,
                      gene_name = gene.id) %>%
        left_join(res0$DE_res_all %>% 
                    dplyr::filter(gene_name == top_genes) %>%
                    dplyr::select(gene_name,stat,pvalue), by = c("gene_name"))
      
      cor.df.time <- data.frame()
      for(j in c("Baseline","Week_8")){
        expResTable_tmp <- expResTable %>% dplyr::filter(timepoint == j)
        cor.df <- NULL
        try(cor.df <- cor.test(x = expResTable_tmp$bmk_trans, y = expResTable_tmp$normExp_trans, method = "pearson"),silent = TRUE)
        if(is.null(cor.df)){
          cor.df <- data.frame(Gene = top_genes, Biomarker = i, timepoint = j, cor.estimate = NA, stat = NA, pvalue = NA)
        }else{
          cor.df <- data.frame(Gene = top_genes, Biomarker = i, timepoint = j, cor.estimate = cor.df$estimate, stat = cor.df$statistic, pvalue = cor.df$p.value)
        }
        cor.df.time <- rbind(cor.df, cor.df.time)
      }
      
      tableResults <- rbind(expResTable,tableResults)
      cor.df.all <- rbind(cor.df.time, cor.df.all)
    }
    
    ## Export table with Wilcoxon test on R (BL vs W8)
    if(!grepl("top",top_genes)){
      expResTable <- clin.data.biomarker2 %>%
        dplyr::mutate(gene.id = top_genes) %>%
        dplyr::rename(bmk_trans = AVAL,
                      normExp_trans = NormCounts,
                      Biomarker = PARAMCD,
                      gene_name = gene.id) 
      
      res_BL <- expResTable %>% dplyr::filter(timepoint == "Baseline")
      res_w8 <- expResTable %>% dplyr::filter(timepoint == "Week_8")
    
      bmk_res <- wilcox.test(as.numeric(res_BL$bmk_trans), 
                  as.numeric(res_w8$bmk_trans), exact = FALSE)
      
      gene_res <- wilcox.test(as.numeric(res_BL$normExp_trans), 
                  as.numeric(res_w8$normExp_trans), exact = FALSE)
      
      wilcox.res <- data.frame(Outcome = res_phenotype, Gene = top_genes, Biomarker = i,
                 Wilcox.stat.bmk =  bmk_res$statistic, pvalue.bmk = bmk_res$p.value,
                 Wilcox.stat.gene = gene_res$statistic, pvalue.gene = gene_res$p.value)
      
      wilcox.res.all <- rbind(wilcox.res,wilcox.res.all)
    }
    
    ## Export table with Wilcoxon test on R vs NR (BL and W8)
    if(!grepl("top",top_genes)){
      expResTable <- clin.data.biomarker2_R_NR %>%
        dplyr::mutate(gene.id = top_genes) %>%
        dplyr::rename(bmk_trans = AVAL,
                      normExp_trans = NormCounts,
                      Biomarker = PARAMCD,
                      gene_name = gene.id) 
      
      res_RBL <- expResTable %>% dplyr::filter(timepoint == "Baseline", get(res_phenotype) == "R")
      res_NRBL <- expResTable %>% dplyr::filter(timepoint == "Baseline", get(res_phenotype) == "NR")
      res_RW8 <- expResTable %>% dplyr::filter(timepoint == "Week_8", get(res_phenotype) == "R")
      res_NRW8 <- expResTable %>% dplyr::filter(timepoint == "Week_8", get(res_phenotype) == "NR")
      
      bmk_BL_res <- wilcox.test(as.numeric(res_RBL$bmk_trans), 
                                as.numeric(res_NRBL$bmk_trans), exact = FALSE)
      
      bmk_W8_res <- wilcox.test(as.numeric(res_RW8$bmk_trans), 
                                as.numeric(res_NRW8$bmk_trans), exact = FALSE)
      
      gene_BL_res <- wilcox.test(as.numeric(res_RBL$normExp_trans), 
                                 as.numeric(res_NRBL$normExp_trans), exact = FALSE)
      
      gene_W8_res <- wilcox.test(as.numeric(res_RW8$normExp_trans), 
                                 as.numeric(res_NRW8$normExp_trans), exact = FALSE)
      
      wilcox.res <- data.frame(Outcome = res_phenotype, Gene = top_genes, Biomarker = i,
                               Wilcox.stat.bmk.BL =  bmk_BL_res$statistic, pvalue.bmk.BL = bmk_BL_res$p.value,
                               Wilcox.stat.gene.BL = gene_BL_res$statistic, pvalue.gene.BL = gene_BL_res$p.value,
                               Wilcox.stat.bmk.W8 =  bmk_W8_res$statistic, pvalue.bmk.W8 = bmk_W8_res$p.value,
                               Wilcox.stat.gene.W8 = gene_W8_res$statistic, pvalue.gene.W8 = gene_W8_res$p.value)
      
      wilcox.res.all_RvsNR <- rbind(wilcox.res,wilcox.res.all_RvsNR)
    }
    
    
    ## Reformat 
    clin.data.biomarker3 <- unique(rbind(clin.data.biomarker2 %>% 
                                           dplyr::select(SUBJID,BIOMARKER=PARAMCD,score=AVAL,timepoint,phenotype=all_of(res_phenotype)) %>%
                                           mutate(BIOMARKER = "BMK"),
                                         
                                         clin.data.biomarker2 %>% 
                                           dplyr::select(SUBJID,sample_name,timepoint,BIOMARKER=gene.id,score=NormCounts,timepoint,phenotype=all_of(res_phenotype)) %>%
                                           group_by(SUBJID,timepoint,phenotype) %>% 
                                           dplyr::summarise(BIOMARKER = biomarkerLab, score = mean(score)) %>% 
                                           dplyr::select(SUBJID,BIOMARKER,score,timepoint,phenotype))) %>%
      arrange(BIOMARKER) %>%
      dplyr::filter(!is.na(score),
                    !is.na(phenotype)) %>%
      mutate(timepoint = if_else(timepoint == "Baseline","BL","W8"))
    
    ## Plot
    pCorr <- ggplot(clin.data.biomarker3) + 
      geom_boxplot(mapping = aes(x = BIOMARKER, y = score, fill = BIOMARKER)) + 
      geom_jitter(mapping = aes(x = BIOMARKER, y = score), colour = "black", size = 0.8) +
      theme_bw() + ggtitle(paste0(top_genes," vs ",i, " (Responders)")) +
      scale_fill_manual(values = c("Top20" = "#0072B2", "Gene" = "#0072B2", "BMK" = "#D55E00")) +
      theme(axis.text.x = element_text(angle = 0, size = 8),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(hjust = 0.5,size = 8),
            strip.text = element_text(size = 7)) + 
      facet_grid(cols = vars(timepoint))
    
    listplots <- c(listplots, list(pCorr))
  }
  
  title <- paste0("5.3.1 Primary Objective: Analysis of Biomarkers associated with the mechanistic effects of MH002 at ",timepointFilt,". Boxplot.",
                  "\nITT population (N = ",length(dds$sample_name),")", 
                  "\nBiomarker Statistical Analysis Plan for MH002 in subjects with mild to moderate Ulcerative Colitis")
  foot1 <- paste0("\nNote: This plot shows the comparison of the z-score transformed of the top DE genes normalized counts and the biomarker scores in the", 
                  "\nITT population at both timepoints: Baseline and Week 8.",
                  "\nAbbreviations: ITT = intent-to-treatment")
  foot2 <- paste0("\nProgram Path: ", file.path(prg.dir, prg.file))
  foot3 <- paste("Status: Under review", sep=" ")
  foot4 <- paste("Date: ",Sys.Date(), sep=" ")
  
  out.resp <- file.path(out.dir,res_phenotype)
  if(!dir.exists(out.resp)){dir.create(out.resp)}
  
  plot.to.tfl(p.ggplot = ggpubr::ggarrange(plotlist = listplots, common.legend = TRUE),
              scale.format.x = 1,
              scale.mar.up = 0.65,
              scale.mar.right = 0.2,
              scale.mar.bottom = 1,
              scale.mar.left = 0.2,
              tfl.header = title, 
              tfl.foot.1 = foot1, 
              tfl.foot.2 = foot2, 
              tfl.foot.3 = foot3, 
              tfl.foot.4 = foot4,
              outDir = out.resp, 
              prefix = paste0(out.prefix,"_", res_phenotype, "_",top_genes,"_Corr"))
}

#######################################################
##               Generate Table TFL                  ##
#######################################################
## RTF output of transformed scores
title <- paste0("5.3.1 Primary Objective: Analysis of Biomarkers associated with the mechanistic effects of MH002. Summary table.",
                "\nITT population (N = ",length(dds$sample_name),")",  
                "\nBiomarker Statistical Analysis Plan for MH002 in subjects with mild to moderate Ulcerative Colitis")
foot1 <- paste0("\nNote: This is a summary table of Zscore-based Biomarker score and gene expressions.")
foot2 <- paste0("\nProgram Path: ", file.path(prg.dir, "01_DE_DEseq2.R"))
foot3 <- paste("Status: Under review", sep=" ")
foot4 <- paste("Date: ",Sys.Date(), sep=" ")

tableResults <- tableResults %>% 
  dplyr::mutate(bmk_zscore = round(bmk_trans,digits = 2),
                normExp_zscore = round(normExp_trans,digits = 2),
                pvalue = round(pvalue,digits = 2))

table.rtf(tab = tableResults,
          rtf.header = title,
          rtf.foot.1 = foot1,
          rtf.foot.2 = foot2,
          rtf.foot.3 = foot3,
          rtf.foot.4 = foot4,
          colW = 0.85,
          outDir = out.resp,
          prefix = paste0(out.prefix,"_", res_phenotype, "_summary"))

## RTF output of correlation results
title <- paste0("5.3.1 Primary Objective: Analysis of Biomarkers associated with the mechanistic effects of MH002. Summary table.",
                "\nITT population (N = ",length(dds$sample_name),")",  
                "\nBiomarker Statistical Analysis Plan for MH002 in subjects with mild to moderate Ulcerative Colitis")
foot1 <- paste0("\nNote: This is a summary table of the correlations tests made between genes and biomarkers using the Responders population (",res_phenotype,").")
foot2 <- paste0("\nProgram Path: ", file.path(prg.dir, "01_DE_DEseq2.R"))
foot3 <- paste("Status: Under review", sep=" ")
foot4 <- paste("Date: ",Sys.Date(), sep=" ")

cor.df.all <- cor.df.all %>% 
  dplyr::mutate(cor.estimate = round(cor.estimate,digits = 2),
                stat = round(stat,digits = 2),
                pvalue = round(pvalue,digits = 2))

table.rtf(tab = cor.df.all,
          rtf.header = title,
          rtf.foot.1 = foot1,
          rtf.foot.2 = foot2,
          rtf.foot.3 = foot3,
          rtf.foot.4 = foot4,
          colW = 1.2,
          outDir = out.resp,
          prefix = paste0(out.prefix,"_", res_phenotype, "_corr.summary"))

## RTF output of wicoxon test results
title <- paste0("5.3.1 Primary Objective: Analysis of Biomarkers associated with the mechanistic effects of MH002. Summary table.",
                "\nITT population (N = ",length(dds$sample_name),")",  
                "\nBiomarker Statistical Analysis Plan for MH002 in subjects with mild to moderate Ulcerative Colitis")
foot1 <- paste0("\nNote: This is a summary table of the correlations tests made between genes and biomarkers using the Responders population (",res_phenotype,").")
foot2 <- paste0("\nProgram Path: ", file.path(prg.dir, "01_DE_DEseq2.R"))
foot3 <- paste("Status: Under review", sep=" ")
foot4 <- paste("Date: ",Sys.Date(), sep=" ")

wilcox.res.all <- wilcox.res.all %>% 
  dplyr::mutate(pvalue.bmk = round(pvalue.bmk,digits = 2),
                pvalue.gene = round(pvalue.gene,digits = 2))

table.rtf(tab = wilcox.res.all,
          rtf.header = title,
          rtf.foot.1 = foot1,
          rtf.foot.2 = foot2,
          rtf.foot.3 = foot3,
          rtf.foot.4 = foot4,
          colW = 1.2,
          outDir = out.resp,
          prefix = paste0(out.prefix,"_", res_phenotype, "_wilcox.summary"))

## RTF output of wicoxon test results (R vs NR)
title <- paste0("5.3.1 Primary Objective: Analysis of Biomarkers associated with the mechanistic effects of MH002. Summary table.",
                "\nITT population (N = ",length(dds$sample_name),")",  
                "\nBiomarker Statistical Analysis Plan for MH002 in subjects with mild to moderate Ulcerative Colitis")
foot1 <- paste0("\nNote: This is a summary table of the correlations tests made between genes and biomarkers using the Responders population (",res_phenotype,").")
foot2 <- paste0("\nProgram Path: ", file.path(prg.dir, "01_DE_DEseq2.R"))
foot3 <- paste("Status: Under review", sep=" ")
foot4 <- paste("Date: ",Sys.Date(), sep=" ")

wilcox.res.all_RvsNR <- wilcox.res.all_RvsNR %>% 
  dplyr::mutate(pvalue.bmk.BL = round(pvalue.bmk.BL,digits = 2),
                pvalue.gene.BL = round(pvalue.gene.BL,digits = 2),
                pvalue.bmk.W8 = round(pvalue.bmk.W8,digits = 2),
                pvalue.gene.W8 = round(pvalue.gene.W8,digits = 2))

table.rtf(tab = wilcox.res.all_RvsNR,
          rtf.header = title,
          rtf.foot.1 = foot1,
          rtf.foot.2 = foot2,
          rtf.foot.3 = foot3,
          rtf.foot.4 = foot4,
          colW = 0.9,
          outDir = out.resp,
          prefix = paste0(out.prefix,"_", res_phenotype, "_wilcox_RvsNR.summary"))

#######################################################
##          Generate heatmap.2x plot in TFL          ##
#######################################################
for(timepointFilt in c("Week_8", "Baseline")){
  title <- paste0("5.3.1 Primary Objective: Analysis of Biomarkers associated with the mechanistic effects of MH002 at ",timepointFilt,". Heatmap plot.",
                "\nITT population (N = ",length(dds$sample_name),")",
                "\nBiomarker Statistical Analysis Plan for MH002 in subjects with mild to moderate Ulcerative Colitis")
  foot1 <- paste0("\nNote: This plot shows the z-score values of the normalized counts for the significantly expressed genes (nominal pvalue) in the ITT population",
                  "\n at timepoint: ",timepointFilt,". The different sample groups are represented at the top of the heatmap as coloured squares",
                  "\n (PBO = yellow, MH002 = grey).",
                  "\nAbbreviations: ITT = intent-to-treatment")
  foot2 <- paste0("\nProgram Path: ", file.path(prg.dir, prg.file))
  foot3 <- paste("Status: Under review", sep=" ")
  foot4 <- paste("Date: ",Sys.Date(), sep=" ")
  
  ## Genes to plot
  top.DE.genes <- res0$DE_res %>%
    dplyr::filter(pvalue < 0.01) %>%
    pull(gene.id)
  
  ## Format sampleTable
  sampleTable_tmp <- sampleTable %>% 
    dplyr::mutate(phenotype = get(res_phenotype))
  
  ## Extract normalized counts
  sampleIDsFilt <- sampleTable_tmp %>% 
    dplyr::filter(timepoint == timepointFilt) %>%
    pull(sample_name)
  
  M <- res0$DE_norm_counts[top.DE.genes,sampleIDsFilt]
  
  ## Plot heatmap.2
  heatmap.2.tfl(exp.matrix = M,
                zscore = TRUE,
                sampleTable = sampleTable_tmp[sampleIDsFilt,],
                title = "Responders vs Non-responders",
                tfl.header = title,
                tfl.foot.1 = foot1,
                tfl.foot.2 = foot2,
                tfl.foot.3 = foot3,
                tfl.foot.4 = foot4,
                outDir = out.resp,
                prefix = paste0(out.prefix,"_",timepointFilt))
}

#######################################################
##               Retrieve R Session info             ##
#######################################################
(session.info <- list(sessionInfo = sessionInfo(), options = options(), RNGkind = RNGkind(), Sys.info = Sys.info(), Date = date()))
