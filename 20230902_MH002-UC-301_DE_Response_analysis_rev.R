
  ###########################################################################
# Project           : MRM Health - Statistical support for MH002-UC-201 - PLX-23-35515 (1545)
# Program name      : 20230902_MH002-UC-301_DE_Response_analysis.R
# Developed in      : R 4.2.0
# Purpose           : Differential expression analysis  (R vs NR) of the 3 response measures included the SAP of MH002-UC-201
##                   Complementary analysis considering nested model contrlloing by SUBJID following Love et al
##                  http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#group-specific-condition-effects-individuals-nested-within-groups
# Inputs            : RSEM Expression values located at Transcriptomics data
#   
# Outputs           : Excel files with Differentially Expressed Genes (DEG) and their enrichment analysis
# Revision History  :
#   Version   Date        Author       Revision 
#   -------   ---------   ----------   ---------------------------------
#   1.0       31082023    Antonio Gomez  Creation
# Reviewed by       : Pablo Villegas , 05/09/2023
# Reference number  : NA                                      
# Validation Level  : Complete: Generation of Results
# Declaration of Confidentiality:
#   - I declare that this program may contain or contain confidential
#     information and that it cannot be shared as it.
###########################################################################



library(tidyverse)
library(data.table)
#####GO enrichment analysis
library(clusterProfiler)
library(org.Hs.eg.db)
library(GO.db)
library(DOSE)
library(annotate)
library(dplyr)
library(enrichplot)
library(EnsDb.Hsapiens.v86)
library(rtf)
library(DESeq2)
library(apeglm)
library(grid)
library(xlsx)
library(ReactomePA)
library(limma)
library(msigdbr)
###paths

base.dir <- "Use_your_path"

#counts
in.counts.file <- file.path(base.dir, "rsem_gcounts_matrix.txt")
#Samplesheet
sample.annot <- file.path(base.dir, "MRM_responders.csv")

##load gene annot from vendor

load(file.path(base.dir,"Annot_Final.RData"))

##load data

count.matrix <- fread(in.counts.file, header = TRUE) %>% 
  column_to_rownames("V1") %>%
 # mutate_all(as.integer) %>%
  as.matrix() 

colnames(count.matrix) <- substr(colnames(count.matrix),24,31)

##load samplesheet
clin.data <- fread(sample.annot, header = TRUE)
rownames(clin.data)<-clin.data$SAMPLE_ID
## Matrix columns and clin.data rows must be in the same order.
## Filter sample IDs in clin.data and count matrix
## Sort colnames in matrix and rownames in clin.data
clin.data <- clin.data %>% 
  dplyr::filter(rownames(.) %in% colnames(count.matrix)) %>%
  arrange(match(rownames(.), colnames(count.matrix)))
clin.data$Time<-ifelse(clin.data$Timepoint=="Baseline","Baseline","Week8")

## Filter count.matrix samples. Get the ones that are in the samplesheet
count.matrix <- count.matrix[,colnames(count.matrix) %in% rownames(clin.data)]


## Create sampleTable object to annotate deseqdataset

sampleTable <- data.frame(C_Response = factor(ifelse(clin.data$cResponse==1,"R","NR")),
                          Remission=factor(ifelse(clin.data$cRemission==1,"R","NR")),
                          E_Response=factor(ifelse(clin.data$eResponse==1,"R","NR")),
                          Treatment=clin.data$Treatment,
                          time=factor(clin.data$Time),
                          SUBJID=factor(clin.data$SUBJ_ID),
                          sample_name = colnames(count.matrix))

rownames(sampleTable) <- sampleTable$sample_name


##Prepare count matrix. RSEM values must be rounded to be used
##Approach from tximport package (from )Soneson et al (2015), using tximport function's part regarding RSEM
count.matrix<-round(count.matrix)
mode(count.matrix)<-"integer"

##remove UDI_0059 (because QC issues) and its neighbor in wk8, UDI_0048

count.matrix<-count.matrix[,!(colnames(count.matrix) %in% c("UDI_0059","UDI_0048"))]
sampleTable<-sampleTable[!(sampleTable$sample_name %in% c("UDI_0059","UDI_0048")),]

##Relevel factors as reference
sampleTable$time<-relevel(sampleTable$time,ref="Baseline")
sampleTable$C_Response<-relevel(sampleTable$C_Response,ref="NR")
sampleTable$Remission<-relevel(sampleTable$Remission,ref="NR")
sampleTable$E_Response<-relevel(sampleTable$E_Response,ref="NR")


##Check genes unexpressed

zero_unexpressed <- (apply(count.matrix, 1, max) == 0)

##remove them from the matrix

count.matrix<-count.matrix[!zero_unexpressed,]

####Clean of pseudogenes. We are not interested in pseudogenes, antisense and/or lincRNA from now
###Build a list of them from ENCODE data an use later to filter the resuls

ts <- transcripts(EnsDb.Hsapiens.v86)
ts.pseudo<-grep("pseudo|antisense|lincRNA", unique(ts$tx_biotype), val=TRUE)
ts<-as.data.frame(ts)


##gene genes id pseudo

genes.pseudo<- ts %>%
  dplyr::filter(tx_biotype %in% ts.pseudo) %>%
  pull(gene_id) %>% unique()

######################Functions###################
########Functions definition######################

#####################


##Enrichment
##get enrichment results for GO, Reactome, WikiPathways and BioCarta
#####################

gse_analysis<-function(genesid){
  
  ##clean genes with na or empty values  
  genesid<-genesid[ genesid != "" ]##remove empty elements from a vector
  genesid <- genesid[!is.na(genesid)]
  
  ##conduct enrichment
  ego2 <- enrichGO(gene         = genesid,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable=T
  )
  
  ##Retrieve pathways in results
  go.res<- ego2 %>% as.data.frame() %>%
    dplyr::select(ID,Description,GeneRatio,BgRatio,pvalue, qvalue,Count) %>%
    mutate_at(c("pvalue","qvalue"),formatC,format="e",digits=3)
  
  
  ##FOr rest of enrichment pathways, we need to convert ENSEMBL to ENTREZ ID
  
  genes_entrez <- mapIds(org.Hs.eg.db,
                         keys=genesid, #Column containing Ensembl gene ids
                         column="ENTREZID",
                         keytype="ENSEMBL",
                         multiVals="first")
  ##clean entrez same way as we did for ENSEMBL
  genes_entrez<-genes_entrez[ genes_entrez != "" ]##remove empty elements from a vector
  genes_entrez <- genes_entrez[!is.na(genes_entrez)]
  
  ##Reactome Enrichment
  Reactome_enrich <- enrichPathway(gene=genes_entrez, pvalueCutoff = 0.05, readable=TRUE)
  ##subset for significant
  Reactome_enrich<-Reactome_enrich@result %>% dplyr::filter(pvalue <0.05) %>%
    mutate_at(c("pvalue","p.adjust","qvalue"),formatC,format="e",digits=2)
  
  ##WikiPathways

  WikiPathways_enrich <- enrichWP(genes_entrez, organism = "Homo sapiens",pvalueCutoff=0.05)
  ##subset for significant
    WikiPathways_enrich <- WikiPathways_enrich@result %>% dplyr::filter(pvalue <0.05)%>%
    mutate_at(c("pvalue","p.adjust","qvalue"),formatC,format="e",digits=2)
  
  ##BioCarta
  ##Prepare category where to enrich and gene_sets
  Biocarta<-msigdbr(species = "human", category = "C2",subcategory = "CP:BIOCARTA")
  gene_sets <- msigdbr(category = "C2", species = "human")
  
  ##Get Gene symbols now for this enrichment
  
  genes_symbol <- mapIds(org.Hs.eg.db,
                         keys=genesid, #Column containing Ensembl gene ids
                         column="SYMBOL",
                         keytype="ENSEMBL",
                         multiVals="first")
  ##DO enrichment
  resBioCarta <- enricher(gene = genes_symbol, 
                          TERM2GENE = gene_sets[, c("gs_name", "gene_symbol")],
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 1)
  ##Subset for significant
  resBioCarta<-as.data.frame(resBioCarta) %>% dplyr::filter(pvalue <0.05)
  
  ##Return all enrichment results
  return(list(GO=go.res,
              Reactome=Reactome_enrich,
              WikiPathways=WikiPathways_enrich,
              Biocarta=resBioCarta))
}   


##get enrichment results for GSEA using all DE results
#####################


##GSEA
gsea_analysis<-function(my.gsea){ 
  ##get all genes, UP and DOWN from all results
  
  ##remove ENCODE 
  my.gsea<-my.gsea %>%
    
    mutate(ENSEMBL=gsub("\\..*","",gene.id))
  
  #columns(org.Hs.eg.db) # returns list of available keytypes
  my.gsea$entrez <- mapIds(org.Hs.eg.db,
                           keys=my.gsea$ENSEMBL, #Column containing Ensembl gene ids
                           column="ENTREZID",
                           keytype="ENSEMBL",
                           multiVals="first")
  
  ##build a gene List
  
  ## feature 1: numeric vector
  geneList<-as.numeric(my.gsea$log2FoldChange)
  
  ## feature 2: named vector
  names(geneList) <-my.gsea$entrez
  
  ##remove NA
  geneList <- geneList[!is.na(names(geneList))]
  
  
  ## feature 3: decreasing order
  geneList<-sort(geneList, decreasing = TRUE)
  
  ego3 <- gseGO(geneList     = geneList,
                OrgDb        = org.Hs.eg.db,
                ont          = "BP",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                verbose      = FALSE)
  
  results<-ego3 %>% as_tibble %>% as.data.frame()
  
  ##Positive
  gsea_up<-results %>% dplyr::filter(enrichmentScore >0) %>% 
    dplyr::select(ID, Description, enrichmentScore,NES,pvalue, qvalues, rank,setSize)%>%
    #dplyr::slice(1:15) %>%
    mutate_at(c("enrichmentScore","NES"),round,digits=2)%>%
    mutate_at(c("pvalue","qvalues"),formatC,format="e",digits=2)
  
  ##Negative
  
  gsea_down<-results %>% dplyr::filter(enrichmentScore <0) %>% 
    dplyr::select(ID, Description, enrichmentScore,NES,pvalue, qvalues, rank,setSize)%>%
    #dplyr::slice(1:15) %>%
    mutate_at(c("enrichmentScore","NES"),round,digits=2)%>%
    mutate_at(c("pvalue","qvalues"),formatC,format="e",digits=2)
  return(list(GSEA_UP=gsea_up,
              GSEA_DOWN=gsea_down
  ))
}


##DE
##include nested model taking into account subject differences within each condition
##Tw options:"Covariate" means adding Treatment as covariate, and "Treat" means that we subset only treated sample (remove PBO ones)
#####################
de_analysis<-function(phenotype,counts,sampletable,test,annot,genes.pseudo,contrast){
  
  ####build dds according with analysis
  
  if(test=="Covariate"){
    
    m1<-model.matrix( formula(paste0(paste("~", phenotype),paste0("+Treatment"),
                                     paste("+", phenotype),paste0(":SUBJID"),
                                     paste("+", phenotype),paste0(":time"))), sampletable)
    
    ##Full Rank design matrix check
    ##Avoid levels without samples
    
    all.zero <- apply(m1, 2, function(x) all(x==0))
    
    ##Remove columns with zero
    
    idx <- which(all.zero)
    m1 <- m1[,-idx]
    m1<-m1[,!(colnames(m1) %in% nonEstimable(m1))]
    
    #check full rank
    
    if(is.fullrank(m1)){
      
      dds <- DESeqDataSetFromMatrix(counts, 
                                    sampletable,design = m1)
      
      
      
      dds<- estimateSizeFactors(dds)
      #### Requires genes to have 2 or more samples with counts of 5 or more.
      
      keep <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= min(table(sampletable[,phenotype])) 
      
      dds <- dds[keep,] 
      
      ddsTC <- DESeq(dds, betaPrior = FALSE)
      ##In contrast, first term gets 1 while the second gets -1, hence, Fold change pos is for R
      res<-results(ddsTC,contrast = list(paste0(phenotype,"R.timeWeek8"),paste0(phenotype,"NR.timeWeek8")))
      
      
      res<-as.data.frame(res)
      res<- res %>% mutate(SIG=case_when(pvalue<0.05 ~ "YES",pvalue>0.05 ~ "NO"))
      res.s<-res %>% dplyr::filter(SIG=="YES")
      
      ##annotate 
      
      
      res.s<-res.s %>% mutate(gene.id=rownames(res.s))
      res.s <-res.s %>%
        left_join(annot %>% dplyr::rename(gene.id=ID),by="gene.id")
      
      ##remove pseudogenes
      res.s<- res.s %>%
        dplyr::filter(!(gsub("\\..*","",gene.id) %in% genes.pseudo))
      res.s<-res.s %>% arrange(desc(log2FoldChange)) %>%
        dplyr::filter(!grepl("^RP|^CTD|^CTB|^AC|^ENSG|^CTC|^KB|^XX|^y|^AL|^AC|^RN|^CTA|^Y_|^AF|^AP",gene)) %>%
        mutate_if(is.numeric,round,digits=5)
      
      sorted.up<-res.s %>% dplyr::filter(log2FoldChange >0)
      sorted.down<-res.s %>% dplyr::filter(log2FoldChange <0)
      ##do the same with all results to conduct GSEA later
      
      res<-res %>% mutate(gene.id=rownames(res))
      res<-res %>%
        left_join(annot %>% dplyr::rename(gene.id=ID),by="gene.id") 
     
      Up_enrich<-gse_analysis(genesid = gsub("\\..*","",sorted.up$gene.id))
      Down_enrich<-gse_analysis(genesid = gsub("\\..*","",sorted.down$gene.id))
      
      
      ##do GSEA 
      
      gsea.enrich<-gsea_analysis(res)
      
      return(list(UP=sorted.up,
                  DOWN=sorted.down,
                  ALL=res,
                  GO_UP=Up_enrich$GO,
                  React_Up=Up_enrich$Reactome,
                  Wiki_Up=Up_enrich$WikiPathways,
                  Biocarta_up=Up_enrich$Biocarta,
                  GO_Down=Down_enrich$GO,
                  React_Down=Down_enrich$Reactome,
                  Wiki_Down=Down_enrich$WikiPathways,
                  Biocarta_Down=Down_enrich$Biocarta,
                  GSEA_UP=gsea.enrich$GSEA_UP,
                  GSEA_DOWN=gsea.enrich$GSEA_DOWN
      ))
      
      
    }else{message("Matrix is not full ranked")}
    
  }else{
    
    ##filter the sampletable and counts matrix to remove PBO
    
    sampletable<-sampletable %>% dplyr::filter(Treatment=="MH002")
    counts<-counts[,colnames(counts) %in% sampletable$sample_name ]
    
    ###get the model matrix 
    
    
    m1<-model.matrix( formula(paste0(paste("~", phenotype),
                                     paste("+", phenotype),paste0(":SUBJID"),
                                     paste("+", phenotype),paste0(":time"))), sampletable)
    ##Full rank design matrix check
    ##Avoid levels without samples
    
    all.zero <- apply(m1, 2, function(x) all(x==0))
    
    ##Remove columns with zero
    
    idx <- which(all.zero)
    m1 <- m1[,-idx]
    m1<-m1[,!(colnames(m1) %in% nonEstimable(m1))]
    
    #check full rank, then do the analysis
    
    if(is.fullrank(m1)){
      
      dds <- DESeqDataSetFromMatrix(counts, 
                                    sampletable,design = m1)
      
      
      
      dds<- estimateSizeFactors(dds)
      #### get genes with norm counts greater than or equal to 5
      
      keep <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= min(table(sampletable[,phenotype])) 
      
      dds <- dds[keep,] 
      
      ddsTC <- DESeq(dds, betaPrior = FALSE)
      res<-results(ddsTC,contrast = list(paste0(phenotype,"R.timeWeek8"),paste0(phenotype,"NR.timeWeek8")))
      
      
      res<-as.data.frame(res)
      res<- res %>% mutate(SIG=case_when(pvalue<0.05 ~ "YES",pvalue>0.05 ~ "NO"))
      res.s<-res %>% dplyr::filter(SIG=="YES")
      
      ##annotate 
      
      
      res.s<-res.s %>% mutate(gene.id=rownames(res.s))
      res.s <-res.s %>%
        left_join(annot %>% dplyr::rename(gene.id=ID),by="gene.id")
      ##remove pseudogenes
      res.s<- res.s %>%
        dplyr::filter(!(gsub("\\..*","",gene.id) %in% genes.pseudo))
      res.s<-res.s %>% arrange(desc(log2FoldChange)) %>%
        dplyr::filter(!grepl("^RP|^CTD|^CTB|^AC|^ENSG|^CTC|^KB|^XX|^y|^AL|^AC|^RN|^CTA|^Y_|^AF|^AP",gene)) %>%
        mutate_if(is.numeric,round,digits=5)
      
      sorted.up<-res.s %>% dplyr::filter(log2FoldChange >0)
      sorted.down<-res.s %>% dplyr::filter(log2FoldChange <0)
      
      Up_enrich<-gse_analysis(genesid =gsub("\\..*","",sorted.up$gene.id))
      Down_enrich<-gse_analysis(genesid =gsub("\\..*","",sorted.down$gene.id))
      
      ##do the same with all results to conduct GSEA later
      
      res<-res %>% mutate(gene.id=rownames(res))
      res<-res %>%
        left_join(annot %>% dplyr::rename(gene.id=ID),by="gene.id") 
      ##do GSEA 
      
      gsea.enrich<-gsea_analysis(res)
      
    
      ###Returns Up, and Down genes, Enrichment for each one (GO, Reactome, Wikipathways and BioCarta)
      ###and GSEA for UP and DOWN
      
      return(list(UP=sorted.up,
                  DOWN=sorted.down,
                  ALL=res,
                  GO_UP=Up_enrich$GO,
                  React_Up=Up_enrich$Reactome,
                  Wiki_Up=Up_enrich$WikiPathways,
                  Biocarta_up=Up_enrich$Biocarta,
                  GO_Down=Down_enrich$GO,
                  React_Down=Down_enrich$Reactome,
                  Wiki_Down=Down_enrich$WikiPathways,
                  Biocarta_Down=Down_enrich$Biocarta,
                  GSEA_UP=gsea.enrich$GSEA_UP,
                  GSEA_DOWN=gsea.enrich$GSEA_DOWN
                  ))
    }else{message("Matrix is not full ranked")}
  }
  
  
  
  
}

##Venn function to get which genes are individually in each model and which ones overlap


venn_function<-function(my.res.model1,my.res.model2,phenotype){
  
  
  overlap_genes<- my.res.model1[my.res.model1 %in% my.res.model2]
  
  
  Table_all_Genes<- data.frame(Gene=unique(c(my.res.model1,my.res.model2)))
  
  ##Specify with ones are in Both models and which ones don't overlap
  Table_all_Genes <-Table_all_Genes %>% 
    mutate(Model=case_when(Gene %in% overlap_genes ~ "Both",
                           Gene %in% my.res.model2~ "M2",
                           Gene %in% my.res.model1~ "M1",))
  
  ###Venn Diagram Res
  
  
  venn.res<-list()
  
  venn.res<-list(
    model_1=my.res.model1,
    model_2=my.res.model2
    
  )
  
  library(ggVennDiagram)
  
  
  venn.plot <- ggVennDiagram(venn.res,
                             category.names = c("Model 1","Model 2"),
                             edge_lty = "dashed", edge_size = 1) + scale_fill_distiller(palette = "RdBu")+
    ggtitle(enquo(phenotype))
  
  
  
  return(list(venn.plot,Table_all_Genes))
}



#Outputs Functions  
###Need to modify a couple of functions in write_xlsx function to cope with empty results in Enrichment if any
##deal with empty data frames if any 

.write_block <- function(wb, sheet, y, rowIndex=seq_len(nrow(y)),
                         colIndex=seq_len(ncol(y)), showNA=TRUE)
{
  rows  <- createRow(sheet, rowIndex)      # create rows
  cells <- createCell(rows, colIndex)      # create cells
  
  for (ic in seq_len(ncol(y)))
    mapply(setCellValue, cells[seq_len(nrow(cells)), colIndex[ic]], y[,ic], FALSE, showNA)
  
  # Date and POSIXct classes need to be formatted
  indDT <- which(sapply(y, function(x) inherits(x, "Date")))
  if (length(indDT) > 0) {
    dateFormat <- CellStyle(wb) + DataFormat(getOption("xlsx.date.format"))
    for (ic in indDT){
      lapply(cells[seq_len(nrow(cells)),colIndex[ic]], setCellStyle, dateFormat)
    }
  }
  
  indDT <- which(sapply(y, function(x) inherits(x, "POSIXct")))
  if (length(indDT) > 0) {
    datetimeFormat <- CellStyle(wb) + DataFormat(getOption("xlsx.datetime.format"))
    for (ic in indDT){
      lapply(cells[seq_len(nrow(cells)),colIndex[ic]], setCellStyle, datetimeFormat)
    }
  }
  
}


write.xlsx.custom <- function(x, file, sheetName="Sheet1",
                              col.names=TRUE, row.names=TRUE, append=FALSE, showNA=TRUE)
{
  if (!is.data.frame(x))
    x <- data.frame(x)    # just because the error message is too ugly
  
  iOffset <- jOffset <- 0
  if (col.names)
    iOffset <- 1
  if (row.names)
    jOffset <- 1
  
  if (append && file.exists(file)){
    wb <- loadWorkbook(file)
  } else {
    ext <- gsub(".*\\.(.*)$", "\\1", basename(file))
    wb  <- createWorkbook(type=ext)
  }  
  sheet <- createSheet(wb, sheetName)
  
  noRows <- nrow(x) + iOffset
  noCols <- ncol(x) + jOffset
  if (col.names){
    rows  <- createRow(sheet, 1)                  # create top row
    cells <- createCell(rows, colIndex=1:noCols)  # create cells
    mapply(setCellValue, cells[1,(1+jOffset):noCols], colnames(x))
  }
  if (row.names)             # add rownames to data x                   
    x <- cbind(rownames=rownames(x), x)
  
  if(nrow(x) > 0) {
    colIndex <- seq_len(ncol(x))
    rowIndex <- seq_len(nrow(x)) + iOffset
    
    .write_block(wb, sheet, x, rowIndex, colIndex, showNA)
  }
  saveWorkbook(wb, file)
  
  invisible()
}

#####Write Results function


write_results<-function(phenotype,results,model2test){
  ##subset results
  Res<-results[[phenotype]]
  
  write.xlsx.custom(Res$UP, file = paste0("/Users/AnGomez/Documents/MRM/Results_NewModel/res_",phenotype,paste0("Model_",model2test),".xlsx"),
                    sheetName = "UP_Genes", append = FALSE,row.names = F)
  # Add a second data set in a new worksheet
  write.xlsx.custom(Res$DOWN, file = paste0("/Users/AnGomez/Documents/MRM/Results_NewModel/res_",phenotype,paste0("Model_",model2test),".xlsx"),
                    sheetName = "DOWN_Genes", append = TRUE,row.names = F)
  
  write.xlsx.custom(Res$GO_UP, file = paste0("/Users/AnGomez/Documents/MRM/Results_NewModel/res_",phenotype,paste0("Model_",model2test),".xlsx"),
                    sheetName = "GO_enrich_UP", append = TRUE,row.names = F)
  
  write.xlsx.custom(Res$GO_Down, file = paste0("/Users/AnGomez/Documents/MRM/Results_NewModel/res_",phenotype,paste0("Model_",model2test),".xlsx"),
                    sheetName = "GO_enrich_DOWN", append = TRUE,row.names = F)
  
  write.xlsx.custom(Res$React_Up, file = paste0("/Users/AnGomez/Documents/MRM/Results_NewModel/res_",phenotype,paste0("Model_",model2test),".xlsx"),
                    sheetName = "Reactome_UP", append = TRUE,row.names = F)
  
  write.xlsx.custom(Res$React_Down, file = paste0("/Users/AnGomez/Documents/MRM/Results_NewModel/res_",phenotype,paste0("Model_",model2test),".xlsx"),
                    sheetName = "Reactome_DOWN", append = TRUE,row.names = F)
  
  write.xlsx.custom(Res$Biocarta_up, file = paste0("/Users/AnGomez/Documents/MRM/Results_NewModel/res_",phenotype,paste0("Model_",model2test),".xlsx"),
                    sheetName = "Biocarta_UP", append = TRUE,row.names = F)
  
  write.xlsx.custom(Res$Biocarta_Down, file = paste0("/Users/AnGomez/Documents/MRM/Results_NewModel/res_",phenotype,paste0("Model_",model2test),".xlsx"),
                    sheetName = "Biocarta_DOWN", append = TRUE,row.names = F)
  
  write.xlsx.custom(Res$Wiki_Up, file = paste0("/Users/AnGomez/Documents/MRM/Results_NewModel/res_",phenotype,paste0("Model_",model2test),".xlsx"),
                    sheetName = "Wiki_UP", append = TRUE,row.names = F)
  
  write.xlsx.custom(Res$Wiki_Down, file = paste0("/Users/AnGomez/Documents/MRM/Results_NewModel/res_",phenotype,paste0("Model_",model2test),".xlsx"),
                    sheetName = "Wiki_DOWN", append = TRUE,row.names = F)
  
  write.xlsx.custom(Res$GSEA_UP, file = paste0("/Users/AnGomez/Documents/MRM/Results_NewModel/res_",phenotype,paste0("Model_",model2test),".xlsx"),
                    sheetName = "GSEA_UP", append = TRUE,row.names = F)
  
  write.xlsx.custom(Res$GSEA_DOWN, file = paste0("/Users/AnGomez/Documents/MRM/Results_NewModel/res_",phenotype,paste0("Model_",model2test),".xlsx"),
                    sheetName = "GSEA_DOWN", append = TRUE,row.names = F)
  
  
}
##############End of Functions definition##########################

######################Analysis###################
########Run analysis on the three IQVIA response definition######################


my.phenotypes<-c("C_Response","E_Response","Remission")

##Model 1

C_response<-list()
E_response<-list()
Remiss<-list()

results_DE_model1<-list(C_response,E_response,Remiss)
names(results_DE_model1)<-c("C_Response","E_Response","Remission")


for( i in 1:length(my.phenotypes)){
  
  results_DE_model1[[i]]<- de_analysis(phenotype =my.phenotypes[[i]] , 
                                       counts = count.matrix,
                                       sampletable = sampleTable,
                                       annot = annot,
                                       genes.pseudo = genes.pseudo,
                                       test="Treat")
}
##Model 2

C_response<-list()
E_response<-list()
Remiss<-list()

results_DE_model2<-list(C_response,E_response,Remiss)
names(results_DE_model2)<-c("C_Response","E_Response","Remission")


for( i in 1:length(my.phenotypes)){
  
  results_DE_model2[[i]]<- de_analysis(phenotype =my.phenotypes[[i]] , 
                                       counts = count.matrix,
                                       sampletable = sampleTable,
                                       annot = annot,
                                       genes.pseudo = genes.pseudo,
                                       test="Covariate")
}


save(results_DE_model1,
     results_DE_model2,
     file="/your_path/Results_New_Model.RData")



##Write functions

sapply(my.phenotypes, function(x) write_results(phenotype = x,
              results = results_DE_model1,
              model2test = "1"
              ))

sapply(my.phenotypes, function(x) write_results(phenotype = x,
                                                results = results_DE_model2,
                                                model2test = "2"
))

##Tables with the comparison between the models


C_response_model1_model2<-venn_function(my.res.model1 = c(results_DE_model1$C_Response$UP$gene,results_DE_model1$C_Response$DOWN$gene),
                                        my.res.model2 = c(results_DE_model2$C_Response$UP$gene,results_DE_model2$C_Response$DOWN$gene,
                                        phenotype="C_Response")
                                        )
E_response_model1_model2<-venn_function(my.res.model1 = c(results_DE_model1$E_Response$UP$gene,results_DE_model1$E_Response$DOWN$gene),
                                        my.res.model2 = c(results_DE_model2$E_Response$UP$gene,results_DE_model2$E_Response$DOWN$gene),
                                        phenotype="E_Response"
)

Rem_response_model1_model2<-venn_function(my.res.model1 = c(results_DE_model1$Remission$UP$gene,results_DE_model1$Remission$DOWN$gene),
                                        my.res.model2 = c(results_DE_model2$Remission$UP$gene,results_DE_model2$Remission$DOWN$gene),
                                        phenotype="Remission"
)

write.xlsx.custom(C_response_model1_model2[[2]], file = paste0("/Users/AnGomez/Documents/MRM/Results_NewModel/Table_C_response_comparison_Models.xlsx"))
write.xlsx.custom(E_response_model1_model2[[2]], file = paste0("/Users/AnGomez/Documents/MRM/Results_NewModel/Table_R_response_comparison_Models.xlsx"))
write.xlsx.custom(Rem_response_model1_model2[[2]], file = paste0("/Users/AnGomez/Documents/MRM/Results_NewModel/Table_Rem_response_comparison_Models.xlsx"))

pdf(file = "/Users/AnGomez/Documents/MRM/Results_NewModel/Table_C_response_comparison_Models.pdf",  
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches
C_response_model1_model2[[1]]

dev.off()


pdf(file = "/Users/AnGomez/Documents/MRM/Results_NewModel/Table_E_response_comparison_Models.pdf",  
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches
E_response_model1_model2[[1]]

dev.off()

pdf(file = "/Users/AnGomez/Documents/MRM/Results_NewModel/Table_Rem_response_comparison_Models.pdf",  
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches
Rem_response_model1_model2[[1]]

dev.off()
