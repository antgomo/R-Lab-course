library(tidyverse)

library(gplots)
library(RColorBrewer)

library(ggsci) # more cool colours
library(heatmap.2x)## approach
library(data.table)
#####GO enrichment analysis
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
#library(GOstats)
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
###paths

base.dir <- "/Users/AnGomez/Documents/MRM"
prg.dir <- file.path(base.dir, qc.status, "PROGRAMS", "transcriptomics")

in.counts.file <- file.path(base.dir, "rsem_gcounts_matrix.txt")
sample.annot <- file.path(base.dir, "MRM_responders.csv")

##load gene annot

load(file.path(base.dir,"Annot_Final.RData"))

##load data

count.matrix <- fread(in.counts.file, header = TRUE) %>% 
  column_to_rownames("V1") %>%
 # mutate_all(as.integer) %>%
  as.matrix() 

colnames(count.matrix) <- substr(colnames(count.matrix),24,31)

clin.data <- fread(sample.annot, header = TRUE)
rownames(clin.data)<-clin.data$SAMPLE_ID
## Matrix columns and clin.data rows must be in the same order.
## Filter sample IDs in clin.data and count matrix
## Sort colnames in matrix and rownames in clin.data
clin.data <- clin.data %>% 
  dplyr::filter(rownames(.) %in% colnames(count.matrix)) %>%
  arrange(match(rownames(.), colnames(count.matrix)))
clin.data$Time<-ifelse(clin.data$Timepoint=="Baseline","Baseline","Week8")

## Filter count.matrix samples
count.matrix <- count.matrix[,colnames(count.matrix) %in% rownames(clin.data)]


## Create sampleTable object to annotate deseqdataset
sampleTable <- data.frame(C_Response = factor(ifelse(clin.data$cResponse==1,"R","NR")),
                          Remission=factor(ifelse(clin.data$cRemission==1,"R","NR")),
                          E_Response=factor(ifelse(clin.data$eResponse==1,"R","NR")),
                          Treatment=clin.data$Treatment,
                          time=factor(clin.data$Time),
                          sample_name = colnames(count.matrix))

rownames(sampleTable) <- sampleTable$sample_name

##Analysis:
##Option 01; subset to treated
##Option 02; consider treatment as covariate


##Prepare count matrix
count.matrix<-round(count.matrix)
mode(count.matrix)<-"integer"

##remove UDI_0059

count.matrix<-count.matrix[,!(colnames(count.matrix) %in% "UDI_0059")]
sampleTable<-sampleTable[!(sampleTable$sample_name %in% "UDI_0059"),]

##Relevel factors as reference
sampleTable$time<-relevel(sampleTable$time,ref="Baseline")
sampleTable$C_Response<-relevel(sampleTable$C_Response,ref="NR")
sampleTable$Remission<-relevel(sampleTable$Remission,ref="NR")
sampleTable$E_Response<-relevel(sampleTable$E_Response,ref="NR")


##Functions

##DE


de_analysis<-function(phenotype,counts,sampletable,test,annot){
  
  ####build dds according with analysis
  
  if(test=="Covariate"){
    
    dds <- DESeqDataSetFromMatrix(counts, 
                                  sampletable, 
                                  formula(paste0(paste("~", phenotype),paste0("+time+Treatment"),paste("+", phenotype),paste0(":time"))))
   
    
    dds<- estimateSizeFactors(dds)
    
    #### filter less than min samples with norm counts greater than or equal to 10
    keep <- rowSums( counts(dds, normalized=TRUE) >= 10 ) >= min(table(sampletable[,phenotype])) 
    
    dds <- dds[keep,] 

    ddsTC <- DESeq(dds, test="LRT", reduced = formula(paste0(paste("~", phenotype),paste0("+time+Treatment")))) 
    
    ##get the LogFold change with R vs NR, pvalues will be the ones from LRT test considering interaction
    res<-results(ddsTC)
    res<-lfcShrink(ddsTC, coef = paste0(phenotype,"_R_vs_NR"),res=res,type = "ashr")
    

    res<-as.data.frame(res)
    res<- res %>% mutate(SIG=case_when(pvalue<0.05 ~ "YES",pvalue>0.05 ~ "NO"))
    res.s<-res %>% dplyr::filter(SIG=="YES")
    
    ##annotate 
    
    
    res.s<-res.s %>% mutate(gene.id=rownames(res.s))
    res.s <-res.s %>%
          left_join(annot %>% dplyr::rename(gene.id=ID),by="gene.id")
    
    res<-res %>% mutate(gene.id=rownames(res))
    res<-res %>%
      left_join(annot %>% dplyr::rename(gene.id=ID),by="gene.id")   
    
    
    ##make heatmap
    
    my.data<-as.data.frame(DESeq2::counts(ddsTC,normalized=T)[rownames(counts(ddsTC)) %in% res.s$gene.id ,])
    ###subset to wk8
    wk8_samples<-sampleTable %>% dplyr::filter(time=="Week8")
    
    my.data<-my.data[,colnames(my.data) %in% wk8_samples$sample_name]


      #####scale
      
    data <- t(scale(t(my.data))) # z-score normalise each row (feature)
    data <- apply(data, 2, function(x) ifelse(x > 4, 4, x)) # compress data within range
    data <- apply(data, 2, function(x) ifelse(x < -4, -4, x)) # compress data within range
      
      
      spcol<-samples
      
    
    return(list(Sig=res.s,ALL=res,data2plot=data,cols2plot=spcol))
    
  }else{
    
    ##filter the sampletable and counts matrix to remove PBO
    
    sampletable<-sampletable %>% dplyr::filter(Treatment=="MH002")
    counts<-counts[,colnames(counts) %in% sampletable$sample_name ]
    
    dds <- DESeqDataSetFromMatrix(counts, 
                                  sampletable, 
                                  formula(paste0(paste("~", phenotype),paste0("+time"),paste("+", phenotype),paste0(":time"))))
    


    dds<- estimateSizeFactors(dds)
    #### filter less than min samples with norm counts greater than or equal to 10
    
    keep <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= min(table(sampletable[,phenotype])) 
    
    dds <- dds[keep,] 

    ddsTC <- DESeq(dds, test="LRT", reduced = formula(paste0(paste("~", phenotype),paste0("+time")) )) 
    
    ##get the LogFold change with R vs NR, pvalues will be the ones from LRT test considering interaction
    

    res<-lfcShrink(ddsTC, coef = paste0(phenotype,"_R_vs_NR"), type = "ashr")
    
    
    res<-as.data.frame(res)
    res<- res %>% mutate(SIG=case_when(pvalue<0.05 ~ "YES",pvalue>0.05 ~ "NO"))
    res.s<-res %>% dplyr::filter(SIG=="YES")
    
    ##annotate 
    
    
    res.s<-res.s %>% mutate(gene.id=rownames(res.s))
    res.s <-res.s %>%
      left_join(annot %>% rename(gene.id=ID),by="gene.id")
    
    res<-res %>% mutate(gene.id=rownames(res))
    res<-res %>%
      left_join(annot %>% rename(gene.id=ID),by="gene.id") 
    ##make heatmap
    
    my.data<-as.data.frame(DESeq2::counts(ddsTC,normalized=T)[rownames(counts(ddsTC)) %in% res.s$gene.id ,])
    ###subset to wk8
    wk8_samples<-sampleTable %>% dplyr::filter(time=="Week8")
    
    my.data<-my.data[,colnames(my.data) %in% wk8_samples$sample_name]
    
    
    #####scale
    
    data <- t(scale(t(my.data))) # z-score normalise each row (feature)
    data <- apply(data, 2, function(x) ifelse(x > 4, 4, x)) # compress data within range
    data <- apply(data, 2, function(x) ifelse(x < -4, -4, x)) # compress data within range
    
    
    spcol<-samples
    
    
    return(list(Sig=res.s,ALL=res,data2plot=data,cols2plot=spcol))
    
  }
  
  
  
  
}


##Heatmap

plot_heatmap<-function(phenotype,data2plot,cols,spcol){
  
  heatmap.2x(as.matrix(data2plot), col=rev(cols), Colv = T,Rowv=TRUE, scale="none", 
             ColSideColors=spcol,
             trace="none", dendrogram="column", 
             cexRow=1, cexCol=1,
             main=paste0("R vs Non NR week8 ", phenotype),
             #labCol=NA,
             labRow=NA, 
             density.info="none",
             hclust=function(x) hclust(x,method="complete"),
             distfun=function(x) as.dist((1-cor(t(x)))/2),srtCol=45)}


##RTF

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
  
  rtf <- RTF(output, width=11, height=8.5, font.size=12, omi=rep(0.75,ncol(tab)))
  addText(rtf, "\\deff {\\fonttbl {\\f5 Times;}}\\n")
  
  startParagraph(rtf)
  
  addText(rtf,bold=F,paste("\\f5\\fs22", rtf.header, "{\\pindtabqr}","", "\n", sep=" "))
  
  endParagraph(rtf)
  
  addTable(rtf, tab, row.names = F, font.size = 11, 
           col.widths=rep(colW,ncol(tab)),
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



##Intersect

intersect_function<-function(res){
  
  
  ##Datasets
  
  # Set 1. Intestinal barrier function
  intestinal_barrier <- c("CLDN1","CLDN2","CLDN3",
                          "CLDN4", "TJP1", "OCLN",
                          "MLCK", "MUC2","MUC3",
                          "MUC5AC", "TFF3")
  
  # Set 2. Epithelial restitution and Wound healing
  epiteh_rest <- c("ANXA2", "CDH1", "CDX2",
                   "GJA1", "MKI-67", "MMP2",
                   "MMP29", "NANOG", "PCNA",
                   "POU5F1", "RND3", "TLR2",
                   "TLR4", "TLR5", "TLR9",
                   "VIL1", "WISP1") 
  
  # Set 3. Immune modulation
  immune_mod <- c("L12p70", "IFNG", "IL4", 
                  "IL1B", "IL17A", "IL21",
                  "IL22", "IL9")
  
  
  
  
  print("Intestinal Barrier")
  
  res %>% dplyr::select(gene,log2FoldChange,pvalue,padj) %>% dplyr::filter(gene %in% intestinal_barrier)
  
  print("Epithelial restitution")
  
  res %>% dplyr::select(gene,log2FoldChange,pvalue,padj) %>% dplyr::filter(gene %in% epiteh_rest)
  
  print("Immune modulation")
  
  res %>% dplyr::select(gene,log2FoldChange,pvalue,padj) %>% dplyr::filter(gene %in% immune_mod)
}

##plot to TFL
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

#########################
##run analysis
###
###
###
#########################
###covariates


analysis_C_response_treatment<-de_analysis(phenotype = "C_Response",
                                 counts = count.matrix,
                                 sampletable = sampleTable,
                                 annot = annot,
                                 test="Covariate"
                                 )

analysis_E_response_treatment<-de_analysis(phenotype = "E_Response",
                                           counts = count.matrix,
                                           sampletable = sampleTable,
                                           annot = annot,
                                           test="Covariate"
)


analysis_Remission_treatment<-de_analysis(phenotype = "Remission",
                                           counts = count.matrix,
                                           sampletable = sampleTable,
                                           annot = annot,
                                           test="Covariate"
)


###Only treatment,not PBO in the model
analysis_C_response<-de_analysis(phenotype = "C_Response",
                                           counts = count.matrix,
                                           sampletable = sampleTable,
                                           annot = annot,
                                           test="Treat"
)

analysis_E_response<-de_analysis(phenotype = "E_Response",
                                           counts = count.matrix,
                                           sampletable = sampleTable,
                                           annot = annot,
                                           test="Treat"
)


analysis_Remission<-de_analysis(phenotype = "Remission",
                                          counts = count.matrix,
                                          sampletable = sampleTable,
                                          annot = annot,
                                          test="Treat"
)

save(analysis_C_response,analysis_C_response_treatment,analysis_E_response,analysis_E_response_treatment, analysis_Remission,analysis_Remission_treatment,
     
     file="/Users/AnGomez/Documents/MRM/Results.RData")
####Filter results
####clean of pseudogenes

ts <- transcripts(EnsDb.Hsapiens.v86)
ts.pseudo<-grep("pseudo|antisense|lincRNA", unique(ts$tx_biotype), val=TRUE)
ts<-as.data.frame(ts)


##gene genes id pseudo

genes.pseudo<- ts %>%
  dplyr::filter(tx_biotype %in% ts.pseudo) %>%
  pull(gene_id) %>% unique()


###Use results

phenotype<-"Remission"

my.res.model1<- analysis_Remission$Sig %>%
  dplyr::filter(!(gsub("\\..*","",gene.id) %in% genes.pseudo))
my.res.model2<- analysis_Remission_treatment$Sig %>%
  dplyr::filter(!(gsub("\\..*","",gene.id) %in% genes.pseudo))


my.res<- my.res %>%
  dplyr::filter(!(gsub("\\..*","",gene.id) %in% genes.pseudo))


my.res.model1<-my.res
my.res.model2<-my.res

###Venn Diagram Res

venn.res<-list()

venn.res<-list(model_1=my.res.model1$gene,
               model_2=my.res.model2$gene
)

library(ggvenn)
ggvenn(
  venn.res, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)



#######################################################
##               Outputs               ##
#######################################################

##UP

sorted.up<-my.res %>% arrange(desc(log2FoldChange)) %>%
  dplyr::filter(!grepl("^RP|^CTD|^CTB|^AC|^ENSG|^CTC|^KB|^XX|^y|^AL|^RN|^CTA|^Y_",gene)) %>%
  mutate_if(is.numeric,round,digits=5)

to_rtf_up<-sorted.up  %>% 
  dplyr::select(gene.id,gene,log2FoldChange,pvalue,padj) %>%
  dplyr::slice(1:20)

##DOWN
sorted.down<-my.res %>% arrange(log2FoldChange) %>%
  dplyr::filter(!grepl("^RP|^CTD|^CTB|^AC|^ENSG|^CTC|^KB|^XX|^y|^AL|^RN|^CTA|^Y_",gene)) %>%
  mutate_if(is.numeric,round,digits=5)

to_rtf_down<-sorted.down %>% 
  dplyr::select(gene.id,gene,log2FoldChange,pvalue,padj) %>%
  dplyr::slice(1:20)

## RTF output
title <- paste0("5.2.6 Secondary Objective: Analysis of Biomarkers associated with response to treatment at week 8 .",
                "\nITT population (N = ",length(sampleTable$sample_name),")",  
                "\nBiomarker Statistical Analysis Plan for MH002 in subjects with mild to moderate Ulcerative Colitis")
foot1 <- paste0("\nNote: This table shows the sample information of the treatment used and the response state assigned in the ITT population.", 
                "\nAbbreviations: R = Responder, NR = Non_Responder")
foot2 <- paste0("\nProgram Path: ", file.path(base.dir))
foot3 <- paste("Status: Under review", sep=" ")
foot4 <- paste("Date: ",Sys.Date(), sep=" ")



table.rtf(tab = to_rtf_up,
          rtf.header = title,
          rtf.foot.1 = foot1,
          rtf.foot.2 = foot2,
          rtf.foot.3 = foot3,
          rtf.foot.4 = foot4,
          colW = 1.75,
          outDir = getwd(),
          prefix = paste0("UP_Genes_",phenotype,".rtf"))


####GO enrichment


UP<- my.res %>%
  dplyr::filter(log2FoldChange > 0) %>%
  pull(gene.id) 
UP<-  gsub("\\..*","",UP)

DOWN<- my.res %>%
  dplyr::filter(log2FoldChange < 0) %>%
  pull(gene.id)
DOWN<-  gsub("\\..*","",DOWN)

genesid<-UP
genesid<-DOWN


genesid<-genesid[ genesid != "" ]##remove empty elements from a vector
genesid <- genesid[!is.na(genesid)]

ego2 <- enrichGO(gene         = genesid,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05,
                 readable=T
)


go.res.up<- ego2 %>% as.data.frame() %>%
            dplyr::select(ID,Description,GeneRatio,BgRatio,pvalue, qvalue,Count) %>%
  mutate_at(c("pvalue","qvalue"),formatC,format="e",digits=2)

go.res.down<- ego2 %>% as.data.frame() %>%
  dplyr::select(ID,Description,GeneRatio,BgRatio,pvalue, qvalue,Count)%>%
  mutate_if(is.numeric,round,digits=3)
mutate_at(c("pvalue","qvalue"),formatC,format="e",digits=2)


my.p<-dotplot(ego2, title=paste0("Responders vs Non Responders DOWN R wk8 Model 2 ",phenotype),showCategory=20)




title <- paste0("5.2.6 Secondary Objective: Analysis of Biomarkers associated with response to treatment at week 8 .",
                "\nGO enrichment Responders DOWN Genes",
                "\nITT population (N = ",length(sampleTable$sample_name),")",  
                "\nBiomarker Statistical Analysis Plan for MH002 in subjects with mild to moderate Ulcerative Colitis")
foot1 <- paste0("\nNote: This plot shows the sample information of the treatment used and the response state assigned in the ITT population.", 
                "\nAbbreviations: R = Responder, NR = Non_Responder")
foot2 <- paste0("\nProgram Path: ", file.path(base.dir))
foot3 <- paste("Status: Under review", sep=" ")
foot4 <- paste("Date: ",Sys.Date(), sep=" ")

plot.to.tfl(p.ggplot = my.p,
          tfl.header = title,
          tfl.foot.1 = foot1,
          tfl.foot.2 = foot2,
          tfl.foot.3 = foot3,
          tfl.foot.4 = foot4,
          outDir = "/Users/AnGomez/Documents/MRM/",
          prefix = paste0("DOWN_Genes_",phenotype))




###GSEA
my.gsea<-my.res
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


## feature 3: decreasing orde
geneList<-sort(geneList, decreasing = TRUE)




ego3 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              pvalueCutoff = 0.1,
              pAdjustMethod = "none",
              verbose      = FALSE)

results<-ego3 %>% as_tibble %>% as.data.frame()

id<-351

my.p2<-gseaplot2(ego3, geneSetID = id, title = ego3$Description[id])

grep("TOR",ego3$Description)

results %>% dplyr::filter(grepl("cell adhesion", Description))

results %>% dplyr::filter(enrichmentScore >0) %>% pull(Description)


##Positive
gsea_up<-results %>% dplyr::filter(enrichmentScore >0) %>% 
  dplyr::select(ID, Description, enrichmentScore,NES,pvalue, qvalues, rank,setSize)%>%
 #dplyr::slice(1:15) %>%
  mutate_at(c("enrichmentScore","NES"),round,digits=2)%>%
  mutate_at(c("pvalue","qvalues"),formatC,format="e",digits=2)




##Negative

gsea_down<-results %>% dplyr::filter(enrichmentScore <0) %>% 
  dplyr::select(ID, Description, enrichmentScore,NES,pvalue, qvalues, rank,setSize)%>%
 dplyr::slice(1:15) %>%
  mutate_at(c("enrichmentScore","NES"),round,digits=2)%>%
  mutate_at(c("pvalue","qvalues"),formatC,format="e",digits=2)





#####Results 

model2test<-2

write.xlsx(sorted.up, file = paste0("/Users/AnGomez/Documents/MRM/Results_",phenotype,paste0("Model_",model2test),".xlsx"),
           sheetName = "UP_Genes", append = FALSE,row.names = F)
# Add a second data set in a new worksheet
write.xlsx(sorted.down, file = paste0("/Users/AnGomez/Documents/MRM/Results_",phenotype,paste0("Model_",model2test),".xlsx"),
           sheetName = "DOWN_Genes", append = TRUE,row.names = F)

write.xlsx(go.res.up, file = paste0("/Users/AnGomez/Documents/MRM/Results_",phenotype,paste0("Model_",model2test),".xlsx"),
           sheetName = "GO_enrich_UP", append = TRUE,row.names = F)

write.xlsx(go.res.down, file = paste0("/Users/AnGomez/Documents/MRM/Results_",phenotype,paste0("Model_",model2test),".xlsx"),
           sheetName = "GO_enrich_DOWN", append = TRUE,row.names = F)

write.xlsx(gsea_up, file = paste0("/Users/AnGomez/Documents/MRM/Results_",phenotype,paste0("Model_",model2test),".xlsx"),
           sheetName = "GSEA_UP", append = TRUE,row.names = F)

write.xlsx(gsea_down, file = paste0("/Users/AnGomez/Documents/MRM/Results_",phenotype,paste0("Model_",model2test),".xlsx"),
           sheetName = "GSEA_DOWN", append = TRUE,row.names = F)














##################
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)

plot_heatmap(data2plot =analysis_C_response$data2plot,
             phenotype="C_response",
             cols=cols,
             spcol = analysis_C_response$cols2plot
             
             
)

plot_heatmap(data2plot =analysis_E_response$data2plot,
             phenotype="E_response",
             cols=cols,
             spcol = analysis_E_response$cols2plot
             
             
)

plot_heatmap(data2plot =analysis_Remission_treatment$data2plot,
             phenotype="Remission",
             cols=cols,
             spcol = analysis_Remission_treatment$cols2plot
             
             
)










