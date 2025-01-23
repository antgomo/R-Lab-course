---
  ###########################################################################
# Project           : MRM Health - Statistical support for MH002-UC-201 - PLX-23-35515 (1545)
# Program name      : 999_Transcriptomic_Analysis_MH002CR_5ASA.R
# Developed in      : R 4.2.0
# Purpose           : Differential expression analysis  (5ASA vs MH002-CR)
##                   Complementary analysis considering nested model controlloing by SUBJID following Love et al
##                  http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#group-specific-condition-effects-individuals-nested-within-groups
# Inputs            : RSEM Expression values located at Transcriptomics data
#   
# Outputs           : Excel files with Differentially Expressed Genes (DEG) and their enrichment analysis
# Revision History  :
#   Version   Date        Author       Revision 
#   -------   ---------   ----------   ---------------------------------
#   1.0       14092023    Antonio Gomez  Creation
#   2.0       02112023    Antonio Gomez  Modification
# Reviewed by       : Pablo Villegas 02112023
# Reference number  : NA                                      
# Validation Level  : NA
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
library(RColorBrewer)
library(heatmap.2x)
library(glue)
library(GSVA)
library(ggrepel)

###paths

base.dir <- "/Users/AnGomez/Documents/MRM"

in.counts.file <- file.path(base.dir, "rsem_gcounts_matrix.txt")
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
                          Treatment=as.factor(clin.data$Treatment),
                          time=factor(clin.data$Time),
                          SUBJID=factor(clin.data$SUBJ_ID),
                          sample_name = colnames(count.matrix))

rownames(sampleTable) <- sampleTable$sample_name



##Prepare count matrix. RSEM values must be rounded to be used
##Approach form tximport package from Soneson et al (2015), using tximport function part dealing with RSEM
count.matrix<-round(count.matrix)
mode(count.matrix)<-"integer"

##remove UDI_0059 (because QC issues) and its neighbor in wk8, UDI_0048

count.matrix<-count.matrix[,!(colnames(count.matrix) %in% c("UDI_0059","UDI_0048"))]
sampleTable<-sampleTable[!(sampleTable$sample_name %in% c("UDI_0059","UDI_0048")),]

##Relevel factors as reference
sampleTable$time<-relevel(sampleTable$time,ref="Baseline")
sampleTable$Treatment<-relevel(sampleTable$Treatment,ref="PBO")



##Check genes unexpressed

zero_unexpressed <- (apply(count.matrix, 1, max) == 0)

##remove them from the matrix

count.matrix<-count.matrix[!zero_unexpressed,]


##use information from the clinet to remove duplicates

sample_client<-readRDS(file.path(base.dir,"LAST/UC201.rnaseq.rds"))

metadata <- sample_client$meta$samples%>%dplyr::filter(`Barcode` != "UDI_0059")

###### exclude samples with high number of PCR duplicates ######
metadata = metadata%>%dplyr::filter(`Pct Dup` <= allowed_duplicates)
allowed_duplicates = 100

metadata<-metadata %>% 
  dplyr::filter(`Pct Dup` <= allowed_duplicates)

metadata<-metadata %>% dplyr::rename(SUBJID='Subject ID')

metadata<-metadata %>% dplyr::rename(sample_name='Barcode')
metadata<-metadata %>% dplyr::rename(PCT_DUP='Pct Dup')



#CR:5ASA contrast is the one of interest. 
#considering all subjects that *only*
# have cResponse (so not eResponse nor cRemission).

##merge sampleTable with client metadata to add Group

sampleTable<-sampleTable %>%
              left_join(metadata %>% dplyr::select(c("sample_name","Group","PCT_DUP")))

####Filter results
####clean of pseudogenes. We are not interested in pseudogenes, antisense and/or lincRNA from now
###Build a list of them from ENCODE data an use later to filter the resuls

ts <- transcripts(EnsDb.Hsapiens.v86)
ts.pseudo<-grep("pseudo|antisense|lincRNA", unique(ts$tx_biotype), val=TRUE)
ts<-as.data.frame(ts)


##gene genes id pseudo


genes.pseudo<- ts %>%
  dplyr::filter(tx_biotype %in% ts.pseudo) %>%
  pull(gene_id) %>% unique()

######################
##gene sets to overlap
#######################

##get pathways and genes object
pathw_genes<- as.list(org.Hs.egGO2ALLEGS)
# Get the pathways we want to check


my.paths<-data.frame("GOID"=c("GO:0090557","GO:0042060","GO:0070160","GO:0098594","GO:0008283","GO:0030154",
                              "GO:0002526","GO:0006955"),
                     "Desc"=c("Barrier_integrity","Wound_Healing","Tight_junction","Mucin","Cell_proliferation",
                              "Cell_differentiation","Inflammation_Response","Immunity_response"
                              )
                     
                     )
gene_paths<-list()

gene_paths<-sapply(my.paths$GOID, function(x) return(as.character(unlist(pathw_genes[x])))
  
  )

Genes2paths<-list()

Genes2paths <- lapply(gene_paths, function(x)
  
                    return(as.character(mapIds(org.Hs.eg.db,
                       keys=x, #Column containing Ensembl gene ids
                       column="ALIAS",
                       keytype="ENTREZID",
                       multiVals="first"))))

##change names to 

names(Genes2paths)<-my.paths$Desc

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
  
  
  ## feature 3: decreasing orde
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
#####################

de_analysis<-function(counts,sampletable,test,annot,genes.pseudo,contrast){
  
  ####build dds according with analysis
    
  ##droplevels
  
  sampletable$Group<-droplevels(sampletable$Group)
  sampletable$SUBJID<-droplevels(sampletable$SUBJID)

  ##order sample table according count matrix
  
  sampletable<-sampletable %>% arrange(factor(sample_name,levels=colnames(counts)))
  
    
    m1<-model.matrix( formula(paste0(paste("~", "Group"),paste0("+SUBJID"))), sampletable)

    ##Avoid levels without samples
    
    all.zero <- apply(m1, 2, function(x) all(x==0))
    
    ##Remove columns with zero
    
    m1<-m1[,!(colnames(m1) %in% nonEstimable(m1))]
    
    #check full rank
    
    if(is.fullrank(m1)){
      
      dds <- DESeqDataSetFromMatrix(counts, 
                                    sampletable,design = m1)
      
      
      
      dds<- estimateSizeFactors(dds)
      #### Requires genes to have 2 or more samples with counts of 5 or more.
      
      keep <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= min(table(sampletable[,"Treatment"])) 
      
      dds <- dds[keep,] 
      
      ddsTC <- DESeq(dds, betaPrior = FALSE)
      ##In contrast, first term gets 1 while the second gets -1, hence, Fold change pos is for R
      res<-results(ddsTC,contrast = list("GroupMH002.CR"))
      
      
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
      
      ##heatmap
      my.data<-as.data.frame(DESeq2::counts(ddsTC,normalized=T)[rownames(counts(ddsTC)) %in% res.s$gene.id ,])
      ###subset to wk8
      wk8_samples<-sampletable %>% dplyr::filter(time=="Week8")
      my.data<-my.data[,colnames(my.data) %in% wk8_samples$sample_name]

      samples<-ifelse(wk8_samples$Group=="MH002-CR","purple","grey")

      #####scale
    
      data <- t(scale(t(my.data))) # z-score normalise each row (feature)
      data <- apply(data, 2, function(x) ifelse(x > 4, 4, x)) # compress data within range
      data <- apply(data, 2, function(x) ifelse(x < -4, -4, x)) # compress data within range
      ##colors
      spcol<-samples
      cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
      
      
      hm2<-heatmap.2x(as.matrix(data), col=rev(cols), Colv = T,Rowv=TRUE, scale="none", ColSideColors=spcol,
                        trace="none", dendrogram="column", 
                        cexRow=1, cexCol=1,
                        main="MH002-CR vs 5ASA week8",
                        #labCol=NA,
                        labRow=NA, 
                        density.info="none",
                        hclust=function(x) hclust(x,method="complete"),
                        distfun=function(x) as.dist((1-cor(t(x)))/2),srtCol=45
                        
      )
      
      
      return(list(UP=sorted.up,
                  DOWN=sorted.down,
                  heatmap=hm2,
                  cols_heatmap=spcol,
                  SIG=res.s,
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
                  GSEA_DOWN=gsea.enrich$GSEA_DOWN,
                  ITT=dim(data)[2]
                  
      ))
      
      
    }else{message("Matrix is not full ranked")}
  
  
}


####################
###ssGSEA
####################

##Single  sample  GSEA  (ssGSEA)  calculates  a  gene  set enrichment  score  per  sample  
#as  the  normalized  difference  in  empirical  cumulative  distribution functions of gene expression ranks 
#inside and outside the gene set.


##test interaction treatment and time using as random SUBJID.Log-likehood maximized

test_ssgsea<-function(counts,pathways.set,sampletable){
  
  ##conduct ssgsea for raw counts on pathways of interest as set of genes
  counts<-as.matrix(counts)
  
  sampletable<-sampletable %>% arrange(factor(sample_name,levels=colnames(counts)))
  
  ssgsea <- gsva (counts,pathways.set,method="ssgsea",kcdf="Poisson", verbose=TRUE,ssgsea.norm=T)
  
  ssgsea<-as.data.frame(t(ssgsea))
  ssgsea<-rownames_to_column(ssgsea,var="sample_name")
  
  res.ssgsea<- ssgsea %>% 
    left_join(sampleTable,by="sample_name")
  
  ####test over a set of columns, in this case, the pathways
  
  pathw<-colnames(res.ssgsea)[2:8]
  
  ##test considering SUBJID as the random variable (repeated measures)
  
  ##list to keep results
  results.path<-list()
  
  for(i in 1:length(pathw))
  { 
    
    ##build formula
    formula.test<-formula(paste0(pathw[i],paste("~", "Treatment"),paste0("*time")))
    ##conduct test
    results.path[[i]]<-anova(lme(formula.test,random = ~1|SUBJID, method = "ML", na.action = na.exclude,data=res.ssgsea))$'p-value'[3]
  }
  
  ###store in a data frame 
  
  names(results.path)<-pathw
  results.path<-as.data.frame(cbind(results.path))
  colnames(results.path)<-"p-value"
  results.path$`p-value`<-as.numeric(results.path$`p-value`)
  results.path<-rownames_to_column(results.path,var= "Pathway")
  
  ##round p-vals
  
  results.path<-results.path %>% mutate(across(where(is.numeric), round, 3))
  return(results.path)
}

      ##Heatmap

pdf_heatmap <- function(res, outDir, outname, title, foot1, foot2, foot3=Sys.Date()) {
  
  cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
  colorBlindGrey8 <- c("#9900FF", "#999999", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  
  filename <- file.path(outDir, paste0(outname, ".pdf"))
  
  # PDF
  pdf(filename, width=8.27, height=11.69)
  par(oma = c(20,0,15,10))  
  par(mar = c(5, 1, 4, 2))
  
  title <- title
  
  foot1 <- foot1
  foot2 <- foot2
  foot3 <- foot3
  
  revRowInd <- match(c(1:length(res$heatmap$rowInd)), res$heatmap$rowInd)
  revColInd <- match(c(1:length(res$heatmap$colInd)), res$heatmap$colInd)
  
  
  heatmap.2x(t(res$heatmap$carpet)[revRowInd, revColInd], 
             col=res$heatmap$col, 
             Rowv=res$heatmap$rowDendrogram, 
             Colv=res$heatmap$colDendrogram, 
             ColSideColors=res$cols_heatmap,
             trace="none", dendrogram="column", 
             cexRow = 0.75, 
             cexCol = 0.75,
             main = "", 
             cex.main = 0.5,
             margins = c(5,5),
             #labCol=NA,
             labRow=NA, 
             density.info="none",
             hclust=function(x) hclust(x,method="complete"),
             distfun=function(x) as.dist((1-cor(t(x)))/2),srtCol=45)
  
  
  legend("bottomleft",c("MH002-CR","5ASA"),pch=20:20,col=c("purple","grey"), bty = "n")
  
  #print header
  grid.text(title, x = .08, y = .93, just = c('left', 'top'), gp = gpar(fontsize = 12, col ='black'))
  
  # print footer
  grid.text(paste0(foot1, '\n', foot2, '\n', foot3), x = .08, y = .02, just = c('left', 'bottom'), gp = gpar(fontsize = 8, col = 'black'))
  
  # Close the PDF file
  dev.off()
}                          

######Volcano plot


pdf_volcano <- function(DE_res, pv.thr,fc.thr, outDir, outname, title, foot1, foot2, foot3=Sys.Date()) {
  
  ##filter DE results from pseudogenes and HAVANA coding genes
  
  DE_res<- DE_res %>%
    dplyr::filter(!(gsub("\\..*","",gene.id) %in% genes.pseudo)) %>%
    dplyr::filter(!grepl("^RP|^CTD|^CTB|^AC|^ENSG|^CTC|^KB|^XX|^y|^AL|^AC|^RN|^CTA|^Y_|^AF|^AP",gene)) 
  
  ##plot
  
  p.volcano.df <- DE_res %>%
    data.frame() %>%
    mutate(THR = ifelse(pvalue >= pv.thr | log2FoldChange > -fc.thr & log2FoldChange < fc.thr, "N", "Y")) 
  # pAdj
  p.volcano <- ggplot(p.volcano.df, mapping = aes(x = log2FoldChange, y = -log10(pvalue)), size=1.5,show.legend = F) +
    geom_point(aes(colour=THR)) +
    scale_color_manual(values = c(Y = "red", N = "grey"))+
    geom_vline(xintercept=fc.thr, linetype="dashed", color = "black") +
    geom_vline(xintercept=-fc.thr, linetype="dashed", color = "black") +
    geom_hline(yintercept=-log10(pv.thr), linetype="dashed", color = "black") +
    geom_text_repel(
      p.volcano.df %>% 
        dplyr::filter(THR == "Y"),
      mapping = aes(x = log2FoldChange, y = -log10(pvalue), label = gene),
      box.padding = 1, max.overlaps = Inf) +
    xlab("log2(FC)") +
    ylab("-log10(pval)") +
    theme_bw()
  
  ##ooutput as pdf
  
  
  # output as a pdf
  output <- file.path(outDir, paste0(outname, "_volcano.pdf"))
  
  
  scale.mar.y = 1.1
  scale.mar.x = 1
  
  
  pdf(output, height=8.5, width=12)
  print(p.volcano +
          theme(axis.text.y = element_text(size=10), 
                axis.text.x = element_text(size=11), 
                axis.title.x = element_text(size=12), 
                axis.title.y = element_text(size=11),
                plot.title = element_text(hjust = 0.5),
                legend.position = "none",
                plot.margin=unit(c(50*scale.mar.y,50*scale.mar.x,50*scale.mar.y,50*scale.mar.x),"mm")))
  # print header
  
  grid.text(title, x = .08, y = .93, just = c('left', 'top'), gp = gpar(fontsize = 12, col = 'black'))
  # print footer
  
  grid.text(paste0(foot1, '\n', foot2, '\n', foot3, '\n'), x = .08, y = 0, just = c('left', 'bottom'), gp = gpar(fontsize = 10, col = 'black'))
  dev.off()
  
}    



######################Analysis###################
##for this analysis, consider only the subjects that are CResponders and rest of 5ASA

my.subjects<-sampleTable %>% dplyr::filter(Group=="MH002-CR") %>% pull(SUBJID) %>% unique() %>% as.character()


##get sample names for these subjects

my.sample.names<- sampleTable %>% dplyr::filter(SUBJID %in% my.subjects ) %>% pull(sample_name)
  
###get all the sample in 5ASA arm
asa.arm.samples<-sampleTable %>% dplyr::filter(Group=="5ASA" ) %>% pull(sample_name)

###

samples2get<-unique(c(my.sample.names,asa.arm.samples))

results_DE_MH002CR_5ASA<- de_analysis( counts = count.matrix[,samples2get],
                                       sampletable = sampleTable %>% dplyr::filter(sample_name %in% samples2get),
                                       annot = annot,
                                       genes.pseudo = genes.pseudo)



    
save(results_DE_MH002CR_5ASA, file = "/Users/AnGomez/Documents/MRM/LAST/Mechanistic/MH002CR_5ASA/Res_matrix_Heatmap_Genes_MH002CR_5ASA.RData")



##SsGSEA results
##subset by different pathways genes

DE_results_by_ptw<-list()

DE_results_by_ptw<- lapply(Genes2paths, function(x){
  
   DEG_genes<-rbind(results_DE_MH002CR_5ASA$UP,results_DE_MH002CR_5ASA$DOWN)
  
   res<-DEG_genes %>% dplyr::filter(DEG_genes$gene %in% x)
   }
  
  )
DE_results_by_ptw<-do.call("rbind",DE_results_by_ptw) 

DE_results_by_ptw<-rownames_to_column(DE_results_by_ptw,var="Pathway")
##transform count matrix ENSEMBL in rows to ALIAS

df.counts<-as.data.frame(count.matrix[,samples2get]) %>%
           mutate(ID=rownames(count.matrix[,samples2get])) %>%
           left_join(annot, by="ID") %>% 
           dplyr::filter(duplicated(gene)==FALSE)

df.counts<- column_to_rownames(df.counts,var="gene") %>%
            dplyr::select(-c(ID))


ssgsea_results<-test_ssgsea(counts = df.counts,
                             pathways.set= Genes2paths,
                             sampletable = sampleTable %>% dplyr::filter(sample_name %in% samples2get)
                             
                             )

###1.	differential expression between baseline and week 8:
#a.	heatmap to see if samples cluster by treatment 


outname <- "5.2.5_heatmap_MH002CR_vs_5ASA_week8"

progDir<-"/datadrive/clients/MRM/MH002_UC_201/DEV/PROGRAMS/transcriptomics/01_DE_DEseq2_MH002CR_5ASA.R"

title <- paste0(glue("5.2.5. Differentially Expressed Genes at week 8 (CFBL) comparing MH002 Clinical Resp (CR) and 5ASA 
                     Heatmap Plot."),
                '\nITT population (N = ',results_DE_MH002CR_5ASA$ITT,")")



foot1 <- glue("\nNote:  This heatmap plot shows the results of the Differentially Expressed Genes (DEG)
at week 8 (Visit 5) performed on the MH002  Clinical Resp and 5ASA populations . 
  In red are shown the genes over expressed. Samples are represented by grey (5ASA) and purple (MH002 CR group).
  Data has been  z-score normalized by row, and constrained within the range [-4, 4]..
                \nAbbreviations: CR = Clinical Response\n")
  
  
  
  foot2 <- paste0("\nProgram Path: ", progDir, "\nStatus: Reviewed")
  foot3 <- Sys.Date()
  
  
  
  pdf_heatmap(results_DE_MH002CR_5ASA, 
              "C:/Users/AnGomez/Documents/MRM/LAST/Mechanistic/MH002CR_5ASA/", 
              outname=outname, 
              title, 
            foot1, 
            foot2, 
            foot3)

##Volcano
  
  outname <- "5.2.5_volcano_MH002CR_vs_5ASA_week8"
  
title <- paste0(glue("5.2.5. Differentially Expressed Genes at week 8 (CFBL) comparing MH002 Clinical Resp (CR) and 5ASA 
                     Volcano Plot."),
                '\nITT population (N = ',results_DE_MH002CR_5ASA$ITT,")")



foot1 <- glue("\nNote:  This Volcano plot shows the results of the Differentially Expressed Genes (DEG)
at week 8 (Visit 5) performed on the MH002  Clinical Resp and 5ASA populations . 
In red are shown the genes significant.
              \nAbbreviations: CR = Clinical Response\n")



pdf_volcano(results_DE_MH002CR_5ASA$ALL, 
            "C:/Users/AnGomez/Documents/MRM/LAST/Mechanistic/MH002CR_5ASA/", 
            pv.thr=0.05,fc.thr=1.5,
            outname=outname, 
            title, 
            foot1, 
            foot2, 
            foot3)

#######################################################
##               Outputs               ##
#######################################################

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


write_results<-function(Res){
  ##subset results

  write.xlsx.custom(Res$UP, file = paste0("/Users/AnGomez/Documents/MRM/LAST/Mechanistic/MH002CR_5ASA/res_MH002CR_5ASA.xlsx"),
             sheetName = "UP_Genes", append = FALSE,row.names = F)
  # Add a second data set in a new worksheet
  write.xlsx.custom(Res$DOWN, file = paste0("/Users/AnGomez/Documents/MRM/LAST/Mechanistic/MH002CR_5ASA/res_MH002CR_5ASA.xlsx"),
             sheetName = "DOWN_Genes", append = TRUE,row.names = F)
  
  write.xlsx.custom(Res$GO_UP, file = paste0("/Users/AnGomez/Documents/MRM/LAST/Mechanistic/MH002CR_5ASA/res_MH002CR_5ASA.xlsx"),
             sheetName = "GO_enrich_UP", append = TRUE,row.names = F)
  
  write.xlsx.custom(Res$GO_Down, file = paste0("/Users/AnGomez/Documents/MRM/LAST/Mechanistic/MH002CR_5ASA/res_MH002CR_5ASA.xlsx"),
             sheetName = "GO_enrich_DOWN", append = TRUE,row.names = F)
  
  write.xlsx.custom(Res$React_Up, file = paste0("/Users/AnGomez/Documents/MRM/LAST/Mechanistic/MH002CR_5ASA/res_MH002CR_5ASA.xlsx"),
             sheetName = "Reactome_UP", append = TRUE,row.names = F)
  
  write.xlsx.custom(Res$React_Down, file = paste0("/Users/AnGomez/Documents/MRM/LAST/Mechanistic/MH002CR_5ASA/res_MH002CR_5ASA.xlsx"),
             sheetName = "Reactome_DOWN", append = TRUE,row.names = F)
  
  write.xlsx.custom(Res$Biocarta_up, file = paste0("/Users/AnGomez/Documents/MRM/LAST/Mechanistic/MH002CR_5ASA/res_MH002CR_5ASA.xlsx"),
             sheetName = "Biocarta_UP", append = TRUE,row.names = F)
  
  write.xlsx.custom(Res$Biocarta_Down, file = paste0("/Users/AnGomez/Documents/MRM/LAST/Mechanistic/MH002CR_5ASA/res_MH002CR_5ASA.xlsx"),
             sheetName = "Biocarta_DOWN", append = TRUE,row.names = F)
  
  write.xlsx.custom(Res$Wiki_Up, file = paste0("/Users/AnGomez/Documents/MRM/LAST/Mechanistic/MH002CR_5ASA/res_MH002CR_5ASA.xlsx"),
             sheetName = "Wiki_UP", append = TRUE,row.names = F)
  
  write.xlsx.custom(Res$Wiki_Down, file = paste0("/Users/AnGomez/Documents/MRM/LAST/Mechanistic/MH002CR_5ASA/res_MH002CR_5ASA.xlsx"),
             sheetName = "Wiki_DOWN", append = TRUE,row.names = F)
  
  write.xlsx.custom(Res$GSEA_UP, file = paste0("/Users/AnGomez/Documents/MRM/LAST/Mechanistic/MH002CR_5ASA/res_MH002CR_5ASA.xlsx"),
             sheetName = "GSEA_UP", append = TRUE,row.names = F)
  
  write.xlsx.custom(Res$GSEA_DOWN, file = paste0("/Users/AnGomez/Documents/MRM/LAST/Mechanistic/MH002CR_5ASA/res_MH002CR_5ASA.xlsx"),
             sheetName = "GSEA_DOWN", append = TRUE,row.names = F)
  write.xlsx.custom(ssgsea_results, file = paste0("/Users/AnGomez/Documents/MRM/LAST/Mechanistic/MH002CR_5ASA/res_MH002CR_5ASA.xlsx"),
                    sheetName = "ssGSEA", append = TRUE,row.names = F)

  write.xlsx.custom(DE_results_by_ptw, file = paste0("/Users/AnGomez/Documents/MRM/LAST/Mechanistic/MH002CR_5ASA/res_MH002CR_5ASA.xlsx"),
                    sheetName = "results_ptw", append = TRUE,row.names = F)
}

write_results(Res = results_DE_MH002CR_5ASA)







