library(tximport)
library(factoextra)
library(FactoMineR)
library(dplyr)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(writexl)
library(biomaRt)
library(heatmap.2x)## approach
library(tidyverse)
library(NMF) # loading for aheatmap plotting function
library(ggsci) # more cool colours
library(ggrepel)
library(ggpubr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(GO.db)
library(DOSE)
library(annotate)
library(ggvenn)
library(readxl)

setwd("/Users/agomez/IDIBELL/Pujana/December24/")

library(RMariaDB)
library(GenomicFeatures)

txdb <- makeTxDbFromGFF("/Users/agomez/kallisto/Genomes/mus_musculus/Mus_musculus.GRCm38.96.gtf")
txdb


k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")


dir<-setwd("/Users/agomez/IDIBELL/Pujana/December24/")

samples <- read.table(file.path(dir, "samples.txt"), header = TRUE)
files <- file.path(dir,as.character(samples$sample), "abundance.h5")


names(files) <- paste0(samples$sample)

txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = TRUE)## for downstream analysis

#We could alternatively generate counts from abundances, using the argument countsFromAbundance, scaled to library size, "scaledTPM", 
#or additionally scaled using the average transcript length, averaged over samples and to library size, "lengthScaledTPM". 
#Using either of these approaches, the counts are not correlated with length, and so the length matrix should not be provided as an offset 
#for downstream analysis packages. 

txi.kallisto.tsv <- tximport(files, type = "kallisto", countsFromAbundance="lengthScaledTPM",tx2gene = tx2gene,ignoreTxVersion = TRUE)## for downstream analysis
tpm<-txi.kallisto.tsv$abundance

ensembl <- useEnsembl(biomart = "ensembl", 
                   dataset = "mmusculus_gene_ensembl",mirror = "www")

annot<-getBM(c("ensembl_gene_id", "mgi_symbol"), mart=ensembl,filters="ensembl_gene_id", values = rownames(tpm))
tpm<-tpm %>% 
  as.data.frame(.) %>%
  mutate(ensembl_gene_id=rownames(tpm))

tpm<- tpm %>%
  left_join(annot)
write_xlsx(as.data.frame(tpm),"TPM_RNAseq_Dec24.xlsx")

####PCA to check the samples

####PCA


pca <- FactoMineR::PCA(t(na.omit(as.matrix(txi.kallisto.tsv$abundance))), scale.unit = T, graph = F, ncp = 20)
factoextra::fviz_pca_ind(pca, axes = c(1,2), habillage=as.factor(samples$condition), repel = T) + coord_fixed() 


###DESEq2 downstream

#Comparisons:
  
#DEG PBC-LNK-WT vs PBC-LNK-HOMO


library(DESeq2)
#The user should make sure the rownames of sampleTable align with the colnames of txi$counts, if there are colnames. The best practice is to read sampleTable from a CSV file, and to construct files from a column of sampleTable, as was shown in the tximport examples above.


sampleTable <- data.frame(condition = factor(samples$condition))

dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, sampleTable, ~condition)

dds<- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

dds$condition<-relevel(dds$condition,ref="PBC_LNK_WT")


dds_diff <- DESeq(dds)


res <- results( dds_diff )  
head(res)##to see order fold change, WT reference
res<-as.data.frame(res)

##not enough statistical power, subset to 0,1
res<- res %>% 
      mutate(SIG=case_when(padj<0.1 ~ "YES",padj>0.1 ~ "NO"))

##remove NA

res<-res[complete.cases(res),]
#ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl") ##run only once
annot<-getBM(c("ensembl_gene_id", "mgi_symbol"), mart=ensembl,filters="ensembl_gene_id", values = rownames(res))
res<-res %>% 
     mutate(ensembl_gene_id=rownames(res))

res<- res %>%
      left_join(annot)


write_xlsx(res,"Res_comp_PBC_LNK_HOMO_PBC_LNK_WT.xlsx")

###significative

res.s<-res %>% 
       filter(SIG=="YES") 
       
###Heatmap


my.data<-as.data.frame(counts(dds_diff)) %>%
         filter(rownames(.) %in% res.s$ensembl_gene_id ) %>% 
         mutate(ensembl_gene_id=rownames(.)) %>%
         left_join(annot) %>%
         filter(!duplicated(mgi_symbol)) %>%
         column_to_rownames(var="ensembl_gene_id")  %>% 
         dplyr::select(-c(mgi_symbol))
         
  
  
color_samples<-ifelse(dds_diff$condition=="PBC_LNK_HOMO","red","blue1")

#####scale

data <- t(scale(t(my.data))) # z-score normalise each row (feature)
data <- apply(data, 2, function(x) ifelse(x > 4, 4, x)) # compress data within range
data <- apply(data, 2, function(x) ifelse(x < -4, -4, x)) # compress data within range


spcol<-color_samples
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)

####merge to get gene names


heatmap.2x(as.matrix(data), col=rev(cols), Colv = T,Rowv=TRUE, scale="none", ColSideColors=spcol,
           trace="none", dendrogram="both", 
           cexRow=1, cexCol=1,
           main="PBC_LNK_HOMO vs PBC_WT_HOMO",
           #labCol=NA,
           labRow=NA, 
           density.info="none",
           hclust=function(x) hclust(x,method="complete"),
           distfun=function(x) as.dist((1-cor(t(x)))/2),srtCol=45
           
)

legend("topright",c("PBC_LNK_HOMO","PBC_WT_HOMO"),pch=20:20,col=c("red","blue1"))






#####Split in UP and DOWN
#"By default, R will choose a reference level for factors based on alphabetical order. 
#Then, if you never tell the DESeq2 functions which level you want to compare against (e.g. which level represents the control group), 
#the comparisons will be based on the alphabetical order of the levels. There are two solutions: you can either explicitly tell 
#results which comparison to make using the contrast argument (this will be shown later), or you can explicitly set the factors levels."

UP<-res.s %>% 
   filter(log2FoldChange>0)

DOWN<-res.s %>% 
  filter(log2FoldChange<0)


####


#####GO

#####GO enrichment analysis


genesid<-UP %>% 
         pull(ensembl_gene_id)
genesid<-DOWN %>% 
      pull(ensembl_gene_id)


genesid<-genesid[ genesid != "" ]##remove empty elements from a vector
genesid <- genesid[!is.na(genesid)]

ego2 <- enrichGO(gene         = genesid,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "none",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05,
                 readable=T
)


dotplot(ego2, title="PBC_LNK_HOMO DOWN",showCategory=15)


write_xlsx(ego2@result,"GO_BP_PBC_LNK_HOMO_UP.xlsx")
write_xlsx(ego2@result,"GO_BP_PBC_LNK_HOMO_DOWN.xlsx")


#simplify <- function(enrichResult, cutoff=0.7, by="p.adjust", select_fun=min) {
  ## GO terms that have semantic similarity higher than `cutoff` are treated as redundant terms
  ## select one representative term by applying `select_fun` to feature specifying by `by`.
  ## user can defined their own `select_fun` function.
  
  ## return an updated `enrichResult` object.
#}
##simplify output from enrichGO by removing redundancy of enriched GO terms
#egobp2 <- clusterProfiler::simplify(ego2, cutoff=0.7, by="p.adjust", select_fun=min,measure="Wang")


####Volcano

res.limma<-res %>%
           mutate(Gene=ifelse(is.na(mgi_symbol) | mgi_symbol ==" ",ensembl_gene_id,mgi_symbol)) %>%
           mutate(SIG=ifelse(padj<0.1,"YES","NO")) %>%
           mutate(pvalue=as.numeric(ifelse(pvalue==0,"2.2e-16",pvalue))) %>%
           mutate(LOG=-log10(pvalue))


vp<-ggplot(data=res.limma,aes(x=log2FoldChange,y=-log10(pvalue)))+
  geom_point(aes(size=-log(pvalue),fill=SIG), colour="black",shape=21)+
  geom_vline(xintercept=0,linetype="dashed",col="red")+
  geom_vline(xintercept=c(-2,2),linetype="dashed",col="blue")+
  #ggtitle("Ulcerative Colitis HI vs LOW") +
  theme(legend.position="none")+
  theme_minimal()





ap<-vp+
  
  geom_text_repel(data=filter(res.limma, res.limma$SIG=="YES"  & log2FoldChange>2), 
                       
                       aes(label=Gene),
                       nudge_x       = 10 - subset(res.limma, res.limma$SIG=="YES"  & log2FoldChange>2)$log2FoldChange,
                       direction    = "y",
                       segment.color = "grey50",
                       hjust        = 0,
                       segment.size = 0.1,
                  max.overlaps=15)+
  
  geom_text_repel(data=filter(res.limma, res.limma$SIG=="YES"  & log2FoldChange< -2), 
                  
                  aes(label=Gene),
                  nudge_x       = -10 - subset(res.limma, res.limma$SIG=="YES"  & log2FoldChange< -2)$log2FoldChangeOG,
                  direction    = "y",
                  segment.color = "grey50",
                  hjust        = 1,
                  segment.size = 0.1,
                  max.overlaps=15)
  
annotate_figure(ap,
                top = text_grob("PBC_LNK_HOMO vs PBC_LNK_WT", color = "black", size = 14)
                # bottom = text_grob("Data source: \n mtcars data set", color = "blue",
                #  hjust = 1, x = 1, face = "italic", size = 10)
                
                # fig.lab = "Figure 1", fig.lab.face = "bold"
)


#Do DEG in 2 overlap with DEG in 4? Is the loss-of-function similar to genetic deletion?

comp_prev<-read_xlsx("../September24/Res_comp_shRNA_PLKO_significative.xlsx") %>%
       filter(SIG=="YES")

actual<-res %>%
  filter(SIG=="YES") %>%
  column_to_rownames(.,var="mgi_symbol")
#%>%
 # mutate(Gene=toupper(mgi_symbol))


gene_df <- orthogene::convert_orthologs(gene_df = actual,
                                        gene_input = "rownames", 
                                        gene_output = "rownames", 
                                        input_species = "mouse",
                                        output_species = "human",
                                        non121_strategy = "drop_both_species",
                                        method = "homologene") %>%
                                       mutate(hgnc_symbol=rownames(.))

my.list<-list(Analysis_sept24=comp_prev$hgnc_symbol,
              Analysis_dec24=rownames(gene_df)) 


###do the venn


venn.plot<-ggvenn(my.list, c("Analysis_sept24", "Analysis_dec24"))       # Venn diagram with three sets

my.genes<-intersect(comp_prev$hgnc_symbol,rownames(gene_df))

my.genes.df<-comp_prev %>% 
             dplyr::select(hgnc_symbol,log2FoldChange,padj) %>%
             filter(hgnc_symbol %in% my.genes & !is.na(hgnc_symbol)) %>%
             dplyr::rename(LFC_sep24=log2FoldChange,padj_sep24=padj) %>% 
             left_join(gene_df %>% dplyr::select(hgnc_symbol,log2FoldChange,padj) %>% 
                        dplyr::rename(LFC_dec24=log2FoldChange,padj_dec24=padj)  %>%
                        filter(hgnc_symbol %in% my.genes & !is.na(hgnc_symbol)),by="hgnc_symbol") %>%
            mutate(padj_sep24=-log10(padj_sep24),padj_dec24=-log10(padj_dec24))
   

write_xlsx(my.genes.df,"Genes_intersection.xlsx")


library("ggpubr")
ggscatter(my.genes.df, x = "LFC_sep24", y = "LFC_dec24", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",alpha = 0.5,
          add.params= list(color="blue",fill="lightgray"),
          xlab = "September 24", ylab = "December 24") +ggtitle("Log Fold Change correlation")

ggscatter(my.genes.df, x = "padj_comp2", y = "padj_comp4", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",alpha = 0.5,
          add.params= list(color="blue",fill="lightgray"),
          xlab = "LFC Loss of Function", ylab = "LFC Genetic deletion") +ggtitle("Adjusted pvals correlation")

####pathways
comp_prev<-read_xlsx("../September24/GO_BP_shRNA_PLKO_DOWN.xlsx") %>%
  filter(qvalue < 0.05)

actual<-ego2@result %>%
  filter(qvalue < 0.05)
#%>%
# mutate(Gene=toupper(mgi_symbol))


my.list<-list(Analysis_sept24=comp_prev$ID,
              Analysis_dec24=actual$ID) 


###do the venn


venn.plot<-ggvenn(my.list, c("Analysis_sept24", "Analysis_dec24"))       # Venn diagram with three sets

my.genes<-intersect(comp_prev$hgnc_symbol,rownames(gene_df))

my.genes.df<-comp_prev %>% 
  dplyr::select(hgnc_symbol,log2FoldChange,padj) %>%
  filter(hgnc_symbol %in% my.genes & !is.na(hgnc_symbol)) %>%
  dplyr::rename(LFC_sep24=log2FoldChange,padj_sep24=padj) %>% 
  left_join(gene_df %>% dplyr::select(hgnc_symbol,log2FoldChange,padj) %>% 
              dplyr::rename(LFC_dec24=log2FoldChange,padj_dec24=padj)  %>%
              filter(hgnc_symbol %in% my.genes & !is.na(hgnc_symbol)),by="hgnc_symbol") %>%
  mutate(padj_sep24=-log10(padj_sep24),padj_dec24=-log10(padj_dec24))


write_xlsx(my.genes.df,"Genes_intersection.xlsx")


library("ggpubr")
ggscatter(my.genes.df, x = "LFC_sep24", y = "LFC_dec24", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",alpha = 0.5,
          add.params= list(color="blue",fill="lightgray"),
          xlab = "September 24", ylab = "December 24") +ggtitle("Log Fold Change correlation")

ggscatter(my.genes.df, x = "padj_comp2", y = "padj_comp4", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",alpha = 0.5,
          add.params= list(color="blue",fill="lightgray"),
          xlab = "LFC Loss of Function", ylab = "LFC Genetic deletion") +ggtitle("Adjusted pvals correlation")
