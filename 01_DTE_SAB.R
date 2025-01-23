##" Differential Transcript Expression (DTE) analysis"
library(dplyr)
library(tidyr)
library(ggplot2)
library(sleuth)
library(annotables)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GO.db)
library(DOSE)
library(annotate)
library(writexl)



###first get the files that have run with kallisto
####aquí els teus outputs de Kallisto
files<-list.dirs('/home/proteome/bioinformatics/Data/Storage/Molecular/GeneSplicing/WholeBlood/2017_IMX_P4RNAseq_NovaSeqK_IA/01_Abundance/Kallisto_Transcriptome_P4',recursive=F)
names(files)<-gsub("/home/proteome/bioinformatics/Data/Storage/Molecular/GeneSplicing/WholeBlood/2017_IMX_P4RNAseq_NovaSeqK_IA/01_Abundance/Kallisto_Transcriptome_P4/","",files)
names(files)<-gsub("_","-",names(files))
summarydata<-data.frame(path=files,sample=names(files))

###aquí el pheno q vols testar, més info extra. En aquest cas consideràvem els cell types, bueno ho treus
###add extra info

pheno<-read.csv("/home/proteome/bioinformatics/Data/Storage/Molecular/GeneExpression/WholeBlood/2017_IMX_P4RNAseq_NovaSeqK_IA/Metadata/Old/Final_Table_GSL.ID.csv")
pheno$sample<-pheno$GSL.ID

####same order as splicing matrix
summarydata<-merge(summarydata,pheno[,c("sample","sex","age","activity","new.plate","celgIdAll","imid")],by="sample")
summarydata$batch<-as.factor(paste0("B",summarydata$new.plate))
##get cells
cells<-readRDS("/home/proteome/bioinformatics/Data/Storage/Molecular/GeneExpression/WholeBlood/2017_IMX_P4RNAseq_NovaSeqK_IA/DerivedFeatures/Quantitative/Cell_Estimates/cell_proportions.rds")
cells$celgIdAll<-rownames(cells)

#summarydata<-merge(summarydata,cells[,c("Granulocytes","celgIdAll")],by="celgIdAll")
summarydata$sex<-as.factor(summarydata$sex)


##get genes to transcripts
##important, aquesta és la DB on passen esl Gens a trànscrits

t2g <- grch37 %>% 
  dplyr::select(symbol, ensgene) %>% 
  dplyr::inner_join(grch37_tx2gene) %>% 
  dplyr::rename(target_id = enstxp, 
                ens_gene = ensgene, 
                ext_gene = symbol)

###estableix el disseny
##establish design
design <- ~ sex + age + batch +imid

#Select samples

my.imid<-c("SLE","RA","UC","PSA","PS","CD")

for( i in 1:length(my.imid)){ 

  
    my.sel<-c("CTRL",my.imid[i])

    summarydf<-summarydata[summarydata$imid %in% my.sel,]
    summarydf$imid<-as.factor(summarydf$imid)
    summarydf$path<-as.character(summarydf$path)


#Import data
#Main point is we have to use sleuth_prep outside Rstudio GUI, otherwise we cannot use parallelization and hence, the loading is very slow
# Create sleuth object for analysis 
###aquí és tà el tema, 
    so <- sleuth_prep(summarydf, 
                   full_model = design, 
                   read_bootstrap_tpm = T,###per anar més ràpid posa F, però en el teu cas són poques mostres, i això no està de més
                   extra_bootstrap_summary = F,##idem q a dalt
                   num_cores=7,##num_cores
                   transformation_function = function(x) log2(x + 0.5))  	
    
    
    
    #full model
    so <- sleuth_fit(so, ~sex + age+batch+imid, 'full')

    ##only covs
    so <- sleuth_fit(so, ~sex + age+batch, 'imid')
   

  ##Results

    so.res <- sleuth_lrt(so, 'imid', 'full')#test here, LRT
    so.betas <- sleuth_wt(so, paste0("imid",my.imid[i]), 'full')##betas here. És el mateix q abans però amb estimació beta 'b', tindràs direcció
    # Sleuth will transform elements in the condition field to 0s and 1s in alphabetical order 
    ##hence CTRL beta values will be relative to the 1 condition for SLE, RA, UC, PS and PSA, except in  CD, wich will be 0
    full_results <- sleuth_results(so.res, 'imid:full', 'lrt', show_all = FALSE)###results LRT
    full_betas <- sleuth_results(so.betas, paste0("imid",my.imid[i]), 'wt', show_all = FALSE)##results Betas
    ###a veure, jo vaig ser molt  "pejigueras" pq diuen q els pvals de LRT molen més, però ambs els de sleuth_wt, ja tens pval i beta, així q pass de LRT
    full_results<-merge(full_results,t2g,by="target_id")
    full_results<-merge(full_results,full_betas[,c("target_id","b")],by="target_id")
    
  
    sleuth_significant <- dplyr::filter(full_results, qval <= 0.05)##subset signifi

        saveRDS(full_results,file=paste0("/home/proteome/bioinformatics/Projects/Internal_Projects/2021_IMX_IAP4_Mining/Results/TranscriptExpression/WholeBlood/01_Features/01_clinicalAssociation/01_caseCtrl/",my.imid[i],"_SAB",".rds"))
    
  #####GO enrichment analysis

   genesid<-sleuth_significant$ext_gene
   genesid<-unique(genesid)

   genesid <- genesid[!is.na(genesid)]

   eg<-bitr(genesid, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

   ego2 <- enrichGO(gene         = eg$ENTREZID,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05,
                 readable=T
   )
  
   
  # saveRDS(ego2,file=paste0("/home/proteome/bioinformatics/Projects/Internal_Projects/2021_IMX_IAP4_Mining/Results/TranscriptExpression/WholeBlood/01_Features/02_Enrichments/01_caseCtrl/",my.imid[i],"_SAB",".rds"))
  
   


   }
  
