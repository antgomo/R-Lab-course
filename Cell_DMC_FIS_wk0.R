library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(doParallel)
library(minfi)
library(DMRcate)
library(EpiDISH)
library(data.table)
library(readxl)


setwd("/media/IMID/Projects/Internal_Projects/2019_GRR_FISMethAnalysis/Paper_Material/Predictor/")
ann850k<-getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
ann850k<-as.data.frame(ann850k)
ann850k<-ann850k[!(ann850k$chr %in% c("chrX","chrY")),]


##FIS original samplesheet
targets <- read.csv("/media/IMID/Projects/Internal_Projects/2019_GRR_FISMethAnalysis/Paper_Material/Targets_all.csv")
targets <- read.csv("//10.4.2.66/binf/Projects/Internal_Projects/2019_GRR_FISMethAnalysis/Paper_Material/Targets_all.csv")
targets$Basename<-paste(idat.folder,targets$xipName, sep="")  # ERROR

###use joint normalized FIS and P6 Data
load("/media/IMID/Projects/Internal_Projects/2019_GRR_FISMethAnalysis/FIS_P6/FIS_P6_FINAL.RData")
load("//10.4.2.66/binf/Projects/Internal_Projects/2019_GRR_FISMethAnalysis/FIS_P6/FIS_P6_FINAL.RData")

###right DAS values coming form DB, in order to remove samples below DAS28<3.2 at wk0
load("/media/IMID/Projects/Internal_Projects/2019_GRR_FISMethAnalysis/Paper_Material/DAS28_BD.RData")
load("//10.4.2.66/binf/Projects/Internal_Projects/2019_GRR_FISMethAnalysis/Paper_Material/DAS28_BD.RData")
pheno_das28<-df
colnames(pheno_das28)[1]<-"Code"

####check days from wk0 to second point
pheno_extended<-readRDS("/media/IMID/Data/Storage/Clinical/2012_INNPACTO_IB/03_Preprocessed/clinicalData_v1.rds")
pheno_extended<-readRDS("//10.4.2.66/binf/Data/Storage/Clinical/2012_INNPACTO_IB/03_Preprocessed/clinicalData_v1.rds")
pheno_extended<-pheno_extended[,c("idCaso","AR_S0_E_0106","AR_S0_E_0104","AR_S0_E_0105")]
colnames(pheno_extended)<-c("Code","EULAR","Point_after_wk0","DAys_to_2nd_point")
pheno_das28<-merge(pheno_das28,pheno_extended,by="Code")

##only betas in blood target
targets<-targets[targets$tissue=="blood",]
targets$Code<-gsub("-S0|-S12|-SRM|-SRTB","",targets$donation)

####remove samples in baseline with DAS28<3.2, for this reason we loaded DAS28 right values
samples2rem<-pheno_das28[pheno_das28$DAS28.S0<3.2,"Code"]### samples IB-30107, IB-30127,IB-30141 and IB-30406
targets<-targets[!(targets$Code %in% samples2rem),]; dim(targets) # 186 samples

##get betas from normalized data
betas<-getBeta(countsEPIC$normalizedData)

##week0
pD.w0<-targets[grep("S0",targets$donation,perl=T),]; dim(pD.w0) # 101

##Load PhenoData with, sex and age
phenoData<-read.csv("/media/IMID/Projects/Internal_Projects/2019_GRR_FISMethAnalysis/Paper_Material/phenodata_FIS_fixed.csv")
phenoData<-read.csv("//10.4.2.66/binf/Projects/Internal_Projects/2019_GRR_FISMethAnalysis/Paper_Material/phenodata_FIS_fixed.csv")
##remove sexual chromosomes and SNPs probes
betas<-rmSNPandCH(betas, dist = 2, mafcut = 0.05, and = TRUE, rmcrosshyb = TRUE, rmXY=TRUE)

####Comparisons
#R0 vs NR0 
library(limma)

##simple limma approach, build design matrix with Cell composition
phenoData3<-phenoData[!(is.na(phenoData$Age)),]
phenoData3<-phenoData3[!(is.na(phenoData3$RESP)),]
gc()

#####merge pheno with pD samplesheet
pD.w0<-merge(pD.w0,phenoData3[,c("antiTNF","RESP","SEX","Age","Code","EULAR2")],by="Code"); dim(pD.w0) # 98
##avoid Tocilizumab, Tofacinitib and others in treatment
pD.w0<-pD.w0[!(pD.w0$antiTNF=="Tocilizumab" | pD.w0$antiTNF=="Tofacitinib" | pD.w0$antiTNF=="Abatacept" | pD.w0$antiTNF=="Other"| pD.w0$antiTNF==""),]
dim(pD.w0) # 62 patients

##ßelect betas from FIS only
betas.fis<-betas[,colnames(betas) %in% pD.w0$xipName]
##same order as in samplesheet
betas.fis<-betas.fis[,match(as.character(pD.w0$xipName),colnames(betas.fis))]#w0
pD.w0$SEX<-droplevels(as.factor(pD.w0$SEX))

##load cell type
##
cell_comp<-countsEPIC$counts
##only FIS samples
cc.df<-cell_comp[rownames(cell_comp) %in% colnames(betas.fis),]
##same order as in samplesheet
cc.df<-cc.df[match(colnames(betas.fis),rownames(cc.df)),]
gc()

##program does not use categorica
pD.w0$R<-ifelse(pD.w0$RESP=="R",1,0)  # Response: 1, Non-Response: 0
ann850k<-ann850k[rownames(ann850k) %in% rownames(betas.fis),]


###design matrix
#A list with the following two items.
#dmct A matrix gives wheter the input CpGs are DMCTs and DMCs. 
#The first column tells whether a CpG is a DMC or not. 
#If the CpG is called as DMC, the value will be 1, otherwise it is 0.
#The following columns give DMCTs for each cell-type.
#If a CpG is a DMCT, the value will be1 (hypermethylated for 
#case compared to control) or -1 (hypomethylated for case compared to
#control).
#Otherwise, the value is 0 (non-DMCT). The rows of this matrix are ordered as the same
#as that of the input beta.m.

#coe This list contains several dataframes, which correspond to each cell-type
#in frac.m. Each dataframe contains all CpGs in input beta.m. 
#All dataframes contain estimated DNAm changes
#(Estimate), standard error (SE), estimated t statistics (t), raw P values (p), 
#and multiple hypothesiscorrected P values (adjP).
Phen<-model.matrix(~pD.w0$SEX+pD.w0$Age)##covariates
rownames(Phen)<-pD.w0$xipName

##run celDMC for each cell type
library(EpiDISH)
results.deconv <- CellDMC(betas.fis, pD.w0$R, cc.df,mc.cores = 6,cov.mod = Phen) 
#
###subset in cell types

##subset
for( i in 1:length(results.deconv$coe)){
  
  my.data<-as.data.frame(results.deconv$coe[i])
  my.data$SIG<-ifelse(my.data[,5]<.05,"TRUE","FALSE")
  #my.data<-my.data[my.data[,5]<.05,]
  
  if(dim(my.data)[1]>0){
    my.data <- merge(my.data,ann850k[,c("chr","pos","UCSC_RefGene_Name","Relation_to_Island","UCSC_RefGene_Group")], by="row.names")
    rownames(my.data)<-my.data$Row.names
    my.data<-my.data[,-1]
    write.csv(my.data,paste0("/media/IMID/Projects/FIS_Meth_analysis/Paper_Material/CellType_Analysis/",names(results.deconv$coe)[i],"_Disc_RvsNR_wk0_cellDMC.csv"),row.names = T)
    
  }else{
    write.csv(my.data,paste0("/media/IMID/Projects/FIS_Meth_analysis/Paper_Material/CellType_Analysis/",names(results.deconv$coe)[i],"_Disc_RvsNR_wk0_cellDMC.csv"),row.names = T)
  }
}

##store for validation, list with all cell types
results.fis<-results.deconv$coe
