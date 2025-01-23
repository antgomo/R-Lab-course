library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(doParallel)

library(minfi)
library(DMRcate)
library(EpiDISH)
library(data.table)
library(readxl)
setwd("/mnt/md127/Analysis_Projects/FIS_Meth/FIS_Dat/")

ann850k<-getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
ann850k<-as.data.frame(ann850k)
idat.folder <- "/mnt/md127/Analysis_Projects/FIS_Meth/FIS_Dat/idats/" ###change for wour path dir

targets <- read.metharray.sheet(base=idat.folder)
targets$Basename<-paste(idat.folder,targets$xipName, sep="")

#load("FIS_norm_cellcomp_data.RData")
###use joint normalized
load("/media/IMID/Projects/FIS_Meth_analysis/FIS_P6/FIS_P6_FINAL.RData")
###begin here after normalized beta pipeline is done


###right DAS values coming form DB
setwd("/mnt/md127/Analysis_Projects/FIS_Meth/")


load("DAS28_BD.RData")

pheno_das28<-df
colnames(pheno_das28)[1]<-"Code"

####check days from wk0 to second point
pheno_extended<-readRDS("/media/IMID/Projects/Project4-Celgene/Data-Storage/Phenotype-Data/2013_longIMID/03_Preprocessed/PHENO_DATA_EXTENDED.rds")
pheno_extended<-pheno_extended[,c("idCaso","AR_S0_E_0106","AR_S0_E_0104","AR_S0_E_0105")]
colnames(pheno_extended)<-c("Code","EULAR","Point_after_wk0","DAys_to_2nd_point")


pheno_das28<-merge(pheno_das28,pheno_extended,by="Code")

samples2rem<-pheno_das28[pheno_das28$DAS28.S0<3.2,"Code"]### samples IB-30107, IB-30127,IB-30141 and IB-30406
##only betas in blood target

targets<-targets[targets$tissue=="blood",]
targets$Code<-gsub("-S0|-S12|-SRM|-SRTB","",targets$donation)

####remove samples in baseline with DAS28<3.2
samples2rem<-pheno_das28[pheno_das28$DAS28.S0<3.2,"Code"]### samples IB-30107, IB-30127,IB-30141 and IB-30406
targets<-targets[!(targets$Code %in% samples2rem),]

betas<-getBeta(countsEPIC$normalizedData)

#betas <-betas[,colnames(betas) %in% targets$xipName]
#targets <-targets[targets$xipName %in% colnames(betas),]

#betas<-betas[,match(targets$xipName,colnames(betas))]

pD.w0<-targets[grep("S0",targets$donation,perl=T),]
phenoData<-read.csv("phenodata_FIS_fixed.csv")
##remove sexual chromosomes and SNPs probes

betas<-rmSNPandCH(betas, dist = 2, mafcut = 0.05, and = TRUE, rmcrosshyb = TRUE, rmXY=TRUE)
####run combat to adjust for sentrix Barcode and Section
#batch<-pheno$batch
#modcombat<-model.matrix(~1, data=pheno)
#combat_mydata<-ComBat(dat=mydata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)


####Comparisons

#R0 vs NR0 

library(limma)

##simple limma approach, build design matrix with Cell composition
##avoid Tocilizumab, Tofacinitib and others in treatmen

phenoData3<-phenoData[!(is.na(phenoData$Age)),]
phenoData3<-phenoData3[!(is.na(phenoData3$RESP)),]

gc()


#####merge pheno with pD samplesheet

pD.w0<-merge(pD.w0,phenoData3[,c("antiTNF","RESP","SEX","Age","Code","EULAR2")],by="Code")
pD.w0<-pD.w0[!(pD.w0$antiTNF=="Tocilizumab" | pD.w0$antiTNF=="Tofacitinib" | pD.w0$antiTNF=="Abatacept" | pD.w0$antiTNF=="Other"| pD.w0$antiTNF==""),]

betas.fis<-betas[,colnames(betas) %in% pD.w0$xipName]
betas.fis<-betas.fis[,match(as.character(pD.w0$xipName),colnames(betas.fis))]#w0
pD.w0$SEX<-droplevels(as.factor(pD.w0$SEX))

##
cell_comp<-countsEPIC$counts
cc.df<-cell_comp[rownames(cell_comp) %in% colnames(betas.fis),]
cc.df<-cc.df[match(colnames(betas.fis),rownames(cc.df)),]

gc()

##program does not use categorical
pD.w0$R<-ifelse(pD.w0$RESP=="R",1,0)


ann850k<-ann850k[!(ann850k$chr %in% c("chrX","chrY")),]
###ensure that we are dealing with same Cpgs in annot and china2
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
##first run combat
## correct for Sentrix
#batch<-as.factor(pD.w0$Sentrix.Barcode)
#modcombat<-model.matrix(~1, data=pD.w0)
#combat_mydata<-ComBat(dat=betas.fis, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

Phen<-model.matrix(~pD.w0$SEX+pD.w0$Age)
rownames(Phen)<-pD.w0$xipName

results.deconv <- CellDMC(betas.fis, pD.w0$R, cc.df,mc.cores = 6,cov.mod = Phen) 
#
###subset in cell types
##Patient UP

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

##store for validation

results.fis<-results.deconv$coe
####BLOCK2
rm(betas.fis)
gc()
####validate in P6
targets<-readRDS("/media/IMID/Projects/Project4-Celgene/Data-Storage/Phenotype-Data/2013_longIMID/03_Preprocessed/PHENO_DATA_EXTENDED.rds")
pheno<-targets[,c("idCaso","AR_S0_E_0106")]
rm(targets)
gc()
AROMA <- read_excel("/media/IMID/Databases/RA_AROMA_collected_June_2018.xls")
AROMA$Code<-gsub("-S0","",AROMA$Code)
targets<-read.csv("/mnt/md127/Analysis_Projects/FIS_Meth/TCA_results/Validation/P6/samplesheets_p6_notinFIS.csv")
colnames(targets)<-gsub(".x|.y","",colnames(targets))
###test
pheno<-pheno[pheno$idCaso %in% targets$Code,]

##get out samples in Pheno data extended
AROMA<-AROMA[!(AROMA$Code %in% pheno$idCaso),] ##
AROMA<-AROMA[AROMA$Code %in% targets$Code,]

AROMA$baseline <- as.numeric(AROMA$DAS28_S0)
AROMA$followup <- as.numeric(AROMA$DAS28_S12)
AROMA$improvement <- - (AROMA$followup-AROMA$baseline)
AROMA$EULAR_S <- "NON_RESP"
AROMA$EULAR_S[is.na(AROMA$improvement)] <- NA
AROMA$EULAR_S[AROMA$improvement > 0.6 & AROMA$improvement <= 1.2 & AROMA$followup < 5.1] <- "MOD"
AROMA$EULAR_S[AROMA$improvement > 1.2] <- "MOD"
AROMA$EULAR_S[AROMA$improvement > 1.2 & AROMA$followup < 3.2] <- "GOOD"

###updated response info
colnames(pheno)<-c("Code","EULAR_S")
response.table<-as.data.frame(rbind(AROMA[,c("Code","EULAR_S")],pheno))

###all
####subset for after comparisons


targets<-merge(targets,response.table,by="Code")
###remove duplicates
targets<-unique(targets)

###establish adequate response, MOD are resp
targets$RESP_F<-ifelse(targets$EULAR_S=="NON_RESP","NR",ifelse(is.na(targets$EULAR_S),NA,"R"))

betas<-betas[,colnames(betas) %in% targets$Name]
betas<-betas[,match(targets$Name,colnames(betas))]


###remove NAs in R
targets<-targets[!(is.na(targets$RESP_F)),]
betas<-betas[,colnames(betas) %in% targets$Name]
betas<-betas[,match(targets$Name,colnames(betas))]
cell_comp<-cell_comp[rownames(cell_comp) %in% targets$Name,]
####subset for after comparisons
pD.w0<-targets[targets$week=="wk0",]
#pD.w12<-targets[targets$week=="wk12",]


library(limma)

####wk0
betas.fis<-betas[,colnames(betas) %in% pD.w0$Name]
betas.fis<-betas.fis[,match(as.character(pD.w0$Name),colnames(betas.fis))]#w0
####wk12
###remove NAs just in case
#pD.w12<-pD.w12[!(is.na(pD.w12$DAS28_S12)),]

#betas.fis<-betas[,colnames(betas) %in% pD.w12$Name]
#betas.fis<-betas.fis[,match(as.character(pD.w12$Name),colnames(betas.fis))]#w12


##

cc.df<-cell_comp[rownames(cell_comp) %in% colnames(betas.fis),]
cc.df<-cc.df[match(colnames(betas.fis),rownames(cc.df)),]


##program does not use categorical

###design matrix

Phen<-model.matrix(~as.factor(pD.w0$S)+pD.w0$Age)
rownames(Phen)<-pD.w0$Name



####fast solution, get all cpgs for all types to test and then substract later

##FIS

results.fis.005<-lapply(results.fis, function(x) return(x[x[5]<.05,]))

fis.cpgs<-unique(unlist(lapply(results.fis.005, function(x) return(unique(rownames(x))))))

pD.w0$R<-ifelse(pD.w0$RESP_F=="R",1,0)


results.deconv <- CellDMC(betas.fis[rownames(betas.fis) %in% fis.cpgs,], pD.w0$R, cc.df,mc.cores = 7,cov.mod = Phen,adjPMethod = "fdr") 

###subset in cell types

for( i in 1:length(results.deconv$coe)){
  
  my.data<-as.data.frame(results.deconv$coe[i])
  
  ###select results in FIS
  index<-grep(names(results.deconv$coe)[i],names(results.fis.005))
  my.file<-results.fis.005[[index]]
  ###check for crossing
  
  my.data<-merge(my.data,my.file,by="row.names")
  rownames(my.data)<-my.data$Row.names
  my.data<-my.data[,-1]
  ###same direction
  my.data$comparison <- sign(my.data[,6]) == sign(my.data[,1])  ####FIS marks direction
  my.data<-my.data[my.data[,5]<.05,]###adjusted p-val from P6
  #my.data<-my.data[my.data[,4]<.05,]###non adjusted p-val from P6
  
  my.data<-my.data[my.data$comparison=="TRUE",]
  #my.data<-my.data[,c(1:5,15)]
  
  if(dim(my.data)[1]>0){
    my.data <- merge(my.data,ann850k[,c("chr","pos","UCSC_RefGene_Name","Relation_to_Island","UCSC_RefGene_Group")], by="row.names")
    rownames(my.data)<-my.data$Row.names
    my.data<-my.data[,-1]
    write.csv(my.data,paste0("/media/IMID/Projects/FIS_Meth_analysis/Paper_Material/CellType_Analysis/",names(results.deconv$coe)[i],"_RvsNR_P6_wk0_cellDMC.csv"),row.names = T)

  }else{
    write.csv(my.data,paste0("/media/IMID/Projects/FIS_Meth_analysis/Paper_Material/CellType_Analysis/",names(results.deconv$coe)[i],"_RvsNR_P6_wk0_cellDMC.csv"),row.names = T)
  }
}



