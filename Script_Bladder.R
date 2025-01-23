library(minfi)
library(DMRcate)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(factoextra)
library(FactoMineR)
library(limma)
library(writexl)
library(readxl)
library(methylGSA)


#- Samples to remove before starting the analyses: only in bladder cancer. Samples to remove are highlighted in red font and annotated as "recurrence".

#- 4th level of analysis: clinicopathological investigations within each cancer type
#4.3. Kidney cancer: Stage I versus others (II, III, IV); Possible covariates: age, gender.

setwd("/Users/agomez/IDIBELL/Joao/alltums/")

idat.folder <- "/Users/agomez/IDIBELL/Joao/alltums/" ###change for wour path dir

targets <- read.metharray.sheet(base=idat.folder)

###load annotation
ann450k<-getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k<-as.data.frame(ann450k)

###samples to remove in Bladder Cancer

torem<-c("TB189","TB202","TB203","TB208","TB204","TB186","TB194","TB207","TB205",
         "TB200","TB209")
targets$ID<-gsub("/Users/agomez/IDIBELL/Joao/alltums/","",targets$Basename)

load("../Res_Jul21/Betas_Norm_ALL.RData")


####Load Specific Phenos
Pheno_Blad<-read.csv("../Pheno_Bladder.csv",sep=";")
Pheno_Kid<-read.csv("../Pheno_Kidney.csv")
###select all bladder

betas2test<-betas[,colnames(betas) %in% Pheno_Blad$ID]
###same order
###merge with general DF, change ID
colnames(targets)[1]<-"ID"
Pheno_Blad<-merge(Pheno_Blad,targets[,c(1,6)],by="ID")
Pheno_Blad<-Pheno_Blad[match(colnames(betas2test),Pheno_Blad$ID),]

##DMP finding
###sex and age as covariates, in this case only sex

####Variable of interest as factor
Pheno_Blad$Type<-as.factor(Pheno_Blad$Type)
##Normal always as referene, hence Cancer is UP
Pheno_Blad$Type<-relevel(Pheno_Blad$Type,ref = "Bladder Normal")
design <- model.matrix(~as.factor(Pheno_Blad$Type)+Pheno_Blad$Age+Pheno_Blad$Gender)
colnames(design)[2]<-"Comp" ##Patient UP

##############

fit <- lmFit(betas2test, design)

# fit the contrasts
fit2 <- eBayes(fit)
results<-topTable(fit2,coef="Comp",n= dim(fit2)[1],adjust.method = "fdr", sort.by="P")

dmp<-results[results$adj.P.Val<.05,]###


dmp <- merge(dmp,ann450k[,c("UCSC_RefGene_Name","Relation_to_Island","UCSC_RefGene_Group")], by="row.names")
dmp<-as.data.frame(dmp)
#rownames(dmp)<-dmp$Row.names
#dmp<-dmp[,-1]
###merge also with betas

write_xlsx(dmp,"../Res_Jul21/Annotated_dmps_Bladder_Canc_Normal.xlsx")

####only prom + islands

proms<-dmp[grep("TSS1500|TSS200|5'UTR|1stExon",dmp$UCSC_RefGene_Group),]
proms.island<-proms[proms$Relation_to_Island=="Island",]
write_xlsx(proms.island,"../Res_Jul21/Annotated_dmps_proms_islands_Bladder_Cancer_Normal.xlsx")

library(M3C)
library(NMF) # loading for aheatmap plotting function
library(ggsci) # more cool colours
library(heatmap.2x)
library(RColorBrewer)

annot<-Pheno_Blad[,c(1,10)]
colnames(annot)[1]<-"ID"
rownames(annot)<-annot$ID


#my.data<-as.data.frame(betas2test[rownames(betas2test) %in% proms.island$Row.names ,])
my.data<-as.data.frame(betas2test[rownames(betas2test) %in% dmp$Row.names ,])


samples<-ifelse(grepl("Cancer",Pheno_Blad$Type),"red","blue1")

#####M3C

data <- t(scale(t(my.data))) # z-score normalise each row (feature)
data <- apply(data, 2, function(x) ifelse(x > 4, 4, x)) # compress data within range
data <- apply(data, 2, function(x) ifelse(x < -4, -4, x)) # compress data within range


spcol<-samples
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)

###For dmp only use first 10000
heatmap.2x(as.matrix(data[1:1000,]), col=rev(cols), Colv = T,Rowv=TRUE, scale="none", ColSideColors=spcol,
           
#heatmap.2x(as.matrix(data), col=rev(cols), Colv = T,Rowv=TRUE, scale="none", ColSideColors=spcol,
           trace="none", dendrogram="column", 
           cexRow=1, cexCol=.7,
           main="Bladder Cancer vs Normal ALL",
           # labCol=NA,
           labRow=NA, 
           density.info="none",
           hclust=function(x) hclust(x,method="complete"),
           distfun=function(x) as.dist((1-cor(t(x)))/2),srtCol=45
           
)

coords <- locator(1)
legend(coords,legend=c("Cancer","Normal"),
       fill=c("red","blue1"), border=FALSE, 
       y.intersp = 0.7, 
       xpd=TRUE, inset=c(0,1), cex=1, bty='n')

####subset first results to prom + islands
results <- merge(results,ann450k[,c("UCSC_RefGene_Name","Relation_to_Island","UCSC_RefGene_Group")], by="row.names")
results<-as.data.frame(results)
results<-results[grep("TSS1500|TSS200|5'UTR|1stExon",results$UCSC_RefGene_Group),]
results<-results[results$Relation_to_Island=="Island",]


results["group"] <- "NotSignificant"
results[which(results['adj.P.Val'] < 0.05 & results['logFC'] > 0.2),"group"] <- "Hypermethylated"
results[which(results['adj.P.Val'] < 0.05 & results['logFC'] < -0.2),"group"] <- "Hypomethylated"

##Volcano
vp<-ggplot(data=results,aes(x=logFC,y=-log10(adj.P.Val),color=group))+
  geom_point(alpha=.6, size=1.2)+
  scale_colour_manual(values = c("green", "red", "lightgrey"))+
  geom_vline(xintercept=-0.2,linetype="dashed",col="black")+
  geom_vline(xintercept=0.2,linetype="dashed",col="black")+
  geom_hline(yintercept=1.3,linetype="dashed",col="black")+
  xlim(-1,1)+
  #geom_vline(xintercept=c(-0.5,0.5),linetype="dashed",col="blue")+
  ggtitle("Bladder Cancer vs Normal Proms + Islands") +
  theme(legend.position="none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Delta Beta")

##add names
#vp<-vp+geom_text_repel(data=dplyr::filter(res, logFC < -.1 & rawP <.05 | logFC > .1 & rawP <.05 ), aes(label=ID),point.padding = 2)

tiff(paste0("VP_Bladder","ALL",".tiff"),units="in", width=5, height=5 ,res=300)
vp
dev.off()






#####GO enrichment analysis

#####GO
library(methylGSA)


###proms
pv.vals<-proms.island$P.Value

names(pv.vals)<-rownames(proms.island)
##all dmps
pv.vals<-dmp$P.Value

names(pv.vals)<-dmp$Row.names


res1 <- methylglm(cpg.pval = pv.vals, minsize = 10, maxsize = 1000, GS.type = "GO",array.type = "450K")

res1<-res1[res1$padj<.05,]

write_xlsx(res1,"../Res_Jul21/GO_enrich_Prostate_Canc_Normal_ALL_DMP.xlsx")

#4.2. Bladder cancer: 
#MIBC versus NMIBC; 
#recurrence versus non-recurrence; 
#papillary low grade versus papillary high grade (this is for NMIBC only); 
#papillary low grade versus papillary high grade versus invasive high grade (here a comparison of 3 groups, if possible). Possible covariates: age, gender.

###
canc<-Pheno_Blad[Pheno_Blad$Type=="Bladder Cancer",]

betascanc<-betas[,colnames(betas) %in% canc$ID]

###same order
canc<-canc[match(colnames(betascanc),canc$ID),]

##DMP finding
###sex and age as covariates, in this case only sex

####Variable of interest as factor
###MIBC versus NMIBC; 

canc$Comp<-as.factor(canc$Muscle.vs.non.muscle)

##recurrence versus non-recurrence
canc$Comp2<-canc$Recurrence.yes.no
##Papillary
canc$Comp3<-canc$Tumor.Classification
####Run depending on the comparison you want to test
design <- model.matrix(~as.factor(canc$Comp2)+canc$Age+as.factor(canc$Gender))##MIBC
#design <- model.matrix(~as.factor(canc$Comp2)+canc$Age+as.numeric(canc$PSA.ao.dx))##TumorStage
###In Papillary, remove Invasive
canc2<-canc[!(grepl("Invasive",canc$Comp3)),]
design <- model.matrix(~as.factor(canc2$Comp3)+canc2$Age+as.factor(canc2$Gender))##Papillary


colnames(design)[2]<-"Comp" ##Patient UP

##############

#fit <- lmFit(betascanc, design)##

fit <- lmFit(betascanc[,colnames(betascanc) %in% canc2$ID] , design)###Paiilary

# fit the contrasts
fit2 <- eBayes(fit)
results<-topTable(fit2,coef="Comp",n= dim(fit2)[1],adjust.method = "fdr", sort.by="P")

#dmp<-results[results$adj.P.Val<.05,]##
dmp<-results[results$P.Value<.01,]##


dmp <- merge(dmp,ann450k[,c("UCSC_RefGene_Name","Relation_to_Island","UCSC_RefGene_Group")], by="row.names")
dmp<-as.data.frame(dmp)


#write_xlsx(dmp,"../Res_Jul21/Annotated_dmps_Bladder_MIBCvsNIMBC.xlsx")
#write_xlsx(dmp,"../Res_Jul21/Annotated_dmps_Bladder_RecvsNon_Rec.xlsx")
write_xlsx(dmp,"../Res_Jul21/Annotated_dmps_Bladder_Papillary.xlsx")

##########comparison of three
#papillary low grade versus papillary high grade versus invasive high grade (here a comparison of 3 groups, if possible). Possible covariates: age, gender.
my.data<-as.data.frame(t(betascanc))

rownames(canc)<-canc$ID

my.data<-merge(my.data,canc[,c(2,3,4)],by="row.names")
my.data<-my.data[,-1]##remove ID
final<-dim(my.data)[2]-3

##store the results

my.res<-matrix(nrow = final,ncol = 3)

for (i in 1:final) {#take into account last columns are to use
  y.df<-my.data[,c(i,432005,432006,432007)]
  result<-aov(y.df[,1]~as.factor(y.df$Tumor.Classification)+y.df$Age+as.factor(y.df$Gender))
  
  res<-summary(result)
  pval<-res[[1]][["Pr(>F)"]][1]
  pval<-signif(pval, digits=3)
  my.res[i,1]<-colnames(my.data)[i]
  my.res[i,2]<-pval
  my.res[i,3]<-p.adjust(pval, method = "fdr", n = final)
  
}

my.res<-as.data.frame(my.res)
my.res<-my.res[my.res$V2<.01,]

####
###heatmap plot
dmp<-read_xlsx("/Users/agomez/IDIBELL/Joao/Res_Dec21/Kidney/NormalCancer/Annotated_dmps_Kidney_Canc_Normal.xlsx")


#annot<-canc[,c(1,12)]##MIBC
#annot<-canc[,c(1,12)]##Recurrence
#annot<-canc2[,c(1,13)]##Pappillary
annot<-canc[,c(1,4)]##Three comp

colnames(annot)[1]<-"ID"
rownames(annot)<-annot$ID

#my.data<-as.data.frame(betascanc[rownames(betascanc) %in% dmp$Row.names ,])
#my.data<-as.data.frame(betascanc[rownames(betascanc) %in% rownames(dmp) ,colnames(betascanc) %in% canc2$ID])##Papillary
my.data<-as.data.frame(betascanc[rownames(betascanc) %in% my.res$V1 ,])


#samples<-ifelse(grepl("NMIBC",canc$Comp),"red","blue1")###MIBC
#samples<-ifelse(grepl("Yes",canc$Comp2),"red","blue1")###Recurrence
#samples<-ifelse(grepl("high",canc2$Comp3),"red","blue1")###Pappillary
samples<-ifelse(grepl("Invasive",canc$Tumor.Classification),"red",ifelse(grepl("high",canc$Tumor.Classification),"green","blue1"))###Pappillary

#####M3C

data <- t(scale(t(my.data))) # z-score normalise each row (feature)
data <- apply(data, 2, function(x) ifelse(x > 4, 4, x)) # compress data within range
data <- apply(data, 2, function(x) ifelse(x < -4, -4, x)) # compress data within range


spcol<-samples
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)

heatmap.2x(as.matrix(data), col=rev(cols), Colv = T,Rowv=TRUE, scale="none", ColSideColors=spcol,
           trace="none", dendrogram="column", 
           cexRow=1, cexCol=.7,
           main="Bladder Cancer Papill Grades",
           # labCol=NA,
           labRow=NA, 
           density.info="none",
           hclust=function(x) hclust(x,method="complete"),
           distfun=function(x) as.dist((1-cor(t(x)))/2),srtCol=45
           
)
rownames(my.res)<-my.res$V1
dmp <- merge(my.res,ann450k[,c("UCSC_RefGene_Name","Relation_to_Island","UCSC_RefGene_Group")], by="row.names")
dmp<-as.data.frame(dmp)

write_xlsx(dmp,"../Res_Jul21/Annotated_dmps_Bladder_Papillary_Grades.xlsx")

#####GO enrichment analysis

#####GO

#pv.vals<-dmp$P.Value
pv.vals<-as.numeric(dmp$V2)#Papillary Grades

names(pv.vals)<-dmp$Row.names


res1 <- methylglm(cpg.pval = pv.vals, minsize = 10, maxsize = 1000, GS.type = "GO",array.type = "450K")

res1<-res1[res1$padj<.05,]



#write_xlsx(res1,"../Res_Jul21/GO_enrich_Bladder_MIBC_Comp.xlsx")
#write_xlsx(res1,"../Res_Jul21/GO_enrich_Bladder_Recurren_Comp.xlsx")
#write_xlsx(res1,"../Res_Jul21/GO_enrich_Bladder_Papillary_Comp.xlsx")
write_xlsx(res1,"../Res_Jul21/GO_enrich_Bladder_Papillary_Grades_Comp.xlsx")




