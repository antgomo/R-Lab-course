library(minfi)
library(DMRcate)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(factoextra)
library(FactoMineR)
library(limma)
library(writexl)
library(readxl)
library(missMethyl)
library(ggpubr)


#- Samples to remove before starting the analyses: only in bladder cancer. Samples to remove are highlighted in red font and annotated as "recurrence".

#- 1st level of analysis: the pancancer strategy - all tumors vs normals
#- 2nd level of analysis: each tumor type (prostate, bladder or kidney) vs its corresponding normal (i.e. prostate cancer vs normal prostate; bladder cancer vs normal bladder; kidney cancer vs normal kidney)
#- 3rd level of analysis: tumor specific signatures (prostate cancer vs bladder cancer vs kidney cancer) - only tumors here

#- 4th level of analysis: clinicopathological investigations within each cancer type
#4.1. Prostate cancer: Gleason 6 versus higher Gleasons; stage 2a/b versus stage 3a/b; no recurrence versus recurrence/no remission. Possible covariates: age; PSA level at diagnosis (continuous variables).
#4.2. Bladder cancer: MIBC versus NMIBC; recurrence versus non-recurrence; papillary low grade versus papillary high grade (this is for NMIBC only); papillary low grade versus papillary high grade versus invasive high grade (here a comparison of 3 groups, if possible). Possible covariates: age, gender.
#4.3. Kidney cancer: Stage I versus others (II, III, IV); Possible covariates: age, gender.

setwd("/Users/agomez/IDIBELL/Joao/alltums/")

idat.folder <- "/Users/agomez/IDIBELL/Joao/alltums/" ###change for our path dir

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
Pheno_Prost<-read.csv("../Pheno_Prostate.csv")
#####Prostate
###select all prostate

betas2test<-betas[,colnames(betas) %in% Pheno_Prost$ID]
###same order
#Pheno_Prost<-Pheno_Prost[match(colnames(betas2test),Pheno_Prost$ID),]
###merge with general DF, change ID
colnames(targets)[1]<-"ID"
Pheno_Prost<-merge(Pheno_Prost,targets[,c(1,6)],by="ID")
Pheno_Prost<-Pheno_Prost[match(colnames(betas2test),Pheno_Prost$ID),]

##DMP finding
###sex and age as covariates, in this case only sex

####Variable of interest as factor
Pheno_Prost$Type<-as.factor(Pheno_Prost$Type)
##Normal always as referene, hence Cancer is UP
Pheno_Prost$Type<-relevel(Pheno_Prost$Type,ref = "Prostate Normal")
design <- model.matrix(~as.factor(Pheno_Prost$Type)+Pheno_Prost$Age)
colnames(design)[2]<-"Comp" ##Patient UP

##############

fit <- lmFit(betas2test, design)

# fit the contrasts
fit2 <- eBayes(fit)
results<-topTable(fit2,coef="Comp",n= dim(fit2)[1],adjust.method = "fdr", sort.by="P")

dmp<-results[results$adj.P.Val<.05,]

dmp <- merge(dmp,ann450k[,c("UCSC_RefGene_Name","Relation_to_Island","UCSC_RefGene_Group")], by="row.names")
dmp<-as.data.frame(dmp)


write_xlsx(dmp,"../Res_Dec21/Prostate/CancNormal/Annotated_dmps_Prostate_Canc_Normal.xlsx")

####only prom + islands

proms<-dmp[grep("TSS1500|TSS200|5'UTR|1stExon",dmp$UCSC_RefGene_Group),]
proms.island<-proms[proms$Relation_to_Island=="Island",]
write_xlsx(proms.island,"../Res_Dec21/Prostate/CancNormal/Annotated_dmps_proms_islands_Prostate_Cancer_Normal.xlsx")

######

###First, get the mean for both groups

##get values

#my.res<-betas2test[rownames(betas2test) %in% dmp$Row.names,]
my.res<-betas2test[rownames(betas2test) %in% proms.island$Row.names,]

my.res<-colMeans(my.res)

my.res<-as.data.frame(t(my.res))
my.res<-t(my.res)
my.res<-as.data.frame(my.res)
my.res$ID<-rownames(my.res)

my.res<-merge(my.res,targets[,c("ID","Type")],by="ID")
###Boxplot

p <- ggboxplot(my.res, x = "Type", y = "V1",
               color = "Type",
               add = "jitter",
               legend="")


p<-p+ ylab("Beta value") + xlab("") + ggtitle("Prostate Normal vs Cancer Proms+Island")+ 
  geom_hline(yintercept = median(my.res[my.res$Type=="Normal","my.res"]), linetype = 2,color="lightblue") # Add horizontal line at bas
# Specify the comparisons you want
my_comparisons <- list( c("Prostate Cancer", "Prostate Normal"))

p<- p + stat_compare_means(method="wilcox",comparisons = my_comparisons,
                           #
                           # label = "p.signif",###only signifficance symbol stars
                           method.args = list(var.equal = TRUE),##make t-test simila to ANOVA 2 groups,
                           symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                               symbols = c("****", "***", "**", "*", "ns"))
)


####


library(M3C)
library(NMF) # loading for aheatmap plotting function
library(ggsci) # more cool colours
library(heatmap.2x)
library(RColorBrewer)

annot<-Pheno_Prost[,c(1,12)]
colnames(annot)[1]<-"ID"
rownames(annot)<-annot$ID


my.data<-as.data.frame(betas2test[rownames(betas2test) %in% proms.island$Row.names ,])
#my.data<-as.data.frame(betas2test[rownames(betas2test) %in% dmp$Row.names ,])



samples<-ifelse(grepl("Cancer",Pheno_Prost$Type),"red","blue1")

#####M3C

data <- t(scale(t(my.data))) # z-score normalise each row (feature)
data <- apply(data, 2, function(x) ifelse(x > 4, 4, x)) # compress data within range
data <- apply(data, 2, function(x) ifelse(x < -4, -4, x)) # compress data within range


spcol<-samples
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)

heatmap.2x(as.matrix(data[sample(1:5000),]), col=rev(cols), Colv = T,Rowv=TRUE, scale="none", ColSideColors=spcol,
           trace="none",
           dendrogram="column", 
          
           cexRow=1, cexCol=.7,
           main="Prostate Cancer vs Normal Proms + Islands",
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

###Volcano

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
  ggtitle("Prostate Cancer vs Normal") +
  theme(legend.position="none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Delta Beta")
#####Volcano proms + islands
###sbset results to proms and islands

proms_islands_annot<- ann450k %>% 
                      filter(Relation_to_Island=="Island") %>%
                      filter(grepl("TSS1500|TSS200|5'UTR|1stExon",UCSC_RefGene_Group))

###subset now
results2 <-results %>%
           filter(rownames(.) %in% proms_islands_annot$Name)

results2["group"] <- "NotSignificant"
results2[which(rownames(results2) %in% proms.island$Row.names & results2['logFC'] > 0.2),"group"] <- "Hypermethylated"
results2[which(rownames(results2)  %in% proms.island$Row.names & results2['logFC'] < -0.2),"group"] <- "Hypomethylated"

##Volcano
vp<-ggplot(data=results2,aes(x=logFC,y=-log10(adj.P.Val),color=group))+
  geom_point(alpha=.6, size=1.2)+
  scale_colour_manual(values = c("green", "red", "lightgrey"))+
  geom_vline(xintercept=-0.2,linetype="dashed",col="black")+
  geom_vline(xintercept=0.2,linetype="dashed",col="black")+
  geom_hline(yintercept=1.3,linetype="dashed",col="black")+
  xlim(-1,1)+
  #geom_vline(xintercept=c(-0.5,0.5),linetype="dashed",col="blue")+
  ggtitle("Prostate Cancer vs Normal Proms + Islands") +
  theme(legend.position="none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Delta Beta")


##########Boxplot


###First, get the mean for both groups

##get values

my.res<-betas2test[rownames(betas2test) %in% dmp$Row.names,]
my.res<-betas2test[rownames(betas2test) %in% proms.island$Row.names,]

my.res<-colMeans(my.res)

my.res<-as.data.frame(my.res)
my.res$ID<-rownames(my.res)

my.res<-merge(my.res,targets[,c("ID","Type")],by="ID")
###Boxplot

p <- ggboxplot(my.res, x = "Type", y = "my.res",
               color = "Type",
               add = "jitter",
               legend="")


p<-p+ ylab("Beta value") + xlab("") + ggtitle("Prostate Normal vs Cancer Gleason6 Proms+Island")+ 
  geom_hline(yintercept = median(my.res[my.res$Type=="Normal","my.res"]), linetype = 2,color="lightblue") # Add horizontal line at bas
# Specify the comparisons you want
my_comparisons <- list( c("Prostate Cancer", "Prostate Normal"))

p<- p + stat_compare_means(method="wilcox",comparisons = my_comparisons,
                           #
                           # label = "p.signif",###only signifficance symbol stars
                           method.args = list(var.equal = TRUE),##make t-test simila to ANOVA 2 groups,
                           symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                               symbols = c("****", "***", "**", "*", "ns"))
)




#####GO enrichment analysis

#####GO
library(methylGSA)
library(missMethyl)

hyper<-dmp[dmp$logFC>0,]
hypo<-dmp[dmp$logFC<0,]


hyper<-proms.island[proms.island$logFC>0,]
hypo<-proms.island[proms.island$logFC<0,]

to.test<-hypo

enrichment_GO <- gometh(as.character(to.test$Row.names),all.cpg = rownames(betas),collection = "GO", 
                        array.type = "450K",plot.bias = F,prior.prob = F,equiv.cpg = T) #prior.prob = T yes, but change
enrichment_GO<-enrichment_GO[enrichment_GO$FDR<0.05,]


write_xlsx(enrichment_GO,"../Res_Dec21/Prostate/Gleason/GO_enrich_Hypo_Promislands_Gleason.xlsx")
#write_xlsx(enrichment_GO,"../Res_Dec21/Prostate/Gleason/GO_enrich_Hypo_ALL_Gleason.xlsx")

####

results$SIG<-ifelse(rownames(results) %in% dmp$Row.names,"YES","NO")###ALL
results$SIG<-ifelse(rownames(results) %in% proms.island$Row.names,"YES","NO")###Proms islands

###
######Get delta beta change

# find genes probes with significantly different methylation statuses in 
# low- and high-grade bladder cancer

Normal <- targets$Sample_Name[targets$Class == "Normal"]
Cancer <- targets$Sample_Name[targets$Class == "Cancer"]
##get dbet matrix using betas

dbet <- data.frame (Normal = rowMeans(betas[, Normal]),
                    Cancer = rowMeans(betas[, Cancer]))
dbet$delta <- abs(dbet$Normal - dbet$Cancer)

###delta to results

results<-merge(results,dbet,by="row.names")


###establish categoty 
#set the decreased and increased like you did:
results["group"] <- "NotSignificant"
results[which(results['adj.P.Val'] < 0.05 & results['logFC'] > 0.2),"group"] <- "Hypermethylated"
results[which(results['adj.P.Val'] < 0.05 & results['logFC'] < -0.2),"group"] <- "Hypomethylated"
####for proms islands
results["group"] <- "NotSignificant"
results[which(results['adj.P.Val'] < 0.05 & results['SIG']=="YES" & results['logFC'] > 0.2),"group"] <- "Hypermethylated"
results[which(results['adj.P.Val'] < 0.05 & results['SIG']=="YES" & results['logFC'] < -0.2),"group"] <- "Hypomethylated"


##Volcano
vp<-ggplot(data=results,aes(x=logFC,y=-log10(adj.P.Val),color=group))+
  geom_point(alpha=.6, size=1.2)+
  scale_colour_manual(values = c("green", "red", "lightgrey"))+
  geom_vline(xintercept=-0.2,linetype="dashed",col="black")+
  geom_vline(xintercept=0.2,linetype="dashed",col="black")+
  geom_hline(yintercept=1.3,linetype="dashed",col="black")+
  xlim(-1,1)+
  #geom_vline(xintercept=c(-0.5,0.5),linetype="dashed",col="blue")+
  ggtitle("Prostate Cancer vs Normal ") +
  theme(legend.position="none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Delta Beta")

##add names
#vp<-vp+geom_text_repel(data=dplyr::filter(res, logFC < -.1 & rawP <.05 | logFC > .1 & rawP <.05 ), aes(label=ID),point.padding = 2)

tiff(paste0("VP_Prostate_","ALL",".tiff"),units="in", width=5, height=5 ,res=300)
vp
dev.off()


####4.1. Prostate cancer: Gleason 6 versus higher Gleasons; 
#stage 2a/b versus stage 3a/b; 
#no recurrence versus recurrence/no remission. Possible covariates: age; PSA level at diagnosis (continuous variables).
##select only prostate cancer
###
canc<-Pheno_Prost[Pheno_Prost$Type=="Prostate Cancer",]

betascanc<-betas[,colnames(betas) %in% canc$ID]

###same order
canc<-canc[match(colnames(betascanc),canc$ID),]

##DMP finding
###sex and age as covariates, in this case only sex

####Variable of interest as factor
##Gleason
canc$Type<-as.factor(ifelse(grepl("6",canc$Gleason),"Gleas6","GleasHigh"))
canc$Type<-relevel(canc$Type,ref = "Gleas6")

##Gleason 3+4 4+3 
canc$Type<-as.factor(ifelse(canc$Gleason=="3+4=7","3+4",ifelse(canc$Gleason=="4+3=7","4+3","6")))
####remvoe 6
canc<-canc[!(canc$Type=="6"),]
canc$Type<-droplevels(canc$Type)

canc$Type<-relevel(canc$Type,ref = "4+3")
##TumorStage
canc$Comp2<-as.factor(ifelse(grepl("2",canc$Tstage),"Tst2","Tst3"))
##Remission
#canc$Comp3<-as.factor(ifelse(grepl("No remission",canc$Recurrence.yes.no),"No_Rem",ifelse(grepl("Yes",canc$Recurrence.yes.no),"Recurren",NA)))

####Run depending on the comparison you want to test
design <- model.matrix(~as.factor(canc$Type)+canc$Age+as.numeric(canc$PSA.ao.dx))##Gleason
#design <- model.matrix(~as.factor(canc$Comp2)+canc$Age+as.numeric(canc$PSA.ao.dx))##TumorStage
###In recurrence vs no remision, remove NAs
#canc2<-canc[!(is.na(canc$Comp3)),]
#design <- model.matrix(~as.factor(canc2$Comp3)+canc2$Age+as.numeric(canc2$PSA.ao.dx))##TumorStage


colnames(design)[2]<-"Comp" ##Patient UP

##############
betascanc<-betascanc[,colnames(betascanc) %in% canc$ID]
fit <- lmFit(betascanc, design)##Gleason and Tumor Stage

#fit <- lmFit(betascanc[,colnames(betascanc) %in% canc2$ID] , design)###Recurrence vs no remission

# fit the contrasts
fit2 <- eBayes(fit)
results<-topTable(fit2,coef="Comp",n= dim(fit2)[1],adjust.method = "fdr", sort.by="P")

dmp<-results[results$P.Value<.01,]###use non adjusted


dmp <- merge(dmp,ann450k[,c("UCSC_RefGene_Name","Relation_to_Island","UCSC_RefGene_Group")], by="row.names")
dmp<-as.data.frame(dmp)

####only prom + islands

proms<-dmp[grep("TSS1500|TSS200|5'UTR|1stExon",dmp$UCSC_RefGene_Group),]
proms.island<-proms[proms$Relation_to_Island=="Island",]
#write_xlsx(proms.island,"../Res_Dec21/Prostate/Gleason/Annotated_dmps_proms_islands_Gleason3_4vsGleason4_3.xlsx")
write_xlsx(proms.island,"../Res_Dec21/Prostate/Gleason/Annotated_dmps_proms_islands_Prostate_TumSt2vsTumSt3.xlsx")


#write_xlsx(dmp,"../Res_Dec21/Prostate/Gleason/Annotated_dmps_Prostate_Gleason3_4vsGleason4_3.xlsx")
write_xlsx(dmp,"../Res_Dec21/Prostate/Gleason/Annotated_dmps_Prostate_TumSt2vsTumSt3.xlsx")
write_xlsx(proms.island,"../Res_Dec21/Prostate/TumorStage/Annotated_dmps_Prostate_PromsIslands_TumSt2vsTumSt3.xlsx")

#write_xlsx(dmp,"../Res_Jul21/Annotated_dmps_Prostate_TumSt2vsTumSt3.xlsx")
#write_xlsx(dmp,"../Res_Jul21/Annotated_dmps_Prostate_TRecvsNoRem.xlsx")


annot<-canc[,c(1,12)]##Gleason
#annot<-canc[,c(1,13)]##TStage
#annot<-canc2[,c(1,14)]##RecvsNoRem

colnames(annot)[1]<-"ID"
rownames(annot)<-annot$ID


my.data<-as.data.frame(betascanc[rownames(betascanc) %in% dmp$Row.names ,])
my.data<-as.data.frame(betascanc[rownames(betascanc) %in% proms.island$Row.names ,])

#my.data<-as.data.frame(betascanc[rownames(betascanc) %in% rownames(dmp) ,colnames(betascanc) %in% canc2$ID])##Rec vs Non Rem


#samples<-ifelse(grepl("High",canc$Type),"red","blue1")###Gleason 6
samples<-ifelse(canc$Type=="3+4","red","blue1")###Gleason 3+4 4+3

#samples<-ifelse(grepl("3",canc$Tstage),"red","blue1")###TumStage
#samples<-ifelse(grepl("Recurren",canc2$Comp3),"red","blue1")###TumStage

#####M3C

data <- t(scale(t(my.data))) # z-score normalise each row (feature)
data <- apply(data, 2, function(x) ifelse(x > 4, 4, x)) # compress data within range
data <- apply(data, 2, function(x) ifelse(x < -4, -4, x)) # compress data within range


spcol<-samples
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)

heatmap.2x(as.matrix(data), col=rev(cols), Colv = T,Rowv=TRUE, scale="none", ColSideColors=spcol,
           trace="none", dendrogram="column", 
           cexRow=1, cexCol=.7,
           main="Prostate Cancer Gleason Prms + Islands ",
           # labCol=NA,
           labRow=NA, 
           density.info="none",
           hclust=function(x) hclust(x,method="complete"),
           distfun=function(x) as.dist((1-cor(t(x)))/2),srtCol=45
           
)

coords <- locator(1)
legend(coords,legend=c("3+4","4+3"),
       fill=c("red","blue1"), border=FALSE, 
       y.intersp = 0.7, 
       xpd=TRUE, inset=c(0,1), cex=1, bty='n')

#####GO enrichment analysis

#####GO
library(methylGSA)
library(missMethyl)

hyper<-dmp[dmp$logFC>0,]
hypo<-dmp[dmp$logFC<0,]


hyper<-proms.island[proms.island$logFC>0,]
hypo<-proms.island[proms.island$logFC<0,]

to.test<-hypo


enrichment_GO <- gometh(as.character(to.test$Row.names),all.cpg = rownames(betas),collection = "GO", 
                        array.type = "450K",plot.bias = F,prior.prob = F,equiv.cpg = T) 
enrichment_GO<-enrichment_GO[enrichment_GO$FDR<0.05,]

write_xlsx(enrichment_GO,"../Res_Dec21/Prostate/TumorStage/GO_enrich_Hypo_All_Tstage2.xlsx")

write_xlsx(enrichment_GO,"../Res_Dec21/Prostate/TumorStage/GO_enrich_Hypo_Promislands_Tstage2.xlsx")

####

results$SIG<-ifelse(rownames(results) %in% dmp$Row.names,"YES","NO")###ALL

###establish category 
#set the decreased and increased like you did:

########Proms + Islands

proms_islands_annot<- ann450k %>% 
  filter(Relation_to_Island=="Island") %>%
  filter(grepl("TSS1500|TSS200|5'UTR|1stExon",UCSC_RefGene_Group))

###subset now
results2 <-results %>%
  filter(rownames(.) %in% proms_islands_annot$Name)

results2["group"] <- "NotSignificant"
results2[which(rownames(results2) %in% proms.island$Row.names & results2['logFC'] > 0.2),"group"] <- "Hypermethylated"
results2[which(rownames(results2)  %in% proms.island$Row.names & results2['logFC'] < -0.2),"group"] <- "Hypomethylated"

##Volcano
vp<-ggplot(data=results2,aes(x=logFC,y=-log10(P.Value),color=group))+
  geom_point(alpha=.6, size=1.2)+
  scale_colour_manual(values = c("green", "red", "lightgrey"))+
  geom_vline(xintercept=-0.2,linetype="dashed",col="black")+
  geom_vline(xintercept=0.2,linetype="dashed",col="black")+
  geom_hline(yintercept=2,linetype="dashed",col="black")+
  xlim(-1,1)+
  #geom_vline(xintercept=c(-0.5,0.5),linetype="dashed",col="blue")+
  ggtitle("Prostate Gleason 3+4 vs 4+3 Proms + Island") +
  theme(legend.position="none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Delta Beta")


###ALL
results["group"] <- "NotSignificant"
results[which(results['P.Value']  < 0.01 & results['logFC'] > 0.2),"group"] <- "Hypermethylated"
results[which(results['P.Value'] < 0.01 & results['logFC'] < -0.2),"group"] <- "Hypomethylated"


##Volcano
vp<-ggplot(data=results,aes(x=logFC,y=-log10(P.Value),color=group))+
  geom_point(alpha=.6, size=1.2)+
  scale_colour_manual(values = c("green", "red", "lightgrey"))+
  geom_vline(xintercept=-0.2,linetype="dashed",col="black")+
  geom_vline(xintercept=0.2,linetype="dashed",col="black")+
  geom_hline(yintercept=2,linetype="dashed",col="black")+
  xlim(-1,1)+
  #geom_vline(xintercept=c(-0.5,0.5),linetype="dashed",col="blue")+
  ggtitle("Prostate Gleason 3+4 vs 4+3 ALL") +
  theme(legend.position="none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Delta Beta")

##add names
#vp<-vp+geom_text_repel(data=dplyr::filter(res, logFC < -.1 & rawP <.05 | logFC > .1 & rawP <.05 ), aes(label=ID),point.padding = 2)

tiff(paste0("VP_Prostate_Gleason_","ALL",".tiff"),units="in", width=5, height=5 ,res=300)
vp
dev.off()


###Boxplot

##get values

my.res<-betascanc[rownames(betascanc) %in% dmp$Row.names,]
my.res<-betascanc[rownames(betascanc) %in% proms.island$Row.names,]

my.res<-colMeans(my.res)

my.res<-as.data.frame(my.res)
my.res$ID<-rownames(my.res)

my.res<-merge(my.res,canc[,c("ID","Type")],by="ID")

my.res$Type<-relevel(my.res$Type,ref="3+4")
###Boxplot

p <- ggboxplot(my.res, x = "Type", y = "my.res",
               color = "Type",
               add = "jitter",
               legend="")


p<-p+ ylab("Beta value") + xlab("") + ggtitle("Prostate Gleason 3+4 vs 4+3")+ 
  geom_hline(yintercept = median(my.res[my.res$Type=="Normal","my.res"]), linetype = 2,color="lightblue") # Add horizontal line at bas
# Specify the comparisons you want
my_comparisons <- list( c("3+4", "4+3"))

p<- p + stat_compare_means(method="wilcox",comparisons = my_comparisons,
                           #
                           # label = "p.signif",###only signifficance symbol stars
                           method.args = list(var.equal = TRUE),##make t-test simila to ANOVA 2 groups,
                           symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                               symbols = c("****", "***", "**", "*", "ns"))
)

