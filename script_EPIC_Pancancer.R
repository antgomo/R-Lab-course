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
library(dplyr)

#- Samples to remove before starting the analyses: only in bladder cancer. Samples to remove are highlighted in red font and annotated as "recurrence".

#- 1st level of analysis: the pancancer strategy - all tumors vs normals
#- 2nd level of analysis: each tumor type (prostate, bladder or kidney) vs its corresponding normal (i.e. prostate cancer vs normal prostate; bladder cancer vs normal bladder; kidney cancer vs normal kidney)
#- 3rd level of analysis: tumor specific signatures (prostate cancer vs bladder cancer vs kidney cancer) - only tumors here

#- 4th level of analysis: clinicopathological investigations within each cancer type
#4.1. Prostate cancer: Gleason 6 versus higher Gleasons; stage 2a/b versus stage 3a/b; no recurrence versus recurrence/no remission. Possible covariates: age; PSA level at diagnosis (continuous variables).
#4.2. Bladder cancer: MIBC versus NMIBC; recurrence versus non-recurrence; papillary low grade versus papillary high grade (this is for NMIBC only); papillary low grade versus papillary high grade versus invasive high grade (here a comparison of 3 groups, if possible). Possible covariates: age, gender.
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

###
load("../Res_Jul21/Betas_Norm_ALL.RData")


betas<-betas[,match(targets$Sample_Name,colnames(betas))]
betas<-betas[,!(colnames(betas) %in% torem)]


####PCA
####
targets<-targets[targets$Sample_Name %in% colnames(betas),]
betas<-betas[,colnames(betas) %in% targets$Sample_Name ]

pca <- FactoMineR::PCA(t(betas), scale.unit = T, graph = F, ncp = 20)
factoextra::fviz_pca_ind(pca, axes = c(1,2), habillage=as.factor(targets$Type), repel = F)+ coord_fixed()

####Fig 1:PCA promoter + islands

###subset betas
betas.pca <- merge(betas,ann450k[,c("UCSC_RefGene_Name","Relation_to_Island","UCSC_RefGene_Group")], by="row.names")

betas.pca<-betas.pca[grep("TSS1500|TSS200|5'UTR|1stExon",betas.pca$UCSC_RefGene_Group),]
betas.pca<-betas.pca[betas.pca$Relation_to_Island=="Island",]
##avoid annotation
rownames(betas.pca)<-betas.pca$Row.names

betas.pca<-betas.pca %>%
      select(-c("UCSC_RefGene_Name","Relation_to_Island","UCSC_RefGene_Group","Row.names"))

pca <- FactoMineR::PCA(t(betas.pca), scale.unit = T, graph = F, ncp = 20)
factoextra::fviz_pca_ind(pca, axes = c(1,2), habillage=as.factor(targets$Type), repel = F)+ coord_fixed()

####Load Specific Phenos
Pheno_Prost<-read.csv("../Pheno_Prostate.csv")
Pheno_Blad<-read.csv("../Pheno_Bladder.csv")
Pheno_Kid<-read.csv("../Pheno_Kidney.csv")
#####Fist analysis, Pan tumoral
##No info related to age and sex in NK, hence no way to use linear model with covariate


##DMP finding
targets$Class<-ifelse(grepl("Normal",targets$Type),"Normal","Cancer")
#design <- model.matrix(~targets$Pheno_G+targets$Batch)
targets$Class<-relevel(as.factor(targets$Class), ref = "Normal")
design <- model.matrix(~targets$Class)

colnames(design)[2]<-"Comp" ##Patient UP

##############

fit <- lmFit(betas, design)

# fit the contrasts
fit2 <- eBayes(fit)
results<-topTable(fit2,coef="Comp",n= dim(fit2)[1],adjust.method = "fdr", sort.by="P")

dmp<-results[results$adj.P.Val<.05,]
dmp <- merge(dmp,ann450k[,c("UCSC_RefGene_Name","Relation_to_Island","UCSC_RefGene_Group")], by="row.names")
dmp<-as.data.frame(dmp)

write_xlsx(dmp,"../Res_Dec21/PanCancer/Annotated_dmps_PanCancer.xlsx")

####only prom + islands

proms<-dmp[grep("TSS1500|TSS200|5'UTR|1stExon",dmp$UCSC_RefGene_Group),]
#proms.island<-merge(proms.island,betas,by="row.names")
proms.island<-proms[proms$Relation_to_Island=="Island",]
#rownames(proms.island)<-proms.island$Row.names
#proms.island<-proms.island[,-1]
write_xlsx(proms.island,"../Res_Dec21/PanCancer/Annotated_dmps_proms_islands_PanCancer.xlsx")

library(M3C)
library(NMF) # loading for aheatmap plotting function
library(ggsci) # more cool colours
library(heatmap.2x)
library(RColorBrewer)

annot<-targets[,c(1,11)]
colnames(annot)[1]<-"ID"
rownames(annot)<-annot$ID


my.data<-as.data.frame(betas[rownames(betas) %in% dmp$Row.names ,])
my.data<-as.data.frame(betas[rownames(betas) %in% proms.island$Row.names ,])

library(heatmap.2x)## approach

class<-ifelse(targets$Class=="Cancer","red","blue1")
###add antoher one
targets$Type<-as.factor(targets$Type)
subtype<-ggsci::pal_jama()(length(levels(targets$Type)))[targets$Type]#5 clusters

spcol <- rbind(subtype, class)

#####M3C
#Add dendrogram to all heatmaps
#Only for PanCancer heatmaps - 
  #Add colour bars indicating tissue type (normal bladder, bladder cancer, normal prostate, prostate cancer, normal kidney, kidney cancer) 
  #and maintain the one that we already have (normal vs cancer)




data <- t(scale(t(my.data))) # z-score normalise each row (feature)
data <- apply(data, 2, function(x) ifelse(x > 4, 4, x)) # compress data within range
data <- apply(data, 2, function(x) ifelse(x < -4, -4, x)) # compress data within range


cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
##as.matrix(data[sample(1:5000),])
heatmap.2x(as.matrix(data[sample(1:9000),]), col=rev(cols), Colv = T,Rowv=TRUE, scale="none", ColSideColors=spcol,
           trace="none", 
           dendrogram="column", 
           cexRow=1, cexCol=.7,
           main="PanCancer Cancer vs Normal ",
          # labCol=NA,
           labRow=NA, 
           density.info="none",
           hclust=function(x) hclust(x,method="complete"),
           distfun=function(x) as.dist((1-cor(t(x)))/2),srtCol=45
           
)

coords <- locator(1)
#legend("bottomleft",c("Cancer","HC"),pch=20:20,col=c("red","blue1"))
legend(coords,legend=c(levels(targets$Type),levels(targets$Class)),
       fill=c(ggsci::pal_jama()(length(levels(targets$Type))),"blue","red"), border=FALSE, 
       y.intersp = 0.7, 
       xpd=TRUE, inset=c(0,1), cex=0.7, bty='n')



#####GO enrichment analysis

#####GO
library(methylGSA)
library(missMethyl)

hyper<-dmp[dmp$logFC>0,]
hypo<-dmp[dmp$logFC<0,]

hyper<-proms.island[proms.island$logFC>0,]
hypo<-proms.island[proms.island$logFC<0,]

to.test<-hypo

pv.vals<-to.test$P.Value

names(pv.vals)<-rownames(to.test)


#res1 <- methylglm(cpg.pval = pv.vals, minsize = 10, maxsize = 1000, GS.type = "GO",array.type = "450K")

#res1<-res1[res1$padj<.05,]

#write_xlsx(res1,"../Res_Dec21/PanCancer/GO_enrich_Hyper_Cancer_Norm.xlsx")

enrichment_GO <- gometh(as.character(to.test$Row.names),all.cpg = rownames(betas),collection = "GO", 
                        array.type = "450K",plot.bias = T,prior.prob = T,equiv.cpg = T) 
enrichment_GO<-enrichment_GO[enrichment_GO$FDR<0.05,]


write_xlsx(enrichment_GO,"../Res_Dec21/PanCancer/GO_enrich_Hypo_Cancer_Norm.xlsx")
write_xlsx(enrichment_GO,"../Res_Dec21/PanCancer/GO_enrich_Hyper_Prom_islands.xlsx")

####

results$SIG<-ifelse(rownames(results) %in% dmp$Row.names,"YES","NO")
results$SIG<-ifelse(rownames(results) %in% proms.island$Row.names,"YES","NO")
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
  ggtitle("PanCancer Cancer vs Normal Proms+ islands") +
  theme(legend.position="none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Delta Beta")

##add names
#vp<-vp+geom_text_repel(data=dplyr::filter(res, logFC < -.1 & rawP <.05 | logFC > .1 & rawP <.05 ), aes(label=ID),point.padding = 2)

tiff(paste0("VP_Prostate_Gleason","ALL",".tiff"),units="in", width=5, height=5 ,res=300)
vp
dev.off()




###First, get the mean for both groups

##get values

my.res<-betas[rownames(betas) %in% dmp$Row.names,]
my.res<-betas[rownames(betas) %in% proms.island$Row.names,]

my.res<-colMeans(my.res)

my.res<-as.data.frame(my.res)
my.res$Sample_Name<-rownames(my.res)

my.res<-merge(my.res,targets[,c("Sample_Name","Tissue","Class")],by="Sample_Name")
###Boxplot


p<-ggboxplot(
  my.res, x = "Class", y = "my.res", 
  facet.by = "Tissue", add = "jitter",color="Class",
  legend=""
) +
  stat_compare_means(
    comparisons = list(c("Normal", "Cancer")), 
    label = "p.signif"
  )


p<-p+ ylab("Beta value") + xlab("") + ggtitle("PanCancer Normal vs Cancer Proms +Islands")+ 
  geom_hline(yintercept = median(my.res[my.res$Class=="Normal","my.res"]), linetype = 2,color="lightblue") # Add horizontal line at bas
# Specify the comparisons you want

p<- p + stat_compare_means(method="wilcox",comparisons = my_comparisons,
                           #
                           # label = "p.signif",###only signifficance symbol stars
                           method.args = list(var.equal = TRUE),##make t-test simila to ANOVA 2 groups,
                           symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                               symbols = c("****", "***", "**", "*", "ns"))
)



###plot
##

##########comparison of three
#select only cancer

###
canc<-targets[grep("Cancer",targets$Type),]

betascanc<-betas[,colnames(betas) %in% canc$Sample_Name]

my.data<-as.data.frame(t(betascanc))

###get extra info

Pheno_Prost<-read.csv("../Pheno_Prostate.csv")
Pheno_Blad<-read.csv("../Pheno_Bladder.csv",sep=";")
Pheno_Kid<-read.csv("../Pheno_Kidney.csv")

extra<-rbind(Pheno_Prost[,1:3],Pheno_Blad[,1:3],Pheno_Kid[,1:3])
canc$ID<-canc$Sample_Name
canc<-merge(canc[,c(10,6)],extra,by="ID")
rownames(canc)<-canc$ID

my.data<-merge(my.data,canc[,c(2,3,4)],by="row.names")
my.data<-my.data[,-1]##remove ID
final<-dim(my.data)[2]-3

##store the results

my.res<-matrix(nrow = final,ncol = 3)

for (i in 1:final) {#take into account last columns are to use
  y.df<-my.data[,c(i,432005,432006,432007)]
  result<-aov(y.df[,1]~as.factor(y.df$Type)+y.df$Age+as.factor(y.df$Gender))
  
  res<-summary(result)
  pval<-res[[1]][["Pr(>F)"]][1]
  pval<-signif(pval, digits=3)
  my.res[i,1]<-colnames(my.data)[i]
  my.res[i,2]<-pval
  my.res[i,3]<-p.adjust(pval, method = "fdr", n = final)
  
}

my.res<-as.data.frame(my.res)
my.res<-my.res[my.res$V3<.05,]

rownames(my.res)<-my.res$V1
dmp <- merge(my.res,ann450k[,c("UCSC_RefGene_Name","Relation_to_Island","UCSC_RefGene_Group")], by="row.names")
dmp<-as.data.frame(dmp)
write_xlsx(dmp,"../Res_Dec21/PanCancer/Annotated_dmps_ALL_Tum_specific signatures.xlsx")

dmp<-read_xlsx("../Res_Dec21/PanCancer/Annotated_dmps_ALL_Tum_specific signatures.xlsx")

annot<-canc[,c(1,2)]##Three comp

colnames(annot)[1]<-"ID"
rownames(annot)<-annot$ID
###depiciting only proms and islands
proms<-dmp[grep("TSS1500|TSS200|5'UTR|1stExon",dmp$UCSC_RefGene_Group),]
proms.island<-proms[proms$Relation_to_Island=="Island",]
write_xlsx(proms.island,"../Res_Dec21/PanCancer/Annotated_dmps_ALL_Tum_specific signatures_proms_islamds.xlsx")


my.data<-as.data.frame(betascanc[rownames(betascanc) %in% proms.island$Row.names ,])
canc<-canc[match(colnames(my.data),canc$ID),]

samples<-ifelse(grepl("Bladder",canc$Type),"red",ifelse(grepl("Kidney",canc$Type),"green","blue1"))###

#####M3C

data <- t(scale(t(my.data))) # z-score normalise each row (feature)
data <- apply(data, 2, function(x) ifelse(x > 4, 4, x)) # compress data within range
data <- apply(data, 2, function(x) ifelse(x < -4, -4, x)) # compress data within range


spcol<-samples
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
##data[sample(1:5000),]
heatmap.2x(as.matrix(data), col=rev(cols), Colv = T,Rowv=TRUE, scale="none", ColSideColors=spcol,
           trace="none", dendrogram="none", 
           cexRow=1, cexCol=.7,
           main="Tumor signature levels Proms + Islands",
           # labCol=NA,
           labRow=NA, 
           density.info="none",
           hclust=function(x) hclust(x,method="complete"),
           distfun=function(x) as.dist((1-cor(t(x)))/2),srtCol=45
           
)

####BOXplot

###First, get the mean for both groups

##get values

my.res<-betascanc[rownames(betascanc) %in% dmp$Row.names,]
my.res<-betascanc[rownames(betascanc) %in% proms.island$Row.names,]

my.res<-colMeans(my.res)

my.res<-as.data.frame(my.res)
my.res$Sample_Name<-rownames(my.res)

my.res<-merge(my.res,targets[,c("Sample_Name","Type")],by="Sample_Name")
###Boxplot

p <- ggboxplot(my.res, x = "Type", y = "my.res",
               color = "Type",
               add = "jitter",
               legend="")


p<-p+ ylab("Beta value") + xlab("") + ggtitle("PanCancer Tumor Signature Proms + Islands")
# Specify the comparisons you want
my_comparisons <- list( c("Bladder Cancer", "Kidney Cancer"),c("Kidney Cancer", "Prostate Cancer"),c("Prostate Cancer", "Bladder Cancer"))

p<- p + stat_compare_means(method="wilcox",comparisons = my_comparisons,
                           #
                           # label = "p.signif",###only signifficance symbol stars
                           method.args = list(var.equal = TRUE),##make t-test simila to ANOVA 2 groups,
                           symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                               symbols = c("****", "***", "**", "*", "ns"))
)





##GO


enrichment_GO <- gometh(proms.island$Row.names,all.cpg = rownames(betas),collection = "GO", 
                        array.type = "450K",plot.bias = T,prior.prob = T,equiv.cpg = T) 
enrichment_GO<-enrichment_GO[enrichment_GO$FDR<0.05,]


write_xlsx(enrichment_GO,"../Res_Dec21/PanCancer/GO_enrich_Signatures_Proms_islands.xlsx")


