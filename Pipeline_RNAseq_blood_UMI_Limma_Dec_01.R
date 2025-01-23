library(tximport)
library(sva)
library(RUVSeq)
library(limma)

#((G_w0+NR_w0)-(G_w12+NR_w12)){ETN/CRT} -  ((G_w0+NR_w0)-(G_w12+NR_w12)){ADA/GOL}

#Certolizumab + Etanercept
#vs
#Adalimumab + Golimumab
setwd("/Users/agomez/IMIDOMICS/P6/")
#https://support.bioconductor.org/p/73872/
#https://support.bioconductor.org/p/92067/
sampleTable<-read.csv("Celgene_Project6_RNAseq_Mastersheet.csv")
rownames(sampleTable)<-sampleTable$GSL.ID

codes<-read.csv("codes.csv")
codes<-codes[,-1]
codes<-unique(codes,MARGIN=1)

sampleTable<-merge(sampleTable,codes,by="donor.code")
###pheno

pheno<-read.csv("Pheno.csv")
sampleTable<-merge(sampleTable,pheno,by="Code",all.x=T)
rownames(sampleTable)<-sampleTable$GSL.ID

library(edgeR)

####PAIRWISE ANALYSIS
###remove samples 5632-RMM-0061, 5632-RMM-0062 -> IB-30203 Non age
sampleTable<-sampleTable[!(is.na(sampleTable$Age)),]

all_cts <- read.delim("C_last.txt",head=T)
colnames(all_cts)<-gsub("X","",colnames(all_cts))
colnames(all_cts)<-gsub(".bam","",colnames(all_cts))
colnames(all_cts)<-gsub("\\.","-",colnames(all_cts))
rownames(all_cts)<-all_cts$Geneid
all_cts<-all_cts[,-1]
all_cts<-all_cts[,-1]

colnames(all_cts)<-substr(colnames(all_cts),1,13)



sampleTable<-sampleTable[rownames(sampleTable) %in% colnames(all_cts),]
sampleTable$sex<-ifelse(sampleTable$Sex=="FEMALE","F","M")

all_cts<-all_cts[,colnames(all_cts) %in% rownames(sampleTable)]
#colnames(all_cts)<-sampleTable$rna.code

###CONDITION 2 : R= R+MOD and NR
sampleTable$COND2<-ifelse(sampleTable$response=="GOOD","G","NR")
samples2test<-sampleTable

###divide by Treatment
samples2test<-samples2test[samples2test$Treat=="Adalimumab",]

samples2test<-samples2test[samples2test$Treat=="Certolizumab" | samples2test$Treat=="Etanercept",]
samples2test<-samples2test[samples2test$Treat=="Adalimumab" | samples2test$Treat=="Golimumab",]

samples2test$Condition<-paste(samples2test$COND2, samples2test$week,sep="_")
samples2test$Condition<-as.factor(samples2test$Condition)
samples2test$Condition<-droplevels(samples2test$Condition)
  ####Analysis

cts<-all_cts[,colnames(all_cts) %in% rownames(samples2test)]
cts<-cts[,match(rownames(samples2test),colnames(cts))]
#y <- DGEList(cts,group=samples2test$Condition)## for COND1
y <- DGEList(cts,group=samples2test$COND2)## for COND2

# filter out very lowly expressed tags, keeping genes that are expressed at a reasonable level in at least
#one treatment condition. Since the smallest group size is 150, we keep genes that achieve at least
#one count per million (cpm) in at least three samples:

keep <- rowSums(cpm(y)>1) >= 15#### change group depending replicates
y <- y[keep, ]

dim(y)

y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y, method="TMM")


####Limma
my.batch<-as.factor(paste("B", samples2test$seq.plate,sep=""))

#design <- model.matrix(~0+samples2test$Code+samples2test$sex + samples2test$Age+my.batch+samples2test$Condition)

##Limma
#There's no point putting in the libType factor, as the batch effect is fully absorbed into the patient-specific blocking factors. 
samples2test$Code<-gsub("-","",samples2test$Code)
samples2test$Time<-ifelse(samples2test$week=="wk0","1","2")
##simple way, due the fact that we don't have many time points
samples2test$S<-paste(samples2test$COND2, samples2test$Time, sep="_")

design <- model.matrix(~0+ samples2test$S+samples2test$sex + samples2test$Age)##with MODs in NR, Cond2
colnames(design)<-gsub("samples2test\\$","",colnames(design))

v <- voom(y, design)

corfit <- duplicateCorrelation(v, design, block =samples2test$Code)

v <- voom(y, design, block =samples2test$Code, correlation =corfit$consensus)

fit <- lmFit(v, design, block = samples2test$Code, correlation =corfit$consensus)

#Which genes respond differently over time in the Resp relative to the NR?

cont.dif <- makeContrasts(
  Dif.w0 =(SG_1-SNR_1),
  Dif.w12 =(SG_2-SNR_2),
  Dif.Time=(SG_1+SNR_1)-(SG_2+SNR_2),
  Int =(SG_1-SNR_1)-(SG_2-SNR_2),
  
  levels=design)


fit2 <- contrasts.fit(fit, cont.dif)
fit2 <- eBayes(fit2)


results <- topTable(fit2,coef=4,num=dim(fit2)[1],sort.by="P",adjust.method = "fdr") ####Comp 1, w0; comp 3 time, comp4 Int
res<-subset(results,results$P.Value<.05)
res<-subset(results,results$adj.P.Val <.01)##Time

#res<-res[!(grepl("^RP1|^RP4|^AC0|^AP0|XXbac|^RP3",rownames(res),perl = T)),]
##subset in range

write.csv(res,"Results2/RvsNR_w0_Golimumab.csv")
write.csv(res,"Results2/w0_w12_Golimumab.csv")
write.csv(res,"Results2/Interaction_Golimumab.csv")




###Heatmap

library(gplots)
library(RColorBrewer)

#Use M3C approach

library(M3C)
library(NMF) # loading for a heatmap plotting function
library(ggsci) # more cool colours
#we run the algorithm using the default settings (100x monte carlo iterations and 100x inner replications).
#We have found the results generally stable using these parameters, 
#although they may be increased or the algorithm may simply be run a few extra times to test stability. 
#Plots from the tool and an .csv file with the numerical outputs may be printed into the working directory 
#(by adding printres = TRUE). We will set the seed in this example, incase you wish to repeat our results exactly (seed = 123). 
#We will add the annotation file for streamlined downstream analyses (des = annot).

library(heatmap.2x)## approach


samples2test$ID<-as.character(samples2test$Library.ID)
annot<-samples2test[samples2test$Time==1,]
#annot<-samples2test###Time diff


annot<-annot[,c(1,37)]

my.data<-as.data.frame(y$counts[rownames(y$counts) %in% rownames(res) ,])
my.data<- cpm(my.data, prior.count=2, log=TRUE)
my.data<-my.data[,colnames(my.data) %in% rownames(annot)]
my.data<-my.data[,match(rownames(annot),colnames(my.data))]

rownames(annot)<-annot$Code
colnames(annot)[1]<-"ID"
colnames(my.data)<-rownames(annot)
#I used the coefficient of variation to remove more less variable genes

data3<-my.data


res.m3c <- M3C(data3, cores=3,des=annot,seed = 123, removeplots = F)

# get the data out of the results list (by using $ - dollar sign), use 2 clusters (see RCSI plot)

data <- res.m3c$realdataresults[[2]]$ordered_data # this is the data
annon <- res.m3c$realdataresults[[2]]$ordered_annotation # this is the annotation
ccmatrix <- res.m3c$realdataresults[[2]]$consensus_matrix # this is the consensus matrix


# normalise and scale the data
data <- t(scale(t(data3))) # z-score normalise each row (feature)
data <- apply(data, 2, function(x) ifelse(x > 4, 4, x)) # compress data within range
data <- apply(data, 2, function(x) ifelse(x < -4, -4, x)) # compress data within range




samples<-ifelse(annot$COND2=="NR","red","black")

cons<-ggsci::pal_futurama()(max(levels(annon$consensuscluster)))[as.factor(annon$consensuscluster)]#5 clusters


#spcol <- rbind(cons, samples)
spcol<-samples
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)


heatmap.2x(data, col=rev(cols), Colv = NA,Rowv=TRUE, scale="none", ColSideColors=spcol,
           trace="none", dendrogram="both", 
           cexRow=1, cexCol=1.4,
           main="G vs NR Adalimumab W0",
           labCol=NA,labRow=NA, density.info="none",
           hclust=function(x) hclust(x,method="complete"),
           distfun=function(x) as.dist((1-cor(t(x)))/2)
           
)

legend("topright",c("GOOD","NR"),pch=20:20,col=c("black","red"))




###NR+MOD
annot<-samples2test

Rw0<-rownames(annot[annot$Condition=="G_wk0",])
Rw12<-rownames(annot[annot$Condition=="G_wk12",])

NRw0<-rownames(annot[annot$Condition=="NR_wk0",])
NRw12<-rownames(annot[annot$Condition=="NR_wk12",])

annot<-annot[c(Rw0,NRw0,Rw12,NRw12),]

my.data<-as.data.frame(y$counts[rownames(y$counts) %in% rownames(res) ,])
my.data<- cpm(my.data, prior.count=2, log=TRUE)
my.data<-my.data



my.data<-my.data[,colnames(my.data) %in% rownames(annot)]
my.data<-my.data[,match(rownames(annot),colnames(my.data))]

##change column names
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)


# normalise and scale the data
data <- t(scale(t(my.data))) # z-score normalise each row (feature)
data <- apply(data, 2, function(x) ifelse(x > 4, 4, x)) # compress data within range
data <- apply(data, 2, function(x) ifelse(x < -4, -4, x)) # compress data within range



##Try to get the genes behaving in some way or another in R and NR
##better use log values not scaled
###R


#######select sample names

data<-as.data.frame(data)

data$uRW0<-rowMeans(data[,Rw0],na.rm=T)
data$uRW12<-rowMeans(data[,Rw12],na.rm=T)

data$uNRW0<-rowMeans(data[,NRw0],na.rm=T)
data$uNRW12<-rowMeans(data[,NRw12],na.rm=T)



# "uRW0"   "uNRW0"  "uRW12"  "uNRW12"


toclust<-data[,c(209,211,210,212)]##with  data NR+MOD
toclust<-data[,c(151,153,152,154)]##with  data NR+MOD
toclust<-data[,c(33,35,34,36)]##with  data NR+MOD


##now define classes

###UP in 0 and DOWN in 12

UP0D12<-toclust[(toclust[,1]>toclust[,2]) & (toclust[,3]<toclust[,4]),]
D0UP12<-toclust[(toclust[,1]<toclust[,2]) & (toclust[,3]>toclust[,4]),]
Trans.up<-toclust[(toclust[,1]>toclust[,2]) & (toclust[,3]>toclust[,4]),,]
Trans.down<-toclust[(toclust[,1]<toclust[,2]) & (toclust[,3]<toclust[,4]),,]


##########plot heatmap

library(ComplexHeatmap)
###first, check number of real clusters

library(RColorBrewer)

data$real_clust<-ifelse(rownames(data) %in% rownames(UP0D12),"1",ifelse(rownames(data) %in% rownames(D0UP12),"2",ifelse(rownames(data) %in% rownames(Trans.up),"3","4")))
cluster_real<-data$real_clust
names(cluster_real)<-rownames(data)

ht<-Heatmap(data[,c(33,35,34,36)], split = cluster_real, row_names_gp = gpar(fontsize = 6)
            ,cluster_columns = FALSE,show_row_names = F,
            
            heatmap_legend_param = list(title = "Gene changes"))
ht<-make_row_cluster(ht)
row_order<-unlist(ht@row_order_list)

png("Results/heatmap_Int_Etanercept_Certo.png",width = 1200, height = 950)

ht

 
draw(ht, column_title = "Gene changes between Response (G vs NR) and Week Golimumab ")

dev.off()
##get clusters
clusters<-as.data.frame(cluster_real)
clusters$Gene<-rownames(clusters)


write.csv(clusters,"Results/Clusters_Genes_G_NR_ALL.csv")


#####Split in UP and DOWN

####


#####GO

#####GO enrichment analysis
library(clusterProfiler)
library(org.Hs.eg.db)
#library(GOstats)
library(GO.db)
library(DOSE)
library(annotate)

genesid<-rownames(res)##Global

#or ( i in unique(clusters$`km$cluster`))
#genesid<-clusters[cluster_real==4,"Gene"]

genesid<-genesid[ genesid != "" ]##remove empty elements from a vector
genesid <- genesid[!is.na(genesid)]

eg<-bitr(genesid, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

ego2 <- enrichGO(gene         = eg$ENTREZID,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "BP",
                 pAdjustMethod = "none",
                 pvalueCutoff  = 0.1,
                 qvalueCutoff  = 0.05,
                 readable=T
)

png("Results2/GO_Time_Response_Etanercept.png",width = 1200, height = 950)
dotplot(ego2, title="G vs NR w0 Golimumab",showCategory=35)
dotplot(ego2, title="w0 vs w12 Golimumab",showCategory=35)
dotplot(ego2, title="Interaction Golimumab",showCategory=35)

dev.off()

dotplot(ego2, title="Cluster RESP DOWN in w0, DOWN in w12 GvsNR",showCategory=35)



