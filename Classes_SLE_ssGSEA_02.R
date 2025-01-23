library(doParallel)
library(dplyr)
library(edgeR)
library(limma)
library(ggplot2)

setwd("/Users/agomez/IMIDOMICS/ssGSEA_classes/")
####load counts

sampleTable<-read.csv("Final_Table_Curated_Sex_FINAL.csv")
rownames(sampleTable)<-sampleTable$celgIdAll

##load cts
load("Raw_Counts_RNAseq_OK.RData")


##load pheno


all_cts<-all_cts[,colnames(all_cts) %in% rownames(sampleTable)]
sampleTable<-sampleTable[rownames(sampleTable) %in% colnames(all_cts),]

all_cts<-all_cts[,match(rownames(sampleTable),colnames(all_cts))]

##select only disease

A <- as.character(sampleTable[sampleTable$imid == "SLE","celgIdAll"])


##Build DG object
cts<-all_cts[,colnames(all_cts) %in% A]

samples2test<-sampleTable[rownames(sampleTable) %in% A,]

group<-samples2test$imid
group<-droplevels.factor(group)

y <- DGEList(cts,group=group)

y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y, method="TMM")
###Normalize using voom

#remove batch
my.batch<-paste("B", samples2test$new.plate,sep="")
samples2test$Sex<-ifelse(samples2test$sex=="FEMALE","F","M")

##same without batch
design <- model.matrix(~1 + samples2test$Sex + samples2test$age)
v1<-cpm(y, log=TRUE, prior.count=3)
#v1 <- voom(y, plot=TRUE)
design <- model.matrix(~samples2test$Sex + samples2test$age)

data <- removeBatchEffect(v1, batch=my.batch, design=design)

##get genes


mv.genes<-read.csv("genes_SLE_classes.csv")
mv.genes<-as.character(mv.genes$Gene)



##get the new beta matrix each time
mv.genes<-rownames(data)[rownames(data) %in% mv.genes]

dat <-data[mv.genes,]

###Heatmap

library(gplots)
library(RColorBrewer)

#Use M3C approach

library(M3C)
library(NMF) # loading for aheatmap plotting function
library(ggsci) # more cool colours
#we run the algorithm using the default settings (100x monte carlo iterations and 100x inner replications).
#We have found the results generally stable using these parameters, 
#although they may be increased or the algorithm may simply be run a few extra times to test stability. 
#Plots from the tool and an .csv file with the numerical outputs may be printed into the working directory 
#(by adding printres = TRUE). We will set the seed in this example, incase you wish to repeat our results exactly (seed = 123). 
#We will add the annotation file for streamlined downstream analyses (des = annot).

res.m3c <- M3C(dat, cores=7,seed = 123, removeplots = F)
stats<-res.m3c$scores
# get the data out of the results list (by using $ - dollar sign), use 2 clusters (see RCSI plot)
data <- res.m3c$realdataresults[[5]]$ordered_data # this is the data
annon <- res.m3c$realdataresults[[5]]$ordered_annotation # this is the annotation
ccmatrix <- res.m3c$realdataresults[[5]]$consensus_matrix # this is the consensus matrix


# normalise and scale the data
data <- t(scale(t(data))) # z-score normalise each row (feature)
data <- apply(data, 2, function(x) ifelse(x > 4, 4, x)) # compress data within range
data <- apply(data, 2, function(x) ifelse(x < -4, -4, x)) # compress data within range



library(heatmap.2x)## approach

library(RColorBrewer)
#samples<-ifelse(annon$cond=="RA_HI","red","black")

cons<-ggsci::pal_futurama()(max(levels(annon$consensuscluster)))[as.factor(annon$consensuscluster)]#5 clusters


spcol<-cons
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)


###order data


hp<-heatmap.2x(data, col=rev(cols), Colv = NA,Rowv=TRUE, scale="none", ColSideColors=spcol,
           trace="none", dendrogram="both", 
           cexRow=.1, cexCol=1.4,
           main="SLE Most variable genes",
           labCol=NA,
         #  labRow=NA, 
           density.info="none",
           hclust=function(x) hclust(x,method="complete"),
           distfun=function(x) as.dist((1-cor(t(x)))/2)
           
)


#####reorder data according heatmap

annon$ID<-rownames(annon)
write.csv(annon, "Classes_SLE_5_GSVA.csv")
##now establish ssGSE


annon<-read.csv("Classes_SLE_5_GSVA.csv")


##Build DG object
cts<-all_cts[,colnames(all_cts) %in% A]

samples2test<-sampleTable[rownames(sampleTable) %in% A,]

samples2test<-samples2test[match(annon$ID,rownames(samples2test)),]
samples2test$Sex<-ifelse(samples2test$sex=="FEMALE","F","M")

cts<-cts[,match(annon$ID,colnames(cts))]
##establish group

my.list.genes<-list()

##another option is to try to get ONLY the genes dirrerent expressed from the ones , 2510 that we have

for ( i in 1:length(unique(annon$consensuscluster))){ 


    group<-ifelse(annon$consensuscluster==i,"CLASS","OTHERS")
    #group<-annon$consensuscluster
    y <- DGEList(cts,group=group)

    keep <- rowSums(cpm(y)>1) >= min(table(group))#### change group depending replicates
    y <- y[keep, ]

    dim(y)

    ##normalize
    y$samples$lib.size <- colSums(y$counts)
    y <- calcNormFactors(y, method="TMM")


    ####Limma
    my.batch<-as.factor(paste("B", samples2test$new.plate,sep=""))


    design <- model.matrix(~samples2test$Sex + my.batch + samples2test$age+ group )

    colnames(design)[16]<-"Comp"
    rownames(design)<-rownames(samples2test)

    ###add design to the object
    y2 <- estimateDisp(y, design)
    #Now estimate gene-specific dispersions:

      #plotBCV(y2)

    ####radical solution
    fit <- glmQLFit(y2, design,robust=T)
    res <- glmQLFTest(fit, coef="Comp")##or  contrast=B.LvsP
    
    results<-topTags(res,n= dim(res$table)[1],adjust.method="BH", sort.by="PValue")
    results<-results$table

##order by FChange, remember that OTHER is UP in the comparison

   # res.sel<-subset(results,results$FDR<.05)
    res.sel<-subset(results,results$FDR<.05)
    res.sel<-res.sel[order(-res.sel$logFC),]

##apply another filter, more than 1 logFC

    res.sel<-res.sel[res.sel$logFC < - .5,]

    my.list.genes[[i]]<-rownames(res.sel)
   names(my.list.genes)[i]<-paste("CLASS",i,sep="_")
}


##Boxplot test

library(ggpubr)

bcells<-y$counts["ISM1",]
#bcells<-y$counts[rownames(y$counts) %in% my.list.genes[[2]] ,]


bcells<-as.data.frame(bcells)

#bcells<-colMeans(bcells)
rownames(annon)<-annon$X
bcells<-merge(bcells,annon,by="row.names")
colnames(bcells)[2]<-"bcells"

p <- ggboxplot(bcells, x = "consensuscluster", y = "bcells" ,##change depending on set of genes
               color = "consensuscluster",
               add = "jitter",
               legend="")



p<-p+ xlab("")+ ylab("Counts per million") + ggtitle("Class I genes")+
  rotate_x_text(angle = 45)

# Specify the comparisons you want
my_comparisons <- list( c("OTHERS", "CLASS"))


p<- p + stat_compare_means(method="t.test",comparisons = my_comparisons,
                           #
                           # label = "p.signif",###only signifficance symbol stars
                           method.args = list(var.equal = TRUE),##make t-test similar to ANOVA 2 groups,
                           
)

p



##get all the genes

my.genes<-unique(unlist(my.list.genes))
###second idea, remove shared genes

library(GSVA)

##Single  sample  GSEA  (ssGSEA)  calculates  a  gene  set enrichment  score  per  sample  
#as  the  normalized  difference  in  empirical  cumulative  distribution functions of gene expression ranks 
#inside and outside the gene set.
#s que funcionava millor al incloure més gens, més que buscar pocs molt diferencials. 

#ssgsea <- gsva (as.matrix(data[rownames(data) %in% my.genes,]),my.list.genes,method="ssgsea", verbose=TRUE,mx.diff=1)
ssgsea <- gsva (y$counts,my.list.genes,method="ssgsea", verbose=TRUE,ssgsea.norm=FALSE)
#ssgsea <- gsva (as.matrix(data[rownames(data) %in% my.genes,]),my.list.genes,method="gsva", kcdf="Gaussian", verbose=TRUE,mx.diff=1)##most probable

ssgsea<-as.data.frame(t(ssgsea))


annon2<-merge(annon,ssgsea, by="row.names")


assignations<-colnames(annon)[5:9][apply(annon[,5:9],1,which.max)]

annon$Assignations<-assignations
##Use Pascual's data

setwd("/mnt/md127/Analysis_Projects/CelGene_2018/Analysis_september/Classes_RNAseq/External_data_validation/")
####load counts

sampleTable<-read.csv("GSE65391_annot.csv")
rownames(sampleTable)<-sampleTable$X
sampleTable<-sampleTable[,-1]
sampleTable$Age<-gsub("age: ","",sampleTable$characteristics_ch1.13)


###samples used 
samples2use<-read.csv("Pascual_matrix_samples.csv",head=F)
samples2use<-as.character(samples2use$V1)

sampleTable$Code<-gsub("subject: ","",sampleTable$characteristics_ch1.2)
samples2test<-sampleTable[sampleTable$Code %in% samples2use,]
##load expression

exprs.matrix<-read.csv("GSE65391_exprs.csv")
rownames(exprs.matrix)<-exprs.matrix$X
exprs.matrix<-exprs.matrix[,-1]

###load matrix annot

annot.matrix<-read.csv("GSE65391_array.annot.csv")
rownames(annot.matrix)<-annot.matrix$X
annot.matrix<-annot.matrix[,-1]

exprs.matrix<-merge(exprs.matrix,annot.matrix,by="row.names")

###remove probe names

exprs.matrix<-exprs.matrix[,-c(1,998)]
exprs.matrix<-aggregate(.~Gene.symbol,data=exprs.matrix,mean)
rownames(exprs.matrix)<-exprs.matrix$Gene.symbol
exprs.matrix<-exprs.matrix[,-1]

###select samples

exprs.matrix<-exprs.matrix[,colnames(exprs.matrix) %in% rownames(samples2test) ]


##calculat gsva

ssgsea_pascual <- gsva (as.matrix(exprs.matrix[rownames(exprs.matrix) %in% my.genes,]),my.list.genes,method="gsva", kcdf="Gaussian", verbose=TRUE,mx.diff=1)##most probable



library(singscore)


# The recommended method for dealing with ties in ranking is 'min', you can
# change by specifying 'tiesMethod' parameter for rankGenes function.
rankData <- rankGenes(tgfb_expr_10_se)

rankData <- rankGenes(data)

# Given the ranked data and gene signature, simpleScore returns the scores and 
# dispersions for each sample
scoredf <- simpleScore(rankData, upSet = tgfb_gs_up, downSet = tgfb_gs_dn)

scoredf <- simpleScore(rankData, upSet = my.list.genes[[1]])

