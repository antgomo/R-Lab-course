library(tximport)
library(sva)
library(RUVSeq)
library(splines)
library(limma)
setwd("/Users/agomez/IMIDOMICS/Proteins/")

####################LOAD proteins

options(stringsAsFactors=F) 

##plasm codes

plasm.codes<-read.csv("Final_plasm_assoc.csv")
colnames(plasm.codes)[2]<-"Assay"
###only basal week experiments

#plasm.codes<-plasm.codes[plasm.codes$W=="0",]

#load proteins
prot1<-read.csv("Olink_inflammation.csv")
prot2<-read.csv("Olink_immuno.csv")
prot3<-read.csv("Olink_cardio.csv")
prot4<-read.csv("sPD-1 sICAM1 CXCL13 anti_PAD4 results.csv")
prot4<-prot4[,-1]
colnames(prot4)[1]<-"Assay"
prot4$Anti.PAD4<-as.numeric(as.character(prot4$Anti.PAD4))
protb<-log2(prot4[,-1])
protb$Assay<-prot4$Assay

###log2


prot.all<-Reduce(function(x,y) merge(x,y,by="Assay",all=TRUE) ,list(prot1,prot2,prot3,protb))

#merge both

prot<-merge(prot.all, plasm.codes,by="Assay")
prot1<-merge(prot1, plasm.codes,by="Assay")
prot2<-merge(prot2, plasm.codes,by="Assay")



###Proteins

prot1<-prot1[,-1]
prot2<-prot2[,-1]


####

###Revisant analisis:
 ## Es varen eliminar per tenir un nombre molt gran de missings:
 ## IXTCD01371 IXTCD01407 IXTCD01433 
##Es varen eliminar per tenir una expresió outlier vs resta:
  ##IXTCD01339 IXTCD01359 IXTCD01393



##select basal 

#prot1<-prot1[prot1$timePoint=="W0",]
#prot2<-prot2[prot2$timePoint=="W0",]


prot3<-prot1
prot3<-prot3[prot3$timePoint=="W0",]


prot3$R<-ifelse(prot3$response_EULAR=="GOOD","R","NR")

prot3<-prot3[!(prot3$donorId %in% c("IXTCD01371","IXTCD01407","IXTCD01433","IXTCD01339","IXTCD01359","IXTCD01393")),]

pca <- prcomp(na.omit(prot3[,1:280]))#avoid because you logged it before
pca$Perc<-round(100*((pca$sdev)^2 / sum(pca$sdev^2) ), digits=2)



##line by line


results.glm<-matrix(nrow=dim(prot3)[2]-16,ncol=5)

for (i in 1:(dim(prot3)[2]-16)){
  prot4<-prot3[!(is.na(prot3[,i])),]
  prot4<-prot3[!(is.na(prot3[,105])),]
  prot4<-prot3[!(is.na(prot3[,106])),]
  
 # fit <- glm(as.factor(prot[,290])~prot[,i]+prot[,279]+prot[,280],family="binomial")##all
  fit <- glm(as.factor(prot4[,108])~prot4[,i]+prot4[,104]+prot4[,105]+prot4[,106],family="binomial")##all
  
  ##extract odss ratio as fold change
  odds<-exp(cbind("Odds ratio" = coef(fit), confint.default(fit, level = 0.95)))
  odds<-log(odds[2,])
  results.glm[i,2]<-coef(summary(fit))[2,4]
  results.glm[i,1]<-colnames(prot4)[i]
  results.glm[i,3:5]<-odds
  
  
}

results.glm<-as.data.frame(results.glm)
colnames(results.glm)<-c("Protein","pval","logFC","CI","CR")

rownames(results.glm)<-results.glm$Protein


write.csv(results.glm,"Prot_log_mod_res_w0_INF.csv")


#You can also use the confint.default function which is based on asymptotic normality.



#results.glm<-rbind()

###plots

results.glm$pval<-as.numeric(as.character(results.glm$pval))
results.glm$SIG<-ifelse(results.glm$pval<.05,"YES","NO")
results.glm$logFC<-as.numeric(as.character(results.glm$logFC))
results.glm$CI<-as.numeric(as.character(results.glm$CI))
results.glm$CR<-as.numeric(as.character(results.glm$CR))

#results.glm<-results.glm[complete.cases(results.glm),]

library(ggplot2)



vp<-ggplot(data=results.glm,aes(x=logFC,y=-log10(pval)))+
  geom_point(aes(size=-log(pval),fill=SIG), colour="black",shape=21)+
  geom_vline(xintercept=0,linetype="dashed",col="red")+
  #geom_vline(xintercept=c(-0.5,0.5),linetype="dashed",col="blue")+
  #ggtitle("Ulcerative Colitis HI vs LOW") +
  theme(legend.position="none")+
  theme_minimal()


library(ggpubr)
library(ggrepel)

vp<-vp+geom_text_repel(data=dplyr::filter(results.glm, pval<0.05), aes(label=Protein),point.padding = 2)

png("Volcano_Cardio_WO.png",width = 800, height = 880,)
vp+ ggtitle("RESP vs NON_RESP W0 IMM") 
dev.off()

results.glm<-results.glm[order(results.glm$logFC,decreasing = T),]

results.glm$Protein <- factor(results.glm$Protein, levels=rev(results.glm$Protein))###ordering from top to last


fp<-ggplot(data=results.glm,aes(x=logFC,y=Protein))+
  geom_point(aes(size=-log(pval),fill=SIG), colour="black",shape=21)+
  geom_errorbarh(aes(xmin=CI,xmax=CR),color="grey",size=.5,alpha=.7)+
  geom_vline(xintercept=0,linetype="dashed")+
  ggtitle("RESP vs NON_RESP W0 Immune Panel")+
  # geom_text(aes(x=2.8,label=SIG),size=4)+
  ylab(NULL)+
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  theme(axis.text.y = element_blank())+
  #scale_y_discrete(expand = c(.05,.05))
  #scale_y_discrete(drop=FALSE)+
  expand_limits(y=-2)



fp<-fp+geom_text_repel(data=dplyr::filter(results.glm, pval<0.1), aes(label=Protein),point.padding = 3)

png("FP_Cardio_WO.png",width = 800, height = 880,)
fp
dev.off()


figure<-ggarrange(fp,vp + rremove("x.text"), 
                 # labels = c("CD vs Control", "HI vs LOW"),
                  ncol = 2, nrow = 1,common.legend = TRUE, legend = "bottom")

annotate_figure(figure,
                top = text_grob("RESP vs NONRESP W0", color = "black", size = 14)
                # bottom = text_grob("Data source: \n mtcars data set", color = "blue",
                #  hjust = 1, x = 1, face = "italic", size = 10)
                
                # fig.lab = "Figure 1", fig.lab.face = "bold"
)


#Boxplot




  p <- ggboxplot(dat2[,c("Anti.PAD4","RESP")], x = "RESP", y = "Anti.PAD4",##change depending on set of genes
                 color = "RESP",
                 add = "jitter",
                 legend="")
  
  y.df<-dat2[,c("Anti.PAD4","RESP")]
  result<-aov(y.df[,1]~as.factor(y.df$RESP))
  
  res<-summary(result)
  pval<-res[[1]][["Pr(>F)"]][1]
  pval=0.03537505
  pval<-signif(pval, digits=3)
  
  p<-p+ ylab("log2(NPX)") + xlab("") + ggtitle(paste("Anti.PAD4",pval,sep=" "))+
    rotate_x_text(angle = 45)+
    geom_hline(yintercept = median(dat2[dat2$RESP=="R","Anti.PAD4"]), linetype = 2,color="lightblue") # Add horizontal line at bas
  # Specify the comparisons you want
  my_comparisons <- list(  c("R", "NR"))
  
  p<- p + stat_compare_means(method="t.test",comparisons = my_comparisons,
                             #
                              label = "p.signif",###only signifficance symbol stars
                            # method.args = list(var.equal = TRUE),##make t-test simila to ANOVA 2 groups,
                             symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                                 symbols = c("****", "***", "**", "*", "ns"))
  )
p



