
library(survival)
library(edgeR)
library(limma)
setwd("/media/IMID/Projects/NIH/")

###get Pheno

library(GEOquery)

GSE<-"GSE157103"
set1<- getGEO(GSE,GSEMatrix=T,destdir="./",AnnotGPL=F)
set1<-set1[[1]]

####get pheno

pheno<-pData(set1)
pheno$Group<-gsub("disease state: ","",pheno$characteristics_ch1)

###clean pheno


pheno<-pheno[,c(1,63:84)]
pheno$Sample_label<-ifelse(grepl("^COVID",pheno$title, perl=T),substring(pheno$title,1,9),substring(pheno$title,1,11))
pheno$Sample_label<-gsub("_$","",pheno$Sample_label)


##load scores

scores<-read.csv("Scores2survival.csv")
rownames(scores)<-scores$X
PHENO<-merge(pheno,scores,by="row.names")
write.csv(PHENO, "PHENO_COVID.csv")
#####split in two groups acording simple pos/neg score

#high.scores<-as.character(scores[scores$score>0,"X"])
#low.scores<-as.character(scores[scores$score<0,"X"])


#####then try to deal with extreme values
percent<-quantile(as.numeric(scores$score), probs = seq(0, 1, by= 0.1))##in 10%

high.scores<-as.character(scores[scores$score>=percent[8],"X"])##8 works for mech intub as censor
low.scores<-as.character(scores[scores$score<=percent[3],"X"])##3 works for mech intub as censor
#med.scores<-as.character(scores[scores$score>=percent[4] & scores$score<=percent[8]  ,"X"])##3 works

###survival analysis

###add group and clean variable names
toplot<-pheno
toplot$G<-ifelse(rownames(toplot) %in% low.scores,"LOW",ifelse(rownames(toplot) %in% high.scores,"HIGH",NA))###

##ICU
colnames(toplot)[14]<-"ICU"

#toplot<-toplot[toplot$ICU=="yes",]
#toplot<-toplot[toplot$ICU=="no",]

###use here survival original data

surv.data<-read.delim("Survival_data.csv")
##order based on Sample_Label
surv.data<-surv.data[order(surv.data$Sample_label),]
###check that to.plot and surv.data are in same direction

###merge

toplot<-merge(toplot,surv.data,by="Sample_label")


toplot$age<-as.numeric(toplot$`age (years):ch1`)
toplot$sex<-as.factor(toplot$`Sex:ch1`)
toplot$mech_vent<-as.factor(toplot$`mechanical ventilation:ch1`)

####remove samples without age or sex

toplot<-toplot[!(is.na(toplot$age) | is.na(toplot$sex)),]
toplot<-toplot[!(toplot$sex=="unknown"),]
toplot$G<-as.factor(toplot$G)
toplot<-toplot[!(is.na(toplot$G)),]
toplot$G<-droplevels(toplot$G)

toplot$sex<-droplevels(toplot$sex)
####subset only to patients

toplot<-toplot[toplot$Group=="COVID-19",]
#toplot<-toplot[toplot$ICU=="yes",]
####censors
toplot$vital_status<-ifelse(toplot$`mechanical ventilation:ch1`=="yes",1,0)###mech vent
toplot$vital_status<-ifelse(toplot$`hospital-free days post 45 day followup (days):ch1`<25,1,0)##HFD45

#toplot$Days<-ifelse(toplot$Free_days<0,5,toplot$Free_days)
####Plotting

library("survminer")


##no covariates
fit<- survfit(Surv(as.numeric(v), as.numeric(vital_status)) ~ G, data = toplot)
surv.plot<-ggsurvplot(fit,pval=T,title = "",
           #xlim=c(0,2000),
           risk.table = TRUE, risk.table.y.text.col = TRUE,
           break.time.by=5
           
)
surv.plot<-surv.plot+ylab("% Hospital days")
#surv.plot<-surv.plot+ylab("% Hospital-free days at day 45")

p<-surv.plot+ ggtitle(paste0("Abatacept Genes COVID "))
p
####with covariates
####sex as mean age

toplot$AGE<-ifelse(toplot$age>59,"1","0")###mean is 59 years old


fit2<- coxph(Surv(as.numeric(Stay), as.numeric(vital_status)) ~ G +AGE+sex, data = toplot)

sum.f<-summary(fit2)
pval<-round(as.numeric(sum.f$coefficients[,"Pr(>|z|)"][1]),digits = 4)
####plot coxph
surv.plot<-ggadjustedcurves(fit2, data = toplot,
                        #    method = "average",
                   variable  = "G",   # Variable of interest
                    #legend.title = "Sex"      # Change legend title
                    palette = "npg",             # nature publishing group color palettes
                    curv.size = 2                # Change line size
)
p<-surv.plot+ylab("% Hospitalization days")
#+ ggtitle(paste0("COXPH Model Censor HFD45 ","pval=",pval))+ylab("% Hospital stay days")
# Drawing survival curves
#pdf("Survival_Results/Survival_COXPH_COVID_ventilation_extreme.pdf",onefile=FALSE)
p+theme(legend.position = "right")

#dev.off()
###HAzard Ration and CI

coxobj_summary <- summary(fit2)
coxobj_summary$conf.int

forest<-ggforest(fit2)

forest<-forest+ggtitle(paste0("Hazard Ratio Censor HFD45 "))




############
###check using several measures if each group is correlated with them

library("ggpubr")

library(reshape)

####recover pheno and this DO NOT remove 0 days

toplot<-merge(toplot,scores,by="row.names")

colnames(toplot)[11]<-"ferritin"
colnames(toplot)[8]<-"ddimer"
colnames(toplot)[12]<-"fibrinogen"
colnames(toplot)[7]<-"crp"
colnames(toplot)[14]<-"ICU"
colnames(toplot)[15]<-"lactate"
colnames(toplot)[3]<-"ApacheII"
colnames(toplot)[6]<-"Charlson"
colnames(toplot)[20]<-"Sofa"
toplot$ferritin<-as.numeric(toplot$fibrinogen)
toplot$ddimer<-as.numeric(toplot$ddimer)
toplot$fibrinogen<-as.numeric(toplot$fibrinogen)
toplot$crp<-as.numeric(toplot$crp)
toplot$lactate<-as.numeric(toplot$lactate)
toplot$ApacheII<-as.numeric(toplot$ApacheII)
toplot$Charlson<-as.numeric(toplot$Charlson)
toplot$Sofa<-as.numeric(toplot$Sofa)
#toplot<-toplot[!(toplot$Vent_days==0),]##recuperation?
###for all samples
toplot$combined_g<-paste0(toplot$G,"_",toplot$ICU)
toplot$combined_g<-gsub("HIGH_no","HIGH_NO_ICU",toplot$combined_g)
toplot$combined_g<-gsub("HIGH_yes","HIGH_ICU",toplot$combined_g)
toplot$combined_g<-gsub("LOW_yes","LOW_ICU",toplot$combined_g)
toplot$combined_g<-gsub("LOW_no","LOW_NO_ICU",toplot$combined_g)
#Corr_Ferritin_Score
png(paste0("Survival/Corr_Score_","lactate",".png"),width = 800, height = 800)

ggscatter(toplot, x = "score", y = "ferritin", color="combined_g",
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "blue", fill = "lightgray"),
    ##      #cor.coeff.args = list(method = "spearman", label.y = 5,label.x = .5, label.sep = "\n"),
          xlab = "Abatacept Score", ylab = "crp")+
  stat_cor(method = "spearman", label.x =10, label.y = 8)  # Add correlation coefficient

dev.off()

##Boxplots better

toplot$combined_g <- factor(toplot$combined_g,levels = c("HIGH_ICU","HIGH_NO_ICU","LOW_ICU", "LOW_NO_ICU"))

p <- ggboxplot(toplot, x = "combined_g", y = "ferritin",##change depending on set of genes
               color = "combined_g",
               add = "jitter",
               legend="",
               outlier.shape = ""
               #ylim=c(-0.4,2)
)


p<-p+ ggtitle(paste0("COVID samples ","ferritin"))+
  rotate_x_text(angle = 45)+xlab(" ")
# Specify the comparisons you want


my_comparisons <- list( c("HIGH_ICU","HIGH_NO_ICU"),c("LOW_ICU", "LOW_NO_ICU"),
                        c("HIGH_ICU","LOW_ICU"),c("HIGH_NO_ICU", "LOW_NO_ICU"),
                        c("HIGH_ICU","LOW_NO_ICU")
                        
                        )

p<- p + stat_compare_means(method="t.test",comparisons = my_comparisons,
                           #
                           # label = "p.signif",###only signifficance symbol stars
                           method.args = list(var.equal = TRUE),##make t-test simila to ANOVA 2 groups,
                           
)

png(paste0("Survival/Boxplot_COVID_","Charlson",".png"),width = 800, height = 800)

p

dev.off()


#########table

totable<-toplot

covid<-totable[totable$Group=="COVID-19",]
noncovid<-totable[totable$Group=="non-COVID-19",]

#crp (mg/l)
#ferritin (ng/ml)

totest<-toplot[,c("crp (mg/l):ch1","Group")]
###categorical
totest<-totest[!(is.na(totest$G)),]
totest$mech_vent<-as.factor(totest$mech_vent)
totest$sex<-droplevels(totest$sex)
chisq.test(table(totest))
###numerical
totest<-totest[!(totest$`crp (mg/l):ch1`=="unknown"),]

mean(as.numeric(totest[totest$Group=="COVID-19",1]))
sd(as.numeric(totest[totest$Group=="COVID-19",1]))

mean(as.numeric(totest[totest$Group=="non-COVID-19",1]))
sd(as.numeric(totest[totest$Group=="non-COVID-19",1]))

t.test(as.numeric(totest$`ferritin (ng/ml):ch1`)~totest$Group)
wilcox.test(as.numeric(totest$`crp (mg/l):ch1`)~totest$Group)

       