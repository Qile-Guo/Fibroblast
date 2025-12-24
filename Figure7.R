#---------Figure7-------#


# Library loading

library(tidyr)
library(data.table)
library(dplyr)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(Seurat)
library(ggsci)
library(RColorBrewer)
library(pheatmap)
library(circlize)
library(grDevices)
library(survminer)
library(survival)
library(ggpubr)
library(corrplot)
library(patchwork)
library(stringr)
library(tibble)
library(clusterProfiler)
library(scales)

######## Figure 7C Survival analysis
# data preparation
mat <- read.csv("EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena",sep="\t")
rownames(mat) <- make.unique(as.character(mat[, 1]))
mat <- mat[, 2:ncol(mat)]
meta1<-read.csv("~/analysis/CAF_Project/survival/TCGA_pancancer/TCGA_phenotype_denseDataOnlyDownload.tsv",sep="\t")
meta2<-read.csv("~/analysis/CAF_Project/survival/TCGA_pancancer/Survival_SupplementalTable_S1_20171025_xena_sp",sep="\t")
sample<-intersect(meta1$sample,meta2$sample)
rownames(meta1)<-meta1$sample
rownames(meta2)<-meta2$sample
meta1<-meta1[sample,]
meta2<-meta2[sample,]
meta<-cbind(meta1,meta2)
newrownames<-rownames(meta)
for (i in (1:length(newrownames))){
    newrownames[i]<-str_replace_all(string =newrownames[i],pattern = "-",replacement = ".")
}
rownames(meta)<-newrownames
adf = t(mat)
meta = meta[,c("sample","X_PATIENT","cancer.type.abbreviation","age_at_initial_pathologic_diagnosis","gender","race","ajcc_pathologic_tumor_stage","OS","OS.time","sample_type")]
use=intersect(rownames(meta),rownames(adf))
meta = meta[use,]
adf = adf[use,]
meta = cbind(meta,adf)
meta = meta[!is.na(meta$OS.time),]
meta = meta[!is.na(meta$age_at_initial_pathologic_diagnosis),]
meta$age = meta$age_at_initial_pathologic_diagnosis
meta$cancerType = meta$cancer.type.abbreviation
meta = meta[!is.na(meta$ajcc_pathologic_tumor_stage),]
meta$stage = meta$ajcc_pathologic_tumor_stage
meta = meta[grepl("Stage I",meta$stage),]
meta$stage = gsub("[ABCD]$","",meta$stage)
meta$cancerType = ifelse(meta$cancerType=="COAD"|meta$cancerType=="READ", "CRC", meta$cancerType)
meta$cancerType = ifelse(meta$cancerType=="LUSC"|meta$cancerType=="LUAD", "NSCLC", meta$cancerType)
meta$cancerType = ifelse(meta$cancerType=="KIRC"|meta$cancerType=="KIRP"|meta$cancerType=="KICH", "RCC", meta$cancerType)
meta$OS.time <- meta$OS.time / 365
meta$OS[meta$OS.time>7] <- 0
meta$OS.time[meta$OS.time>7] <- 7
meta$OS.use=meta$OS
meta$OS.time.use=meta$OS.time
# survival analysis
do.cox = function(df, type, cancer){
  df$group = factor(df$group,levels = c('low','high'))
  if (length(unique(na.omit(df$group))) < 2) {
    stop("分组变量必须包含至少两个有效水平")
  }
  if (type == "reg") {
    formula = if (length(unique(df$gender)) > 1) {
      Surv(OS.time.use, OS.use) ~ group + gender + age + stage
    } else {
      Surv(OS.time.use, OS.use) ~ group + age + stage
    }
  } else {
    formula = Surv(OS.time.use, OS.use) ~ group
  }
  
  cox = coxph(formula, data = df, na.action = na.omit)
  cox.summary = summary(cox)
  group_rows = grep("^group", rownames(cox.summary$conf.int))
  if (length(group_rows) == 0) {
    stop("模型中未检测到分组变量的结果，请检查变量命名")
  }
  HR = round(cox.summary$conf.int[group_rows, 1], 2)
  HR.lower = round(cox.summary$conf.int[group_rows, 3], 2)
  HR.upper = round(cox.summary$conf.int[group_rows, 4], 2)
  HR.range = sprintf("(%.1f-%.1f)", HR.lower, HR.upper)
  coef = cox.summary$coefficients[group_rows, 1]
  coef.se = cox.summary$coefficients[group_rows, 3]
  Pval = round(cox.summary$coefficients[group_rows, 5], 4)
  group = sub("^group", "", rownames(cox.summary$conf.int)[group_rows])
  return(data.frame(
    cancerType = cancer,
    group = group,
    HR = HR,
    HR.range = HR.range,
    coef = coef,
    coef.se = coef.se,
    Pval = Pval
  ))
}
mat<-mat[,rownames(meta)]
rn<-rownames(mat)
mat<-apply(mat,2,as.numeric)
rownames(mat)<-rn
# survival analysis
meta$ADAM12<-as.numeric(mat["ADAM12",])/as.numeric(mat["LUM",])
dir="ADAM12_LUM.out"
for (i in c("ADAM12")) {
  dir.create(sprintf("%s/HR1/%s",dir,i), F, T)
  types = c("reg")
  for (type in types){
    cdat = data.frame(cancerType = 'new', group='new', HR=1, HR.range='(0-0)', coef=1, coef.se=1, Pval=1)
    for (cancer in c("panCan", unique(meta$cancerType))){
      pdat = meta
      if(cancer!="panCan"){
        pdat = meta[meta$cancerType==cancer,]
      }
      res.cut <- surv_cutpoint(pdat, time = "OS.time.use", 
                               event = "OS.use", 
                               variables = i
                               )
      res.cat <- surv_categorize(res.cut)
      pdat$group <- res.cat[,i] 
      pdat$group = factor(pdat$group)
      
      test = do.cox(pdat,type,cancer)
      cdat = rbind(cdat,test)
    }
    cdat = cdat[cdat$cancerType!='new',]
    cdat1 = cdat %>% group_by(group) %>% mutate(adj.Pval=p.adjust(Pval, method="BH"))  
    write.table(cdat1, file=sprintf("%s/HR1/%s/HR.sepCan_%s.txt",dir,i,type), sep="\t", quote=F, col.names=T, row.names=F)
    cdat2 = read.table(sprintf("%s/HR1/%s/HR.sepCan_%s.txt",dir,i,type), sep="\t", header=T, stringsAsFactors=F)
    cdat2$HR.p = cdat2$HR
    cdat2$HR.p = ifelse(cdat2$HR.p>2, 2, cdat2$HR.p)
    cdat2$HR.range.p = sprintf("%.1f\n%s", cdat2$HR, cdat2$HR.range)
    cdat2[!is.na(cdat2$HR) & abs(cdat2$HR)>1000, "HR"] = NA
    cdat2[is.na(cdat2$HR), "HR.range.p"] = "Inf"
    cdat2[cdat2$HR.range.p=="Inf","HR.p"] = 1
    cdat2$sig = ifelse(cdat2$adj.Pval <0.05, "adj.P<0.05", NA)
    cdat2$sig = factor(cdat2$sig, levels=c("adj.P<0.05"))
    sig.color = c("#EE2C2C")
    names(sig.color) = c("adj.P<0.05")
    p = ggplot(cdat2, aes(x=cancerType, y=group, fill=HR.p)) +
      geom_tile(size=3) +
      geom_tile(aes(color=sig),size=1.2, alpha=0) +
      geom_text(aes(label=HR.range.p), angle=0, color="black", size=3) +
      scale_fill_gradient2( high="#EEAEEE", mid="#FFFFFF", low="#1C86EE", midpoint=1, limits=c(0,2), name="HR", na.value="#FFFFFF") +
      scale_colour_manual(name="significance",values=sig.color) +
      theme_classic() + 
      theme(axis.text.x=element_text(vjust=1, hjust=0, angle=315)) +
      xlab('cancerType') + ylab("group") + ggtitle(sprintf("baseline: %s", levels(meta$immuneType)[1]))
    
    width = 5 + length(unique(cdat2$cancerType))*0.55
    height = 1.75 + length(unique(cdat2$cellType))*0.55
    ggsave(p, file=sprintf("%s/HR1/%s/COX.sepCan_%s.pdf",dir,i,type), width=width, height=height, limitsize=F)
    ###meta
    ssdir = sprintf("%s/HR1/%s/meta/",dir,i)
    dir.create(ssdir, F, T)
    cdat3 = read.table(sprintf("%s/HR1/%s/HR.sepCan_%s.txt",dir,i,type), sep="\t", header=T, stringsAsFactors=F)
    cdat3 = cdat3[cdat3$cancerType!='panCan',]
    cdat3 = cdat3[cdat3$cancerType %in% c("BRCA","ESCA","NSCLC","BLCA","STAD","HNSC","SKCM","RCC","CRC","LIHC","THCA","PAAD","MESO","CHOL"),]
    cdat3 = cdat3[order(cdat3[,3]),]
    meta = meta::metagen( TE=cdat3$coef, seTE=cdat3$coef.se, studlab=cdat3$cancerType,
                          comb.fixed=F, comb.random=T, prediction=F, sm="HR")
    pdf(sprintf("%s/%s.pdf",ssdir,type), width=14, height=10)
    meta::forest(meta, layout="JAMA", test.overall.random=T, digits.pval=4,
                 colgap.forest.left="5cm", zero.pval=T)
    dev.off()
    print(meta)
  } 
}

######## Figure 7D Immune-infiltration analysis
# data preparation
mat <- read.csv("./TCGA_pancancer/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena", sep = "\t")
rownames(mat) <- make.unique(mat[[1]])
mat <- as.matrix(mat[, -1, drop = FALSE])
meta <- read.csv("~/analysis/CAF_Project/survival/TCGA_pancancer/Survival_SupplementalTable_S1_20171025_xena_sp", sep = "\t")
rownames(meta) <- meta$sample
rownames(meta) <- gsub("-", ".", rownames(meta))
meta <- meta %>%
  mutate(
    cancer.type.abbreviation = case_when(
      cancer.type.abbreviation %in% c("READ", "COAD") ~ "CRC",
      cancer.type.abbreviation %in% c("LUAD", "LUSC") ~ "NSCLC",
      cancer.type.abbreviation %in% c("KIRC", "KIRP", "KICH") ~ "RCC",
      TRUE ~ as.character(cancer.type.abbreviation)
    )
  )

# align samples
common_samples <- intersect(colnames(mat), rownames(meta))
mat <- mat[, common_samples, drop = FALSE]
meta <- meta[common_samples, , drop = FALSE]
if (!is.numeric(mat)) {
  mat <- apply(mat, 2, as.numeric)
  rownames(mat) <- rownames(mat)  # 保持行名
}
# read infiltration data
infil <- read.csv("infiltration_estimation_for_tcga.csv")
rownames(infil) <- gsub("-", ".", infil$cell_type)
infil$cell_type <- NULL  # 移除重复列
# align samples
common_samples2 <- intersect(colnames(mat), rownames(infil))
mat <- mat[, common_samples2, drop = FALSE]
meta <- meta[common_samples2, , drop = FALSE]
infil <- infil[common_samples2, , drop = FALSE]
# add gene expr
infil$ADAM12 <- mat["ADAM12", ]
infil$LUM <- mat["LUM", ]
infil$cancer <- meta[rownames(infil), "cancer.type.abbreviation"]
# calculate correlation
method_list <- data.frame(
  method = sapply(strsplit(colnames(infil)[2:120], "_"), `[`, 2),
  celltype = colnames(infil)[2:120],
  stringsAsFactors = FALSE
)
method_select <- c("CIBERSORT", "EPIC")
method_list <- method_list[method_list$method %in% method_select, ]
gene <- "ADAM12"
results_list <- list()
for (c in unique(infil$cancer)) {
  sub_infil <- infil[infil$cancer == c, , drop = FALSE]
  gene_ratio <- sub_infil[, gene] / sub_infil$LUM
  for (i in method_list$celltype) {
    cor_test <- cor.test(sub_infil[, i], gene_ratio, method = "spearman")
    results_list[[length(results_list) + 1]] <- data.frame(
      celltype = i,
      rho = cor_test$estimate,
      p = cor_test$p.value,
      cancer = c,
      method = method_list$method[method_list$celltype == i],
      stringsAsFactors = FALSE
    )
  }
}
data <- do.call(rbind, results_list)
rownames(data) <- NULL
data <- data[data$celltype %in% c("T.cell.CD8._CIBERSORT", "T.cell.CD8._EPIC"), ]
data$fdr <- p.adjust(data$p, method = "BH")
data$cell <- sapply(strsplit(data$celltype, "_"), `[`, 1)
data$pvalue <- ifelse(data$fdr < 0.05, "sig", "ns")
selected_cancers <- c('BRCA', 'ESCA', 'NSCLC', 'BLCA', 'STAD', 'HNSC', 'SKCM', 'RCC', 'CRC', 'LIHC', 'THCA', 'PAAD', 'MESO', 'CHOL')
ADAM12_CD8 <- data %>% filter(cancer %in% selected_cancers) %>% mutate(celltype = paste0("ADAM12_", celltype))
# visualization
p <- ggplot(ADAM12_CD8, aes(x = cancer, y = celltype, shape = pvalue, color = rho)) +
  geom_point(size = 5) +
  scale_color_gradient2(
    breaks = c(-1, -0.5, 0, 0.5, 1),
    low = "#001d73", mid = "#e8edff", high = "#ff4c4c", midpoint = 0
  ) +
  scale_shape_manual(values = c(7, 15)) +
  labs(title = "LUM normalization") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
print(p)


######## Figure 7E Immune therapy response

CRC_p_meta<-read.csv("CRC_p_meta.csv")
colnames(CRC_p_meta)[1]<-"Patient"

load("CRC_ICB.Rdata")
CRC_ICB_CAF<-subset(CRC_ICB,SubCellType %in% c("c72_Fibro_SOX6", "c73_Fibro_ADAMDEC1",
                                               "c74_Fibro_C7","c75_Fibro_CCL19","c76_Fibro_PI16",
                                               "c77_Fibro_ACTA2","c78_Fibro_DES","c79_Fibro_FAP","c80_Fibro_MKI67"))

CRC_ICB_CAF_T<-subset(CRC_ICB_CAF,Tissue =="Tumor")


data_sc<-CRC_ICB_CAF_T@meta.data

ADAM12<-FetchData(CRC_ICB_CAF_T,vars = c(c("ADAM12","LUM")))
ADAM12<-as.data.frame(ADAM12)
ADAM12<-ADAM12[rownames(data_sc),]
data_sc<-cbind(data_sc,ADAM12)



data_l<-melt(data_sc,id.vars = c("Patient" ,  "Treatment" ,"Cell","Response"),measure.vars = "ADAM12")

data_l$point<-data_l$Treatment
data_l$point<-gsub("IV","Post",data_l$point)
data_l$point<-gsub("III","On",data_l$point)
data_l$point<-gsub("II","On",data_l$point)
data_l$point<-gsub("I","Baseline",data_l$point)

data_l$Condition<-paste0(data_l$point,"_",data_l$Response)
data_l$responder<-data_l$Response
data_l$responder<-gsub("CR","R",data_l$responder)
data_l$responder<-gsub("PR","R",data_l$responder)
data_l$responder<-gsub("SD","NR",data_l$responder)


data_l_mean <- data_l %>%
  group_by(Patient,variable,point,Condition,Response,responder) %>%
  summarise(value=mean(value))

data_l_mean_m1<-dcast(data_l_mean,Patient+variable+Response+responder~point ,value.var = "value")

data_l_mean_m1$Change<-data_l_mean_m1$On-data_l_mean_m1$Baseline


p<-ggplot(data= data_l_mean_m1[!is.na(data_l_mean_m1$Change),],
          aes(x = responder, y = Change,color=responder))+
  theme_classic()+
  geom_point(size=2)+
  geom_boxplot(outliers = F,fill=NA)+
  theme_classic()+
  scale_color_manual(values = c("#3F89B9","#D7832D"))+ 
  theme(plot.title = element_text(hjust = 0.5,size=10,face = "bold"), 
        legend.position="right", 
        axis.title.x =element_text(size=10), 
        axis.title.y=element_text(size=10),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.text=element_text(size=10),
        legend.title = element_text(size=10),
        panel.grid = element_blank())+
  labs(x="",y="")+NoLegend()+ 
  theme(strip.background = element_blank(),strip.text = element_text(size = 10,color = "black") )+
  theme(
    strip.background = element_rect(
      color = NA, fill = NA),
    panel.grid = element_blank())+
  geom_signif(comparisons = list(c("R", "NR")),
              test = t.test,
              test.args = list(alternative="two.side"),
              y_position = c(0.5),
              parse = T,tip_length = 0,color="black")

p

############## Tumor regression

data_l_mean_m2<-merge(data_l_mean_m1,CRC_p_meta,by="Patient")

p2<-ggplot(data =data_l_mean_m2,aes(x=Change,y=Tumor.Regression.Ratio))+
  geom_point(color="dimgray",size=2.8)+
  theme_classic()+
  theme(panel.grid = element_blank(),legend.position = "none")+
  scale_color_manual(values = c("#4678A0","#D6A859","#832B20"))+ 
  labs(x="ADAM12 dynamics",y="Tumor regression")+
  theme(plot.title = element_text(hjust = 0.5,size=24), 
        legend.position="right", 
        axis.title.x =element_text(size=12,color="black"), 
        axis.title.y=element_text(size=12,color="black"),
        axis.text=element_text(size=12,color="black")) +
  geom_smooth(method = 'lm', formula = y ~ x, se = T,color="dimgray",fill="dimgray",alpha = 0.1,linewidth=0.5)+ 
stat_cor(method = "pearson")

p2

