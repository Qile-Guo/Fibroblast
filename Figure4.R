#---------Figure4-------#

# Library loading
library(clusterProfiler)
library(org.Mm.eg.db)
library(homologene)
library(enrichplot)
library(DOSE)
library(pheatmap)
library(grDevices)
library(tidyr)
library(data.table)
library(dplyr)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(Seurat)
library(ggsci)
library(RColorBrewer)
library(ggrepel)
library(fmsb)
library(GSVA)

obj<-readRDS("coinjection_fibroblast.rds")
color<-c("#CAB2D6","#0055AA","#FFACAA","#EAC862","#7FD2FF","#00C19B","#B2DF8A","#007ED3","#FFACAA","#C40003","#FF9D1E","#C3EF00","#894FC6" )

######### Figure E Radar plot for GSVA
# GSVA
ref<-read.csv("Ex_vivo_program.csv")
Idents(obj) <- "Condition" 
for (i in (1:length(colnames(ref)))){
  if (i==1){
    gene_name=na.omit(ref[,i])
    gene_name<-homologene(gene_name, inTax = 9606, outTax = 10090)[,2]
    path_ID=rep(colnames(ref)[i],length(gene_name))
  }
  if (i>1){
    gene_name<-c(gene_name,(homologene(na.omit(ref[,i]), inTax = 9606, outTax = 10090)[,2]))
    path_ID=c(path_ID,rep(colnames(ref)[i],length((homologene(na.omit(ref[,i]), inTax = 9606, outTax = 10090)[,2]))))
  }
}
path_name<-as.data.frame(cbind(path_ID,gene_name))
colnames(path_name)<-c("term","gene")
expr <- AverageExpression(obj, assays = "RNA", slot = "data")[[1]]
expr <- as.matrix(expr[rowSums(expr)>0,])
GSVA_Set = path_name %>% split(x = .$gene, f = .$term)
gsva_res <- gsva(expr, gset.idx.list = GSVA_Set, kcdf="Gaussian",method = "gsva",parallel.sz=1)
# Visualization
gsva_res_sub<-gsva_res[c("Progenitor.like","ECM.remodeling","Inteferon.response","Inflammation","Collagen","Myofibro.like"),]
data<-reshape2::melt(gsva_res_sub)
colnames(data)<-c("Module","Condition","Value")
data$Value<-as.numeric(data$Value)
gsva_res_sub<-cbind(gsva_res_sub,as.data.frame(rep(0.5,6)))#apply(gsva_res_sub,1,max)
gsva_res_sub<-cbind(gsva_res_sub,as.data.frame(rep(-0.5,6)))
colnames(gsva_res_sub)<-c("Control","KO","Max","Min")
gsva_res_sub<-t(gsva_res_sub)
gsva_res_sub<-as.data.frame(gsva_res_sub[c(3,4,1,2),])
radarchart(
  gsva_res_sub, axistype = 1,
  pcol = c("#BAB9B9", "#518D9F", "#FC4E07"), pfcol = scales::alpha(c("#BAB9B9", "#518D9F", "#FC4E07"),0.5), plwd = 2, plty = 1,
  cglcol = "grey", cglty = 1, cglwd = 0.8,
  axislabcol = "grey", 
  vlcex = 1,
  caxislabels = c(-0.5, -0.25, 0, 0.25, 0.5))
legend(
 x = "bottom", legend = rownames(gsva_res_sub[-c(1,2),]), horiz = TRUE,
  bty = "n", pch = 20 , col = c("#BAB9B9", "#518D9F", "#FC4E07"),
  text.col = "black", cex = 1, pt.cex = 1.5
)

######## Figure F Volcano Plot for DEGs
# calculate DEGs
Idents(obj)<-"Condition"
marker.list<-FindMarkers(obj,ident.1="KO",ident.2="Control",only.pos = FALSE,logfc.threshold = 0)
up_gene<-rownames(subset(marker.list,(avg_log2FC>0.25)&(p_val_adj<0.01)))
down_gene<-rownames(subset(marker.list,(avg_log2FC<(-0.25))&(p_val_adj<0.01)))
# Visualization
cut_off_pvalue = 0.01 
cut_off_logFC = 0.25   
marker.list$change = ifelse(marker.list$p_val_adj < cut_off_pvalue & abs(marker.list$avg_log2FC) >= cut_off_logFC, ifelse(marker.list$avg_log2FC> cut_off_logFC ,'Up','Down'),'Stable')
marker.list$gene<-rownames(marker.list)
p <- ggplot(
  marker.list, aes(x = avg_log2FC, y = -log10(p_val_adj), colour=change)) +
  geom_point(alpha=1, size=1) +
  scale_color_manual(values=c("#70bce5", "#d2dae2","#ee908e"))+
  geom_vline(xintercept=c((-1*cut_off_logFC),cut_off_logFC),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())
label_gene=c("Ly6c1","Ifit1","Lum","Pi16","Ly6a","Irf7","Isg15","Ifi205","Mfap2","Cd34","Mfap5","Ifi204","Cd55","Mmp3","Cxcl1","Col5a3","Col5a2","Lrrc15","Tagln","Spp1","Inhba","Ltbp2","Tnc","Timp1","Tpm2")
marker.list$label <- ifelse(marker.list$gene %in% label_gene,as.character(marker.list$gene), "")
p2<-p + geom_text_repel(aes(x=marker.list$avg_log2FC,y=-log10(marker.list$p_val_adj),label = marker.list$label),point.padding =1,force=1,nudge_x=0.5,nudge_y=0,max.overlaps = 10000)
ggsave("Fibro_DEG_volcano.pdf",width=5,height=4)

######## Figure G UMAP of fibroblast
Idents(obj)<-"Annotation"
levels(obj)<-c("Pi16_CAF","Apod_CAF","Lrrc15_CAF","Mki67_CAF","Mt1_CAF")
p<-DimPlot(obj,reduction="umap",cols=color,pt.size = 3)
ggsave("Dimplot_Fibro.pdf",width=7,height=5)

####### Figure H Heatmap of markers
df <- AverageExpression(obj, verbose=F)$RNA
marker<-c("Pi16","Dpt","Ifit1","Apod","C3","Ifi205","Acta2","Lrrc15","Col12a1","Mki67","Top2a","Ube2c","Mt1","Mt2")
df<-df[marker,]
pheatmap(df[,c("Pi16_CAF","Apod_CAF","Lrrc15_CAF","Mki67_CAF","Mt1_CAF")],scale="row",cluster_rows = FALSE,cluster_cols = FALSE,border="white",gaps_row=c(3,6,8,11),filename="marker_heatmap.pdf")

####### Figure I Density plot
data=obj[["umap"]]@cell.embeddings
data <- as.data.frame(data)
data$Condition <- obj@meta.data$Condition
p<-ggplot()+geom_point() +
  stat_density_2d(data[data$Condition=="Control",],mapping=aes(x=UMAP_1, y=UMAP_2,fill=..level..) , geom = "polygon")+
  scale_fill_distiller(palette=7, direction=1)+
  theme_void()+
  scale_x_continuous(limits=c(min(data$UMAP_1) - 0.1*diff(range(data$UMAP_1)),
                              max(data$UMAP_1) + 0.1*diff(range(data$UMAP_1))))+
  scale_y_continuous(limits=c(min(data$UMAP_2) - 0.1*diff(range(data$UMAP_2)),
                              max(data$UMAP_2) + 0.1*diff(range(data$UMAP_2))))
ggsave("Density_Fibro_Ctrl.pdf",width=6,height=4.5)
p<-ggplot()+geom_point() +
  stat_density_2d(data[data$Condition=="KO",],mapping=aes(x=UMAP_1, y=UMAP_2,fill=..level..) , geom = "polygon")+
  scale_fill_distiller(palette=7, direction=1)+
  theme_void()+
  scale_x_continuous(limits=c(min(data$UMAP_1) - 0.1*diff(range(data$UMAP_1)),
                              max(data$UMAP_1) + 0.1*diff(range(data$UMAP_1))))+
  scale_y_continuous(limits=c(min(data$UMAP_2) - 0.1*diff(range(data$UMAP_2)),
                              max(data$UMAP_2) + 0.1*diff(range(data$UMAP_2))))
ggsave("Density_Fibro_KO.pdf",width=6,height=4.5)