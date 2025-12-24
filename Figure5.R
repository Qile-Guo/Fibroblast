#---------Figure5-------#


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
library(AUCell)

######## Figure I
obj<-readRDS("fibro_cell_B16_KPC_MC38.rds")
selected_genes<-read.csv("ProgB_ProgC_genes.csv")
# calculate gene expression
Idents(obj)<-'condition_model'
avg_exp <- AverageExpression(obj,assays = "RNA",features = selected_genes,group.by = "condition_model")$RNA
dat <- as.matrix(avg_exp) 
models <- unique(gsub("(.*)_.*", "\\1", colnames(dat)))
ratio_matrix <- sapply(models, function(model){
  ko_col <- paste0(model, "_KO")
  ctrl_col <- paste0(model, "_Ctrl")
  if(!all(c(ko_col, ctrl_col) %in% colnames(dat))){
    message(paste("缺少", model, "的配对数据"))
    return(rep(NA, nrow(dat)))
  }
  log2( (dat[, ko_col]) / (dat[, ctrl_col]) )
})
# visualization
rownames(ratio_matrix) <- rownames(dat)
ratio_matrix <- na.omit(ratio_matrix)
norm_data<-t(ratio_matrix )
pheatmap(
  mat = norm_data[c("B16", "KPC", "MC38"),],
  scale = "none",          # 行标准化（Z-score）
  cluster_rows = FALSE,  # 保持基因顺序
  cluster_cols = FALSE,  # 保持分组顺序
  border_color = NA,
  show_colnames = TRUE,
  main = "ECM Genes Expression",
  display_numbers = F
)

######## Figure J UMAP with density contour
#### Example: B16
library(ggnewscale)
B16_CAF<-readRDS("./B16_CAF_cleaned2_0613.RDS")

data<-B16_CAF@reductions$umap@cell.embeddings
data<-cbind(data,B16_CAF@meta.data)

dd<-subset(data,Condition=="Ctrl" )
dd<-subset(data,Condition=="KO"  )

col<-c("#86B8AE", "#C5D2C0" ,"#D8D8D8" ,"#DDB5A5" ,"#ACD16A" ,"#B1B4D4", "#7EAAC9")
p<-ggplot(dd, aes(x=umap_1, y=umap_2))+
  geom_point(data = subset(dd,cell_pro2 !=c("Apod+_Pi16+_Fib")),
             aes(color=cell_pro2),size=3)+
  geom_point(data = subset(dd,cell_pro2 %in% c("Apod+_Pi16+_Fib")),
             color=alpha(col[6],alpha=0.2),size=3)+
  scale_color_manual(values = alpha(col[c(1:5,7)],alpha=0.2))+ 
  new_scale_color()+
  stat_density_2d(aes(color=..level..))+
  scale_color_gradientn(colors = colorRampPalette(brewer.pal(name = "Greys",n=9)[c(2:8)])(200))+
theme_void()+
  scale_x_continuous(limits=c(min(data$umap_1) - 
                                0.1*diff(range(data$umap_1)),max(data$umap_1) + 0.1*diff(range(data$umap_1))))+
  scale_y_continuous(limits=c(min(data$umap_2) - 
                                0.1*diff(range(data$umap_2)),max(data$umap_2) + 0.1*diff(range(data$umap_2))))

p