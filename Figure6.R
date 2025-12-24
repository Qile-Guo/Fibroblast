#---------Figure6-------#


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

######## GO barplot

GO_up_plot<-read.csv("CAF_up_100.csv")
GO_down_plot<-read.csv("CAF_down_100.csv")

GO_plot<-rbind(GO_up_plot,GO_down_plot)

p<-ggplot(GO_plot, aes(gene_ratio, term_name, 
                       fill = negative_log10_of_adjusted_p_value)) +
  geom_bar(stat="identity",alpha=0.8) +
  labs(y = "", x = "Gene Ratio",fill="-Log10(Padj)")+
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0,size=12,face = "bold"), 
        legend.position="bottom", 
        axis.title.x =element_text(size=12), 
        axis.title.y=element_text(size=12),
        axis.ticks.y = element_blank(),
        axis.text=element_text(size=12),
        legend.title = element_text(size=12),
        panel.grid = element_blank(),
        axis.text.y = element_blank())+  
  scale_fill_gradientn(colours = colorRampPalette(brewer.pal(n=9,"Blues")[2:9])(50))+  
  geom_text(aes(label = term_name), x=0,size = 4.5,hjust = 0)

p

################### Umap plot + density plot

load("Tcell_cleaned.Rdata")

T_color<-c("#85BBB1","#F8A979","#4891C2","#D4C1DE","#88638B","#8663AC")

DimPlot(Tcell_cleaned,group.by = "Anno",cols = T_color,alpha = 0.5,label = F,pt.size = 0.5)

data<-Tcell_cleaned@reductions$umap@cell.embeddings
data<-cbind(data,Tcell_cleaned$Condition)
data<-as.data.frame(data)

data$umap_1<-as.numeric(data$umap_1)
data$umap_2<-as.numeric(data$umap_2)

p<-ggplot(data[data$V3=="KO",], aes(x=umap_1, y=umap_2) ) +
  stat_density_2d(aes(fill =  after_stat(level)), geom = "polygon")+
  scale_fill_distiller(palette=7, direction=1)+
  theme_void()+
  scale_x_continuous(limits=c(min(data$umap_1) - 0.1*diff(range(data$umap_1)),
                              max(data$umap_1) + 0.1*diff(range(data$umap_1))))+
  scale_y_continuous(limits=c(min(data$umap_2) - 0.1*diff(range(data$umap_2)),
                              max(data$umap_2) + 0.1*diff(range(data$umap_2))))+
  theme( legend.position="right")


p

p<-ggplot(data[data$V3=="WT",], aes(x=umap_1, y=umap_2) ) +
  stat_density_2d(aes(fill =  after_stat(level)), geom = "polygon")+
  scale_fill_distiller(palette=7, direction=1)+
  theme_void()+
  scale_x_continuous(limits=c(min(data$umap_1) - 0.1*diff(range(data$umap_1)),
                              max(data$umap_1) + 0.1*diff(range(data$umap_1))))+
  scale_y_continuous(limits=c(min(data$umap_2) - 0.1*diff(range(data$umap_2)),
                              max(data$umap_2) + 0.1*diff(range(data$umap_2))))+
  theme( legend.position="right")


p


################# Gene Expression

CD8T<-subset(Tcell_cleaned,Anno %in% c("T3_CD8T_Mki67", "T4_CD8T_Ccl4", "T5_CD8T_Mcm3"))

CD8T$Condition<-factor(CD8T$Condition,c("WT","KO"))
VlnPlot(CD8T,features = "Ifng",group.by = "Condition",cols = c("#7491B8","#D7605A"))
VlnPlot(CD8T,features = "Cxcr6",group.by = "Condition",cols = c("#7491B8","#D7605A"))


################# CellChat
library(CellChat)
load("Immune_integrated_cleaned.Rdata")

Immune_WT<-subset(Immune_integrated_cleaned, Condition =="WT")

cellchat_immune_WT <- createCellChat(object = Immune_WT,
                                     meta = Immune_WT@meta.data,
                                     group.by = "assign")

CellChatDB <- CellChatDB.mouse  
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB,search = "Secreted Signaling") 
cellchat_immune_WT@DB <- CellChatDB.use

cellchat_immune_WT <- subsetData(cellchat_immune_WT) 

cellchat_immune_WT <- identifyOverExpressedGenes(cellchat_immune_WT)
cellchat_immune_WT <- identifyOverExpressedInteractions(cellchat_immune_WT)
cellchat_immune_WT <- computeCommunProb(cellchat_immune_WT,type =  "truncatedMean", raw.use = TRUE, population.size = TRUE) 

cellchat_immune_WT <- filterCommunication(cellchat_immune_WT, min.cells = 10)
cellchat_immune_WT <- computeCommunProbPathway(cellchat_immune_WT)

cellchat_immune_WT <- aggregateNet(cellchat_immune_WT)


Immune_KO<-subset(Immune_integrated_cleaned, Condition =="KO" )

cellchat_immune_KO <- createCellChat(object = Immune_KO,
                                     meta = Immune_KO@meta.data,
                                     group.by = "assign")

CellChatDB <- CellChatDB.mouse  
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB,search = "Secreted Signaling") 
cellchat_immune_KO@DB <- CellChatDB.use

cellchat_immune_KO <- subsetData(cellchat_immune_KO) 

cellchat_immune_KO <- identifyOverExpressedGenes(cellchat_immune_KO)
cellchat_immune_KO <- identifyOverExpressedInteractions(cellchat_immune_KO)
cellchat_immune_KO <- computeCommunProb(cellchat_immune_KO, type =  "truncatedMean",raw.use = T, population.size = TRUE) 

cellchat_immune_KO <- filterCommunication(cellchat_immune_KO, min.cells = 10)
cellchat_immune_KO <- computeCommunProbPathway(cellchat_immune_KO)

cellchat_immune_KO <- aggregateNet(cellchat_immune_KO)


cellchat_immune_WT<-netAnalysis_computeCentrality(cellchat_immune_WT)
cellchat_immune_KO<-netAnalysis_computeCentrality(cellchat_immune_KO)

###### Merge
object.list <- list( WT = cellchat_immune_WT,KO = cellchat_immune_KO)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))


p<-plotGeneExpression(cellchat, signaling = "IFN-II", split.by = "datasets", 
                      colors.ggplot = T,color.use = c("#7491B8","#D7605A"),group.by = "assign")

p

p1<-plotGeneExpression(cellchat, signaling = "CXCL", split.by = "datasets", 
                       colors.ggplot = T,color.use = c("#7491B8","#D7605A"),group.by = "assign",features = c("Cxcr6","Cxcl16"))
p1



pathways.show="IFN-II"
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) 

netVisual_aggregate(cellchat_immune_KO, 
                    signaling = "IFN-II",
                    layout = "hierarchy" ,
                    edge.weight.max = weight.max,
                    small.gap = 0.1,big.gap = 0.1,
                    vertex.receiver = c(1,6),
                    color.use = c("#D6604D",pal_npg()(10),"gray"),
                    thresh = 0.05)

netVisual_aggregate(cellchat_immune_WT, 
                    signaling = "IFN-II",
                    layout = "hierarchy" ,
                    edge.weight.max = weight.max,
                    small.gap = 0.1,big.gap = 0.1,
                    vertex.receiver = c(1,6),
                    color.use = c("#4393C3",pal_npg()(10),"gray"),
                    thresh = 0.05)




