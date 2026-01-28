#---------Figure3-------#


# Library loading

library(tidyr)
library(data.table)
library(dplyr)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(Seurat)
library(ggsci)
library(ggrepel)
library(rstatix)
library(ggpubr) 
library(clusterProfiler)

######## top in silico genes performance in ex vivo screenings

screen_effect <- read.csv("./Screen_effect_i.csv")
effect_on_prog <- read.csv("./Effect_on_programme.csv")

effect_on_prog_crispri <- effect_on_prog %>%
  filter(Library == "CRISPRi") %>%
  select(Perturbation = 1, Program = 2, Effect = 3)

top20_genes <- c('ADAM12', 'NTM', 'CTHRC1', 'FAP', 'TNFAIP6', 'LOX')

# Reshape for heatmap
top20_matrix <- effect_on_prog_crispri %>%
  filter(Perturbation %in% top20_genes) %>%
  pivot_wider(names_from = Program, values_from = Effect) %>%
  column_to_rownames("Perturbation")

# Reorder programs and transpose
top20_matrix <- t(top20_matrix[, c(1, 3, 5, 4, 2)])
rownames(top20_matrix) <- c("P-A", "P-B", "P-C", "P-D", "P-E")

color_palette <- brewer.pal(n = 11, name = "RdBu")
pheatmap(top20_matrix[, top20_genes], 
         cluster_cols = FALSE, 
         cluster_rows = FALSE,
         color = rev(colorRampPalette(c(color_palette[1:6], "white", color_palette[7:11]))(256)),
         angle_col = "90",
         border_color = "gray",
         main = "CRISPRi Effect on Programs")

################## Figure 3C Expression correlation
# read genes
genelist1<-read.csv("Top_PBPC_genes.csv")
genelist2<-genelist1$gene
igenes<-subset(genelist1,library=="CRISPRi")$gene
# read in vivo fibroblast
obj<-readRDS("Invivo_fibro_cluster.rds")
obj<-subset(obj,(annotation!='6 Pericyte')&(annotation!='7 SMC')) # remove fibro-like cells
mp<-read.csv("Program_gene_list.csv")
n<-ncol(obj@meta.data)
obj<-AddModuleScore(obj,list(na.omit(mp[,"ProgramB"])),name = "ProgramB")
obj<-AddModuleScore(obj,list(na.omit(mp[,"ProgramC"])),name = "ProgramC")
colnames(obj@meta.data)<-c(colnames(obj@meta.data[,1:n]),"ProgramB","ProgramC")
# read DEG files
combine<-read.csv("In_vivo_dataset_tumor_vs_normal_DEG.csv")
genelist<-rownames(combine)
# correlation
data = t(as.data.frame(obj@assays$RNA@data[genelist,]))
data_cor<-matrix(,nrow=length(genelist),ncol=2)
data_cor<-as.data.frame(data_cor)
rownames(data_cor)<-genelist
colnames(data_cor)<-c('PB_cor','PC_cor')
for (i in genelist){
  c1<-cor.test(data[,i],obj@meta.data[,'ProgramB'])$estimate 
  data_cor[i,'PB_cor']<-c1
  c2<-cor.test(data[,i],obj@meta.data[,'ProgramC'])$estimate
  data_cor[i,'PC_cor']<-c2
}
data_cor_sub<-data_cor[genelist2,]
combine_sub<-combine[genelist2,]
data_new<-merge(data_cor_sub,combine_sub,by="row.names")
rownames(data_new)<-data_new$Row.names
# visualization
data_new$Gene<-rownames(data_new)
data_new$library<-ifelse(data_new$Gene %in% igenes,"CRISPRi","CRISPRa")
label_genes<-c("ADAM12","COL10A1","LTBP1","MFAP2","SPOCK1","CXCL8","NOX4","SGK1","PLAC9","TCF21","LGALS3","PAMR1","BMP4","FOSB","ECM1","ANXA1","PBX1")
data_new$Label<-ifelse(data_new$Gene %in% label_genes,data_new$Gene,"")
p1<-ggplot(data_new)+
  geom_point(aes(x=PB_cor,y=PC_cor,color=as.numeric(avg_log2FC),shape=library,size=abs(as.numeric(avg_log2FC))))+
  scale_color_gradient2("Log2FC\n[Tumor versus Normal]",high='#e16e08',mid='#ececec',low='#0a6801',midpoint=0)+
  geom_hline(aes(yintercept=0),colour="black", linetype="dashed") +
  geom_vline(aes(xintercept=0), colour="black", linetype="dashed")+
  geom_text_repel(aes(x=PB_cor,y=PC_cor,label=Label),max.overlaps=300)+
  scale_size(range=c(3,8))+
  theme_bw()
ggsave("PB_PC_cor_top_genes.pdf",width=8,height=5)


####################### ADAM12 isoform TPM
isoform_data <- read.csv("./isoform_data.csv")
sample_info <- read.csv("./human_bulk_0313.csv") %>% rename(sample = 1)

plot_data <- isoform_data %>%
  left_join(sample_info, by = "sample") %>%
  filter(Condition == "CRC-0111-T", 
         gene_label %in% c("ADAM12-L", "ADAM12-S"),
         Treatment %in% c("DMSO", "TGFb1")) %>%
  mutate(log_tpm = log10(tpm + 1),
         gene_label = factor(gene_label, levels = c("ADAM12-S", "ADAM12-L")))

# Statistics
stat_test <- plot_data %>%
  group_by(Treatment) %>%
  t_test(log_tpm ~ gene_label) %>%
  add_xy_position(x = "Treatment", dodge = 0.8)

# Visualization
ggplot(plot_data, aes(x = Treatment, y = log_tpm, fill = gene_label)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.7, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(shape = 21, color = "black", 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              size = 2.5) +
  stat_pvalue_manual(stat_test, label = "p = {p}", tip.length = 0.01) +
  scale_fill_manual(values = alpha(rev(c("#99C7E1", "#FAD2AA")), 0.8)) +
  theme_bw(base_size = 14) +
  labs(x="",y="Log10(TPM+1)") +
  theme(plot.title = element_text(hjust = 0.5,size=24), 
        legend.text = element_text(size=12,colour = "black"), 
        axis.title.x =element_text(size=13,colour = "black"), 
        axis.title.y=element_text(size=13,colour = "black"),
        axis.text=element_text(size=13,colour = "black"))

##################### Volcano plot

hvg_gene<-read.csv("Bulk_DGE.csv")

hvg_gene[which(hvg_gene$log2FoldChange >= 0.5 & hvg_gene$padj < 0.05),'sig'] <- 'up'
hvg_gene[which(hvg_gene$log2FoldChange <= -0.5 & hvg_gene$padj < 0.05),'sig'] <- 'down'
hvg_gene[which(abs(hvg_gene$log2FoldChange) <= 0.5 | hvg_gene$padj >= 0.05),'sig'] <- 'none'


hvg_gene<-as.data.frame(hvg_gene)
hvg_gene<-cbind(hvg_gene,-log10(hvg_gene$padj))
colnames(hvg_gene)[ncol(hvg_gene)]<-"log_p_val_adj"

color<-brewer.pal(n=11,"RdBu")

p <- ggplot(
  hvg_gene, aes(x = log2FoldChange, y =log_p_val_adj, colour=sig)) +
  geom_point(size=0.6) +
  scale_color_manual(values=c(color[9],"gray",color[3]))+
  geom_vline(xintercept=c(-0.5,0.5),lty=3,col="#ADADAD",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=3,col="#ADADAD",lwd=0.8) +
  labs(x="Log2(Fold Change)",y="-Log10(adj.P-value)")+   
  theme_bw()+   
  theme(panel.grid = element_blank())+ 
  theme(plot.title = element_text(hjust = 0.5,size=24), 
        legend.position="none", 
        axis.title.x =element_text(size=11), 
        axis.title.y=element_text(size=11),
        axis.text=element_text(size=11))

p


gene<-c("ADAM12","IFIT1","IFIT2","IFIT3","CXCL10","COL1A1","COL12A1","LOXL2","TAGLN","COL3A1")

p1<-p+ geom_text_repel(data = hvg_gene[gene,], 
                       aes(x = log2FoldChange, y = log_p_val_adj, label = Gene),
                       size = 3.5,color="black",
                       box.padding = unit(0.1, "lines"), 
                       segment.color = "black",   #连线的颜色
                       segment.size = 0.1,  #连线的粗细
)
p1


################## GSEA
library(GSEA)
library(enrichplot)

# HALLMARK_INTERFERON_ALPHA_RESPONSE and the pan-fibroblast TGFβ response signature (Pan-F-TBRS)
geneset_human_bulk<-read.csv("./geneset_human_bulk.csv") 

down<-subset(hvg_gene,log2FoldChange < 0)
up<-subset(hvg_gene,log2FoldChange> 0 )
down<-down[order(down$log2FoldChange,decreasing = F),]
up<-up[order(up$log2FoldChange,decreasing = T),]

geneList<-c(up$log2FoldChange,rev(down$log2FoldChange))
names(geneList)<-c(up$Gene,rev(down$Gene))


egmt <- GSEA(geneList, TERM2GENE=geneset_human_bulk, 
             verbose=FALSE,  pvalueCutoff=1)
gsea_results <- egmt@result
gsea_results

p<-gseaplot2(egmt, geneSetID = c("HALLMARK_INTERFERON_ALPHA_RESPONSE",
                                 "F_TBRS"), 
             base_size = 12,pvalue_table = F,
             color = c("#853021","#7DB0AB"))

p





