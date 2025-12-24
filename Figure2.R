#---------Figure2-------#


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


###################### Figure 2B Neighborhood analysis
# data preparation
mtx <- read.csv("merfish_cc1.csv") #adata.obs from merfish data
mtx<-subset(mtx,mtx[,"cell_type"]=="Fibro")
mtx<-mtx[!(is.na(mtx$X)),]
calculate_neighbor_score <- function(data, sort_col, top_prop, distance_threshold, 
                                     score_col = "ISG_MMP", decreasing = TRUE) {
  if (decreasing) {
    sorted_data <- data[order(data[[sort_col]], decreasing = TRUE), ]
  } else {
    sorted_data <- data[order(data[[sort_col]], decreasing = FALSE), ]
  }
  n_cells <- nrow(sorted_data)
  selected_cells <- sorted_data[1:floor(n_cells * top_prop), ]
  neighbor_scores <- sapply(1:nrow(selected_cells), function(i) {
    distances <- sqrt(
      (data$x - selected_cells$x[i])^2 + 
        (data$y - selected_cells$y[i])^2
    )
    neighbors <- data[distances <= distance_threshold, ]
    if (nrow(neighbors) > 0) {
      return(mean(neighbors[[score_col]], na.rm = TRUE))
    } else {
      return(NA)
    }
  })
  data.frame(
    ScoreValue = neighbor_scores,
    Distance = distance_threshold,
    ScoreType = ifelse(decreasing, "Program C high", "Program C low"),
    stringsAsFactors = FALSE
  )
}
# set parameter
top_prop <- 0.5
distance_thresholds <- c(0, 10, 20, 30, 40, 50, 60)
results <- lapply(distance_thresholds, function(dist) {
  high_scores <- calculate_neighbor_score(
    mtx, "Myo_Co", top_prop, dist, 
    score_col = "ISG_MMP", decreasing = TRUE
  )
  
  low_scores <- calculate_neighbor_score(
    mtx, "Myo_Co", top_prop, dist, 
    score_col = "ISG_MMP", decreasing = FALSE
  )
  rbind(high_scores, low_scores)
})
df <- do.call(rbind, results)
df$ScoreValue <- as.numeric(df$ScoreValue)
df <- df[!is.na(df$ScoreValue), ]  # 移除NA值
# summary
summary_data <- df %>%
  group_by(Distance, ScoreType) %>%
  summarise(
    Mean = mean(ScoreValue, na.rm = TRUE),
    SE = sd(ScoreValue, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )
t_test_results <- df %>%
  group_by(Distance) %>%
  summarise(
    p_value = tryCatch({
      t.test(ScoreValue ~ ScoreType)$p.value
    }, error = function(e) NA),
    .groups = "drop"
  )
combined_data <- summary_data %>%
  left_join(t_test_results, by = "Distance")
combined_data$Distance <- factor(
  combined_data$Distance, 
  levels = sort(unique(combined_data$Distance))
)
# visualization
p <- ggplot(combined_data, aes(x = Distance, y = Mean, color = ScoreType, group = ScoreType)) +
  geom_point(position = position_dodge(width = 0.3), size = 4) +
  geom_errorbar(
    aes(ymin = Mean - SE, ymax = Mean + SE),
    position = position_dodge(width = 0.3),
    width = 0.2
  ) +
  geom_line(position = position_dodge(width = 0.3), size = 1) +
  labs(
    x = "Distance",
    y = "Neighbourhood Program_B Score",
    title = "Program C High vs Low: Neighbourhood Analysis"
  ) +
  scale_color_manual(values = c('#b10e12', '#a2d8e4')) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_blank(),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, size = 16)
  )
print(p)


######## Regulation correlation between programs

load("effect_prog_top_CRISPRi.Rdata")
load("effect_prog_top_CRISPRa.Rdata")

color<-pal_npg()(10)

p1<-ggplot(data =effect_prog_top_CRISPRi,
           aes(x=Effect_myo,y=Effect_isg))+
  geom_point(aes(color=color))+
  theme_linedraw()+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  scale_color_manual(values=c(color[c(1,4)],"gray"))+
  labs(x="Perturbation effects on ProgC",
       y="Perturbation effects on ProgB")+
  theme(plot.title = element_text(hjust = 0.5,size=24), 
        legend.position="none", 
        axis.title.x =element_text(size=10), 
        axis.title.y=element_text(size=10),
        axis.text=element_text(size=10)) +
  geom_smooth(method = 'lm', formula = y ~ x, 
              se = T,color="dimgray",fill="dimgray",alpha = 0.1,linewidth=0.5)+ 
  stat_cor(data=effect_prog_top_CRISPRi, method = "pearson")

p1

p2<-ggplot(data =effect_prog_top_CRISPRa,
           aes(x=Effect_myo,y=Effect_isg))+
  geom_point(aes(color=color))+
  theme_linedraw()+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  scale_color_manual(values=c(color[c(1,4)],"gray"))+
  labs(x="Perturbation effects on ProgC",
       y="Perturbation effects on ProgB")+
  theme(plot.title = element_text(hjust = 0.5,size=24), 
        legend.position="none", 
        axis.title.x =element_text(size=10), 
        axis.title.y=element_text(size=10),
        axis.text=element_text(size=10)) +
  geom_smooth(method = 'lm', formula = y ~ x, 
              se = T,color="dimgray",fill="dimgray",alpha = 0.1,linewidth=0.5)+ 
  stat_cor(data=effect_prog_top_CRISPRa, method = "pearson")

p2

######## Bulk RNA-seq analysis

########## read STAR results & Preprocess

ct_file =  "countMatrix_featureCounts.tsv"
ct_summary = "countMatrix_featureCounts.tsv.summary"


ct = readr::read_tsv(ct_file, comment = "#")
ct_qc = readr::read_tsv(ct_summary)

colnames(ct_qc)[2:ncol(ct_qc)] <- sapply(strsplit(colnames(ct_qc)[2:ncol(ct_qc)], split = "\\/"), `[[`,9)
features <- ct_qc[[1]]
ct_qc <- as.matrix(ct_qc[, 2:ncol(ct_qc)])
rownames(ct_qc) <- features

ct_qc <- as.data.frame(t(ct_qc)) %>% 
  rownames_to_column("SampleName")

ct_qc <- ct_qc %>%
  rowwise() %>% 
  mutate(Unassigned = sum(c_across(Unassigned_Unmapped:Unassigned_Ambiguity))) %>%
  mutate(mapping_rate = Assigned/(Assigned+Unassigned))

genemeta <- ct[, 1:6]
ct <- ct[, 7:ncol(ct)]
colnames(ct) <- sapply(strsplit(colnames(ct), split = "\\/"), `[[`, 9)
ct <- as.matrix(ct)
rownames(ct) <- genemeta$Geneid


genemap <- readr::read_tsv("gene-map-cellranger-hg38.tsv")
genemap <- genemap %>% 
  as.data.frame() %>% 
  column_to_rownames('GENEID')

ov_genes <- intersect(rownames(genemap), rownames(ct))
genemap <- genemap[ov_genes,]

rownames(ct_qc)<-ct_qc$SampleName
ct_qc<-as.data.frame(ct_qc)
rownames(ct_qc)<-ct_qc$SampleName
ct <- ct[rownames(genemap),]

rownames(genemeta)<-genemeta$Geneid
genemeta<-genemeta[rownames(genemap),]

sampleinfo<-read.csv("sample_info.csv")
rownames(sampleinfo)<-sampleinfo$Sample
sample_analysed<-subset(sampleinfo,Condition=="CRC-0111-T")

ov_samples<-intersect(rownames(sample_analysed),colnames(ct))
ct_remo = ct[, ov_samples]


sample_analysed<-sample_analysed[ov_samples,]
ct_qc<-ct_qc[ov_samples,]

keep <- rowSums(ct_remo>0) >= floor(0.75*ncol(ct_remo))
ct_remo <- ct_remo[keep,] 

genemap<-genemap[rownames(ct_remo),]
genemeta<-as.data.frame(genemeta)
rownames(genemeta)<-genemeta$Geneid
genemeta<-genemeta[rownames(ct_remo),]


########### Deseq2
library(DESeq2)
library(dplyr)

rownames(ct_remo)<-genemap$SYMBOL

dds <-DESeqDataSetFromMatrix(countData = ct_remo,
                             colData = sample_analysed,
                             design = ~ Treatment)

dds <- dds[rowSums(counts(dds)) > 1,]  
dep <- DESeq(dds)
res <- results(dep, contrast = c('Treatment', 'Tgfb', 'DMSO'))


diff = res
diff <- na.omit(diff)  
diff <- data.frame(diff, stringsAsFactors = FALSE, check.names = FALSE)
diff$Gene<-rownames(diff)

res1<-diff
res1[which(res1$log2FoldChange >= 1 & res1$padj < 0.05),'sig'] <- 'up'
res1[which(res1$log2FoldChange <= -1 & res1$padj < 0.05),'sig'] <- 'down'
res1[which(abs(res1$log2FoldChange) <= 1 | res1$padj >= 0.05),'sig'] <- 'none'
table(res1$sig)

################## GSVA
library(GSVA)

load("NMF_module.Rdata")

dat = ct_remo
GSVA_result <- gsva(expr=dat, 
                    gset.idx.list=NMF_module, 
                    kcdf="Poisson" ,
                    verbose=T, 
                    parallel.sz = 10)


################### Radar plot
library(fmsb)

load("score_module_df.Rdata")

color<-c("#4089B6", "#CC5C4A")

score_module_df[4,]<-1
score_module_df[5,]<-0
rownames(score_module_df)[c(4,5)]<-c("Max","Min")


radarchart(
  score_module_df[c(4,5,1:3),], 
  axistype = 1,
  pcol =  c("darkgray",color), 
  pfcol = scales::alpha(c("darkgray",color),0.5), 
  plwd = 1, plty = 1,
  cglcol = "grey", 
  cglty = 1, 
  cglwd = 0.8,
  axislabcol = "grey", 
  vlcex = 1, 
  vlabels = colnames(score_module_df),
  caxislabels = c(NA,NA,NA,NA,NA)
)


###################### Figure 2K & L

#### Example: CRISPRi, ProgB

effect_on_prog<-read.csv("./Effect_on_programme.csv")
effect_on_prog_CRISPRi<-subset(effect_on_prog,Library =="CRISPRi")
colnames(effect_on_prog_CRISPRi)[1:3]<-c("Perturbation","Program","Effect")

effect_on_prog_CRISPRi_m<-reshape2::dcast(effect_on_prog_CRISPRi, 
                                          Program~Perturbation,
                                          value.var = 'Effect')

rownames(effect_on_prog_CRISPRi_m)<-effect_on_prog_CRISPRi_m[,1]
effect_on_prog_CRISPRi_m<-effect_on_prog_CRISPRi_m[,-1]
effect_on_prog_CRISPRi_m<-as.data.frame(t(effect_on_prog_CRISPRi_m))


y<-effect_on_prog_CRISPRi_m[order(effect_on_prog_CRISPRi_m$`Interferon-responsive`,decreasing = T),]
y$Effect<-y$`Interferon-responsive`
y$target<-rownames(y)

y$Rank<-1:432
y$color="No"
y$color[1:20]<-"Top20 upregulate targets"
y$color[413:432]<-"Top20 downregulate targets"


color<-pal_npg()(10)

gg<-c("ADAM12","COL10A1","SPOCK1","CXCL8","NOX4","SGK1","TGFBR2")

ggplot(data = y,aes(x=Rank,y=Effect))+geom_point(aes(color=color),size=1)+theme_linedraw()+
  theme(panel.grid = element_blank(),legend.position = "none")+
  scale_color_manual(values=c("gray",color[c(4)],color[c(1)]))+
  labs(y="Perturbation Effects")+
  geom_hline(yintercept= 0,lty=2,col="black",lwd=0.4)+
  theme(plot.title = element_text(hjust = 0.5,size=15), 
        legend.position="right", 
        axis.title.x =element_text(size=12), 
        axis.title.y=element_text(size=12),
        axis.text=element_text(size=12),
        legend.title=element_blank(),
        legend.text = element_text(siz=12)
        # legend.key.size = unit(20,"mm")
  )+guides(color = guide_legend(override.aes = list(size = 4)))+
  geom_text_repel(data = y[gg,],
                  aes(x=Rank,y=Effect,label=target),
                  size = 3,color="black",
                  box.padding = unit(0.2, "lines"), 
                  segment.color = "black",   
                  segment.size = 0.2 
  )


