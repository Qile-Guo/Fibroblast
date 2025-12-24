#---------Figure1-------#


# Library loading

library(tidyr)
library(data.table)
library(pheatmap)
library(dplyr)
library(reshape2)
library(ggplot2)
library(Seurat)
library(RColorBrewer)
library(scatterplot3d)
library(ggpubr)


########### Figure 1B; Expression correlation
process_dataset <- function(path, signature, test_genes) {
    # file path
    info <- unlist(strsplit(path, "/"))
    organ <- info[2]
    dataset <- info[3]
    folder_path <- paste0("/insilico_datasets/", organ, "/", dataset)
    path1 <- list.files(path = folder_path, pattern = "\\.mtx$", full.names = TRUE)[1]
    path2 <- list.files(path = folder_path, pattern = "\\.txt$", full.names = TRUE)[1]
    path3 <- list.files(path = folder_path, pattern = "Cells.csv", full.names = TRUE)[1]
    # read data
    mtx <- Matrix::readMM(path1)
    genes <- fread(path2, header = FALSE)
    cells <- fread(path3)
    rownames(mtx) <- genes$V1
    colnames(mtx) <- cells$cell_name
    obj <- CreateSeuratObject(counts = mtx, project = dataset, assay = "RNA")
    obj <- AddMetaData(obj, cells)
    obj <- subset(obj, cell_type == "Fibroblast")
    # preprocessing
    obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
    return(obj)
}
calculate_correlations <- function(obj, test_genes, signature, signature_name, dataset_name) {
    # get expression data
    gene_expression_matrix <- GetAssayData(obj, slot = "data")
    valid_genes <- test_genes[test_genes %in% rownames(gene_expression_matrix)]
    gene_expression_data <- as.data.frame(t(as.matrix(gene_expression_matrix[valid_genes, ])))
    obj@meta.data <- cbind(obj@meta.data, gene_expression_data)
    # signature genes
    sig_genes <- unlist(strsplit(signature[signature_name, 2], ", "))
    # signature scores
    obj <- AddModuleScore(obj, features = list(sig_genes), name = "sig_score")
    # correlation
    meta_data <- obj@meta.data
    columns_to_correlate <- meta_data[, valid_genes]
    sig_scores <- meta_data[, "sig_score1"]
    rho_df <- data.frame(
        rho = numeric(length(valid_genes)),
        p_value = numeric(length(valid_genes)),
        gene = character(length(valid_genes)),
        dataset = character(length(valid_genes)),
        signature = character(length(valid_genes)),
        stringsAsFactors = FALSE
    )
    for (i in seq_along(valid_genes)) {
        gene <- valid_genes[i]
        gene_expression <- columns_to_correlate[, gene]
        result <- tryCatch({
            cor.test(gene_expression, sig_scores, method = "spearman", exact = FALSE)
        }, error = function(e) {
            list(estimate = NA, p.value = NA)
        })
        rho_df[i, ] <- c(
            result$estimate,
            result$p.value,
            gene,
            dataset_name,
            signature_name
        )
    }
    return(rho_df)
}
combine_results <- function(rho_df_all, signature_filter) {
    filtered_df <- subset(rho_df_all, signature == signature_filter)
    genes <- unique(na.omit(filtered_df$gene))
    combine_df <- data.frame(
        gene = genes,
        combine_rho = numeric(length(genes)),
        combine_p = numeric(length(genes)),
        stringsAsFactors = FALSE
    )
    for (i in seq_along(genes)) {
        gene <- genes[i]
        sub_table <- subset(filtered_df, gene == gene)
        combine_df$combine_rho[i] <- median(na.omit(as.numeric(sub_table$rho)))
        p_values <- as.numeric(na.omit(sub_table$p_value))
        if (length(p_values) > 0) {
            p_values[p_values < 1e-300] <- 1e-300
            X_squared <- -2 * sum(log(p_values))
            combine_df$combine_p[i] <- pchisq(X_squared, df = 2 * length(p_values), lower.tail = FALSE)
        } else {
            combine_df$combine_p[i] <- NA
        }
    }
    return(combine_df)
}
# read data
dataset <- read.csv("./insilico_screen_dataset.csv")
paths <- dataset$Path
signature <- read.csv("fibroblast_related_gene_signature.csv", row.names = 1)
test_genes <- read.csv("priority_target_list.csv")$gene
test_genes_df <- data.frame(gene = test_genes)
rho_df_all <- data.frame()
# for each dataset
for (path in paths) {
      cat("Processing:", path, "\n")
      obj <- process_dataset(path, signature, test_genes)
      sig_names <- rownames(signature)
      dataset_name <- basename(dirname(path))
      for (sig_name in sig_names) {
      rho_df <- calculate_correlations(obj, test_genes, signature, sig_name, dataset_name)
      if (!is.null(rho_df)) {
            rho_df_all <- rbind(rho_df_all, rho_df)
      }
      }
}
combine1 <- combine_results(rho_df_all, 'LRRC15+ Fib [Buechler et al.]')
test_genes_df$LRRC15_Fib <- combine1$combine_rho[match(test_genes_df$gene, combine1$gene)]
combine2 <- combine_results(rho_df_all, 'Fibroblast activation')
test_genes_df$Fibroblast_activation <- combine2$combine_rho[match(test_genes_df$gene, combine2$gene)]
# plot
highlight_genes <- c("CTHRC1", "INHBA", "COL11A1", "ADAM12", "NTM", "FAP", "HAS2", "EDNRA", "LOXL2")
test_genes_df$label <- ifelse(test_genes_df$gene %in% highlight_genes, test_genes_df$gene, "")
p <- ggscatterhist(
      test_genes_df, x = 'LRRC15_Fib', y = 'Fibroblast_activation', 
      color = '#014055', size = 3, alpha = 0.7,
      margin.params = list(fill = "grey", color = "black", size = 0.2),
      label = 'label',
      repel = TRUE,
      label.select = highlight_genes,
      font.label = c(12, "plain", "black"),
      ggtheme = theme_minimal()
)
p <- p + stat_cor(method = "spearman") + stat_smooth(method = "lm", se = TRUE, color = "red", alpha = 0.2)
ggsave("priority_myCAFsignature_cor.pdf", plot = p, width = 7, height = 6)

########### Figure 1C; Target rank
pdf("priority_3Dscatter.pdf", width=8, height=8)
itg_rank<-read.csv("priority_list.csv")
subcellular_shapes <- c(17,16,18)
names(subcellular_shapes) <- unique(itg_rank$subcellular_location)
itg_rank$pch <- subcellular_shapes[itg_rank$subcellular_location]
rank_colors <- colorRampPalette(c("#440000","#e9c17f", "#c7c7c7"))(54) 
itg_rank$pcolor <- rank_colors[rank(itg_rank$overall_rank)]
with(itg_rank, {
   s3d <- scatterplot3d(norm_expr, -auc_score, HR_score,     
                        color=pcolor,                 
                        pch=pch,                    
                        type="h",                    
                        main="CAF targets",
                        xlab="normal_expr",
                        ylab="HR",
                        zlab="AUC score",
                        cex.symbols = 2,
                        xlim=c(0.5, 3),
                        zlim=c(1.1, 1.5),
                        ylim=c(-0.9, -0.5))
    s3d.coords <- s3d$xyz.convert(norm_expr, -auc_score, HR_score)
    text(s3d.coords$x, s3d.coords$y, 
         labels=itg_rank$label, 
         cex=.9, pos=4)
})
dev.off()

########### Figure 1F; Coregulated programme

effect_on_gene_CRISPRi<-read.csv("effect_on_gene_CRISPRi.csv")
effect_on_gene_CRISPRa<-read.csv("effect_on_gene_CRISPRa.csv")
load("NMF_gene_filter.Rdata")

effect_on_gene_CRISPRi<-reshape2::dcast(effect_on_gene_CRISPRi, 
                                        Gene~Perturbation,
                                        value.var = 'Effect')

rownames(effect_on_gene_CRISPRi)<-effect_on_gene_CRISPRi$Gene
effect_on_gene_CRISPRi<-effect_on_gene_CRISPRi[,-1]

NMF_gene_cor_i <- cor(t(effect_on_gene_CRISPRi), method = 'spearman')
NMF_gene_cor_i<-NMF_gene_cor_i[NMF_gene_filter,NMF_gene_filter]
NMF_gene_cor_i[NMF_gene_cor_i == 1]=NA

color=brewer.pal(n=11,"RdYlBu")[c(2:11)]


pheatmap(NMF_gene_cor_i,
         cluster_cols = F,
         cluster_rows = F,
         scale = "none",color = rev(colorRampPalette(color)(200)))


effect_on_gene_CRISPRa<-reshape2::dcast(effect_on_gene_CRISPRa, 
                              Gene~Perturbation,
                              value.var = 'Effect')

rownames(effect_on_gene_CRISPRa)<-effect_on_gene_CRISPRa$Gene
effect_on_gene_CRISPRa<-effect_on_gene_CRISPRa[,-1]

NMF_gene_cor_a <- cor(t(effect_on_gene_CRISPRa), method = 'spearman')
NMF_gene_cor_a<-NMF_gene_cor_a[NMF_gene_filter,NMF_gene_filter]
NMF_gene_cor_a[NMF_gene_cor_a == 1]=NA

pheatmap(NMF_gene_cor_a,
         cluster_cols = F,
         cluster_rows = F,
         scale = "none",color = rev(colorRampPalette(color)(200)))


############# Figure 1H; GO barplot

Prog_GO<-read.csv("Programme_GO.csv")

p<-ggplot(data=Prog_GO,
          aes(x= Gene_ratio, y= term_name))+
  geom_bar(stat="identity",fill="gray") +
  labs(y = "", x = "Gene Ratio")+
  theme_linedraw()+
  theme(plot.title = element_text(hjust = 0.5,size=14,face = "bold"), 
        legend.position="none", 
        axis.title.x =element_text(size=10), 
        axis.title.y=element_text(size=10),
        axis.text.x = element_text(angle=0, hjust=1, vjust=1,color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.text=element_text(size=10),
        legend.title = element_text(size=10),
        panel.grid = element_blank())+
  facet_grid(vars(Prog),scales="free",switch = "y")+
  scale_y_discrete(position = "right")

p

###################### cofunctional sets


effect_on_prog<-read.csv("./Effect_on_programme.csv")
effect_on_prog_CRISPRi<-subset(effect_on_prog,Library =="CRISPRi")
colnames(effect_on_prog_CRISPRi)[1:3]<-c("Perturbation","Program","Effect")


effect_on_prog_CRISPRi_m<-reshape2::dcast(effect_on_prog_CRISPRi, Program~Perturbation,
                                          value.var = 'Effect')
rownames(effect_on_prog_CRISPRi_m)<-effect_on_prog_CRISPRi_m$Program
effect_on_prog_CRISPRi_m<-effect_on_prog_CRISPRi_m[,-1]

effect_on_prog_CRISPRi_m<-as.data.frame(t(effect_on_prog_CRISPRi_m))

effect_on_prog_CRISPRi_sort<-subset(effect_on_prog_CRISPRi, Program =="myCAF-like")
effect_on_prog_CRISPRi_sort<-effect_on_prog_CRISPRi_sort[order(effect_on_prog_CRISPRi_sort$Effect),]


target_cor_f <- cor(t(effect_on_prog_CRISPRi_m[effect_on_prog_CRISPRi_sort$Perturbation[c(1:100)],1:4]), method = 'pearson')
p1<-pheatmap(target_cor_f,cluster_rows = T,cluster_cols = T)
p1


effect_on_prog_CRISPRa<-subset(effect_on_prog,Library =="CRISPRa")
colnames(effect_on_prog_CRISPRa)[1:3]<-c("Perturbation","Program","Effect")


effect_on_prog_CRISPRa_m<-reshape2::dcast(effect_on_prog_CRISPRa, Program~Perturbation,
                                          value.var = 'Effect')
rownames(effect_on_prog_CRISPRa_m)<-effect_on_prog_CRISPRa_m$Program
effect_on_prog_CRISPRa_m<-effect_on_prog_CRISPRa_m[,-1]

effect_on_prog_CRISPRa_m<-as.data.frame(t(effect_on_prog_CRISPRa_m))

effect_on_prog_CRISPRa_sort<-subset(effect_on_prog_CRISPRa, Program =="myCAF-like")
effect_on_prog_CRISPRa_sort<-effect_on_prog_CRISPRa_sort[order(effect_on_prog_CRISPRa_sort$Effect),]


target_cor_f <- cor(t(effect_on_prog_CRISPRa_m[effect_on_prog_CRISPRa_sort$Perturbation[c(1:100)],1:4]), method = 'pearson')
p1<-pheatmap(target_cor_f,cluster_rows = T,cluster_cols = T)
p1


