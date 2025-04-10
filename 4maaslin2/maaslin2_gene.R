rm(list = ls())
setwd('/Users/yangwanyi/Downloads/ENIGMA_STUDY_DATA/4maaslin2/')

# 加载MaAsLin2
library(Maaslin2)
library(dplyr)


# 输入数据
gene_data <- read.csv("../source_data/ARG68_gene_abundance_20250224.csv", row.names = 1)
rownames(gene_data) <- gene_data$gene_id
gene_name <- gene_data[,c(1,2)]
gene_data <- gene_data[,-c(1,2)]

gene_data <- t(gene_data)

keep <- apply(gene_data, 2, mean) > 1E-4 &
  apply(gene_data > 0, 2, sum) / nrow(gene_data) > 0.1 # mean abundance > 1E-4、Prevalence > 0.1
gene_data <- gene_data[, keep]

save(gene_data,file = 'gene_data.rds')

metadata <- read.csv('../3diversity/new_metadata.csv',row.names = 1)
metadata$Group <- ifelse(metadata$Group != 'Control','CD','Control')

fit_data <- Maaslin2(
  gene_data, metadata, 
  output = 'output_Group_gene',
  min_abundance = 0, 
  min_prevalence = 0, 
  normalization = "NONE",
  transform = "LOG",
  fixed_effects = c('Group'),
  reference = c("Group,Control"),
  standardize = FALSE, 
  cores=4)

metadata1 <- metadata
metadata1$Region <- ifelse(metadata$Region == 'AUS','AUS','other')

fit_data <- Maaslin2(
  gene_data, metadata1, 
  output = 'output_Region_gene_AUS',
  min_abundance = 0, 
  min_prevalence = 0, 
  normalization = "NONE",
  transform = "LOG",
  fixed_effects = c('Region'),
  reference = c("Region,other"),
  standardize = FALSE, 
  plot_heatmap = FALSE,
  plot_scatter = FALSE,
  cores=4)

metadata2 <- metadata
metadata2$Region <- ifelse(metadata$Region == 'HK','HK','other')

fit_data <- Maaslin2(
  gene_data, metadata2, 
  output = 'output_Region_gene_HK',
  min_abundance = 0, 
  min_prevalence = 0, 
  normalization = "NONE",
  transform = "LOG",
  fixed_effects = c('Region'),
  reference = c("Region,other"),
  standardize = FALSE, 
  plot_heatmap = FALSE,
  plot_scatter = FALSE,
  cores=4)


metadata3 <- metadata
metadata3$Region <- ifelse(metadata$Region == 'MC','MC','other')

fit_data <- Maaslin2(
  gene_data, metadata3, 
  output = 'output_Region_gene_MC',
  min_abundance = 0, 
  min_prevalence = 0, 
  normalization = "NONE",
  transform = "LOG",
  fixed_effects = c('Region'),
  reference = c("Region,other"),
  standardize = FALSE, 
  plot_heatmap = FALSE,
  plot_scatter = FALSE,
  cores=4)
# 
# 
# sig_res <- read.csv('output_Group_gene/significant_results.tsv',sep = '\t',header = T)
# sig_res$CD_status <- ifelse(sig_res$coef>0,'CD enriched','CD depleted')
# rownames(sig_res) <- sig_res$feature
# sig_res <- sig_res[,'CD_status',drop=FALSE]
# sig_res <- sig_res %>%
#   arrange(CD_status)
# 
# microbiome_data_df <- read.csv('output_Group_gene/features/filtered_data_norm_transformed.tsv',sep = '\t',header = T,row.names = 1)
# pheatmap_df <- microbiome_data_df[,rownames(sig_res)]
# 
# p1 <- pheatmap::pheatmap(t(pheatmap_df),
#                          annotation_row = sig_res,
#                          cutree_cols = 3,scale = 'column',
#                          cluster_rows = F,
#                          # cluster_rows = FALSE,
#                          annotation_col = metadata[,c('Group','Gender','Region','Smoke')],
#                          show_rownames = FALSE,
#                          show_colnames = FALSE)
# 
# heatmap_clusters <- cutree(p1$tree_col, k = 3)  # 把列聚成7类
# 
# metadata$Cluster <- as.factor(heatmap_clusters[rownames(metadata)])  # 匹配样本顺序
# 
# # 再画一个更新后的热图
# p2 <- pheatmap::pheatmap(t(pheatmap_df),
#                          annotation_row = sig_res,
#                          # scale = 'row',
#                          cutree_cols = 2,
#                          cluster_rows = FALSE,
#                          annotation_col = metadata[,c('Cluster','Group','Gender','Region','Smoke')],
#                          show_rownames = FALSE,
#                          show_colnames = FALSE,
#                          filename = 'gene_heatmap.pdf',
#                          width = 10,height = 12)
