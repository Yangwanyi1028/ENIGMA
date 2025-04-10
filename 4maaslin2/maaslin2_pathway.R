rm(list = ls())
setwd('/Users/yangwanyi/Downloads/ENIGMA_STUDY_DATA/4maaslin2/')

# 加载MaAsLin2
library(Maaslin2)
library(dplyr)


# 输入数据
path_data <- read.csv("../source_data/humann_pathabundance_relab_AUS_HK_KM_479_unstratified.csv", row.names = 1)
path_data <- t(path_data) 
keep <- apply(path_data, 2, mean) > 1E-4 &
  apply(path_data > 0, 2, sum) / nrow(path_data) > 0.1 # mean abundance > 1E-4、Prevalence > 0.1
path_data <- path_data[, keep]

save(path_data,file = 'path_data.rds')
metadata <- read.csv('../3diversity/new_metadata.csv',row.names = 1)
metadata$Group <- ifelse(metadata$Group != 'Control','CD','Control')

# 替换列名中的非法字符（如冒号、空格等）为下划线
# 更新数据框的列名
new_col_names <- paste("Pathway", 1:ncol(path_data), sep="")
# 设置新的行名
original_col_names <- colnames(path_data)

# 创建映射表
mapping_table <- data.frame(
  OriginalRowName = original_col_names,
  NewRowName = new_col_names,
  stringsAsFactors = FALSE  # 避免将字符列转换为因子
)
write.csv(mapping_table,'mapping_table_pathway_name.csv')
colnames(path_data) <- new_col_names


fit_data <- Maaslin2(
  path_data, metadata, 
  output = 'output_Group_pathway_filter',
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
  path_data, metadata1, 
  output = 'output_Region_pathway_AUS',
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
  path_data, metadata2, 
  output = 'output_Region_pathway_HK',
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
  path_data, metadata3, 
  output = 'output_Region_pathway_MC',
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

# sig_res <- read.csv('output_Group_pathway_filter/significant_results.tsv',sep = '\t',header = T)
# dim(sig_res)
# sig_res <- read.csv('output_Group_pathway_unfilter/significant_results.tsv',sep = '\t',header = T)
# sig_res$CD_status <- ifelse(sig_res$coef>0,'CD enriched','CD depleted')
# rownames(sig_res) <- sig_res$feature
# sig_res <- sig_res[,'CD_status',drop=FALSE]
# sig_res <- sig_res %>%
#   arrange(CD_status)

# microbiome_data_df <- read.csv('output_Group_pathway/features/filtered_data_norm_transformed.tsv',sep = '\t',header = T,row.names = 1)
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
# heatmap_clusters <- cutree(p1$tree_col, k = 7)  # 把列聚成7类
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
#                          filename = 'pathway_heatmap.pdf',
#                          width = 10,height = 12)
