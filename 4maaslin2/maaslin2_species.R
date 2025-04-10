rm(list = ls())
setwd('/Users/yangwanyi/Downloads/ENIGMA_STUDY_DATA/4maaslin2/')

# 加载MaAsLin2
library(Maaslin2)
library(dplyr)


# 输入数据
microbiome_data <- read.csv("../3diversity/new_species_meta_df.csv", row.names = 1)
microbiome_data$Group <- ifelse(microbiome_data$Group != "Control","CD","Control")
# keep <- apply(microbiome_data[,10:ncol(microbiome_data)], 2, mean) > 1E-4 & 
#   apply(microbiome_data[,10:ncol(microbiome_data)] > 0, 2, sum) / nrow(microbiome_data[,10:ncol(microbiome_data)]) > 0.1 # mean abundance > 1E-4、Prevalence > 0.1
# microbiome_data <- microbiome_data[, keep]
microbiome_data <- microbiome_data[,10:ncol(microbiome_data)]
dim(microbiome_data)
saveRDS(microbiome_data,file = 'microbiome_data.rds')

metadata <- read.csv('../3diversity/new_metadata.csv',row.names = 1)

# 运行MaAsLin2
fit_data <- Maaslin2(
  microbiome_data, metadata, 
  output = 'output_Group_species',
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
  microbiome_data, metadata1, 
  output = 'output_Region_species_AUS_unfilter',
  min_abundance = 0,       # 过滤低丰度物种（>1%平均丰度）
  min_prevalence = 0,       # 过滤低存在率物种（>10%样本存在）
  normalization = "CLR",      # 使用CLR处理组成型数据
  transform = "NONE",         # CLR后无需再次转换
  analysis_method = "LM",     # 保持线性模型（CLR后数据更接近正态）
  fixed_effects = c('Region'),
  reference = c("Region,other"),
  correction = "BH",          # 默认FDR校正
  standardize = FALSE, 
  plot_heatmap = FALSE,
  plot_scatter = FALSE,
  cores=4)

metadata2 <- metadata
metadata2$Region <- ifelse(metadata$Region == 'HK','HK','other')

fit_data <- Maaslin2(
  microbiome_data, metadata2, 
  output = 'output_Region_species_HK_unfilter',
  min_abundance = 0,       # 过滤低丰度物种（>1%平均丰度）
  min_prevalence = 0,       # 过滤低存在率物种（>10%样本存在）
  normalization = "CLR",      # 使用CLR处理组成型数据
  transform = "NONE",         # CLR后无需再次转换
  analysis_method = "LM",     # 保持线性模型（CLR后数据更接近正态）
  fixed_effects = c('Region'),
  reference = c("Region,other"),
  correction = "BH",          # 默认FDR校正
  standardize = FALSE, 
  plot_heatmap = FALSE,
  plot_scatter = FALSE,
  cores=4)


metadata3 <- metadata
metadata3$Region <- ifelse(metadata$Region == 'MC','MC','other')

fit_data <- Maaslin2(
  microbiome_data, metadata3, 
  output = 'output_Region_species_MC_unfilter',
  min_abundance = 0,       # 过滤低丰度物种（>1%平均丰度）
  min_prevalence = 0,       # 过滤低存在率物种（>10%样本存在）
  normalization = "CLR",      # 使用CLR处理组成型数据
  transform = "NONE",         # CLR后无需再次转换
  analysis_method = "LM",     # 保持线性模型（CLR后数据更接近正态）
  fixed_effects = c('Region'),
  reference = c("Region,other"),
  correction = "BH",          # 默认FDR校正
  standardize = FALSE, 
  plot_heatmap = FALSE,
  plot_scatter = FALSE,
  cores=4)


metadata1 <- metadata
metadata1$Region <- ifelse(metadata$Region == 'AUS','AUS','other')

fit_data <- Maaslin2(
  microbiome_data, metadata1, 
  output = 'output_Region_species_AUS_filter',
  min_abundance = 0.01,       # 过滤低丰度物种（>1%平均丰度）
  min_prevalence = 0.1,       # 过滤低存在率物种（>10%样本存在）
  normalization = "CLR",      # 使用CLR处理组成型数据
  transform = "NONE",         # CLR后无需再次转换
  analysis_method = "LM",     # 保持线性模型（CLR后数据更接近正态）
  fixed_effects = c('Region'),
  reference = c("Region,other"),
  correction = "BH",          # 默认FDR校正
  standardize = FALSE, 
  plot_heatmap = FALSE,
  plot_scatter = FALSE,
  cores=4)

metadata2 <- metadata
metadata2$Region <- ifelse(metadata$Region == 'HK','HK','other')

fit_data <- Maaslin2(
  microbiome_data, metadata2, 
  output = 'output_Region_species_HK_filter',
  min_abundance = 0.01,       # 过滤低丰度物种（>1%平均丰度）
  min_prevalence = 0.1,       # 过滤低存在率物种（>10%样本存在）
  normalization = "CLR",      # 使用CLR处理组成型数据
  transform = "NONE",         # CLR后无需再次转换
  analysis_method = "LM",     # 保持线性模型（CLR后数据更接近正态）
  fixed_effects = c('Region'),
  reference = c("Region,other"),
  correction = "BH",          # 默认FDR校正
  standardize = FALSE, 
  plot_heatmap = FALSE,
  plot_scatter = FALSE,
  cores=4)


metadata3 <- metadata
metadata3$Region <- ifelse(metadata$Region == 'MC','MC','other')

fit_data <- Maaslin2(
  microbiome_data, metadata3, 
  output = 'output_Region_species_MC_filter',
  min_abundance = 0.01,       # 过滤低丰度物种（>1%平均丰度）
  min_prevalence = 0.1,       # 过滤低存在率物种（>10%样本存在）
  normalization = "CLR",      # 使用CLR处理组成型数据
  transform = "NONE",         # CLR后无需再次转换
  analysis_method = "LM",     # 保持线性模型（CLR后数据更接近正态）
  fixed_effects = c('Region'),
  reference = c("Region,other"),
  correction = "BH",          # 默认FDR校正
  standardize = FALSE, 
  plot_heatmap = FALSE,
  plot_scatter = FALSE,
  cores=4)




# sig_res <- read.csv('output_Group/significant_results.tsv',sep = '\t',header = T)
# sig_res$CD_status <- ifelse(sig_res$coef>0,'CD enriched','CD depleted')
# rownames(sig_res) <- sig_res$feature
# sig_res <- sig_res[,'CD_status',drop=FALSE]
# sig_res <- sig_res %>%
#   arrange(CD_status)
# 
# microbiome_data_df <- read.csv('output_Group/features/filtered_data_norm_transformed.tsv',sep = '\t',header = T,row.names = 1)
# pheatmap_df <- microbiome_data_df[,rownames(sig_res)]
# 
# p1 <- pheatmap::pheatmap(t(pheatmap_df),
#                    annotation_row = sig_res,
#                    cutree_cols = 7,cluster_rows = F,
#                    # cluster_rows = FALSE,
#                    annotation_col = metadata[,c('Group','Gender','Region','Smoke')],
#                    show_rownames = FALSE,
#                    show_colnames = FALSE)
# 
# heatmap_clusters <- cutree(p1$tree_col, k = 7)  # 把列聚成7类
# 
# metadata$Cluster <- as.factor(heatmap_clusters[rownames(metadata)])  # 匹配样本顺序
# 
# # 再画一个更新后的热图
# p2 <- pheatmap::pheatmap(t(pheatmap_df),
#                          annotation_row = sig_res,
#                          # scale = 'row',
#                          cutree_cols = 7,
#                          cluster_rows = FALSE,
#                          annotation_col = metadata[,c('Cluster','Group','Gender','Region','Smoke')],
#                          show_rownames = FALSE,
#                          show_colnames = FALSE,
#                          filename = 'species_heatmap.pdf',
#                          width = 10,height = 12)
