rm(list = ls())
setwd('/Users/yangwanyi/Downloads/ENIGMA_STUDY_DATA/4maaslin2/')

# 加载MaAsLin2
library(Maaslin2)
library(dplyr)
library(pheatmap)
library(stringr)


# 输入metadata
metadata <- read.csv('../3diversity/new_metadata.csv',row.names = 1)
metadata$Group <- ifelse(metadata$Group != 'Control','CD','Control')
metadata$study_id <- rownames(metadata)
metadta_with_diseas_status <- readxl::read_excel('../source_data/ENIGMA_metadata.xlsx')
length(intersect(rownames(metadata),metadta_with_diseas_status$study_id))
metadta_with_diseas_status <- metadta_with_diseas_status[,c('study_id','cd_act')]
metadata <- merge(metadata,metadta_with_diseas_status,by='study_id')
dim(metadata)

metadata <- rename(metadata,'Disease_activity'='cd_act')
# 把元数据中"disease_activity"列的空值填为"Unknown"
# metadata$Disease_activity <- ifelse(is.na(metadata$Disease_activity),'Unknown',metadata$Disease_activity)


# 数据预处理函数
process_disease_activity <- function(df) {
  df %>%
    mutate(
      # 提取study_id第一个数字并转换为患者类型
      patient_type = ifelse(str_extract(study_id, "\\d") == "1", 
                            "CD", "Control"),
      
      # 动态填充策略
      Disease_activity = case_when(
        # 健康样本处理
        patient_type == "Control" & is.na(Disease_activity) ~ "Control",
        
        # CD患者处理
        patient_type == "CD" & is.na(Disease_activity) ~ "Unknown",
        
        # 非NA值直接继承
        TRUE ~ Disease_activity
      )
    ) %>%
    select(-patient_type)  # 移除临时列
}

metadata <- process_disease_activity(metadata)
rownames(metadata) <- metadata$study_id
metadata <- metadata[,-1]

# 读取元数据
colnames(metadata)
metadata <- metadata[,c('Group','Region','Disease_activity')]

# 输入pathway数据
path_data <- read.csv("../source_data/new_humann_pathabundance_relab_AUS_HK_KM_479_unstratified.csv", row.names = 1)
# path_data <- t(path_data)
# keep <- apply(path_data, 2, mean) > 1E-5 &
#   apply(path_data > 0, 2, sum) / nrow(path_data) > 0.1 # mean abundance > 1E-4、Prevalence > 0.1
# path_data <- path_data[, keep]
# path_data <- t(path_data)
path_data <- log(path_data+1)
range(path_data)
# 输入物种数据
microbiome_data <- read.csv("../3diversity/new_species_meta_df.csv", row.names = 1)
microbiome_data$Group <- ifelse(microbiome_data$Group != "Control","CD","Control")
keep <- apply(microbiome_data[,10:ncol(microbiome_data)], 2, mean) > 1E-4 & 
  apply(microbiome_data[,10:ncol(microbiome_data)] > 0, 2, sum) / nrow(microbiome_data[,10:ncol(microbiome_data)]) > 0.1 # mean abundance > 1E-4、Prevalence > 0.1
microbiome_data <- microbiome_data[, keep]
# metadata <- microbiome_data[,1:9]
microbiome_data <- microbiome_data[,10:ncol(microbiome_data)]
microbiome_data <- t(microbiome_data)
microbiome_data <- log(microbiome_data+1)
range(path_data)
# 输入基因数据
gene_data <- read.csv("../source_data/new_ARG68_gene_abundance_20250224.csv", row.names = 1)

gene_name <- gene_data[,1,drop=FALSE]
gene_data <- gene_data[,-1]
gene_data <- log(log(gene_data+1)+1)
range(gene_data)
# 把三个特征数据合并到一起
combined_data <- rbind(gene_data, microbiome_data, path_data)

# 定义一个Z-score标准化的函数
zscore <- function(x) {
  if (sd(x, na.rm = TRUE) == 0) {  # 如果标准差为0，说明这个列的值全是一样的
    return(rep(0, length(x)))  # 返回全是0的列
  } else {
    return((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))  # 否则，标准化为Z-score
  }
}

# 对所有数据做Z-score标准化
combined_data <- apply(combined_data, 2, zscore)  # 按列应用Z-score函数
combined_data <- t(scale(t(combined_data)))  # 对行做进一步标准化

# 限制Z-score的值在[-3, 3]范围内，太大的值就设为3，太小的设为-3
zlim = c(-3, 3)
combined_data <- pmin(pmax(combined_data, zlim[1]), zlim[2])
range(combined_data)  # 检查数据范围

# 读取行注释信息
# row_annotation = read.csv('Combined_Annotations_with_CD_Status.csv',row.names = 1)
path_maaslin_sig_res <- read.csv('output_Group_pathway_filter/significant_results.tsv',
                                 sep = '\t',header = T)
mapping_pathway_name <- read.csv('mapping_table_pathway_name.csv')
path_maaslin_sig_res_new <- merge(path_maaslin_sig_res,mapping_pathway_name,by.x='feature',by.y="NewRowName")
path_maaslin_sig_res_new$feature_type <- "Pathway"
path_maaslin_sig_res_new$feature <- path_maaslin_sig_res_new$OriginalRowName
path_maaslin_sig_res_new <- path_maaslin_sig_res_new[,-c(10,11)]

gene_maaslin_sig_res <- read.csv('output_Group_gene/significant_results.tsv',
                                 sep = '\t',header = T)
gene_maaslin_sig_res$feature_type <- "Gene"
sp_maaslin_sig_res <- read.csv('output_Group_species/significant_results.tsv',
                                 sep = '\t',header = T)
sp_maaslin_sig_res$feature_type <- "Species"
sig_res <- rbind(path_maaslin_sig_res_new,gene_maaslin_sig_res,sp_maaslin_sig_res)
sig_res$CD_status <- ifelse(sig_res$coef>0,"CD-enriched","CD-depleted")

rownames(sig_res) <- sig_res$feature
row_annotation <- sig_res[,c('feature_type','CD_status')]

combined_data <- combined_data[sig_res$feature,]

# 定义列注释颜色
# 定义地区的颜色
region_colors <- c("#85A894", "#3E6B7E", "#FCE4A8")
names(region_colors) <- c("Australia", "Hong Kong", "China")

# 定义疾病活动状态的颜色
disease_colors <- c("#Ecc03f", "#FC8002", "#B9181A", "grey","#ADDB88")
names(disease_colors) <- c("Remission", "Mild Act", "Moderate", "Unknown", "Control")

group_colors <- c('#AB6355',"#90a5a6")
names(group_colors) <- c("CD", "Control")

# 定义行注释颜色
Feature_Type_colors <- c("#F89FA8","#66BC98","#E26844")
names(Feature_Type_colors) <- c("Gene", "Species", "Pathway")
CD_Status_colors = c('#AB6355',"#90a5a6")
names(CD_Status_colors) <- c("CD-enriched", "CD-depleted")


# 将所有注释颜色整合到一个列表中
annotation_colors <- list(
  Region = region_colors,
  Disease_activity = disease_colors,
  Group = group_colors,
  feature_type = Feature_Type_colors,
  CD_status = CD_Status_colors
)



# 画一个热图
p1 = pheatmap(
  combined_data,
  cutree_cols = 4,  # 把列分成7个类
  clustering_distance_rows = "correlation",  # 行的聚类用相关性
  scale = "none",  # 不再额外标准化
  clustering_distance_cols = "correlation",  # 列的聚类用相关性
  clustering_method = "ward.D2",  # 聚类算法用ward.D2
  annotation_col = metadata,  # 列注释来源于元数据
  annotation_colors = annotation_colors,  # 列注释的颜色  
  annotation_row = row_annotation,  # 行注释来源于特征类型
  show_rownames = FALSE,  # 不显示行名字
  show_colnames = FALSE,  # 显示列名字
  # main = "Heatmap of VFG, Species, Path Features with Sample Clustering",  # 热图标题
  filename = 'sample_features.png',
  width = 10,height = 8
)

# 获取热图中列的聚类结果
heatmap_clusters <- cutree(p1$tree_col, k = 4)  # 把列聚成7类

# 把聚类信息加入元数据
metadata$Cluster <- as.factor(heatmap_clusters[rownames(metadata)])  # 匹配样本顺序


# region_colors <- c("#85A894", "#3E6B7E", "#FCE4A8")
# names(region_colors) <- c("AUS", "HK", "KM")
region_colors <- c("#85A894", "#3E6B7E", "#FCE4A8")
names(region_colors) <- c("Australia", "Hong Kong", "China")

# 定义疾病活动状态的颜色
disease_colors <- c("#Ecc03f", "#FC8002", "#B9181A", "grey","#ADDB88")
names(disease_colors) <- c("Remission", "Mild Act", "Moderate", "Unknown", "Control")


group_colors <- c('#AB6355',"#90a5a6")
names(group_colors) <- c("CD", "Control")


cluster_colors = c("#A1A9D0","#F0988C","#B883D4","#CFEAF1")# ,"#C4A5DE","#F6CA5E","#9E99EE"
names(cluster_colors) <- levels(metadata$Cluster)  # 给颜色命名

annotation_colors <- list(
  Region = region_colors,
  Disease_activity = disease_colors,
  Group = group_colors,
  feature_type = Feature_Type_colors,
  CD_status = CD_Status_colors,
  Cluster = cluster_colors
)

# 再画一个更新后的热图
p2 <- pheatmap(
  combined_data,
  cutree_cols = 4,
  clustering_distance_rows = "correlation",
  scale = "none",
  clustering_distance_cols = "correlation",
  clustering_method = "ward.D2",
  annotation_col = metadata,  # 这次包含了Cluster信息
  annotation_colors = annotation_colors,
  annotation_row = row_annotation,
  show_rownames = FALSE,
  show_colnames = FALSE,  # 这次不显示列名字
  # main = "Heatmap of VFG, Species, Path Features with Sample Clustering"  # 热图标题
  filename = 'sample_features_new.png',
  width = 10,height = 8
)

# 打印每个Cluster对应的样本
for (cluster in levels(metadata$Cluster)) {
  cat(sprintf("Cluster %s samples:\n", cluster))  # 打印Cluster编号
  print(rownames(metadata[metadata$Cluster == cluster, ]))  # 打印样本名字
}

# 把结果保存到CSV文件
write.csv(metadata, file = "clustered_samples.csv", row.names = TRUE)  # 保存为CSV

# 打印metadata看一看
print(metadata)


# 计算每个Cluster的中心（平均值）
cluster_centers <- do.call(cbind, lapply(levels(metadata$Cluster), function(cluster) {
  colMeans(combined_data[, metadata$Cluster == cluster], na.rm = TRUE)  # 计算每个Cluster的均值
}))
colnames(cluster_centers) <- levels(metadata$Cluster)  # 给列名加上Cluster名字

# 计算每个样本到各Cluster中心的欧几里得距离
distance_matrix <- sapply(1:ncol(cluster_centers), function(cluster_idx) {
  apply(combined_data, 2, function(sample) {
    sqrt(sum((sample - cluster_centers[, cluster_idx])^2))  # 欧几里得距离
  })
})
colnames(distance_matrix) <- levels(metadata$Cluster)  # 给距离矩阵列名加上Cluster名字

# 针对两两 Cluster 计算 Wilcoxon 检验
clusters <- levels(metadata$Cluster)  # 获取所有 Cluster 的名称
pairwise_wilcoxon_results <- list()  # 用于保存所有两两比较结果

for (i in 1:(length(clusters) - 1)) {
  for (j in (i + 1):length(clusters)) {
    cluster1 <- clusters[i]
    cluster2 <- clusters[j]
    
    # 获取 cluster1 和 cluster2 中样本的距离
    samples_in_cluster1 <- rownames(metadata[metadata$Cluster == cluster1, ])
    samples_in_cluster2 <- rownames(metadata[metadata$Cluster == cluster2, ])
    
    # 提取对应的样本距离
    distances_cluster1 <- distance_matrix[samples_in_cluster1, cluster2]
    distances_cluster2 <- distance_matrix[samples_in_cluster2, cluster1]
    
    # 进行 Wilcoxon 检验
    wilcox_test <- wilcox.test(distances_cluster1, distances_cluster2, paired = FALSE)
    
    # 保存检验结果
    pairwise_wilcoxon_results[[paste(cluster1, "vs", cluster2)]] <- list(
      p.value = wilcox_test$p.value,
      statistic = wilcox_test$statistic
    )
  }
}

# 打印两两比较的结果
cat("Pairwise Wilcoxon Test Results:\n")
for (comparison in names(pairwise_wilcoxon_results)) {
  result <- pairwise_wilcoxon_results[[comparison]]
  cat(sprintf(
    "%s: p-value = %.5f, statistic = %.2f\n",
    comparison, result$p.value, result$statistic
  ))
}

# 如果需要，将结果保存为CSV文件
pairwise_results_df <- do.call(rbind, lapply(names(pairwise_wilcoxon_results), function(comparison) {
  result <- pairwise_wilcoxon_results[[comparison]]
  data.frame(
    Comparison = comparison,
    P_Value = result$p.value,
    Statistic = result$statistic
  )
}))
write.csv(pairwise_results_df, file = "pairwise_wilcoxon_results.csv", row.names = FALSE)




# 加载绘图库
library(ggplot2)

# 定义一个函数来统计每个 Cluster 的分布情况并绘制饼图
plot_cluster_pie <- function(metadata, cluster_col = "Cluster", group_col, plot_title_prefix, color_palette) {
  clusters <- levels(metadata[[cluster_col]])  # 获取所有 Cluster
  
  for (cluster in clusters) {
    # 筛选当前 cluster 的样本
    samples_in_cluster <- metadata[metadata[[cluster_col]] == cluster, ]
    # 统计 group_col 的分布
    group_counts <- table(samples_in_cluster[[group_col]])
    group_counts_df <- as.data.frame(group_counts)
    colnames(group_counts_df) <- c(group_col, "Count")
    
    # 计算百分比
    group_counts_df$Percentage <- group_counts_df$Count / sum(group_counts_df$Count) * 100
    
    # 绘制饼图
    p <- ggplot(group_counts_df, aes(x = "", y = Percentage, fill = !!sym(group_col))) +
      geom_bar(stat = "identity", width = 1, color = "white") +
      coord_polar("y", start = 0) +
      scale_fill_manual(values = color_palette) +  # 使用指定的颜色
      # labs(
      #   title = paste(plot_title_prefix, " - Cluster", cluster),
      #   x = NULL,
      #   y = NULL,
      #   fill = group_col
      # ) +
      theme_minimal() +
      theme(
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14)
      )
    
    # 显示饼图
    print(p)
    
    # 保存饼图为 PDF 或 PNG 文件
    ggsave(
      filename = paste0("Pie_", plot_title_prefix, "_Cluster_", cluster, ".png"),
      plot = p,
      width = 6,
      height = 6
    )
  }
}

# 绘制每个 Cluster 的 `disease_activity` 饼图，使用指定颜色
plot_cluster_pie(
  metadata, 
  group_col = "Disease_activity", 
  plot_title_prefix = "Disease_activity", 
  color_palette = disease_colors
)

# 绘制每个 Cluster 的 `region` 饼图，使用指定颜色
plot_cluster_pie(
  metadata, 
  group_col = "Region", 
  plot_title_prefix = "region", 
  color_palette = region_colors
)
