rm(list = ls())
# 加载我们需要用到的一些R库
library(made4)  # 用于数据分析
library(dynamicTreeCut)  # 用于动态剪切树
library(pheatmap)  # 用于画热图

# 设置工作目录，这里是你放数据的文件夹路径
setwd('/Users/wanyiddl/Downloads/ENIGMA_STUDY_DATA/5pheatmap_CD/')

# source data file path
file_path <- '/Users/wanyiddl/Downloads/ENIGMA_STUDY_DATA/newsource_data20250319/'
# 定义我们将用到的几个数据文件的路径
metadata_path <- paste0(file_path,'副本metadata_match20250108_towanyi.xlsx')  # 元数据文件路径
vfg_path <- paste0(file_path,"ARG_68_479samples.csv")  # VFG特征数据
species_path <- paste0(file_path,"species_101_479samples.csv")  # 物种特征数据
path_path <- paste0(file_path,"path_68_479samples.csv")  # 路径特征数据

# 读取数据文件，按行名读取
vfg_data <- read.csv(vfg_path, row.names = 1)  # VFG特征数据
species_data <- read.csv(species_path, row.names = 1)  # 物种特征数据
path_data <- read.csv(path_path, row.names = 1)  # 路径特征数据

# 把三个特征数据合并到一起
combined_data <- rbind(vfg_data, species_data, path_data)
saveRDS(combined_data,file = 'raw_combined_data.rds')
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

# 读取元数据
library(readxl)  # 用来读Excel文件
metadata <- read_excel(metadata_path)  # 读取元数据

# 把元数据中"disease_activity"列的空值填为"Unknown"
metadata = as.data.frame(metadata)  # 转为数据框
rownames(metadata) = metadata$sampleID  # 设置行名为sampleID
metadata = metadata[,-1]  # 去掉第一列

# 定义列注释颜色
# 定义地区的颜色
region_colors <- c("#85A894", "#3E6B7E", "#FCE4A8")
names(region_colors) <- c("AUS", "HK", "MC")

# 定义疾病活动状态的颜色
disease_colors <- c("#7e9bb7", "#FC8002", "#B9181A", "grey","#ADDB88")
names(disease_colors) <- c("Remission", "Mild Act", "Moderate", "CD_Unknown", "Control")

group_colors <- c('#AB6355',"#FC9148")
names(group_colors) <- c("CD", "Controls")

# 读取行注释信息
row_annotation_file <- paste0(file_path,'Combined_Annotations_with_CD_Status.csv')
row_annotation = read.csv(row_annotation_file,row.names = 1)
# 定义行注释颜色
Feature_Type_colors <- c("#F89FA8","#66BC98","#E26844","#5Fc9E2")
names(Feature_Type_colors) <- c("VFG", "ARG","Species", "Pathway")
CD_Status_colors = c("#ab6355", "#90a5a7")
names(CD_Status_colors) <- c("CD_enriched", "CD_depleted")


# 将所有注释颜色整合到一个列表中
annotation_colors <- list(
  region = region_colors,
  disease_activity = disease_colors,
  group = group_colors,
  Feature_Type = Feature_Type_colors,
  CD_Status = CD_Status_colors
)



# 画一个热图
p1 = pheatmap(
  combined_data,
  cutree_cols = 6,  # 把列分成7个类
  clustering_distance_rows = "correlation",  # 行的聚类用相关性
  scale = "none",  # 不再额外标准化
  clustering_distance_cols = "correlation",  # 列的聚类用相关性
  clustering_method = "ward.D2",  # 聚类算法用ward.D2
  annotation_col = metadata,  # 列注释来源于元数据
  annotation_colors = annotation_colors,  # 列注释的颜色  
  annotation_row = row_annotation,  # 行注释来源于特征类型
  show_rownames = FALSE,  # 不显示行名字
  show_colnames = FALSE,  # 显示列名字
  main = "Heatmap of VFG, Species, Path Features with Sample Clustering"  # 热图标题
)

# 获取热图中列的聚类结果
heatmap_clusters <- cutree(p1$tree_col, k = 6)  # 把列聚成7类

# 把聚类信息加入元数据
metadata$Cluster <- as.factor(heatmap_clusters[rownames(metadata)])  # 匹配样本顺序


cluster_colors = c("#B1A1A0","#F0988C","#B983D4","#CFEAF1","#A4A5DE","#F6CA5E") #"#9E99EE"
names(cluster_colors) <- levels(metadata$Cluster)  # 给颜色命名

annotation_colors <- list(
  region = region_colors,
  disease_activity = disease_colors,
  group = group_colors,
  Feature_Type = Feature_Type_colors,
  CD_Status = CD_Status_colors,
  Cluster = cluster_colors
)

# 再画一个更新后的热图
p2 <- pheatmap(
  combined_data,
  cutree_cols = 6,
  clustering_distance_rows = "correlation",
  scale = "none",
  clustering_distance_cols = "correlation",
  clustering_method = "ward.D2",
  annotation_col = metadata,  # 这次包含了Cluster信息
  annotation_colors = annotation_colors,
  annotation_row = row_annotation,
  show_rownames = FALSE,
  show_colnames = FALSE,  # 这次不显示列名字
  filename = 'heatmap.png',width = 10,height = 10
  # main = "Heatmap of VFG, Species, Path Features with Sample Clustering"  # 热图标题
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
saveRDS(metadata,file='metadata.rds')


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


plot_cluster_pie <- function(metadata, cluster_col = "Cluster", group_col, plot_title_prefix, color_palette) {
  clusters <- levels(metadata[[cluster_col]])
  plot_list <- list()  # 用于存储所有饼图对象
  
  for (i in seq_along(clusters)) {
    cluster <- clusters[i]
    samples_in_cluster <- metadata[metadata[[cluster_col]] == cluster, ]
    group_counts <- table(samples_in_cluster[[group_col]])
    group_counts_df <- as.data.frame(group_counts)
    colnames(group_counts_df) <- c(group_col, "Count")
    
    group_counts_df$Percentage <- group_counts_df$Count / sum(group_counts_df$Count) * 100
    
    p <- ggplot(group_counts_df, aes(x = "", y = Percentage, fill = !!sym(group_col))) +
      geom_bar(stat = "identity", width = 1, color = "white") +
      coord_polar("y", start = 0) +
      scale_fill_manual(values = color_palette) +
      labs(
        title = paste("Cluster", cluster),
        x = NULL,
        y = NULL,
        fill = group_col
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = "none"  # 隐藏单个图例
      )
    
    # 单独保存每个饼图
    ggsave(
      filename = paste0("Pie_", plot_title_prefix, "_Cluster_", cluster, ".png"),
      plot = p,
      width = 4,
      height = 4
    )
    
    plot_list[[i]] <- p
  }
  
  # 组合所有饼图并添加共享图例
  combined_plot <- cowplot::plot_grid(
    plotlist = plot_list,
    nrow = 1,
    align = "h"
  ) + theme(plot.margin = margin(20, 20, 20, 20))
  
  # 添加共享图例
  legend <- cowplot::get_legend(plot_list[[1]] + theme(legend.position = "bottom"))
  final_plot <- cowplot::plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(1, 0.15))
  
  # 保存组合图
  ggsave(
    filename = paste0("Combined_Pie_", plot_title_prefix, ".png"),
    plot = final_plot,
    width = length(clusters) * 4,
    height = 5
  )
  
  print(final_plot)
}

# 绘制每个 Cluster 的 `disease_activity` 饼图，使用指定颜色
plot_cluster_pie(
  metadata, 
  group_col = "disease_activity", 
  plot_title_prefix = "Disease Activity Distribution", 
  color_palette = disease_colors
)



# 绘制每个 Cluster 的 `region` 饼图，使用指定颜色
plot_cluster_pie(
  metadata, 
  group_col = "region", 
  plot_title_prefix = "Region Distribution", 
  color_palette = region_colors
)
