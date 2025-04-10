rm(list = ls())
setwd('/Users/wanyiddl/Downloads/ENIGMA_STUDY_DATA/6network_cls6/')
library(compositions)
library(statsExpressions)
library(caret)
library(Hmisc)
library(dplyr)

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
saveRDS(combined_data, file = 'raw_combined_data.rds')

# 读取节点类型数据
row_annotation_file <- paste0(file_path,'Combined_Annotations_with_CD_Status.csv')
row_annotation = read.csv(row_annotation_file,row.names = 1)

head(row_annotation)

row_annotation$id <- rownames(row_annotation)

# 加载数据
microbiome_data <- readRDS('raw_combined_data.rds')
microbiome_data <- as.data.frame(t(microbiome_data))

range(microbiome_data)
dim(microbiome_data)
class(microbiome_data)

microbiome_data$study_id <- rownames(microbiome_data)
phenotype_data <- readRDS("/Users/wanyiddl/Downloads/ENIGMA_STUDY_DATA/5pheatmap_CD/cls6/metadata.rds")
phenotype_data$study_id <- rownames(phenotype_data)
phenotype_data <- phenotype_data[, c('study_id', 'Cluster')]
merge_data <- merge(microbiome_data, phenotype_data, by = 'study_id')

# 定义节点类型
node_types <- c(
  rep("Gene", nrow(vfg_data)),      # VFG 特征标记为 Gene
  rep("Species", nrow(species_data)), # 物种特征标记为 Species
  rep("Pathway", nrow(path_data))    # 路径特征标记为 Pathway
)


library(igraph)

filter_low_degree <- function(nodes, edges, min_degree=5) {
  # 创建图对象
  g <- graph_from_data_frame(d=edges, vertices=nodes, directed=FALSE)
  
  # 迭代过滤直到所有节点满足度要求
  repeat {
    # 计算当前度
    current_degree <- degree(g)
    
    # 标记需要保留的节点（度 >= min_degree）
    keep_nodes <- names(current_degree[current_degree >= min_degree])
    
    # 如果没有需要过滤的节点则退出
    if (length(keep_nodes) == vcount(g)) break
    
    # 重新创建子图
    g <- induced_subgraph(g, keep_nodes)
  }
  
  # 返回过滤后的边和节点
  list(
    nodes = data.frame(id = V(g)$name, type = V(g)$type),
    edges = get.data.frame(g)
  )
}

# 定义检测函数
check_bidirectional_edges <- function(edge_df) {
  # 创建排序后的节点对标识列
  edge_df$sorted_pair <- apply(edge_df[, c("source", "target")], 1, function(x) {
    paste(sort(x), collapse = "--")
  })
  
  # 查找重复的排序对
  duplicates <- edge_df[duplicated(edge_df$sorted_pair) | 
                          duplicated(edge_df$sorted_pair, fromLast = TRUE), ]
  
  if (nrow(duplicates) > 0) {
    cat("发现双向边：\n")
    print(duplicates[order(duplicates$sorted_pair), ])
    cat("\n共计", nrow(duplicates)/2, "对双向边\n")
  } else {
    cat("未发现双向边\n")
  }
}



for (i in unique(merge_data$Cluster)) {
  print(paste('Cluster:', i))
  
  # 提取当前 Cluster 的数据
  corr_df <- merge_data[merge_data$Cluster == i, ]
  
  # 移除不需要的列（Cluster 和 study_id）
  corr_df$Cluster <- NULL
  corr_df$study_id <- NULL
  # 转换为相对丰度（行和为1）
  corr_df_rel <- corr_df / rowSums(corr_df)
  # CLR 转换
  clr_corr_df <- clr(corr_df + 1e-6)
  
  # 计算相关系数矩阵
  cor_result <- rcorr(as.matrix(clr_corr_df), type = "pearson")# spearman
  
  # 提取微生物特征之间的相关性
  r_matrix <- cor_result$r  # 全部微生物特征的相关性矩阵
  p_matrix <- cor_result$P  # 全部微生物特征的 p 值矩阵
  print(range(r_matrix))
  # 提取下三角（不含对角线）的p值
  lower_tri <- lower.tri(p_matrix, diag = FALSE)
  p_flat <- p_matrix[lower_tri]
  
  # FDR校正
  fdr_p_flat <- p.adjust(p_flat, method = "BH")
  # 重建校正后的p值矩阵
  fdr_p <- matrix(NA, nrow = nrow(p_matrix), ncol = ncol(p_matrix))
  fdr_p[lower_tri] <- fdr_p_flat
  fdr_p[upper.tri(fdr_p)] <- t(fdr_p)[upper.tri(fdr_p)]
  # FDR 校正
  # fdr_p <- matrix(p.adjust(p_matrix, method = "BH"), nrow = nrow(p_matrix))
  rownames(fdr_p) <- rownames(p_matrix)
  colnames(fdr_p) <- colnames(p_matrix)
  
  # 筛选显著相关性（|r| > 0.1 且 FDR < 0.05）
  # significant <- which(r_matrix > 0.7 | r_matrix < -0.7 & fdr_p < 0.01, arr.ind = TRUE)
  significant <- which(abs(r_matrix) > 0.5 & fdr_p < 0.01, arr.ind = TRUE)
  edges <- data.frame(
    source = rownames(r_matrix)[significant[, 1]],
    target = colnames(r_matrix)[significant[, 2]],
    weight = r_matrix[significant],
    p = fdr_p[significant]
  )
  print(range(edges$weight))
  print(dim(edges[edges$weight<0,]))
  # 节点列表（微生物特征）
  nodes <- data.frame(
    id = colnames(clr_corr_df),
    type = node_types  # 根据特征来源分配节点类型
  )
  # # 保存节点和边列表
  write.csv(nodes, paste0("cluster_", i, "_nodes.csv"), row.names = FALSE)
  write.csv(edges, paste0("cluster_", i, "_edges.csv"), row.names = FALSE)
  # 移除孤立节点
  connected_nodes <- unique(c(edges$source, edges$source))
  nodes <- nodes[nodes$id %in% connected_nodes, ]
  
  # 过滤低相关性边
  edges <- edges[abs(edges$weight) > 0.5, ]
  # ==== 新增代码开始 ==== <<<<
  # 获取节点类型映射字典
  type_dict <- setNames(nodes$type, nodes$id)
  
  # 给边添加类型属性
  edges$source_type <- type_dict[edges$source]
  edges$target_type <- type_dict[edges$target]
  
  # 过滤跨类型边（保留不同类型节点间的边）
  edges <- subset(edges, source_type != target_type)
  edges$source_type <- edges$target_type <- NULL  # 删除临时列
  # ==== 新增代码结束 ==== <<<<
  # 在保存节点和边之前添加以下过滤逻辑
  filtered <- filter_low_degree(nodes, edges, min_degree=5)

  # 更新过滤后的节点和边
  nodes <- filtered$nodes
  edges <- filtered$edges
  print(names(nodes))
  names(nodes) <- c('Label','type')
  names(edges) <- c("source","target","weight","p")
  # edges <- edges[,c('from','to')]
  # 移除孤立节点（二次保险）
  connected_nodes <- unique(c(edges$source, edges$target))
  nodes <- nodes[nodes$Label %in% connected_nodes, ]
  nodes <- merge(nodes, row_annotation, by.x = 'Label', by.y = 'id')
  nodes$Feature_Type <- NULL
  nodes <- nodes %>% mutate(Polygon = case_when(
    type == "Gene" ~ 1,
    type == "Pathway" ~ 3,
    type == "Species" ~ 4
  ))
  nodes$ID <- nodes$Label
  # nodes <- nodes %>%
  #   mutate(Label = case_when(
  #     type == "Gene" ~ "",
  #     type == "Pathway" ~ "",
  #     TRUE ~ Label
  #   ))
  nodes$Color <- NULL
  nodes <- nodes %>%
    mutate(Color = case_when(
      CD_Status == "CD_enriched" ~ "#ab6355",
      CD_Status == "CD_depleted" ~ "#90a5a7"))
  
  edges$cor_type <- ifelse(edges$weight>0,'positive','negative')
  print(range(edges$weight))
  edges$weight <- round(abs(edges$weight),2)
  edges$p <- NULL
  edges$r <- NULL

  edges <- edges %>%
    distinct()
  print(range(edges$weight))
  check_bidirectional_edges(edges)

  # print(range(edges$r))
  # print(dim(edges[edges$r<0,]))
  # 保存过滤后的文件（新增_filt后缀以示区别）
  write.csv(nodes, paste0("cluster_", i, "_nodes_filt.csv"), row.names=FALSE)
  write.csv(edges, paste0("cluster_", i, "_edges_filt.csv"), row.names=FALSE)

}



