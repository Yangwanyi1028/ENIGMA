rm(list = ls())
setwd('/Users/yangwanyi/Downloads/ENIGMA_STUDY_DATA/3diversity/')
library(vegan)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(dplyr)


group_info <- read.csv('../1maketable/metadata_filter.csv',row.names = 1)

species_abundance <- readxl::read_excel('../source_data/merged_species_abundance_mp4_20240614_605samples.xlsx')
species_abundance <- species_abundance[species_abundance$clade_name %in% group_info$study_id,]
ordered_indices <- match(group_info$study_id, species_abundance$clade_name)

# 根据排序后的索引重新排列 species_abundance
species_abundance_ordered <- species_abundance[ordered_indices, ]
table(species_abundance_ordered$clade_name == group_info$study_id)

species_abundance_ordered <- species_abundance_ordered[,4:ncol(species_abundance)]
species_abundance_ordered <- as.matrix(species_abundance_ordered)
rownames(species_abundance_ordered) <- group_info$study_id
species_abundance_int <- round(species_abundance_ordered * 1000)

species_abundance_ordered[1:4,1:4]
species_abundance_int[1:4,1:4]
group_info[1:4,1:4]

# 计算Bray-Curtis距离，样本与样本之间的差异, 0-1，数值越大表示差异越大
bray_curtis_dist <- vegdist(species_abundance_ordered, method = "bray")


# 1. 将Bray-Curtis距离矩阵转换为数据框
dist_df <- as.data.frame(as.matrix(bray_curtis_dist))
dist_df$Sample1 <- rownames(dist_df)
dist_df_long <- dist_df %>%
  tidyr::gather(key = "Sample2", value = "Distance", -Sample1) %>%
  filter(Sample1 < Sample2)  # 避免重复计算（下三角矩阵）

# 2. 添加分组信息
dist_df_long <- dist_df_long %>%
  left_join(group_info %>% select(study_id, Group), by = c("Sample1" = "study_id")) %>%
  left_join(group_info %>% select(study_id, Group), by = c("Sample2" = "study_id")) %>%
  rename(Group1 = Group.x, Group2 = Group.y) %>%
  mutate(ComparisonType = case_when(
    Group1 == Group2 ~ paste0("Within-", Group1),  # 组内比较
    TRUE ~ paste(Group1, Group2, sep = "-vs-")     # 组间比较
  ))

# 3. 筛选有意义的比较（排除重复的组间比较，如A-vs-B和B-vs-A）
unique_comparisons <- dist_df_long %>%
  mutate(ComparisonType = ifelse(Group1 < Group2,
                                 paste(Group1, Group2, sep = "-vs-"),
                                 paste(Group2, Group1, sep = "-vs-"))) %>%
  filter(Group1 != Group2) %>%
  distinct(ComparisonType) %>%
  pull(ComparisonType)

# 合并组内和组间数据
# 合并组内和组间数据（修复命名不一致问题）
plot_data <- dist_df_long %>%
  # 统一组间比较的命名顺序（按字母顺序排列组名）
  mutate(
    ComparisonType = case_when(
      # 组内比较
      Group1 == Group2 ~ paste0("Within-", Group1),
      # 组间比较：按字母顺序排列组名（确保命名唯一性）
      Group1 < Group2 ~ paste(Group1, Group2, sep = "-vs-"),
      TRUE ~ paste(Group2, Group1, sep = "-vs-")
    )
  ) %>%
  # 定义因子水平（包含所有可能的比较类型）
  mutate(
    ComparisonType = factor(
      ComparisonType,
      levels = c(
        "Within-CD active", "Within-CD inactive", "Within-Control",
        "CD active-vs-CD inactive", 
        "CD active-vs-Control", 
        "CD inactive-vs-Control"
      )
    )
  )

# 检查是否仍有NA
sum(is.na(plot_data$ComparisonType))  # 应为0


# 1. 定义需要比较的组间配对
comparisons <- list(
  c("Within-CD active", "Within-CD inactive"),
  c("Within-CD active", "Within-Control"),
  c("Within-CD inactive", "Within-Control"),
  c("CD active-vs-CD inactive", "CD active-vs-Control"),
  c("CD active-vs-CD inactive", "CD inactive-vs-Control"),
  c("CD active-vs-Control", "CD inactive-vs-Control")
)
# 绘制箱线图
p1 <- ggplot(plot_data, aes(x = ComparisonType, y = Distance, fill = ComparisonType)) +
      geom_boxplot(alpha = 0.7) +
      labs(
        title = "Bray-Curtis Distance Comparisons",
        x = "Comparison Type",
        y = "Bray-Curtis Distance"
      ) +
      stat_compare_means(
        comparisons = comparisons,
        method = "wilcox.test",
        label = "p.signif",  # 显示星号（*表示p < 0.05，** < 0.01等）
        step.increase = 0.1,  # 调整标签的垂直位置避免重叠
        symnum.args = list(
          cutpoints = c(0, 0.001, 0.01, 0.05, 1),
          symbols = c("***", "**", "*", "ns")
        )
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1,size = 12),
        legend.position = "none",
        plot.title = element_text(size = 16),  # 设置标题字体大小
        axis.title.x = element_text(size = 14),  # 设置x轴标题字体大小
        axis.title.y = element_text(size = 14),  # 设置y轴标题字体大小
        axis.text.y = element_text(size = 12),   # 设置y轴刻度标签字体大小
        legend.title = element_text(size = 12),  # 设置图例标题字体大小
        legend.text = element_text(size = 10)    # 设置图例标签字体大小
      ) +
      scale_fill_manual(values = c(
        "Within-CD active" = "#BCBD22",
        "Within-CD inactive" = "#377EB8",
        "Within-Control" = "#4DAF4A",
        "CD active-vs-CD inactive" = "#9467BD",
        "CD active-vs-Control" = "#FF7F00",
        "CD inactive-vs-Control" = "#AEC7E8"
      ))

p2 <- p1 + 
  coord_cartesian(ylim = c(0, 1))
ggsave('bray-curtis_boxplot.pdf',width = 8,height = 6)