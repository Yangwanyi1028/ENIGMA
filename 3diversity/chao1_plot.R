rm(list = ls())
setwd('/Users/yangwanyi/Downloads/ENIGMA_STUDY_DATA/3diversity/')
library(vegan)
library(tidyverse)
library(ggplot2)
library(ggpubr)

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
# 计算Chao1 Alpha多样性
chao1_diversity <- estimateR(species_abundance_int)

# 将Chao1多样性与组信息合并
chao1_df <- data.frame(
  study_id = rownames(species_abundance_int),
  chao1 = chao1_diversity["S.chao1", ]
)

# 合并组信息
chao1_df <- merge(chao1_df, group_info, by = "study_id")

# 比较CD患者与健康对照组的Chao1多样性
chao1_df <- chao1_df %>%
  mutate(Region = case_when(
    Region == "Australia" ~ "AUS", 
    Region == "Hong Kong" ~ "HK", 
    Region == "Mainland China" ~ "MC", 
  ))

region_color <- list("AUS" = "#85a894",
                     "HK" = "#3f6b7e",
                     "MC" = "#fce4a8")

p1 <- ggplot(chao1_df, aes(x = Group, y = chao1, fill = Region)) +  # 设置 fill 为 newregion
  geom_boxplot(outlier.shape = NA) +
  # geom_jitter(width = 0.2, alpha = 0.6) +
  stat_compare_means(
    method = "kruskal.test",  # 非参数检验
    label = "p.format",
    label.y = 500
  ) +
  labs(title = "Chao1 Alpha Diversity Across Regions",x = 'Group',y = 'Chao1 Diversity') +
  scale_fill_manual(values = unlist(region_color))  +# 使用自定义颜色 
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16),  # 设置标题字体大小
    axis.title.x = element_text(size = 14),  # 设置x轴标题字体大小
    axis.title.y = element_text(size = 14),  # 设置y轴标题字体大小
    axis.text.x = element_text(size = 12),   # 设置x轴刻度标签字体大小
    axis.text.y = element_text(size = 12),   # 设置y轴刻度标签字体大小
    legend.title = element_text(size = 12),  # 设置图例标题字体大小
    legend.text = element_text(size = 10)    # 设置图例标签字体大小
  ) +
  scale_y_continuous(limits = c(0, 550))


ggsave('chao1_group.png',width = 6,height = 6)

group_color <- list("CD active" = "firebrick",# "#9467BD"
                    "CD inactive" = "#f96b0a",# "#AEC7E8"
                    "Control" = "#90a5a6")
# list("CD active" = "#9467BD",
#                     "CD inactive" = "#AEC7E8",
#                     "Control" = "#BCBD22")

p2 <- ggplot(chao1_df, aes(x = Region  , y = chao1, fill = Group )) +  # 设置 fill 为 newregion
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  # geom_jitter(width = 0.2, alpha = 0.6) +
  stat_compare_means(
    method = "kruskal.test",  # 非参数检验
    label = "p.format",
    label.y = 500
  ) +
  labs(title = "Chao1 Alpha Diversity Across Diease Groups",x = 'Region',y = 'Chao1 Diversity') +
  scale_fill_manual(values = unlist(group_color))  +# 使用自定义颜色 
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16),  # 设置标题字体大小
    axis.title.x = element_text(size = 14),  # 设置x轴标题字体大小
    axis.title.y = element_text(size = 14),  # 设置y轴标题字体大小
    axis.text.x = element_text(size = 12),   # 设置x轴刻度标签字体大小
    axis.text.y = element_text(size = 12),   # 设置y轴刻度标签字体大小
    legend.title = element_text(size = 12),  # 设置图例标题字体大小
    legend.text = element_text(size = 10)    # 设置图例标签字体大小
  ) +
  scale_y_continuous(limits = c(0,550))


# 保存图形
ggsave(p2, file = 'chao1_region.png', width = 6, height = 6)
