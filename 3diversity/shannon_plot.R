rm(list = ls())
setwd('/Users/yangwanyi/Downloads/ENIGMA_STUDY_DATA/3diversity/')
library(dplyr)
library(plyr)
library(vegan)    # 生态学分析常用包
library(ggplot2)  # 绘图包
library(ggpubr) # PCoA分析包

class_colors <- c("#1F77B4","#FF7F0E","#E377C2","#2CA02C",
                  "#9467BD","#17BECF","#AEC7E8","#D62728",
                  "#7F7F7F","#BCBD22","pink")

# metadata
metadata_479 = read.csv('../1maketable/metadata_filter.csv',row.names = 1)
sample_479 <- metadata_479$study_id


# 物种丰度表
species_df = readxl::read_excel('../source_data/merged_species_abundance_mp4_20240614_605samples.xlsx')
new_species_df <- species_df[species_df$clade_name %in% sample_479,]
dim(new_species_df)
new_species_df <- new_species_df[,-c(2,3)]

metadata1 <- metadata_479[,c('study_id','Group','Gender','Smoke',
                             'Region','History.of.bowel.resection.surgery',
                             "C.Reactive.Protein",'Inflammatory','Stricturing',
                             'Penetrating')]

new_species_meta_df <- merge(metadata1,new_species_df,by.x = 'study_id',by.y = 'clade_name')
new_species_meta_df[1:4,1:4]

df = new_species_meta_df
write.csv(new_species_df,'new_species_df.csv')
write.csv(metadata1,'new_metadata.csv',row.names = F)

# 1. 过滤低丰度物种（保留在>10%样本中出现的物种）
meta_num = dim(metadata1)[2]+1
keep <- colSums(df[,meta_num:ncol(df)] > 0) >= nrow(df)*0.1
filtered_df <- df[, c(rep(TRUE,meta_num-1), keep)]  # 保留前3列元数据

# 2. 添加伪计数处理零值（适用于Shannon计算）
pseudo_df <- filtered_df
pseudo_df[,meta_num:ncol(pseudo_df)] <- pseudo_df[,meta_num:ncol(pseudo_df)] + 0.001

# 3. 计算Shannon多样性指数
shannon <- diversity(pseudo_df[,meta_num:ncol(pseudo_df)], index = "shannon")
analysis_df <- cbind(filtered_df[,1:meta_num-1], Shannon = shannon)

# 定义颜色
region_color <- list("AUS" = "#85a894",
                     "HK" = "#3f6b7e",
                     "MC" = "#fce4a8")

# 绘制分组箱线图
p1 <- ggplot(analysis_df, aes(x = Group, y = Shannon, fill = Region)) +  # 设置 fill 为 newregion
  geom_boxplot(outlier.shape = NA) +
  # geom_jitter(width = 0.2, alpha = 0.6) +
  stat_compare_means(
    method = "kruskal.test",  # 非参数检验
    label = "p.format",
    label.y = max(analysis_df$Shannon) * 1.1
  ) +
  labs(title = "Shannon Alpha Diversity Across Regions",x = 'Group',y = 'Shannon Diversity') +
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
  )

# 保存图形
ggsave(p1, file = 'shannon_group.pdf', width = 6, height = 6)

group_color <- list("CD active" = "firebrick",# "#9467BD"
                     "CD inactive" = "#f96b0a",# "#AEC7E8"
                     "Control" = "#90a5a6") #"firebrick","#AB6355","#90a5a6"

# "firebrick","#1F77B4","#90a5a6"


p2 <- ggplot(analysis_df, aes(x = Region , y = Shannon, fill = Group )) +  # 设置 fill 为 newregion
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  # geom_jitter(width = 0.2, alpha = 0.6) +
  stat_compare_means(
    method = "kruskal.test",  # 非参数检验
    label = "p.format",
    label.y = max(analysis_df$Shannon) * 1.1
  ) +
  labs(title = "Shannon Alpha Diversity Across Diease Groups",x = 'Region',y = 'Shannon Diversity') +
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
  )

# 保存图形
ggsave(p2, file = 'shannon_region.pdf', width = 6, height = 6)

# 
# # 定义分组变量
# group_vars <- c('Group', 'Region', 
#                 # 'Gender', 'Smoke', 
#                 # 'History.of.bowel.resection.surgery',
#                 # "C.Reactive.Protein",'Inflammatory','Stricturing',
#                 # 'Penetrating'
#                 )
# 
# # 定义每个分组变量的颜色
# group_colors <- list(
#   # Group = c("Control" = "#85a894", "Case" = "#3f6b7e"),
#   # Gender = c("Male" = "#d6604d", "Female" = "#4393c3"),
#   # Smoke = c("Current smoker" = "#ace4a8", "Never smoker" = "#4393c3", "Past/Ex smoker" = "#fce4a8"),
#   Region = c("Australia" = "#85a894", "Hong Kong" = "#3f6b7e", "China" = "#fce4a8"),
#   # History.of.bowel.resection.surgery = c("Yes" = "#d6604d", "No" = "#4393c3"),
#   # C.Reactive.Protein = c("<10" = "#4393c3", '>=10' = "#d6604d"),
#   # Inflammatory = c("Yes" = "#d6604d", "No" = "#4393c3"),
#   # Stricturing = c("Yes" = "#d6604d", "No" = "#4393c3"),
#   # Penetrating = c("Yes" = "#d6604d", "No" = "#4393c3")
# )
# 
# # 循环绘制箱线图
# for (var in group_vars[-1]) {
#   p <- ggplot(analysis_df, aes(x = Group, y = Shannon, fill = !!sym(var))) +
#     geom_boxplot(outlier.shape = NA) +
#     geom_jitter(width = 0.2, alpha = 0.6) +
#     stat_compare_means(
#       method = "kruskal.test",  # 非参数检验
#       label = "p.format",
#       label.y = max(analysis_df$Shannon) * 1.1
#     ) +
#     labs(title = paste("Alpha Diversity by", var)) +
#     scale_fill_manual(values = group_colors[[var]]) +  # 使用自定义颜色
#     theme_minimal()
#   
#   # 保存图形
#   ggsave(p, file = paste0(var, '_shannon.png'), width = 6, height = 6)
# }
# 
# # # 分面显示各地区的疾病组对比
# # p2 <- ggplot(analysis_df, aes(x=cd_act, y=Shannon, fill=cd_act)) +
# #   geom_violin(alpha=0.7) +
# #   geom_boxplot(width=0.2, fill="white") +
# #   facet_wrap(~newregion, scales="free_x") +
# #   stat_compare_means(
# #     comparisons = list(c("CD", "HC")),  # 假设有健康对照(HC)
# #     method = "wilcox.test"
# #   )
# # ggsave(p2,file='region_diseasestatus_shannon.pdf',width = 10,height = 10)
# # 
# # 
# # 
# # ## Observed spiecies
# # ## HMM? FDR？
