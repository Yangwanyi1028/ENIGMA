rm(list = ls())
setwd('/Users/yangwanyi/Downloads/ENIGMA_STUDY_DATA/2pcoa/')
library(dplyr)
library(plyr)
library(vegan)    # 生态学分析常用包
library(ggplot2)  # 绘图包
library(ape)      # PCoA分析包

sample_479 = read.csv('../00sample_name_mapping/sample_479.csv')
sample_479 <- sample_479$clade_name
species_df = readxl::read_excel('../source_data/merged_species_abundance_mp4_20240614_605samples.xlsx')

new_species_df <- species_df[species_df$clade_name %in% sample_479,]
dim(new_species_df)


group_df = new_species_df[,c(1,2,3)]
species_abun_df = new_species_df[,-c(1,2,3)]

meta <- read.csv('../1maketable/ywy_metadata.csv')

group_df1 = merge(meta[,c('study_id','cdai')],group_df,by.x = 'study_id',by.y = 'clade_name')

group_df1 <- group_df1 %>%
  mutate(Group = case_when(
    Group %in% c("HHM", "FDR_HHM", "FDR") ~ "Controls",  # 将"HHM", "FDR_HHM", "FDR"替换为"Controls"
    TRUE ~ Group  # 其他值保持不变
  )) %>%
  mutate(Group = case_when(
    cdai < 150 ~ "CD inactive",  # 如果 cdai < 150，设置为 "CD inactive"
    cdai >= 150 ~ "CD active" ,  # 如果 cdai >= 150，设置为 "CD active"
    TRUE ~ Group  # 其他值保持不变
  ))
group_df1 = group_df1[,-2]
# 查看修改后的数据
print(group_df1)

write.csv(group_df1,'group_df.csv',row.names = F)
rownames(species_abun_df) <- new_species_df$clade_name
write.csv(species_abun_df,'species_abun_df.csv',row.names = T)

# 步骤1：读取数据（示例数据结构）
# 物种丰度矩阵（行：样本，列：物种）
species_matrix <- as.matrix(species_abun_df)
# 分组信息（与样本顺序一致）
group_info <- group_df

# 步骤2：计算距离矩阵（Bray-Curtis距离）
dist_matrix <- vegdist(species_matrix, method="bray")

# 步骤3：执行PCoA分析
pcoa_result <- pcoa(dist_matrix, correction="none")

# 步骤4：提取坐标和解释度
coordinates <- pcoa_result$vectors[,1:2]  # 取前两轴
explained_var <- round(pcoa_result$values$Relative_eig[1:2]*100, 2)

# 步骤5：合并分组信息
df_plot <- data.frame(
  Sample = rownames(coordinates),
  PCo1 = coordinates[,1],
  PCo2 = coordinates[,2],
  Group = group_info$Group
)

# 步骤6：按照示例代码风格绘制图形
level_order <- unique(df_plot$Group)  # 确定 Group 的顺序
level_order <- factor(1:length(level_order), labels = level_order)
df_plot$Group <- factor(df_plot$Group, levels = levels(level_order))

# 计算每个组的凸包
find_hull <- function(df_plot) df_plot[chull(df_plot$PCo1, df_plot$PCo2),]
hulls <- ddply(df_plot, "Group", find_hull)

# 计算每个组的中心点
PCo1_mean <- tapply(df_plot$PCo1, df_plot$Group, mean)
PCo2_mean <- tapply(df_plot$PCo2, df_plot$Group, mean)
mean_point <- data.frame(PCo1_mean, PCo2_mean)

# 绘制图形
pcoa_plot <- ggplot(df_plot, aes(PCo1, PCo2, color = Group, shape = Group)) +
  theme_classic() +  # 去掉背景框
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
  geom_point(size = 4) +  # 设置点的透明度、大小
  geom_polygon(data = hulls, alpha = 0.2, aes(fill = Group, color = Group), size = 0.2, show.legend = FALSE) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", 'red', 'blue')) +  # 自定义颜色
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", 'red', 'blue')) +  # 自定义填充颜色
  scale_shape_manual(values = c(20, 20, 20, 20, 20)) +  # 设置点的形状
  scale_x_continuous(limits = c(min(df_plot$PCo1) - 0.1, max(df_plot$PCo1) + 0.1)) +
  scale_y_continuous(limits = c(min(df_plot$PCo2) - 0.1, max(df_plot$PCo2) + 0.1)) +
  theme(panel.grid = element_line(color = 'black', linetype = 1, size = 0.1),
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.title = element_blank()
  ) +
  labs(x = paste0('PCo1 (', explained_var[1], '%)'), y = paste0('PCo2 (', explained_var[2], '%)'))

# 显示图形
pcoa_plot

# 保存图形
png(filename = 'pcoa_plot_group.png', width = 2500, height = 2000, res = 300)
pcoa_plot
dev.off()








############################
          #Region#
############################
colors <- c(AUS = "#85a894", HK = "#3f6b7e", MC = "#fce4a8") ##?? kunming??

# 步骤5：合并分组信息
df_plot <- data.frame(
  Sample = rownames(coordinates),
  PCo1 = coordinates[,1],
  PCo2 = coordinates[,2],
  Region = group_info$Region
)

# 步骤6：绘制ggplot2图形（按照示例代码风格）
level_order <- unique(df_plot$Region)  # 确定level顺序
level_order <- factor(1:length(level_order), labels = level_order)
df_plot$Region <- factor(df_plot$Region, levels = levels(level_order))


# 计算每个组的凸包
find_hull <- function(df_plot) df_plot[chull(df_plot$PCo1, df_plot$PCo2),]
hulls <- ddply(df_plot, "Region", find_hull)

# 计算每个组的中心点
PCo1_mean <- tapply(df_plot$PCo1, df_plot$Region, mean)
PCo2_mean <- tapply(df_plot$PCo2, df_plot$Region, mean)
mean_point <- data.frame(PCo1_mean, PCo2_mean)

# 绘制图形
pcoa_plot <- ggplot(df_plot, aes(PCo1, PCo2, color=Region, shape=Region)) +
  theme_classic() +  # 去掉背景框
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
  geom_polygon(data = hulls, alpha = 0.2, aes(fill=Region, color=Region), size=0.2, show.legend = FALSE) +
  geom_point(size = 4) +  # 设置点的透明度、大小
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  scale_shape_manual(values = c(20, 20, 20)) +
  scale_x_continuous(limits = c(min(df_plot$PCo1) - 0.1, max(df_plot$PCo1) + 0.1)) +
  scale_y_continuous(limits = c(min(df_plot$PCo2) - 0.1, max(df_plot$PCo2) + 0.1)) +
  theme(panel.grid = element_line(color = 'black', linetype = 1, size = 0.1),
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.title = element_blank()
  ) +
  labs(x = paste0('PCo1 (', explained_var[1], '%)'), y = paste0('PCo2 (', explained_var[2], '%)'))

# 显示图形
pcoa_plot

# 保存图形
png(filename = 'pcoa_plot_region.png', width = 2500, height = 2000, res = 300)
pcoa_plot
dev.off()