rm(list = ls())
setwd('/Users/yangwanyi/Downloads/ENIGMA_STUDY_DATA/4maaslin2/')
cd <- read.csv('output_Group_species/significant_results.tsv',sep = '\t')
cd_enriched <- cd[cd$coef>0,]
cd_depleted <- cd[cd$coef<0,]
venn_aus_sp <- read.csv('output_Region_species_AUS_filter/significant_results.tsv',sep='\t',header = T)
venn_aus_hk <- read.csv('output_Region_species_HK_filter/significant_results.tsv',sep='\t',header = T)
venn_aus_mc <- read.csv('output_Region_species_MC_filter/significant_results.tsv',sep='\t',header = T)

# 假设 venn_aus_sp$feature, venn_aus_hk$feature, venn_aus_mc$feature 是你的三个集合
set1_enriched  <- intersect(venn_aus_sp$feature,cd_enriched$feature)
set2_enriched  <- intersect(venn_aus_hk$feature,cd_enriched$feature)
set3_enriched  <- intersect(venn_aus_mc$feature,cd_enriched$feature)


# 创建数据框
venn_data_enriched <- data.frame(
  feature = unique(c(set1_enriched, set2_enriched, set3_enriched)),
  AUS = as.numeric(unique(c(set1_enriched, set2_enriched, set3_enriched)) %in% set1_enriched),
  HK = as.numeric(unique(c(set1_enriched, set2_enriched, set3_enriched)) %in% set2_enriched),
  MC = as.numeric(unique(c(set1_enriched, set2_enriched, set3_enriched)) %in% set3_enriched)
)

# 绘制 UpSet 图
pdf("upset_plot_enriched.pdf", width = 6, height = 6)  # 设置 PDF 文件名和尺寸
grid.newpage()
upset(venn_data_enriched, 
      sets = c("AUS", "HK", "MC"),
      main.bar.color = "#ab6354",    # 设置主柱状图的颜色
      sets.bar.color = c("#85A894", "#3E6B7E", "#FCE4A8"),
      show.numbers = "yes",          # 显示数字
      nintersects = NA,
      text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5)  # 调整文本大小
      )
# 关闭 PDF 设备
dev.off()

set1_depleted <- intersect(venn_aus_sp$feature,cd_depleted$feature)
set2_depleted <- intersect(venn_aus_hk$feature,cd_depleted$feature)
set3_depleted <- intersect(venn_aus_mc$feature,cd_depleted$feature)



# 创建数据框
venn_data_depleted <- data.frame(
  feature = unique(c(set1_depleted, set2_depleted, set3_depleted)),
  AUS = as.numeric(unique(c(set1_depleted, set2_depleted, set3_depleted)) %in% set1_depleted),
  HK = as.numeric(unique(c(set1_depleted, set2_depleted, set3_depleted)) %in% set2_depleted),
  MC = as.numeric(unique(c(set1_depleted, set2_depleted, set3_depleted)) %in% set3_depleted)
)

# 打开 PDF 设备
pdf("upset_plot_depleted.pdf", width = 6, height = 6)  # 设置 PDF 文件名和尺寸
grid.newpage()
# 绘制 UpSet 图
upset(
  venn_data_depleted,
  sets = c("AUS", "HK", "MC"),  
  main.bar.color = "#90a5a6",    # 设置主柱状图的颜色
  sets.bar.color = c("#85A894", "#3E6B7E", "#FCE4A8"),  # 设置集合柱状图的颜色
  show.numbers = "yes",          # 显示数字
  text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5)  # 调整文本大小
)

# 关闭 PDF 设备
dev.off()

