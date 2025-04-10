rm(list = ls())
library(VennDiagram)
setwd('/Users/yangwanyi/Downloads/ENIGMA_STUDY_DATA/4maaslin2/')
venn_aus_sp <- read.csv('output_Region_species_AUS/significant_results.tsv',sep='\t',header = T)
venn_aus_hk <- read.csv('output_Region_species_HK/significant_results.tsv',sep='\t',header = T)
venn_aus_mc <- read.csv('output_Region_species_MC/significant_results.tsv',sep='\t',header = T)

dim(venn_aus_sp)
dim(venn_aus_hk)
dim(venn_aus_mc)
# 假设 venn_aus_sp$feature, venn_aus_hk$feature, venn_aus_mc$feature 是你的三个集合
set1 <- venn_aus_sp$feature
set2 <- venn_aus_hk$feature
set3 <- venn_aus_mc$feature

# 绘制韦恩图
venn.diagram(
  x = list(set1, set2, set3),
  category.names = c("AUS", "HK", "MC"),
  filename = "venn_diagram_unfilter.png",  # 保存为PNG文件
  output = TRUE,
  compression = "lzw",
  lwd = 2,width = 1000,height = 1000,resolution = 300,
  lty = 'blank',
  fill = c("#85a894", "#3f6b7e", "#fce4a8"),
  cex = 1,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  disable.logging = T
)


venn_aus_sp <- read.csv('output_Region_species_AUS_filter/significant_results.tsv',sep='\t',header = T)
venn_aus_hk <- read.csv('output_Region_species_HK_filter/significant_results.tsv',sep='\t',header = T)
venn_aus_mc <- read.csv('output_Region_species_MC_filter/significant_results.tsv',sep='\t',header = T)

dim(venn_aus_sp)
dim(venn_aus_hk)
dim(venn_aus_mc)
# 假设 venn_aus_sp$feature, venn_aus_hk$feature, venn_aus_mc$feature 是你的三个集合
set1 <- venn_aus_sp$feature
set2 <- venn_aus_hk$feature
set3 <- venn_aus_mc$feature

# 绘制韦恩图
venn.diagram(
  x = list(set1, set2, set3),
  category.names = c("AUS", "HK", "MC"),
  filename = "venn_diagram_filter.png",  # 保存为PNG文件
  output = TRUE,
  compression = "lzw",
  lwd = 2,width = 1000,height = 1000,resolution = 300,
  lty = 'blank',
  fill = c("#85a894", "#3f6b7e", "#fce4a8"),
  cex = 1,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  disable.logging = T
)


