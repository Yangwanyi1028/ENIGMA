setwd('/home/wanyiyang/ENIGMA_STUDY_DATA/7association/vfg_asso/')
library(ggplot2)
library(pheatmap)
library(dplyr)
library(readr)
library(tidyr)
library(grid)
library(tibble)

# 自定义函数：保存 pheatmap 为 PNG
save_pheatmap_png <- function(x, filename, width = 1800, height = 3000, res = 300) {
  png(filename, width = width, height = height, res = res)
  grid::grid.draw(x)
  dev.off()
}

# 自定义函数：绘制热图
plotAssociationsDag3HeatmapV2 <- function(inData, phenosToPlot, statToPlot, featuresToPlot, 
                                          nrFeaturesToPlot = 50, nrColors = 13, 
                                          sortPhenos = FALSE, retData = FALSE, fdrCutoff = 0.05) {
  # 筛选需要绘制的表型和特征
  inData <- inData %>%
    filter(phenotype %in% phenosToPlot & vfg_gene %in% featuresToPlot)
  
  # 创建矩阵：行是特征，列是表型
  heatmap_data <- inData %>%
    select(phenotype, vfg_gene, !!sym(statToPlot)) %>%
    tidyr::pivot_wider(
      names_from = phenotype, 
      values_from = !!sym(statToPlot),
      values_fn = mean,  # 如果有重复值，取平均值
      values_fill = NA   # 填充缺失值为 NA
    ) %>%
    column_to_rownames("vfg_gene")
  
  # 计算 Z-score 标准化
  heatmap_data <- scale(heatmap_data, center = TRUE, scale = TRUE)
  
  # 排序
  if (sortPhenos) {
    heatmap_data <- heatmap_data[, order(colnames(heatmap_data))]
  }
  
  # 设置颜色
  paletteLength <- nrColors
  myColor <- colorRampPalette(c("Orange", "white", "Darkblue"))(paletteLength)
  myBreaks <- c(seq(min(heatmap_data, na.rm = TRUE), 0, length.out = ceiling(paletteLength / 2) + 1), 
                seq(max(heatmap_data, na.rm = TRUE) / paletteLength, max(heatmap_data, na.rm = TRUE), length.out = floor(paletteLength / 2)))

  # 添加显著性标注
  significance_data <- inData %>%
    select(phenotype, vfg_gene, fdr) %>%
    tidyr::pivot_wider(
      names_from = phenotype, 
      values_from = fdr,
      values_fn = mean,  # 如果有重复值，取平均值
      values_fill = NA   # 填充缺失值为 NA
    ) %>%
    column_to_rownames("vfg_gene")
  
  # 根据 FDR 添加 +, - 符号
  significance_data <- ifelse(significance_data < fdrCutoff, "+", "")
  significance_data[heatmap_data < 0 & significance_data == "+"] <- "-"  # 如果值为负数，标注为 "-"
  significance_data[is.na(heatmap_data)] <- ""  # 缺失值不显示标注

  # 绘制热图
  phm <- pheatmap(heatmap_data, 
                  color = myColor,
                  breaks = myBreaks,
                  cluster_rows = FALSE,  # 不聚类行
                  cluster_cols = FALSE,  # 不聚类列
                  # main = paste("Heatmap of", statToPlot, "(Z-score)"),
                  silent = TRUE,
                  fontsize = 10,  # 字体大小
                  border_color = "#EEEEEE",  # 边框颜色
                  na_col = "white",  # 缺失值颜色
                  display_numbers = significance_data,  # 显著性标注
                  fontsize_number = 8)  # 标注字体大小
  
  if (retData) {
    return(list(heatmap_data, phm))
  } else {
    return(phm)
  }
}

# 定义一个函数来绘制热图
generate_heatmap <- function(data_file, output_file, title, stat_to_plot = "estimate") {
  # 读取数据
  data <- read_csv(data_file, show_col_types = FALSE)
  
  # 确保数据格式符合要求
  if (!"metadata_var" %in% colnames(data) || !"vfg_gene" %in% colnames(data)) {
    stop("数据文件缺少必要的列：metadata_var 或 vfg_gene")
  }
  
  # 准备数据
  data$phenotype <- data$metadata_var  # 将 metadata_var 作为 phenotype
  data$effect.size <- data$estimate   # 将 estimate 作为 effect size
  data$FDR <- data$fdr                # 将 fdr 作为 FDR
  
  # 提取需要绘制的变量和特征
  phenos_to_plot <- unique(data$phenotype)
  features_to_plot <- unique(data$vfg_gene)
  
  # 调用自定义热图函数
  phm <- plotAssociationsDag3HeatmapV2(
    inData = data,
    phenosToPlot = phenos_to_plot,
    statToPlot = stat_to_plot,
    featuresToPlot = features_to_plot,
    nrFeaturesToPlot = length(features_to_plot),  # 绘制所有特征
    nrColors = 13,  # 配色梯度
    sortPhenos = FALSE,  # 不排序
    retData = TRUE
  )
  
  # 保存热图
  save_pheatmap_png(phm[[2]], filename = output_file, width = 1800, height = 3000, res = 300)
}

# 定义输入文件、输出文件和标题
heatmap_files <- list(
  list(file = "biomarkers_pathology_associations.csv", output = "biomarkers_pathology_heatmap.png", title = "Biomarkers Pathology Heatmap"),
  list(file = "disease_associations.csv", output = "disease_heatmap.png", title = "Disease Heatmap"),
  list(file = "earlylife_diet_associations.csv", output = "earlylife_diet_heatmap.png", title = "Early Life Diet Heatmap"),
  list(file = "medication_associations.csv", output = "medication_heatmap.png", title = "Medication Heatmap"),
  list(file = "recent_additives_associations.csv", output = "recent_additives_heatmap.png", title = "Recent Additives Heatmap"),
  list(file = "threeday_nutrients_associations.csv", output = "threeday_nutrients_heatmap.png", title = "Three-Day Nutrients Heatmap")
)

# 遍历文件并生成热图
for (heatmap in heatmap_files) {
  generate_heatmap(heatmap$file, heatmap$output, heatmap$title)
}
