setwd('/home/wanyiyang/ENIGMA_STUDY_DATA/7association/cb_all/')
library(dplyr)       # For data manipulation (mutate, left_join, etc.)
library(purrr)       # For functional programming (map_dfr)
library(ggplot2)     # For data visualization (ggplot)
library(readr)       # For reading/writing CSV files (read_csv, write_csv)
library(tidyverse)
library(broom)
library(compositions) # For CLR transformation
library(survminer)    # ggforest 函数所在包
library(pheatmap)
library(readr)
library(tidyr)
library(grid)
library(tibble)

# 1. Load and prepare data ---------------------------------------------------
# Load clinical data
clinical_data <- read.csv("../category_data/cb_all.csv", row.names = 1) %>%
  mutate(across(c(
                  "Age.4.12.months.Home.Grown.Fruit.Veg", "Age.1.5.years.Home.Grown.Fruit.Veg",
                  "Age.5.10.years.Home.Grown.Fruit.Veg", "Age.10.18.years.Home.Grown.Fruit.Veg",
                  "Age.4.12.months.Milk.Beverages", "Age.1.5.years.Milk.Beverages",
                  "Age.5.10.years.Milk.Beverages", "Age.10.18.years.Milk.Beverages",
                  "Age.4.12.months.Processed.Meat.Seafood", "Age.1.5.years.Processed.Meat.Seafood",
                  "Age.5.10.years.Processed.Meat.Seafood", "Age.10.18.years.Processed.Meat.Seafood",
                  "Age.4.12.months.Processed.Carbs","Age.4.12.months.Processed.Fruit", "Age.1.5.years.Processed.Fruit",
                  "Age.5.10.years.Processed.Fruit",
                  "Age.10.18.years.Processed.Fruit", "Age.4.12.months.Processed.Veg", "Age.1.5.years.Processed.Veg",
                  "Age.5.10.years.Processed.Veg","Age.10.18.years.Processed.Veg", "Age.4.12.months.Fast.food", "Age.1.5.years.Fast.food",
                  "Age.5.10.years.Fast.food", "Age.10.18.years.Fast.food", 'Age.4.12.months.Soft.Drink', "Age.1.5.years.Soft.Drink",
                  "Age.5.10.years.Soft.Drink", "Age.10.18.years.Soft.Drink", 'Age.4.12.months.Packaged.Snacks',
                  "Age.1.5.years.Packaged.Snacks", "Age.5.10.years.Packaged.Snacks", 'Age.10.18.years.Packaged.Snacks'
                  ), as.factor))

# Load demographic data
demographic_data <- read.csv("../category_data/Demographics_Baseline.csv", row.names = 1) %>%
  select(sample_id, Age.at.enrolment, Gender, Smoking.status, Region, BMI, Alcohol.Consumption)

# Define analysis variables for clinical data
clinical_vars <- colnames(clinical_data)[-1] # 除去 sample_id 的所有列
continuous_vars <- c("Age.at.enrolment", "BMI", "P80", "CMC","CRN","AlSiO","SO32",
                     "TiO2","Asp","Suc","Sac","TotalAd","TotalAs","TotalEM",
                     "P80.Annual.Intake",   "CMC.Annual.Intake",  'CRN.Annual.Intake',
                     'AlSiO.Annual.Intake', "SO32.Annual.Intake", "TiO2.Annual.Intake", "Asp.Annual.Intake", "Suc.Annual.Intake",  
                     "Sac.Annual.Intake"   ,"TotalAd.Annual.Intake" ,"TotalAs.Annual.Intake", "TotalEM.Annual.Intake", "Protein.intake",
                     "Total.Fat.intake","Saturated.Fat.intake", "Monounsaturated.Fat.intake", "Polyunsaturated.Fat.intake",
                     "Omega.3.Fatty.Acid.intake", "Omega.6.Fatty.Acid.intake", "Cholesterol.Intake", "Alcohol.Intake",   "Fibre.Intake",
                     "Vitamin.B1","Vitamin.B2","Vitamin.B3","Vitamin.B6","Vitamin.B12","Vitamin.C",  
                     "Vitamin.E.Tocopherol.Acetate", "Beta.Carotene","Folic.Acid","Folate","Iron",       
                     "Zinc","Selenium","Magnesium","Caffeine","Vitamin.D")
categorical_vars <- setdiff(clinical_vars, continuous_vars)

# 定义列注释信息
clinical_annotations <- list(
  `Recent Addictives` = names(clinical_data)[8:31],
  `Threeday Nutrients` = names(clinical_data)[32:57],
  `Early Life Diet` = names(clinical_data)[58:96]
)

# 定义行注释信息
feature_status_annotations = read.csv('../../newsource_data20250319/Combined_Annotations_with_CD_Status.csv')
feature_status_annotations <- feature_status_annotations[,c(1,3)]
names(feature_status_annotations) <- c('features','CD_status')

# 自定义函数：保存 pheatmap 为 PNG，动态调整大小
save_pheatmap_png <- function(x, filename, nrow, ncol, res = 300) {
  # 根据行列数动态调整宽度和高度
  width <- max(3000, ncol * 50) + 1000  # 每列 50 像素，最小宽度 3000
  height <- max(3000, nrow * 50) # 每行 50 像素，最小高度 3000
  
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
    filter(phenotype %in% phenosToPlot & feature_name %in% featuresToPlot)
  
  # 创建矩阵：行是特征，列是表型
  heatmap_data <- inData %>%
    select(phenotype, feature_name, !!sym(statToPlot)) %>%
    tidyr::pivot_wider(
      names_from = phenotype, 
      values_from = !!sym(statToPlot),
      values_fn = mean,  # 如果有重复值，取平均值
      values_fill = NA   # 填充缺失值为 NA
    ) %>%
    column_to_rownames("feature_name")
  
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
    select(phenotype, feature_name, fdr) %>%
    tidyr::pivot_wider(
      names_from = phenotype, 
      values_from = fdr,
      values_fn = mean,  # 如果有重复值，取平均值
      values_fill = NA   # 填充缺失值为 NA
    ) %>%
    column_to_rownames("feature_name")
  
  # 根据 FDR 添加 +, - 符号
  significance_data <- ifelse(significance_data < fdrCutoff, "+", "")
  significance_data[heatmap_data < 0 & significance_data == "+"] <- "-"  # 如果值为负数，标注为 "-"
  significance_data[is.na(heatmap_data)] <- ""  # 缺失值不显示标注

  # 创建列注释数据框
  annotation_col <- data.frame(
    Category = rep("Other", length(phenosToPlot)),
    row.names = phenosToPlot
  )
  
  # 填充注释类别
  annotation_col$Category[phenosToPlot %in% clinical_annotations$`Recent Addictives`] <- "Recent Addictives"
  annotation_col$Category[phenosToPlot %in% clinical_annotations$`Threeday Nutrients`] <- "Threeday Nutrients"
  annotation_col$Category[phenosToPlot %in% clinical_annotations$`Early Life Diet`] <- "Early Life Diet"
  
  # 创建行注释数据框
  annotation_row <- feature_status_annotations %>%
    filter(features %in% rownames(heatmap_data)) %>%
    column_to_rownames("features")
  
  # 定义行注释颜色
  annotation_colors <- list(
    CD_status = c(
      "CD_enriched" = "#ec9191",  # 红色
      "CD_depleted" = "#65b6f8"   # 蓝色
    ),
    Category = c(
      "Recent Addictives"  = "#df7a7ba4",    # 红色
      "Threeday Nutrients" = "#5e7f9c",  # 蓝色
      "Early Life Diet" = "#abc281",        # 绿色
      "Other" = "grey"
    )
  )
  # 绘制热图
  phm <- pheatmap(heatmap_data, 
                  color = myColor,
                  breaks = myBreaks,
                  cluster_rows = FALSE,  # 不聚类行
                  cluster_cols = FALSE,  # 不聚类列
                  annotation_col = annotation_col,  # 添加列注释
                  annotation_row = annotation_row,  # 添加列注释
                  annotation_colors = annotation_colors,  # 注释颜色
                  fontsize = 10,  # 字体大小
                  border_color = "#EEEEEE",  # 边框颜色
                  na_col = "white",  # 缺失值颜色
                  display_numbers = significance_data,  # 显著性标注
                  fontsize_number = 8,  # 标注字体大小
                  angle_col = 315,
                  annotation_names_col = TRUE,  # 显示注释名称
                  annotation_legend = TRUE)     # 显示图例
  
  if (retData) {
    return(list(heatmap_data, phm))
  } else {
    return(phm)
  }
}

# 修改 generate_heatmap 函数，动态调整热图大小
generate_heatmap <- function(data_file, output_file, title, stat_to_plot = "estimate") {
  # 读取数据
  data <- read_csv(data_file, show_col_types = FALSE)
  
  # 确保数据格式符合要求
  if (!"metadata_var" %in% colnames(data) || !"feature_name" %in% colnames(data)) {
    stop("数据文件缺少必要的列：metadata_var 或 feature_name")
  }
  
  # 准备数据
  data$phenotype <- data$metadata_var  # 将 metadata_var 作为 phenotype
  data$effect.size <- data$estimate   # 将 estimate 作为 effect size
  data$FDR <- data$fdr                # 将 fdr 作为 FDR
  
  # 提取需要绘制的变量和特征
  phenos_to_plot <- unique(data$phenotype)
  features_to_plot <- unique(data$feature_name)
  
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
  
  # 动态调整图像大小
  nrow <- length(features_to_plot)  # 行数为特征数
  ncol <- length(phenos_to_plot)    # 列数为表型数
  
  # 保存热图
  save_pheatmap_png(phm[[2]], filename = output_file, nrow = nrow, ncol = ncol, res = 300)
}

# Define file paths and their identifiers
data_files <- list(
  arg = list(
    path = "../../newsource_data20250319/ARG_68_479samples.csv",
    prefix = "arg"
  ),
  path = list(
    path = "../../newsource_data20250319/path_68_479samples.csv",
    prefix = "path"
  ),
  species = list(
    path = "../../newsource_data20250319/species_101_479samples.csv",
    prefix = "species"
  ),
  vfg = list(
    path = "../../newsource_data20250319/VFG_173_479samples.csv",
    prefix = "vfg"
  )
)

# Function to process each dataset
process_dataset <- function(file_path, prefix) {
  # Load abundance data
  abundance_data <- read.csv(file_path, row.names = 1) %>%
    t() %>% 
    as.data.frame() %>%
    rownames_to_column("sample_id")
  
  feature_columns <- colnames(abundance_data)[-1]
  
  # Merge with clinical and demographic data
  merged_data <- clinical_data %>%
    left_join(demographic_data, by = "sample_id") %>%
    left_join(abundance_data, by = "sample_id")
  
  # Perform association analysis
  association_results <- map_dfr(feature_columns, function(feature) {
    map_dfr(clinical_vars, function(clinical_var) {
      tryCatch({
        if (clinical_var %in% continuous_vars) {
          # 连续变量：直接进行线性回归
          lm(reformulate(c(clinical_var), feature), 
             data = merged_data) %>%
            tidy() %>%
            filter(term == clinical_var) %>%
            mutate(
              metadata_var = clinical_var,
              feature_name = feature,
              fdr = p.adjust(p.value, method = "fdr")
            )
        } else {
          # 分类型变量：因子化后进行线性回归
          lm(reformulate(c(sprintf("factor(%s)", clinical_var)), feature), 
             data = merged_data) %>%
            tidy() %>%
            filter(grepl(clinical_var, term)) %>%
            mutate(
              metadata_var = clinical_var,
              feature_name = feature,
              fdr = p.adjust(p.value, method = "fdr")
            )
        }
      }, error = function(e) NULL)
    })
  })
  
  # 检查结果是否包含所有生物标志物变量
  print(unique(association_results$metadata_var)) # 调试：确保所有变量都被包含
  
  # 确保 association_results 包含置信区间列
  if (!("conf.low" %in% colnames(association_results)) || !("conf.high" %in% colnames(association_results))) {
    association_results <- association_results %>%
      mutate(
        conf.low = estimate - 1.96 * std.error, # 计算下置信区间
        conf.high = estimate + 1.96 * std.error # 计算上置信区间
      )
  }
  
  # 生成可视化
  association_plot <- ggplot(association_results, aes(estimate, -log10(p.value), 
                                                      color = fdr < 0.05, shape = metadata_var)) +
    geom_point(aes(size = abs(estimate)), alpha = 0.8) +
    geom_vline(xintercept = 0, linetype = 2) +
    scale_shape_manual(values = seq(16, 16 + length(unique(association_results$metadata_var)) - 1)) + # 动态生成形状
    labs(title = "ARG-Clinical Data Associations",
         x = "Effect Size (Estimate)",
         y = "-log10(p-value)",
         color = "FDR < 0.05",
         shape = "Clinical Variable")
  
  # 绘制森林图
  forest_plot <- ggplot(association_results, aes(x = estimate, y = metadata_var)) +
    geom_point(aes(color = fdr < 0.05), size = 3) + # 效应值点
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, color = fdr < 0.05), height = 0.2) + # 置信区间
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) + # 显著性颜色
    labs(
      title = "Forest Plot of ARG-Clinical Associations",
      x = "Effect Size (Estimate)",
      y = "Clinical Variable",
      color = "FDR < 0.05"
    ) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10)) # 调整 y 轴文本大小
  
  # 保存结果
  write_csv(association_results, paste0(prefix, "_clinical_associations.csv"))
  ggsave(paste0(prefix, "_clinical_volcano.png"), association_plot, width = 12, height = 8)
  ggsave(paste0(prefix, "_clinical_forest.png"), forest_plot, width = 12, height = 8)
  
  # Generate and save heatmap
  generate_heatmap(
    paste0(prefix, "_clinical_associations.csv"),
    paste0(prefix, "_clinical_associations_annotated.png"),
    paste0(toupper(prefix), " Clinical Associations")
  )
  
  return(association_results)
}

# 新增函数：处理 demographic_data 和 feature 文件的关联分析
process_demographic_association <- function(file_path, prefix) {
  # Load abundance data
  abundance_data <- read.csv(file_path, row.names = 1) %>%
    t() %>% 
    as.data.frame() %>%
    rownames_to_column("sample_id")
  
  feature_columns <- colnames(abundance_data)[-1]
  
  # Merge with demographic data
  merged_data <- demographic_data %>%
    left_join(abundance_data, by = "sample_id")
  
  # Perform association analysis
  association_results <- map_dfr(feature_columns, function(feature) {
    map_dfr(colnames(demographic_data)[-1], function(demographic_var) {
      tryCatch({
        if (demographic_var %in% continuous_vars) {
          # 连续变量：直接进行线性回归
          lm(reformulate(c(demographic_var), feature), 
             data = merged_data) %>%
            tidy() %>%
            filter(term == demographic_var) %>%
            mutate(
              metadata_var = demographic_var,
              feature_name = feature,
              fdr = p.adjust(p.value, method = "fdr")
            )
        } else {
          # 分类型变量：因子化后进行线性回归
          lm(reformulate(c(sprintf("factor(%s)", demographic_var)), feature), 
             data = merged_data) %>%
            tidy() %>%
            filter(grepl(demographic_var, term)) %>%
            mutate(
              metadata_var = demographic_var,
              feature_name = feature,
              fdr = p.adjust(p.value, method = "fdr")
            )
        }
      }, error = function(e) NULL)
    })
  })
  
  # 确保 association_results 包含置信区间列
  if (!("conf.low" %in% colnames(association_results)) || !("conf.high" %in% colnames(association_results))) {
    association_results <- association_results %>%
      mutate(
        conf.low = estimate - 1.96 * std.error, # 计算下置信区间
        conf.high = estimate + 1.96 * std.error # 计算上置信区间
      )
  }
  
  # 生成可视化
  association_plot <- ggplot(association_results, aes(estimate, -log10(p.value), 
                                                      color = fdr < 0.05, shape = metadata_var)) +
    geom_point(aes(size = abs(estimate)), alpha = 0.8) +
    geom_vline(xintercept = 0, linetype = 2) +
    scale_shape_manual(values = seq(16, 16 + length(unique(association_results$metadata_var)) - 1)) + # 动态生成形状
    labs(title = "Demographic-Feature Associations",
         x = "Effect Size (Estimate)",
         y = "-log10(p-value)",
         color = "FDR < 0.05",
         shape = "Demographic Variable")
  
  # 绘制森林图
  forest_plot <- ggplot(association_results, aes(x = estimate, y = metadata_var)) +
    geom_point(aes(color = fdr < 0.05), size = 3) + # 效应值点
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, color = fdr < 0.05), height = 0.2) + # 置信区间
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) + # 显著性颜色
    labs(
      title = "Forest Plot of Demographic-Feature Associations",
      x = "Effect Size (Estimate)",
      y = "Demographic Variable",
      color = "FDR < 0.05"
    ) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10)) # 调整 y 轴文本大小
  
  # 保存结果
  write_csv(association_results, paste0(prefix, "_demographic_associations.csv"))
  ggsave(paste0(prefix, "_demographic_volcano.png"), association_plot, width = 12, height = 8)
  ggsave(paste0(prefix, "_demographic_forest.png"), forest_plot, width = 12, height = 8)
  
  # Generate and save heatmap
  generate_heatmap(
    paste0(prefix, "_demographic_associations.csv"),
    paste0(prefix, "_demographic_associations_annotated.png"),
    paste0(toupper(prefix), " Demographic Associations")
  )
  
  return(association_results)
}

# Process all datasets
all_results <- map(data_files, ~process_dataset(.x$path, .x$prefix))

# 在主流程中调用新函数
demographic_results <- map(data_files, ~process_demographic_association(.x$path, .x$prefix))