setwd('/home/wanyiyang/ENIGMA_STUDY_DATA/7association/cd_only/')
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

# Define data files and their corresponding names
data_files <- list(
  ARG = list(
    path = "../../newsource_data20250319/ARG_68_479samples.csv",
    prefix = "CD_only_ARG"
  ),
  PATH = list(
    path = "../../newsource_data20250319/path_68_479samples.csv",
    prefix = "CD_only_PATH"
  ),
  SPECIES = list(
    path = "../../newsource_data20250319/species_101_479samples.csv",
    prefix = "CD_only_SPECIES"
  ),
  VFG = list(
    path = "../../newsource_data20250319/VFG_173_479samples.csv",
    prefix = "CD_only_VFG"
  )
)

# Load clinical data
clinical_data <- read.csv("../category_data/cb_cd_only.csv", row.names = 1) %>%
  mutate(across(c("Disease.location", "Isolated.upper.disease", "Inflammatory", 
                  "Stricturing", "Penetrating", "Perianal.disease.modifier", "C.Reactive.Protein",
                  "Any.treatment.for.CD", "Any.Topical.Treatment", "Any.ORAL.aminosalicilate",
                  "Any.systemic.steroids", "Any.ORAL.immunosuppression", "Any.Biological.therapy",
                  "CRP.Elevated.over.reference.range", "History.of.bowel.resection.surgery",
                  "History.of.Cancer.in.last.5.years"
                  ), as.factor))

# Load demographic data
demographic_data <- read.csv("../category_data/Demographics_Baseline.csv", row.names = 1) %>%
  select(sample_id, Age.at.enrolment, Gender, Smoking.status, Region, BMI, Alcohol.Consumption)

# Define analysis variables for clinical data
clinical_vars <- colnames(clinical_data)[-1] # 除去 sample_id 的所有列
continuous_vars <- c("CDAI", "Haemoglobin", "White.Blood.Cell.Count", "Haematocrit.revised",
                     "Erythrocyte.Sedimentation.Rate")
categorical_vars <- setdiff(clinical_vars, continuous_vars)

# Define column annotations
clinical_annotations <- list(
  disease_status = c("Disease.location", "Isolated.upper.disease", "Inflammatory",
                     "Stricturing", "Penetrating", "Perianal.disease.modifier", "CDAI"),
  medication_status = c("Any.treatment.for.CD", "Any.Topical.Treatment", 
                        "Any.ORAL.aminosalicilate", "Any.systemic.steroids",
                        "Any.ORAL.immunosuppression", "Any.Biological.therapy"),
  biomarkers = c("Haemoglobin", "White.Blood.Cell.Count", "C.Reactive.Protein",
                 "CRP.Elevated.over.reference.range", "Haematocrit.revised",
                 "Erythrocyte.Sedimentation.Rate", "History.of.bowel.resection.surgery",
                 "History.of.Cancer.in.last.5.years")
)

# Define row annotations
feature_status_annotations <- read.csv('../../newsource_data20250319/Combined_Annotations_with_CD_Status.csv')
feature_status_annotations <- feature_status_annotations[, c(1, 3)]
names(feature_status_annotations) <- c('features', 'CD_status')

# Custom function: Save pheatmap as PNG
save_pheatmap_png <- function(x, filename, nrow, ncol, res = 300) {
  width <- max(3000, ncol * 50) + 1000  # Minimum width 3000
  height <- max(3000, nrow * 50)        # Minimum height 3000
  png(filename, width = width, height = height, res = res)
  grid::grid.draw(x)
  dev.off()
}

# Custom function: Plot heatmap
plotAssociationsDag3HeatmapV2 <- function(inData, phenosToPlot, statToPlot, featuresToPlot, 
                                          nrFeaturesToPlot = 50, nrColors = 13, 
                                          sortPhenos = FALSE, retData = FALSE, fdrCutoff = 0.05,
                                          addColAnnotation = TRUE) {
  # Filter data
  inData <- inData %>%
    filter(phenotype %in% phenosToPlot & feature_name %in% featuresToPlot)
  
  # Create matrix
  heatmap_data <- inData %>%
    select(phenotype, feature_name, !!sym(statToPlot)) %>%
    tidyr::pivot_wider(
      names_from = phenotype, 
      values_from = !!sym(statToPlot),
      values_fn = mean,
      values_fill = NA
    ) %>%
    column_to_rownames("feature_name")
  
  # Z-score normalization
  heatmap_data <- scale(heatmap_data, center = TRUE, scale = TRUE)
  
  # Sort phenotypes
  if (sortPhenos) {
    heatmap_data <- heatmap_data[, order(colnames(heatmap_data))]
  }
  
  # Define colors
  paletteLength <- nrColors
  myColor <- colorRampPalette(c("Orange", "white", "Darkblue"))(paletteLength)
  myBreaks <- c(seq(min(heatmap_data, na.rm = TRUE), 0, length.out = ceiling(paletteLength / 2) + 1), 
                seq(max(heatmap_data, na.rm = TRUE) / paletteLength, max(heatmap_data, na.rm = TRUE), length.out = floor(paletteLength / 2)))

  # Add significance annotations
  significance_data <- inData %>%
    select(phenotype, feature_name, fdr) %>%
    tidyr::pivot_wider(
      names_from = phenotype, 
      values_from = fdr,
      values_fn = mean,
      values_fill = NA
    ) %>%
    column_to_rownames("feature_name")
  
  significance_data <- ifelse(significance_data < fdrCutoff, "+", "")
  significance_data[heatmap_data < 0 & significance_data == "+"] <- "-"
  significance_data[is.na(heatmap_data)] <- ""

  # Column annotations
  annotation_col <- data.frame(
    Category = rep("Other", length(phenosToPlot)),
    row.names = phenosToPlot
  )
  if (addColAnnotation) {
    annotation_col$Category[phenosToPlot %in% clinical_annotations$disease_status] <- "Disease Status"
    annotation_col$Category[phenosToPlot %in% clinical_annotations$medication_status] <- "Medication Status"
    annotation_col$Category[phenosToPlot %in% clinical_annotations$biomarkers] <- "Biomarkers"
  }

  # Row annotations
  annotation_row <- feature_status_annotations %>%
    filter(features %in% rownames(heatmap_data)) %>%
    column_to_rownames("features")
  
  # Define annotation colors
  annotation_colors <- list(
    CD_status = c(
      "CD_enriched" = "#ec9191",  # 红色
      "CD_depleted" = "#65b6f8"   # 蓝色
    ),
    Category = c(
      "Disease Status" = "#df7a7ba4",
      "Medication Status" = "#5e7f9c",
      "Biomarkers" = "#abc281",
      "Other" = "grey"
    )
  )
  
  # Plot heatmap
  phm <- pheatmap(heatmap_data, 
                  color = myColor,
                  breaks = myBreaks,
                  cluster_rows = FALSE,
                  cluster_cols = FALSE,
                  annotation_col = if (addColAnnotation) annotation_col else NULL,
                  annotation_row = annotation_row,
                  annotation_colors = annotation_colors,
                  fontsize = 10,
                  border_color = "#EEEEEE",
                  na_col = "white",
                  display_numbers = significance_data,
                  fontsize_number = 8,
                  angle_col = 315,
                  annotation_names_col = TRUE,
                  annotation_legend = TRUE)
  
  if (retData) {
    return(list(heatmap_data, phm))
  } else {
    return(phm)
  }
}

# Generate heatmap
generate_heatmap <- function(data_file, output_file, title, stat_to_plot = "estimate", addColAnnotation = TRUE) {
  data <- read_csv(data_file, show_col_types = FALSE)
  
  if (!"metadata_var" %in% colnames(data) || !"feature_name" %in% colnames(data)) {
    stop("数据文件缺少必要的列：metadata_var 或 feature_name")
  }
  
  data$phenotype <- data$metadata_var
  data$effect.size <- data$estimate
  data$FDR <- data$fdr
  
  phenos_to_plot <- unique(data$phenotype)
  features_to_plot <- unique(data$feature_name)
  
  phm <- plotAssociationsDag3HeatmapV2(
    inData = data,
    phenosToPlot = phenos_to_plot,
    statToPlot = stat_to_plot,
    featuresToPlot = features_to_plot,
    nrFeaturesToPlot = length(features_to_plot),
    nrColors = 13,
    sortPhenos = FALSE,
    retData = TRUE,
    addColAnnotation = addColAnnotation
  )
  
  nrow <- length(features_to_plot)
  ncol <- length(phenos_to_plot)
  
  save_pheatmap_png(phm[[2]], filename = output_file, nrow = nrow, ncol = ncol, res = 300)
}

# Process datasets
process_dataset <- function(file_path, prefix) {
  abundance_data <- read.csv(file_path, row.names = 1) %>%
    t() %>% 
    as.data.frame() %>%
    rownames_to_column("sample_id")
  
  feature_columns <- colnames(abundance_data)[-1]
  
  merged_data <- clinical_data %>%
    left_join(demographic_data, by = "sample_id") %>%
    left_join(abundance_data, by = "sample_id")
  
  association_results <- map_dfr(feature_columns, function(feature) {
    map_dfr(clinical_vars, function(clinical_var) {
      tryCatch({
        if (clinical_var %in% continuous_vars) {
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
  
  if (!("conf.low" %in% colnames(association_results)) || !("conf.high" %in% colnames(association_results))) {
    association_results <- association_results %>%
      mutate(
        conf.low = estimate - 1.96 * std.error,
        conf.high = estimate + 1.96 * std.error
      )
  }
  
  write_csv(association_results, paste0(prefix, "_clinical_associations.csv"))
  generate_heatmap(
    paste0(prefix, "_clinical_associations.csv"),
    paste0(prefix, "_clinical_associations_annotated.png"),
    paste0(toupper(prefix), " Clinical Associations")
  )
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
  
  # 保存结果文件
  write_csv(association_results, paste0(prefix, "_demographic_associations.csv"))
  
  # 生成热图
  generate_heatmap(
    paste0(prefix, "_demographic_associations.csv"),
    paste0(prefix, "_demographic_associations_annotated.png"),
    paste0(toupper(prefix), " Demographic Associations"),
    addColAnnotation = FALSE
  )
  
  return(association_results)
}

# 在主流程中调用新函数
demographic_results <- map(data_files, ~process_demographic_association(.x$path, .x$prefix))
# Process all datasets
all_results <- map(data_files, ~process_dataset(.x$path, .x$prefix))