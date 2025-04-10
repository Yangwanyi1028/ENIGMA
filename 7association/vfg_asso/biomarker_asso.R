setwd("/home/wanyiyang/ENIGMA_STUDY_DATA/7association/vfg_asso/")
library(dplyr)       # For data manipulation (mutate, left_join, etc.)
library(purrr)       # For functional programming (map_dfr)
library(broom)       # For tidying model outputs (tidy)
library(ggplot2)     # For data visualization (ggplot)
library(readr)       # For reading/writing CSV files (read_csv, write_csv)

# Clear environment, keeping only specific variables
library(tidyverse)
library(broom)
library(compositions) # For CLR transformation

# 1. Load and prepare data ---------------------------------------------------
# Load ARG abundance data (samples in columns)
arg_data <- read.csv("../../newsource_data20250319/VFG_173_479samples.csv",row.names = 1) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id")


# Load Biomarkers_Pathology data
biomarkers_pathology <- read.csv("../category_data/Biomarkers_Pathology.csv") %>%
  mutate(across(c("crp", "crp_yn", "bowel_resection", "malignancy"), as.factor)) # 将分类型变量转换为因子

arg_columns <- colnames(arg_data)[-1]

# 2. Merge with demographic data ---------------------------------------------
# Load demographic data
demo_data <- read.csv("../category_data/Demographics_Baseline.csv") %>%
  select(sample_id, age_final, gender, region)

# Merge full dataset
full_data_biomarkers <- biomarkers_pathology %>%
  left_join(demo_data, by = "sample_id") %>%
  left_join(arg_data, by = "sample_id")

# Define analysis variables for Biomarkers_Pathology
biomarkers_vars <- colnames(biomarkers_pathology)[-1] # 除去 sample_id 的所有列
continuous_vars <- c("hb", "wbc", "cd_hct_r", "esr") # 连续变量
categorical_vars <- setdiff(biomarkers_vars, continuous_vars) # 分类型变量

# Perform association analysis for all Biomarkers_Pathology variables
biomarkers_results <- map_dfr(arg_columns, function(arg) {
  map_dfr(biomarkers_vars, function(biomarker_var) {
    tryCatch({
      if (biomarker_var %in% continuous_vars) {
        # 连续变量：直接进行线性回归
        lm(reformulate(c(biomarker_var, "age_final", "gender", "region"), arg), 
           data = full_data_biomarkers) %>%
          tidy() %>%
          filter(term == biomarker_var) %>% # 筛选当前变量的结果
          mutate(
            metadata_var = biomarker_var,  # 添加变量名称
            biomarker = arg,                # 添加ARG基因名称
            fdr = p.adjust(p.value, method = "fdr") # 多重检验校正
          )
      } else {
        # 分类型变量：因子化后进行线性回归
        lm(reformulate(c(sprintf("factor(%s)", biomarker_var), "age_final", "gender", "region"), arg), 
           data = full_data_biomarkers) %>%
          tidy() %>%
          filter(grepl(biomarker_var, term)) %>% # 筛选当前变量的结果
          mutate(
            metadata_var = biomarker_var,  # 添加变量名称
            biomarker = arg,                # 添加ARG基因名称
            fdr = p.adjust(p.value, method = "fdr") # 多重检验校正
          )
      }
    }, error = function(e) NULL) # 捕获错误，避免中断
  })
})

# 检查结果是否包含所有生物标志物变量
print(unique(biomarkers_results$metadata_var)) # 调试：确保所有变量都被包含

# 生成可视化
biomarkers_plot <- ggplot(biomarkers_results, aes(estimate, -log10(p.value), 
                                                  color = fdr < 0.05, shape = metadata_var)) +
  geom_point(aes(size = abs(estimate)), alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_shape_manual(values = seq(16, 16 + length(unique(biomarkers_results$metadata_var)) - 1)) + # 动态生成形状
  labs(title = "ARG-Biomarkers Pathology Associations",
       x = "Effect Size (Estimate)",
       y = "-log10(p-value)",
       color = "FDR < 0.05",
       shape = "Biomarker Variable")

# 保存结果
write_csv(biomarkers_results, "biomarkers_pathology_associations.csv")
ggsave("biomarkers_pathology_volcano.png", biomarkers_plot, width = 12, height = 8)
