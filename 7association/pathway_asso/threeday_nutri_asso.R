setwd('/home/wanyiyang/ENIGMA_STUDY_DATA/7association/pathway_asso/')

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
arg_data <- read.csv("../../newsource_data20250319/path_68_479samples.csv", row.names = 1) %>%
  t() %>% # Transpose to samples-as-rows format
  as.data.frame() %>%
  rownames_to_column("sample_id")


arg_columns <- colnames(arg_data)[-1]



# 2. Merge with demographic data (from previous categories) -------------------
# Load demographic data (assuming you saved it previously)
demo_data <- read.csv("../category_data/Demographics_Baseline.csv") %>%
  select(sample_id, age_final, gender, region)

# Load ThreeDay_Nutrients data
threeday_nutrients <- read.csv("../category_data/ThreeDay_Nutrients.csv") %>%
  mutate(across(where(is.character), as.factor)) # 将字符型变量转换为因子

# Merge full dataset
full_data_nutrients <- threeday_nutrients %>%
  left_join(demo_data, by = "sample_id") %>%
  left_join(arg_data, by = "sample_id")

# Define analysis variables for ThreeDay_Nutrients
nutrients_vars <- colnames(threeday_nutrients)[-1] # 除去 sample_id 的所有列

# Perform association analysis for all ThreeDay_Nutrients variables
nutrients_results <- map_dfr(arg_columns, function(arg) {
  map_dfr(nutrients_vars, function(nutrient_var) {
    tryCatch({
      # 检查变量类型：连续变量 vs 分类型变量
      if (is.numeric(full_data_nutrients[[nutrient_var]])) {
        # 连续变量：直接进行线性回归
        lm(reformulate(c(nutrient_var, "age_final", "gender", "region"), arg), 
           data = full_data_nutrients) %>%
          tidy() %>%
          filter(term == nutrient_var) %>% # 筛选当前变量的结果
          mutate(
            metadata_var = nutrient_var,  # 添加变量名称
            pathway = arg,               # 添加ARG基因名称
            fdr = p.adjust(p.value, method = "fdr") # 多重检验校正
          )
      } else {
        # 分类型变量：因子化后进行线性回归
        lm(reformulate(c(sprintf("factor(%s)", nutrient_var), "age_final", "gender", "region"), arg), 
           data = full_data_nutrients) %>%
          tidy() %>%
          filter(grepl(nutrient_var, term)) %>% # 筛选当前变量的结果
          mutate(
            metadata_var = nutrient_var,  # 添加变量名称
            pathway = arg,               # 添加ARG基因名称
            fdr = p.adjust(p.value, method = "fdr") # 多重检验校正
          )
      }
    }, error = function(e) NULL) # 捕获错误，避免中断
  })
})

# 检查结果是否包含所有营养变量
print(unique(nutrients_results$metadata_var)) # 调试：确保所有变量都被包含

# 生成可视化
nutrients_plot <- ggplot(nutrients_results, aes(estimate, -log10(p.value), 
                                                color = fdr < 0.05, shape = metadata_var)) +
  geom_point(aes(size = abs(estimate)), alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_shape_manual(values = seq(16, 16 + length(unique(nutrients_results$metadata_var)) - 1)) + # 动态生成形状
  labs(title = "ARG-ThreeDay Nutrients Associations",
       x = "Effect Size (Estimate)",
       y = "-log10(p-value)",
       color = "FDR < 0.05",
       shape = "Nutrient Variable")

# 保存结果
write_csv(nutrients_results, "threeday_nutrients_associations.csv")
ggsave("threeday_nutrients_volcano.png", nutrients_plot, width = 12, height = 8)