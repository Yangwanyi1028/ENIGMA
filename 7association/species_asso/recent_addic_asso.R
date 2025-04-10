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
arg_data <- read.csv("../../newsource_data20250319/species_101_479samples.csv",row.names = 1) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id")


arg_columns <- colnames(arg_data)[-1]

# 2. Merge with demographic data (from previous categories) -------------------
# Load demographic data (assuming you saved it previously)
demo_data <- read.csv("../category_data/Demographics_Baseline.csv") %>%
  select(sample_id, age_final, gender, region)

# Load recent additives data
recent_additives <- read.csv("../category_data/Recent_Additives.csv") %>%
  mutate(across(where(is.character), as.factor)) # 将字符型变量转换为因子

# Merge full dataset
full_data_additives <- recent_additives %>%
  left_join(demo_data, by = "sample_id") %>%
  left_join(arg_data, by = "sample_id")

# Define analysis variables for recent additives
additives_vars <- colnames(recent_additives)[-1] # 除去 sample_id 的所有列

# Perform association analysis for all recent additives variables
additives_results <- map_dfr(arg_columns, function(arg) {
  map_dfr(additives_vars, function(additive_var) {
    tryCatch({
      # 检查变量类型：连续变量 vs 分类型变量
      if (is.numeric(full_data_additives[[additive_var]])) {
        # 连续变量：直接进行线性回归
        lm(reformulate(c(additive_var, "age_final", "gender", "region"), arg), 
           data = full_data_additives) %>%
          tidy() %>%
          filter(term == additive_var) %>% # 筛选当前变量的结果
          mutate(
            metadata_var = additive_var,  # 添加变量名称
            species = arg,               # 添加ARG基因名称
            fdr = p.adjust(p.value, method = "fdr") # 多重检验校正
          )
      } else {
        # 分类型变量：因子化后进行线性回归
        lm(reformulate(c(sprintf("factor(%s)", additive_var), "age_final", "gender", "region"), arg), 
           data = full_data_additives) %>%
          tidy() %>%
          filter(grepl(additive_var, term)) %>% # 筛选当前变量的结果
          mutate(
            metadata_var = additive_var,  # 添加变量名称
            species = arg,               # 添加ARG基因名称
            fdr = p.adjust(p.value, method = "fdr") # 多重检验校正
          )
      }
    }, error = function(e) NULL) # 捕获错误，避免中断
  })
})

# 检查结果是否包含所有添加剂变量
print(unique(additives_results$metadata_var)) # 调试：确保所有变量都被包含

# 生成可视化
additives_plot <- ggplot(additives_results, aes(estimate, -log10(p.value), 
                                                color = fdr < 0.05, shape = metadata_var)) +
  geom_point(aes(size = abs(estimate)), alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_shape_manual(values = seq(16, 16 + length(unique(additives_results$metadata_var)) - 1)) + # 动态生成形状
  labs(title = "ARG-Recent Additives Associations",
       x = "Effect Size (Estimate)",
       y = "-log10(p-value)",
       color = "FDR < 0.05",
       shape = "Additive Variable")

# 保存结果
write_csv(additives_results, "recent_additives_associations.csv")
ggsave("recent_additives_volcano.png", additives_plot, width = 12, height = 8)