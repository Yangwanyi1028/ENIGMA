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
# Load early life diet data
earlylife_diet <- read.csv("../category_data/EarlyLife_Diet.csv") %>%
  mutate(across(everything(), as.factor)) # 将所有列转换为因子

arg_columns <- colnames(arg_data)[-1]

# 2. Merge with demographic data ---------------------------------------------
# Load demographic data
demo_data <- read.csv("../category_data/Demographics_Baseline.csv") %>%
  select(sample_id, age_final, gender, region)

# Merge full dataset
full_data_earlylife <- earlylife_diet %>%
  left_join(demo_data, by = "sample_id") %>%
  left_join(arg_data, by = "sample_id")

# Define analysis variables for early life diet
earlylife_vars <- colnames(earlylife_diet)[-1] # 除去 sample_id 的所有列

# Perform association analysis for all early life diet variables
earlylife_results <- map_dfr(arg_columns, function(arg) {
  map_dfr(earlylife_vars, function(earlylife_var) {
    tryCatch({
      # 分类型变量：因子化后进行线性回归
      lm(reformulate(c(sprintf("factor(%s)", earlylife_var), "age_final", "gender", "region"), arg), 
         data = full_data_earlylife) %>%
        tidy() %>%
        filter(grepl(earlylife_var, term)) %>% # 筛选当前变量的结果
        mutate(
          metadata_var = earlylife_var,  # 添加变量名称
          earlylife_diet = arg,                # 添加ARG基因名称
          fdr = p.adjust(p.value, method = "fdr") # 多重检验校正
        )
    }, error = function(e) NULL) # 捕获错误，避免中断
  })
})

# 检查结果是否包含所有早期饮食变量
print(unique(earlylife_results$metadata_var)) # 调试：确保所有变量都被包含

# 生成可视化
earlylife_plot <- ggplot(earlylife_results, aes(estimate, -log10(p.value), 
                                                color = fdr < 0.05, shape = metadata_var)) +
  geom_point(aes(size = abs(estimate)), alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_shape_manual(values = seq(16, 16 + length(unique(earlylife_results$metadata_var)) - 1)) + # 动态生成形状
  labs(title = "ARG-Early Life Diet Associations",
       x = "Effect Size (Estimate)",
       y = "-log10(p-value)",
       color = "FDR < 0.05",
       shape = "Early Life Diet Variable")

# 保存结果
write_csv(earlylife_results, "earlylife_diet_associations.csv")
ggsave("earlylife_diet_volcano.png", earlylife_plot, width = 12, height = 8)