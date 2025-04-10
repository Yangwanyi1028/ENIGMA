
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
# Load ARG abundance data (samples in columns)
arg_data <- read.csv("../../newsource_data20250319/species_101_479samples.csv",row.names = 1) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id")
# Load early life metadata (assuming same sample order)

arg_columns <- colnames(arg_data)[-1]


# 2. Merge with demographic data (from previous categories) -------------------
# Load demographic data (assuming you saved it previously)
demo_data <- read.csv("../category_data/Demographics_Baseline.csv") %>%
  select(sample_id, age_final, gender, region)

# Load disease characteristics data
disease <- read.csv("../category_data/Disease_Characteristics.csv") %>%
  mutate(across(starts_with("cd_"), factor))

# Merge full dataset
full_data <- disease %>%
  left_join(demo_data, by = "sample_id") %>%
  left_join(arg_data, by = "sample_id")

# Define analysis variables
disease_vars <- c("cd_site1", "cd_site_l4", "cd_perianal", "cdai", "cd_act",
                  "cd_behaviour___1", "cd_behaviour___2", "cd_behaviour___3")

# Perform association analysis for all disease variables
disease_results <- map_dfr(arg_columns, function(arg) {
  map_dfr(disease_vars, function(disease_var) {
    tryCatch({
      # 检查变量类型：连续变量 vs 分类型变量
      if (is.numeric(full_data[[disease_var]])) {
        # 连续变量：直接进行线性回归
        lm(reformulate(c(disease_var, "age_final", "gender", "region"), arg), 
           data = full_data) %>%
          tidy() %>%
          filter(term == disease_var) %>% # 筛选当前变量的结果
          mutate(
            metadata_var = disease_var,  # 添加变量名称
            species = arg,              # 添加ARG基因名称
            fdr = p.adjust(p.value, method = "fdr") # 多重检验校正
          )
      } else {
        # 分类型变量：因子化后进行线性回归
        lm(reformulate(c(sprintf("factor(%s)", disease_var), "age_final", "gender", "region"), arg), 
           data = full_data) %>%
          tidy() %>%
          filter(grepl(disease_var, term)) %>% # 筛选当前变量的结果
          mutate(
            metadata_var = disease_var,  # 添加变量名称
            species = arg,              # 添加ARG基因名称
            fdr = p.adjust(p.value, method = "fdr") # 多重检验校正
          )
      }
    }, error = function(e) NULL) # 捕获错误，避免中断
  })
})

# 检查结果是否包含所有疾病变量
print(unique(disease_results$metadata_var)) # 调试：确保所有变量都被包含

# 生成可视化
disease_plot <- ggplot(disease_results, aes(estimate, -log10(p.value), 
                                            color = fdr < 0.05, shape = metadata_var)) +
  geom_point(aes(size = abs(estimate)), alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_shape_manual(values = c(16, 17, 15, 18, 19, 20, 21, 22)) + # 提供足够的形状
  labs(title = "ARG-Disease Characteristics Associations",
       x = "Effect Size (Estimate)",
       y = "-log10(p-value)",
       color = "FDR < 0.05",
       shape = "Disease Variable")

# 保存结果
write_csv(disease_results, "disease_associations.csv")
ggsave("disease_volcano.png", disease_plot, width = 12, height = 8)
