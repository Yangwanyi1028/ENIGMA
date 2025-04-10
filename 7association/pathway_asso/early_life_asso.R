# rm(list = ls())
# setwd('/home/wanyiyang/ENIGMA_STUDY_DATA/7association/')
# library(tidyverse)
# library(broom)
# library(compositions) # For CLR transformation

# # 1. Load and prepare data ---------------------------------------------------
# # Load ARG abundance data (samples in columns)
# arg_data <- read.csv("../source_data/new_ARG68_gene_abundance_20250224.csv", row.names = 1) %>%
#   t() %>% # Transpose to samples-as-rows format
#   as.data.frame() %>%
#   rownames_to_column("sample_id")
# arg_data <- arg_data[-1,]
# # Load early life metadata (assuming same sample order)
# early_meta <- read.csv("category_data/EarlyLife_Diet.csv") %>% # Save your metadata to this file
#   mutate(across(c(breast_fed, baby_food, family_food), as.factor))



# # 2. Merge with demographic data (from previous categories) -------------------
# # Load demographic data (assuming you saved it previously)
# demo_data <- read.csv("category_data/Demographics_Baseline.csv") %>%
#   select(sample_id, age_final, gender, region)

# # Merge all data
# full_data <- early_meta %>%
#   left_join(demo_data, by = "sample_id") %>%
#   left_join(arg_data, by = "sample_id") %>%
#   mutate(
#     gender = factor(gender),
#     region = factor(region)
#   )

# # 3. Compositional data transformation ---------------------------------------
# # CLR transform ARG abundances
# arg_columns <- colnames(arg_data)[-1]
# clr_transform <- full_data %>%
#   select(all_of(arg_columns)) %>%
#   clr() %>%
#   as_tibble()

# full_data <- full_data %>%
#   select(-all_of(arg_columns)) %>%
#   bind_cols(clr_transform)

# # 4. Association analysis function -------------------------------------------
# run_early_life_association <- function(metadata_var, arg_var) {
#   formula <- as.formula(
#     paste(arg_var, "~", metadata_var, "+ age_final + gender + region")
#   )
  
#   lm(formula, data = full_data) %>%
#     tidy() %>%
#     filter(term == metadata_var) %>%
#     mutate(
#       metadata_variable = metadata_var,
#       arg_gene = arg_var,
#       fdr = p.adjust(p.value, method = "fdr")
#     )
# }

# # 5. Batch analysis ----------------------------------------------------------
# # Define variables to analyze
# early_life_vars <- c(
#   "breast_fed", "baby_food", "family_food",
#   grep("^age_", colnames(full_data), value = TRUE)
# )

# # Run all combinations
# results <- map_dfr(arg_columns, function(arg) {
#   map_dfr(early_life_vars, ~ {
#     tryCatch(
#       run_early_life_association(.x, arg),
#       error = function(e) NULL
#     )
#   })
# })

# # 6. Result visualization ----------------------------------------------------
# volcano_plot <- results %>%
#   ggplot(aes(x = estimate, y = -log10(p.value), 
#              color = fdr < 0.05, size = abs(estimate))) +
#   geom_point(alpha = 0.7) +
#   geom_hline(yintercept = -log10(0.05), linetype = 2) +
#   scale_color_manual(values = c("gray70", "red")) +
#   facet_wrap(~metadata_variable, scales = "free_x") +
#   labs(title = "ARG-Early Life Associations with Batch Correction",
#        x = "Effect Size", y = "-log10(p-value)") +
#   theme_minimal()
setwd('/home/wanyiyang/ENIGMA_STUDY_DATA/7association/pathway_asso/')

# # 7. Save results ------------------------------------------------------------
# write_csv(results, "early_life_arg_associations.csv")
# ggsave("association_volcano_plots.png", volcano_plot, 
#        width = 16, height = 12, dpi = 300)
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
          pathway = arg,                # 添加ARG基因名称
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