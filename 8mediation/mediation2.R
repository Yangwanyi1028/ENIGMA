# Clear environment and set working directory
rm(list = ls())
setwd('/home/wanyiyang/ENIGMA_STUDY_DATA/8mediation/')

# Load required packages
packages <- c("vegan", "lme4", "randomForest", "mediation", "ggplot2",
              "dplyr", "parallel", "foreach", "doParallel", "data.table","magrittr")
invisible(lapply(packages, library, character.only = TRUE))

# Function to load and clean data
load_clean_data <- function(diet_path, cdai_path, microbiome_path) {
  # Load datasets
  diet <- fread(diet_path)
  cdai <- fread(cdai_path)[, .(sample_id, CDAI)]
  microbiome <- as.data.frame(t(read.csv(microbiome_path,row.names = 1)))
  
  # Clean column names
  setnames(diet, make.names(names(diet)))
  microbiome$sample_id <- rownames(microbiome)
  
  # Clean and merge datasets
  diet_clean <- na.omit(diet)
  microbiome_clean <- microbiome[rowSums(microbiome[, -ncol(microbiome)]) > 0, ]
  
  # Merge all datasets using dplyr
  merged_data <- diet_clean %>%
    inner_join(cdai, by = "sample_id") %>%
    inner_join(microbiome_clean, by = "sample_id") %>%
    na.omit()
  
  return(merged_data)
}


# 修改后的并行中介分析函数
perform_parallel_mediation <- function(merged_data, diet_factors, species_columns, 
                                       n_cores = detectCores() - 1) {
  
  registerDoParallel(cores = n_cores)
  
  results <- foreach(treat = diet_factors, .combine = rbind) %:%
    foreach(species = species_columns, .combine = rbind) %dopar% {
      
      data <- data.frame(
        X = merged_data[[treat]],
        Y = merged_data$CDAI,
        M = merged_data[[species]]
      ) %>% na.omit()
      
      # 常规中介模型 -----------------------------------------------------------
      model.m <- lm(M ~ X, data)
      model.y <- lm(Y ~ X + M, data)
      med_result <- mediate(model.m, model.y, treat = "X", mediator = "M", 
                            boot = TRUE, sims = 1000)  # 修改点1：增加到1000次
      
      # 逆向中介模型 -----------------------------------------------------------
      model.m_inv <- lm(Y ~ X, data)
      model.y_inv <- lm(M ~ X + Y, data)
      med_result_inv <- mediate(model.m_inv, model.y_inv, treat = "X", 
                                mediator = "Y", boot = TRUE, sims = 1000)  # 修改点2
      
      # 提取完整结果 -----------------------------------------------------------
      data.frame(
        diet = treat,
        species = species,
        # 常规模型指标
        ACME_est = med_result$d.avg,            # 平均因果中介效应
        ACME_CI_low = med_result$d.avg.ci[1],   # 95% CI下限
        ACME_CI_high = med_result$d.avg.ci[2],  # 95% CI上限
        ADE_est = med_result$z.avg,             # 平均直接效应
        ADE_CI_low = med_result$z.avg.ci[1],
        ADE_CI_high = med_result$z.avg.ci[2],
        Total_est = med_result$tau.coef,        # 总效应
        Total_CI_low = med_result$tau.ci[1],
        Total_CI_high = med_result$tau.ci[2],
        Pval_mediate = med_result$d.avg.p,
        Pval_direct = med_result$z.avg.p,
        
        # 逆向模型指标
        ACME_inv_est = med_result_inv$d.avg,
        ACME_inv_CI_low = med_result_inv$d.avg.ci[1],
        ACME_inv_CI_high = med_result_inv$d.avg.ci[2],
        ADE_inv_est = med_result_inv$z.avg,
        ADE_inv_CI_low = med_result_inv$z.avg.ci[1],
        ADE_inv_CI_high = med_result_inv$z.avg.ci[2],
        Total_inv_est = med_result_inv$tau.coef,
        Total_inv_CI_low = med_result_inv$tau.ci[1],
        Total_inv_CI_high = med_result_inv$tau.ci[2],
        Pval_mediate_inverse = med_result_inv$d.avg.p,
        Pval_direct_inverse = med_result_inv$z.avg.p
      )
    }
  
  stopImplicitCluster()
  return(results)
}





# Add correlation analysis function
perform_correlation_analysis <- function(merged_data, diet_factors, species_columns) {
  # Initialize results dataframe
  correlation_results <- data.frame()
  
  for(diet in diet_factors) {
    for(species in species_columns) {
      # Diet-Species correlation
      cor_diet_sp <- cor.test(merged_data[[diet]], 
                              merged_data[[species]], 
                              method = "spearman")
      
      # Species-CDAI correlation
      cor_sp_cdai <- cor.test(merged_data[[species]], 
                              merged_data$CDAI, 
                              method = "spearman")
      
      # Store results
      result <- data.frame(
        diet = diet,
        species = species,
        Cor_diet_sp = cor_diet_sp$estimate,
        Pval_diet_sp = cor_diet_sp$p.value,
        Cor_sp_cdai = cor_sp_cdai$estimate,
        Pval_sp_cdai = cor_sp_cdai$p.value
      )
      
      correlation_results <- rbind(correlation_results, result)
    }
  }
  
  # FDR correction
  correlation_results$Qval_diet_sp <- p.adjust(correlation_results$Pval_diet_sp, 
                                               method = "BH")
  correlation_results$Qval_sp_cdai <- p.adjust(correlation_results$Pval_sp_cdai, 
                                               method = "BH")
  
  return(correlation_results)
}

# Modify main analysis function
main_analysis <- function(merged_data) {
  # Extract relevant columns
  diet_factors <- names(merged_data)[2:39]
  species_columns <- names(merged_data)[42:142]
  
  # Step 1: Correlation analysis
  correlation_results <- perform_correlation_analysis(merged_data, 
                                                      diet_factors, 
                                                      species_columns)
  
  # Save all correlation results
  fwrite(correlation_results, 
         "mediation_results/all_correlation_results.txt",
         quote = FALSE, sep = "\t", row.names = FALSE)
  
  # Save significant correlation results
  significant_pairs <- correlation_results %>%
    filter(Qval_diet_sp < 0.05 & Qval_sp_cdai < 0.05)
  
  fwrite(significant_pairs, 
         "mediation_results/significant_correlation_pairs.txt",
         quote = FALSE, sep = "\t", row.names = FALSE)
  
  # Step 2: Perform mediation analysis only on significant pairs
  mediation_results <- perform_parallel_mediation(merged_data, 
                                                  unique(significant_pairs$diet),
                                                  unique(significant_pairs$species))
  print(mediation_results)
  # Merge correlation and mediation results
  final_results <- merge(significant_pairs, mediation_results, 
                         by = c("diet", "species"))
  
  # Filter for significant mediation effects
  # 修改后的结果筛选逻辑
  final_results <- final_results %>%
    mutate(
      Qval_mediate = p.adjust(Pval_mediate, method = "BH"),
      Qval_mediate_inverse = p.adjust(Pval_mediate_inverse, method = "BH")
    ) %T>% 
    # 使用tee运算符写入中间结果（%T>%保持管道继续）
    { write.csv(., "mediation_results/mediation_analysis_results.csv", row.names = FALSE, quote = FALSE) } %>% 
    filter(
      (Qval_mediate < 0.05 & Pval_mediate_inverse > 0.05) |
        (Pval_mediate > 0.05 & Qval_mediate_inverse < 0.05)
    ) %>%
    # 添加效应方向判断
    mutate(
      Mediation_Direction = case_when(
        Qval_mediate < 0.05 & ACME_est > 0 ~ "Positive Mediation",
        Qval_mediate < 0.05 & ACME_est < 0 ~ "Negative Mediation",
        Qval_mediate_inverse < 0.05 & ACME_inv_est > 0 ~ "Inverse Positive",
        Qval_mediate_inverse < 0.05 & ACME_inv_est < 0 ~ "Inverse Negative"
      )
    )
  
  # Save final mediation results
  fwrite(final_results, 
         "mediation_results/diet_species_cdai_mediation_results.txt",
         quote = FALSE, sep = "\t", row.names = FALSE)
  final_report <- final_results %>%
    # 先过滤有效结果
    filter(!is.na(Mediation_Direction)) %>%
    # 格式化报告语句
    mutate(
      # 提取方向关键词
      Direction_Type = case_when(
        str_detect(Mediation_Direction, "Positive") ~ "正向",
        str_detect(Mediation_Direction, "Negative") ~ "负向"
      ),
      
      # 统一置信区间格式
      ACME_CI = sprintf("[%.2f, %.2f]", 
                        round(ACME_CI_low, 2),
                        round(ACME_CI_high, 2)),
      
      # 构建完整报告语句
      Report_Statement = glue::glue(
        "*[{species}]通过[{Direction_Type}]中介作用连接[{diet}]与CDAI",
        "（ACME={round(ACME_est, 2)}，95%CI{ACME_CI}）*"
      )
    ) %>%
    # 按效应量绝对值排序
    arrange(desc(abs(ACME_est)))
  
  # 输出报告文件
  writeLines(final_report$Report_Statement, "mediation_results/mediation_findings.md")
  
  
  return(final_results)
}

# Load and prepare data
merged_data <- load_clean_data(
  '../7association/category_data/EarlyLife_Diet.csv',
  '../7association/category_data/Disease_Characteristics.csv',
  # '../newsource_data20250319/species_101_479samples1.csv'
  '../source_data/merged_species_abundance_mp4_20240614_605samples.csv'
)

# Run analysis
results <- main_analysis(merged_data)

