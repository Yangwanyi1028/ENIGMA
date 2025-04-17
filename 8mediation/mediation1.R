# Clear environment and set working directory
rm(list = ls())
setwd('/home/wanyiyang/ENIGMA_STUDY_DATA/8mediation/')

# Load required packages
packages <- c("vegan", "lme4", "randomForest", "mediation", "ggplot2",
             "dplyr", "parallel", "foreach", "doParallel", "data.table")
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

# Load and prepare data
merged_data <- load_clean_data(
  '../7association/category_data/EarlyLife_Diet.csv',
  '../7association/category_data/Disease_Characteristics.csv',
  '../newsource_data20250319/species_101_479samples1.csv'
)

# Parallel Mediation Analysis
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
      
      # Regular mediation
      model.m <- lm(M ~ X, data)
      model.y <- lm(Y ~ X + M, data)
      med_result <- mediate(model.m, model.y, treat = "X", mediator = "M", 
                           boot = TRUE, sims = 500)
      
      # Inverse mediation
      model.m_inv <- lm(Y ~ X, data)
      model.y_inv <- lm(M ~ X + Y, data)
      med_result_inv <- mediate(model.m_inv, model.y_inv, treat = "X", 
                              mediator = "Y", boot = TRUE, sims = 500)
      
      data.frame(
        diet = treat,
        species = species,
        Pval_mediate = summary(med_result)$d.avg.p,
        Pval_direct = summary(med_result)$z.avg.p,
        Pval_mediate_inverse = summary(med_result_inv)$d.avg.p,
        Pval_direct_inverse = summary(med_result_inv)$z.avg.p
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
  
  # Merge correlation and mediation results
  final_results <- merge(significant_pairs, mediation_results, 
                         by = c("diet", "species"))
  
  # Filter for significant mediation effects
  final_results <- final_results %>%
    mutate(
      Qval_mediate = p.adjust(Pval_mediate, method = "BH"),
      Qval_mediate_inverse = p.adjust(Pval_mediate_inverse, method = "BH")
    ) %>%
    filter(
      (Qval_mediate < 0.05 & Pval_mediate_inverse > 0.05) |
        (Pval_mediate > 0.05 & Qval_mediate_inverse < 0.05)
    )
  
  # Save final mediation results
  fwrite(final_results, 
         "mediation_results/diet_species_cdai_mediation_results.txt",
         quote = FALSE, sep = "\t", row.names = FALSE)
  
  return(final_results)
}


# Run analysis
results <- main_analysis(merged_data)
