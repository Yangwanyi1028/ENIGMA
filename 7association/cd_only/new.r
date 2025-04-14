setwd('/home/wanyiyang/ENIGMA_STUDY_DATA/7association/cd_only')

# Load required libraries
library(dplyr)       # For data manipulation
library(purrr)       # For functional programming
library(ggplot2)     # For data visualization
library(readr)       # For reading/writing CSV files
library(tidyverse)
library(broom)
library(compositions) # For CLR transformation
library(survminer)    # For ggforest function
library(pheatmap)
library(readr)
library(tidyr)
library(grid)
library(tibble)

# ==================== DATA LOADING FUNCTIONS ====================

load_clinical_data <- function() {
    clinical_data <- read.csv("../category_data/cb_cd_only.csv", row.names = 1) %>%
    mutate(across(c("Disease.location", "Isolated.upper.disease", "Inflammatory", 
                  "Stricturing", "Penetrating", "Perianal.disease.modifier", "C.Reactive.Protein",
                  "Any.treatment.for.CD", "Any.Topical.Treatment", "Any.ORAL.aminosalicilate",
                  "Any.systemic.steroids", "Any.ORAL.immunosuppression", "Any.Biological.therapy",
                  "CRP.Elevated.over.reference.range", "History.of.bowel.resection.surgery",
                  "History.of.Cancer.in.last.5.years"
                  ), as.factor))
}

load_demographic_data <- function() {
  read.csv("../category_data/Demographics_Baseline.csv", row.names = 1) %>%
    select(sample_id, Age.at.enrolment, Gender, Smoking.status, Region, BMI, Alcohol.Consumption)
}

load_feature_annotations <- function() {
  read.csv('../../newsource_data20250319/Combined_Annotations_with_CD_Status.csv')[,c(1,3)] %>%
    setNames(c('features','CD_status'))
}

load_abundance_data <- function(file_path) {
  read.csv(file_path, row.names = 1) %>%
    t() %>% 
    as.data.frame() %>%
    rownames_to_column("sample_id")
}

# ==================== VARIABLE DEFINITIONS ====================

# Define clinical variable categories
define_clinical_variables <- function(clinical_data) {
  clinical_vars <- colnames(clinical_data)[-1] # Exclude sample_id
  
  continuous_vars <- c("Age.at.enrolment", "BMI", "P80", "CMC","CRN","AlSiO","SO32",
                       "TiO2","Asp","Suc","Sac","TotalAd","TotalAs","TotalEM",
                       "P80.Annual.Intake", "CMC.Annual.Intake", 'CRN.Annual.Intake',
                       'AlSiO.Annual.Intake', "SO32.Annual.Intake", "TiO2.Annual.Intake", 
                       "Asp.Annual.Intake", "Suc.Annual.Intake", "Sac.Annual.Intake",
                       "TotalAd.Annual.Intake", "TotalAs.Annual.Intake", "TotalEM.Annual.Intake", 
                       "Protein.intake", "Total.Fat.intake", "Saturated.Fat.intake", 
                       "Monounsaturated.Fat.intake", "Polyunsaturated.Fat.intake",
                       "Omega.3.Fatty.Acid.intake", "Omega.6.Fatty.Acid.intake", 
                       "Cholesterol.Intake", "Alcohol.Intake", "Fibre.Intake",
                       "Vitamin.B1", "Vitamin.B2", "Vitamin.B3", "Vitamin.B6", "Vitamin.B12", "Vitamin.C",  
                       "Vitamin.E.Tocopherol.Acetate", "Beta.Carotene", "Folic.Acid", "Folate", "Iron",       
                       "Zinc", "Selenium", "Magnesium", "Caffeine", "Vitamin.D")
  
  categorical_vars <- setdiff(clinical_vars, continuous_vars)
  
  list(
    clinical_vars = clinical_vars,
    continuous_vars = continuous_vars,
    categorical_vars = categorical_vars
  )
}

# Define clinical annotations
define_clinical_annotations <- function(clinical_data) {
  list(
    `Recent Addictives` = names(clinical_data)[8:31],
    `Threeday Nutrients` = names(clinical_data)[32:57],
    `Early Life Diet` = names(clinical_data)[58:96]
  )
}

# ==================== HEATMAP VISUALIZATION FUNCTIONS ====================

# Function to save pheatmap as PNG with square cells
save_pheatmap_png <- function(pheatmap_object, filename, nrow, ncol, cell_size = 50, res = 300, extra_width = 0) {
  # Calculate dimensions to ensure square cells
  width <- ncol * cell_size + extra_width + 1000
  height <- nrow * cell_size
  
  # Save the plot
  png(filename, width = width, height = height, res = res)
  grid::grid.draw(pheatmap_object)
  dev.off()
}

# Function to prepare data for heatmap
prepare_heatmap_data <- function(inData, phenosToPlot, statToPlot, featuresToPlot) {
  # Filter data
  filtered_data <- inData %>%
    filter(phenotype %in% phenosToPlot & feature_name %in% featuresToPlot)
  
  # Create matrix for heatmap
  heatmap_data <- filtered_data %>%
    select(phenotype, feature_name, !!sym(statToPlot)) %>%
    tidyr::pivot_wider(
      names_from = phenotype, 
      values_from = !!sym(statToPlot),
      values_fn = mean,  # Average for duplicate values
      values_fill = NA   # Fill missing values with NA
    ) %>%
    column_to_rownames("feature_name")
  
  # Z-score standardization
  scaled_data <- scale(heatmap_data, center = TRUE, scale = TRUE)
  
  return(scaled_data)
}

# Function to prepare significance data
prepare_significance_data <- function(inData, phenosToPlot, featuresToPlot, heatmap_data, fdrCutoff = 0.05) {
  # Create significance matrix
  significance_data <- inData %>%
    filter(phenotype %in% phenosToPlot & feature_name %in% featuresToPlot) %>%
    select(phenotype, feature_name, fdr) %>%
    tidyr::pivot_wider(
      names_from = phenotype, 
      values_from = fdr,
      values_fn = mean,
      values_fill = NA
    ) %>%
    column_to_rownames("feature_name")
  
  # Add +/- symbols based on significance and direction
  significance_symbols <- ifelse(significance_data < fdrCutoff, "+", "")
  significance_symbols[heatmap_data < 0 & significance_symbols == "+"] <- "-"
  significance_symbols[is.na(heatmap_data)] <- ""
  
  return(significance_symbols)
}

# Function to prepare annotations
prepare_annotations <- function(phenosToPlot, featuresToPlot, clinical_annotations, feature_status_annotations) {
  # Column annotations
  annotation_col <- data.frame(
    Category = rep("Other", length(phenosToPlot)),
    row.names = phenosToPlot
  )
  
  # Fill categories
  annotation_col$Category[phenosToPlot %in% clinical_annotations$`Recent Addictives`] <- "Recent Addictives"
  annotation_col$Category[phenosToPlot %in% clinical_annotations$`Threeday Nutrients`] <- "Threeday Nutrients"
  annotation_col$Category[phenosToPlot %in% clinical_annotations$`Early Life Diet`] <- "Early Life Diet"
  
  # Row annotations
  annotation_row <- feature_status_annotations %>%
    filter(features %in% featuresToPlot) %>%
    column_to_rownames("features")
  
  # Calculate gaps for visual separation
  gaps_col <- which(diff(as.numeric(factor(annotation_col$Category))) != 0)
  
  # Define annotation colors
  annotation_colors <- list(
    CD_status = c(
      "CD_enriched" = "#ab6355",  # Red
      "CD_depleted" = "#90a5a7"   # Blue
    ),
    Category = c(
      "Recent Addictives"  = "#df7a7ba4",  # Red
      "Threeday Nutrients" = "#5e7f9c",    # Blue
      "Early Life Diet" = "#abc281",       # Green
      "Other" = "grey"
    )
  )
  
  return(list(
    annotation_col = annotation_col,
    annotation_row = annotation_row,
    gaps_col = gaps_col,
    annotation_colors = annotation_colors
  ))
}

# Main heatmap plotting function
plot_associations_heatmap <- function(inData, phenosToPlot, statToPlot, featuresToPlot, 
                                     fdrCutoff = 0.05, addColAnnotation = TRUE,
                                     clinical_annotations, feature_status_annotations) {
  
  # Prepare heatmap data
  heatmap_data <- prepare_heatmap_data(inData, phenosToPlot, statToPlot, featuresToPlot)
  
  # Prepare significance data
  significance_symbols <- prepare_significance_data(inData, phenosToPlot, featuresToPlot, heatmap_data, fdrCutoff)
  
  # Prepare annotations
  annotations <- prepare_annotations(phenosToPlot, featuresToPlot, clinical_annotations, feature_status_annotations)
  
  # Set color palette
  paletteLength <- 13
  myColor <- colorRampPalette(c("Orange", "white", "Darkblue"))(paletteLength)
  myBreaks <- c(
    seq(min(heatmap_data, na.rm = TRUE), 0, length.out = ceiling(paletteLength / 2) + 1),
    seq(max(heatmap_data, na.rm = TRUE) / paletteLength, max(heatmap_data, na.rm = TRUE),
        length.out = floor(paletteLength / 2))
  )
  
  # Create pheatmap
  phm <- pheatmap(
    heatmap_data, 
    color = myColor,
    breaks = myBreaks,
    cluster_rows = FALSE,  # No row clustering
    cluster_cols = FALSE,  # No column clustering
    annotation_col = if (addColAnnotation) annotations$annotation_col else NULL,
    annotation_row = annotations$annotation_row,
    annotation_colors = annotations$annotation_colors,
    fontsize = 10,
    border_color = NA,
    display_numbers = significance_symbols,
    fontsize_number = 8,
    angle_col = 315,
    annotation_names_col = TRUE,
    annotation_legend = TRUE,
    gaps_col = annotations$gaps_col
  )
  
  return(list(data = heatmap_data, plot = phm))
}

# Function to generate and save heatmap
generate_heatmap <- function(data_file, output_file, title, stat_to_plot = "estimate", 
                             addColAnnotation = TRUE, clinical_annotations, feature_status_annotations) {
  # Read data
  data <- read_csv(data_file, show_col_types = FALSE)
  
  # Validate data
  if (!"metadata_var" %in% colnames(data) || !"feature_name" %in% colnames(data)) {
    stop("Data file is missing required columns: metadata_var or feature_name")
  }
  
  # Prepare data
  data$phenotype <- data$metadata_var
  data$effect.size <- data$estimate
  data$FDR <- data$fdr
  
  # Extract variables and features
  phenos_to_plot <- unique(data$phenotype)
  features_to_plot <- unique(data$feature_name)
  
  # Calculate the maximum length of feature names
  max_feature_name_length <- max(nchar(features_to_plot))
  
  # Create heatmap
  heatmap_result <- plot_associations_heatmap(
    inData = data,
    phenosToPlot = phenos_to_plot,
    statToPlot = stat_to_plot,
    featuresToPlot = features_to_plot,
    addColAnnotation = addColAnnotation,
    clinical_annotations = clinical_annotations,
    feature_status_annotations = feature_status_annotations
  )
  
  # Dynamically calculate cell size and adjust heatmap dimensions
  cell_size <- 50  # Default cell size in pixels
  nrow <- length(features_to_plot)
  ncol <- length(phenos_to_plot)
  
  # Adjust extra width specifically for ARG demographic associations
  if (grepl("arg_demographic_associations", output_file)) {
    extra_width <- 1500  # Custom width for ARG heatmap
  } else {
    extra_width <- max(0, (max_feature_name_length - 10) * 40)  # Default width adjustment
  }
  
  # Save heatmap with square cells
  save_pheatmap_png(
    heatmap_result$plot, 
    filename = output_file, 
    nrow = nrow, 
    ncol = ncol, 
    cell_size = cell_size, 
    res = 300,
    extra_width = extra_width
  )
  
  return(heatmap_result)
}

# ==================== STATISTICAL ANALYSIS FUNCTIONS ====================

# Function to perform association analysis
perform_association_analysis <- function(feature, variables, data, continuous_vars) {
  map_dfr(variables, function(var) {
    tryCatch({
      if (var %in% continuous_vars) {
        # Continuous variable: direct linear regression
        lm(reformulate(c(var), feature), data = data) %>%
          tidy() %>%
          filter(term == var) %>%
          mutate(
            metadata_var = var,
            feature_name = feature,
            fdr = p.adjust(p.value, method = "fdr")
          )
      } else {
        # Categorical variable: factor in linear regression
        lm(reformulate(c(sprintf("factor(%s)", var)), feature), data = data) %>%
          tidy() %>%
          filter(grepl(var, term)) %>%
          mutate(
            metadata_var = var,
            feature_name = feature,
            fdr = p.adjust(p.value, method = "fdr")
          )
      }
    }, error = function(e) NULL)
  })
}

# Add confidence intervals to results
add_confidence_intervals <- function(results) {
  if (!("conf.low" %in% colnames(results)) || !("conf.high" %in% colnames(results))) {
    results <- results %>%
      mutate(
        conf.low = estimate - 1.96 * std.error,
        conf.high = estimate + 1.96 * std.error
      )
  }
  return(results)
}

# ==================== MAIN PROCESSING FUNCTIONS ====================

# Process clinical associations
process_dataset <- function(file_path, prefix, clinical_data, demographic_data, 
                           variable_definitions, clinical_annotations, feature_status_annotations) {
  # Load abundance data
  abundance_data <- load_abundance_data(file_path)
  feature_columns <- colnames(abundance_data)[-1]
  
  # Merge data
  merged_data <- clinical_data %>%
    left_join(demographic_data, by = "sample_id") %>%
    left_join(abundance_data, by = "sample_id")
  
  # Perform association analysis
  association_results <- map_dfr(feature_columns, function(feature) {
    perform_association_analysis(
      feature, 
      variable_definitions$clinical_vars, 
      merged_data, 
      variable_definitions$continuous_vars
    )
  })
  
  # Add confidence intervals
  association_results <- add_confidence_intervals(association_results)
  
  # Save results
  write_csv(association_results, paste0(prefix, "_clinical_associations.csv"))
  
  # Generate heatmap
  generate_heatmap(
    paste0(prefix, "_clinical_associations.csv"),
    paste0(prefix, "_clinical_associations_annotated.png"),
    paste0(toupper(prefix), " Clinical Associations"),
    clinical_annotations = clinical_annotations,
    feature_status_annotations = feature_status_annotations
  )
  
  return(association_results)
}

# Process demographic associations
process_demographic_association <- function(file_path, prefix, demographic_data, 
                                          variable_definitions, clinical_annotations, 
                                          feature_status_annotations) {
  # Load abundance data
  abundance_data <- load_abundance_data(file_path)
  feature_columns <- colnames(abundance_data)[-1]
  
  # Merge with demographic data
  merged_data <- demographic_data %>%
    left_join(abundance_data, by = "sample_id")
  
  # Perform association analysis
  association_results <- map_dfr(feature_columns, function(feature) {
    perform_association_analysis(
      feature, 
      colnames(demographic_data)[-1], 
      merged_data, 
      variable_definitions$continuous_vars
    )
  })
  
  # Add confidence intervals
  association_results <- add_confidence_intervals(association_results)
  
  # Save results
  write_csv(association_results, paste0(prefix, "_demographic_associations.csv"))
  
  # Generate heatmap
  generate_heatmap(
    paste0(prefix, "_demographic_associations.csv"),
    paste0(prefix, "_demographic_associations_annotated.png"),
    paste0(toupper(prefix), " Demographic Associations"),
    addColAnnotation = FALSE,
    clinical_annotations = clinical_annotations,
    feature_status_annotations = feature_status_annotations
  )
  
  return(association_results)
}

# ==================== MAIN EXECUTION ====================

main <- function() {
  # Load data
  clinical_data <- load_clinical_data()
  demographic_data <- load_demographic_data()
  feature_status_annotations <- load_feature_annotations()
  
  # Define variables and annotations
  variable_definitions <- define_clinical_variables(clinical_data)
  clinical_annotations <- define_clinical_annotations(clinical_data)
  
  # Define file paths
  data_files <- list(
    arg = list(
      path = "../../newsource_data20250319/ARG_68_479samples1.csv",
      prefix = "arg"
    ),
    path = list(
      path = "../../newsource_data20250319/path_68_479samples1.csv",
      prefix = "path"
    ),
    species = list(
      path = "../../newsource_data20250319/species_101_479samples1.csv",
      prefix = "species"
    ),
    vfg = list(
      path = "../../newsource_data20250319/VFG_173_479samples1.csv",
      prefix = "vfg"
    )
  )
  
  # Process all datasets
  all_results <- map(data_files, ~process_dataset(
    .x$path, 
    .x$prefix, 
    clinical_data, 
    demographic_data, 
    variable_definitions,
    clinical_annotations,
    feature_status_annotations
  ))
  
  # Process demographic associations
  demographic_results <- map(data_files, ~process_demographic_association(
    .x$path, 
    .x$prefix, 
    demographic_data, 
    variable_definitions,
    clinical_annotations,
    feature_status_annotations
  ))
  
  # Generate heatmap for ARG demographic associations
  generate_heatmap(
    data_file = "arg_demographic_associations.csv",
    output_file = "arg_demographic_associations_annotated.png",
    title = "ARG Demographic Associations",
    addColAnnotation = FALSE,
    clinical_annotations = clinical_annotations,
    feature_status_annotations = feature_status_annotations
  )
  
  return(list(clinical = all_results, demographic = demographic_results))
}

# Execute main function
results <- main()
