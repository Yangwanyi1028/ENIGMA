rm(list = ls())
setwd('/home/wanyiyang/ENIGMA_STUDY_DATA/7association/')

ibd_data <- read.csv('../source_data/new_meta_data.csv')


if("study_id" %in% names(ibd_data)) {
  names(ibd_data)[names(ibd_data) == "study_id"] <- "sample_id"
}
# Define variable groups based on your classification
variable_groups <- list(
  Demographics_Baseline = c(
    "sample_id", "record_id", "age_final", "gender", "smoke",
    "region",  "bmi", "alcohol"
  ),
  
  Disease_Characteristics = c("sample_id" ,
    "cd_site1", "cd_site_l4", "cd_behaviour___1",
    "cd_behaviour___2", "cd_behaviour___3", "cd_perianal",
    "cdai", "cd_act"
  ),
  
  EarlyLife_Diet = c("sample_id" ,
    "breast_fed", "baby_food", "family_food",
    grep("^age_4_12_months_3$|^age_1_5_years_3$|^age5_10_years_3$|^age_10_18_years_3$|
         ^age_4_12_months_5$|^age_1_5_years_5$|^age_5_10_years_5$|^age_10_18_years_5$|
         ^age_4_12_months_6$|^age_1_5_years_6$|^age_5_10_years_6$|^age_10_18_years_6$|
         ^age_4_12_months_7$|^age_1_5_years_7$|^age_5_10_years_7$|^age_10_18_years_7$|
         ^age_4_12_months_8$|^age_1_5_years_8$|^age_5_10_years_8$|^age_10_18_years_8$|
         ^age_4_12_months_9$|^age_1_5_yers_9$|^age_5_10_years_9$|^age_10_18_years_9$|
         ^age_4_12_months_10$|^age_1_5_years_10$|^age_5_10_years_10$|^age_10_18_years_10$|
         ^age_4_12_months_11$|^age_1_5_years_11$|^age_5_10_years_11$|^age_10_18_years_11$|
         ^age_4_12_months_12$|^age_1_5_years_12$|^age_5_10_years_12$|^age_10_18_years_12$", 
         names(ibd_data), value = TRUE)
  ),
  
  Recent_Additives = c("sample_id" ,
    "p80", "cmc", "crn", "alsio", "so32", "tio2", "asp", "suc", "sac",
    "totalad", "totalas", "totalem", "p80mgkgday", "cmcmgkgday",
    "crnmgkgday", "alsiomgkgday", "so32mgkgday", "tio2mgkgday",
    "aspmgkgday", "sucmgkgday", "sacmgkgday", "totaladmgkgday",
    "totalasmgkgday", "totalemmgkgday"
  ),
  
  ThreeDay_Nutrients = c("sample_id" ,
    "protein_g", "total_fat_g", "sfa_g", "mufa_g", "pufa_g",
    "omega3fa_g", "omega6fa_g", "cholesterol_mg", "alcohol_g",
    "fiber_g", "b1_mg", "b2_mg", "b3_mg", "b6_mg", "b12_mcg",
    "vitaminc_mg", "vitaminetocoa_mg", "bcarotene_mcg",
    "folicacid_mcg_dfe", "folate_mcg", "iron_mg", "zinc_mg",
    "selenium_mcg", "magnesium_mg", "caffeine_mg", "vitd_iu"
  ),
  
  Medications = c("sample_id" ,
    "cd_tx_noTx", "cd_tx_topical", "cd_tx_oralASA",
    "cd_tx_SysSteroids", "cd_tx_immunosup",
    grep("^cd_therapy_|^drug_", names(ibd_data), value = TRUE)
  ),
  
  Biomarkers_Pathology = c("sample_id" ,
    "hb", "wbc", "crp", "crp_yn", "cd_hct_r", "esr",
    "bowel_resection", "malignancy"
  )
)

# Create output directory if needed
if(!dir.exists("category_data")) dir.create("category_data")

# Save each group as CSV
lapply(names(variable_groups), function(group_name) {
  subset_data <- ibd_data[, variable_groups[[group_name]]]
  filename <- paste0("category_data/", group_name, ".csv")
  write.csv(subset_data, file = filename, row.names = FALSE, na = "")
})

message("Files saved in: ", normalizePath("category_data"))
