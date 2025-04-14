rm(list = ls())
setwd("/home/wanyiyang/ENIGMA_STUDY_DATA/7association/")

ibd_data <- read.csv('../source_data/new_meta_data.csv') 


if("study_id" %in% names(ibd_data)) {
  names(ibd_data)[names(ibd_data) == "study_id"] <- "sample_id"
}
# Define variable groups based on your classification
variable_groups <- list(
  Demographics_Baseline = c(
    "sample_id", "age_final", "gender", "smoke",
    "region",  "bmi", "alcohol"
  ),
  
  Disease_Characteristics = c("sample_id" ,
    "cd_site1", "cd_site_l4", "cd_behaviour___1",
    "cd_behaviour___2", "cd_behaviour___3", "cd_perianal",
    "cdai"
  ),
  
  EarlyLife_Diet = c("sample_id" ,
    "breast_fed", "baby_food", "family_food",
    "age_4_12_months_3","age_1_5_years_3","age5_10_years_3",                
    "age_10_18_years_3","age_4_12_months_5",              
    "age_1_5_years_5","age_5_10_years_5","age_10_18_years_5",              
    "age_4_12_months_6","age_1_5_years_6","age_5_10_years_6",               
    "age_10_18_years_6","age_4_12_months_7","age_1_5_years_7",                
    "age_5_10_years_7","age_10_18_years_7","age_4_12_months_8",              
    "age_1_5_years_8","age_5_10_years_8","age_10_18_years_8",              
    "age_4_12_months_9","age_1_5_yers_9","age_5_10_years_9",               
    "age_10_18_years_9","age_4_12_months_10","age_1_5_years_10",               
    "age_5_10_years_10","age_10_18_years_10","age_4_12_months_11",             
    "age_1_5_years_11","age_5_10_years_11","age_10_18_years_11",             
    "age_4_12_months_12","age_1_5_years_12","age_5_10_years_12","age_10_18_years_12"      
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

patient_status <- read.csv('/home/wanyiyang/ENIGMA_STUDY_DATA/1maketable/metadata_filter.csv')
patient_status <- patient_status[,c('study_id','Group')]
names(patient_status) <- c('sample_id','group')
patient_status$group <- ifelse(patient_status$group == 'Control', 'Control','CD' )
patient_status_cd <- patient_status[patient_status$group == "CD",]$sample_id                                 
                                 
disease_dt = read.csv('category_data/Disease_Characteristics.csv')
# disease_dt <- na.omit(disease_dt)
names(disease_dt) <- c('sample_id','Disease location','Isolated upper disease','Inflammatory','Stricturing',
                       'Penetrating','Perianal disease modifier','CDAI')
disease_dt <- disease_dt[disease_dt$sample_id %in% patient_status_cd,]
write.csv(disease_dt,'category_data/Disease_Characteristics.csv')


medication <- read.csv('category_data/Medications.csv')
medication <- medication[,1:7]
names(medication) <- c('sample_id',
                       "Any treatment for CD","Any Topical Treatment","Any ORAL aminosalicilate",
                       "Any systemic steroids","Any ORAL immunosuppression","Any Biological therapy"
)
medication <- medication[medication$sample_id %in% patient_status_cd,]
write.csv(medication,'category_data/Medications.csv')


biomarker_dt <- read.csv('category_data/Biomarkers_Pathology.csv',stringsAsFactors = FALSE)
biomarker_dt$crp <- sapply(biomarker_dt$crp, function(x) {
  x <- trimws(x) # 去除空格
  if (grepl("^<", x)) {
    as.numeric(sub("<", "", x)) / 2 # 将 "<5" 转换为 2.5
  } else {
    as.numeric(x) # 转换为数值
  }
})

biomarker_dt <- biomarker_dt %>%
  mutate(crp = case_when(
    is.na(crp) ~ "Unknown",
    crp < 5 ~ "<5",
    crp >= 5 & crp < 10 ~ "5-10",
    crp >= 10 ~ ">=10"
  ))

summary(biomarker_dt)
names(biomarker_dt) <- c('sample_id',
                         "Haemoglobin",
                         "White Blood Cell Count",
                         "C-Reactive Protein",
                         "CRP Elevated over reference range",
                         "Haematocrit revised",
                         "Erythrocyte Sedimentation Rate",
                         "History of bowel resection/surgery",
                         "History of Cancer in last 5 years"
)
biomarker_dt <- biomarker_dt[biomarker_dt$sample_id %in% patient_status_cd,]
biomarker_dt$`C-Reactive Protein`
write.csv(biomarker_dt,'category_data/Biomarkers_Pathology.csv')




demo_dt <- read.csv('category_data/Demographics_Baseline.csv')
demo_dt$age_final <- round(demo_dt$age_final,0)
names(demo_dt) <- c('sample_id','Age at enrolment','Gender','Smoking status','Region','BMI','Alcohol Consumption')
write.csv(demo_dt,'category_data/Demographics_Baseline.csv')

early_dt <- read.csv('category_data/EarlyLife_Diet.csv')
names(early_dt) <- c('sample_id','Breastfeeding Status','Baby Foof Origin', 'Child Food Origin', 
                     'Age 4 12 months Home Grown Fruit Veg','Age 1 5 years Home Grown Fruit Veg','Age 5 10 years Home Grown Fruit Veg', 
                     'Age 10 18 years Home Grown Fruit Veg',
                     "Age 4 12 months Milk Beverages", "Age 1 5 years Milk Beverages",
                     "Age 5 10 years Milk Beverages", "Age 10 18 years Milk Beverages", "Age 4 12 months Processed Meat Seafood",
                     "Age 1 5 years Processed Meat Seafood","Age 5 10 years Processed Meat Seafood",
                     "Age 10 18 years Processed Meat Seafood","Age 4 12 months Processed Carbs",
                     "Age 1 5 years Processed Carbs", "Age 5 10 years Processed Carbs",
                     "Age 10 18 years Processed Carbs", "Age 4 12 months Processed Fruit",
                     "Age 1 5 years Processed Fruit", "Age 5 10 years Processed Fruit",
                     "Age 10 18 years Processed Fruit", "Age 4 12 months Processed Veg",
                     "Age 1 5 years Processed Veg", "Age 5 10 years Processed Veg",
                     "Age 10 18 years Processed Veg", "Age 4 12 months Fast food",
                     "Age 1 5 years Fast food", "Age 5 10 years Fast food",
                     "Age 10 18 years Fast food", "Age 4 12 months Soft Drink","Age 1 5 years Soft Drink",
                     "Age 5 10 years Soft Drink", "Age 10 18 years Soft Drink", "Age 4 12 months Packaged Snacks",
                     "Age 1 5 years Packaged Snacks", "Age 5 10 years Packaged Snacks", "Age 10 18 years Packaged Snacks"
)

write.csv(early_dt,'category_data/EarlyLife_Diet.csv')


recent_intake <- read.csv('category_data/Recent_Additives.csv')
names(recent_intake) <- c('sample_id',"P80","CMC","CRN","AlSiO","SO32","TiO2","Asp",
                          "Suc","Sac","TotalAd","TotalAs",
                          "TotalEM","P80 Annual Intake",# (mg/kg body weight/day)
                          "CMC Annual Intake","CRN Annual Intake",
                          "AlSiO Annual Intake","SO32 Annual Intake",
                          "TiO2 Annual Intake","Asp Annual Intake",
                          "Suc Annual Intake","Sac Annual Intake",
                          "TotalAd Annual Intake","TotalAs Annual Intake",
                          "TotalEM Annual Intake"
)
write.csv(recent_intake,'category_data/Recent_Additives.csv')



threeday_dt <- read.csv('category_data/ThreeDay_Nutrients.csv')
names(threeday_dt) <- c('sample_id',"Protein intake" ,
                       "Total Fat intake","Saturated Fat intake",
                       "Monounsaturated Fat intake","Polyunsaturated Fat intake",
                       "Omega 3 Fatty Acid intake","Omega 6 Fatty Acid intake",
                       "Cholesterol Intake","Alcohol Intake",
                       "Fibre Intake","Vitamin B1",
                       "Vitamin B2","Vitamin B3","Vitamin B6",
                       "Vitamin B12","Vitamin C","Vitamin E/Tocopherol Acetate",
                       "Beta Carotene","Folic Acid",
                       "Folate" ,"Iron" ,"Zinc","Selenium","Magnesium","Caffeine","Vitamin D"
)
write.csv(threeday_dt,'category_data/ThreeDay_Nutrients.csv')


# 将数据框放入列表中
data_frames <- list(disease_dt, medication, biomarker_dt)
# 按 sample_id 合并
cb_cd_only <- reduce(data_frames, full_join, by = "sample_id")
print(dim(cb_cd_only))
print(names(cb_cd_only))
data_frames <- list(demo_dt, recent_intake, threeday_dt, early_dt)

cb_all <- reduce(data_frames, full_join, by = "sample_id")
write.csv(cb_cd_only,'category_data/cb_cd_only.csv')
write.csv(cb_all,'category_data/cb_all.csv')
