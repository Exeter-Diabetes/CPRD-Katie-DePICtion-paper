
# Define incident cohort of insulin-treated T1 and T2 diagnosed aged 18+ years at diagnosis:

## All IDs in download diagnosed within registration and between 01/03/2003-01/03/2013 - max 10 years follow up
## With linkage data

# Pull in patient features from previous script
## Diagnosis date based on earliest diabetes code: exclude those diagnosed between -30 and +90 days of registration and those diagnosed <18 years
## Diabetes type based on (latest) type codes, and only keep if current diagnosis is Type 1 or 2
## Current treatment, and exclude those not on insulin
## Earliest insulin/OHA and time to insulin
## Current BMI, HDL, total cholesterol, triglyceride and HbA1c - as far back as diagnosis
## Antibodies, C-peptide, family history, features at diagnosis, autoimmune conditions
## Hospitalisation for hypoglycaemia and DKA outcomes

# Run lipid model


############################################################################################

# Setup
library(tidyverse)
library(aurum)
library(EHRBiomarkr)
library(lubridate)
rm(list=ls())

cprd = CPRDData$new(cprdEnv = "diabetes-jun2024", cprdConf = "~/.aurum.yaml")

analysis = cprd$analysis("dpctn_paper_inci")


############################################################################################

# Set index dates

first_index_date <- as.Date("2013-03-01")
second_index_date <- as.Date("2023-03-01")


# Point to patient variables coded up in previous script

analysis = cprd$analysis("dpctn_paper")
diagnosis_dates <- diagnosis_dates %>% analysis$cached("diagnosis_dates")
earliest_type_code <- earliest_type_code %>% analysis$cached("earliest_type_code")
diabetes_type <- diabetes_type %>% analysis$cached("diabetes_type")
diabetes_meds <- diabetes_meds %>% analysis$cached("diabetes_meds")
biomarkers <- biomarkers %>% analysis$cached("biomarkers")
gad <- gad %>% analysis$cached("gad")
ia2 <- ia2 %>% analysis$cached("ia2")
c_pep <- c_pep %>% analysis$cached("c_pep")
fh_positive_negative <- fh_positive_negative %>% analysis$cached("fh_positive_negative")
at_diag_features <- at_diag_features %>% analysis$cached("at_diag_features")
autoimmune <- autoimmune %>% analysis$cached("autoimmune")
dk_hypo_contacts <- dk_hypo_contacts %>% analysis$cached("dk_hypo_contacts")


############################################################################################

# Find IDs of those registered at some point between 01/03/2003-01/03/2013 and pull in static patient features

# Include general patient variables as per all_diabetes_cohort:
## gender
## DOB (derived as per: https://github.com/Exeter-Diabetes/CPRD-Codelists/blob/main/readme.md#general-notes-on-implementation, cached for all IDs in 'diagnosis_date_dob' table from all_patid_dob script)
## pracid
## prac_region
## ethnicity 5-category and 16-category (derived as per: https://github.com/Exeter-Diabetes/CPRD-Codelists#ethnicity)
## regstartdate	
## gp_record_end (earliest of last collection date from practice, deregistration and 30/06/2024 (latest date in records))	
## death_date	(from ONS)
## whether patient has HES data available
## IMD decile

analysis = cprd$analysis("diabetes_cohort")
dob <- dob %>% analysis$cached("dob")

analysis = cprd$analysis("all_patid")
ethnicity <- ethnicity %>% analysis$cached("ethnicity")


analysis = cprd$analysis("dpctn_paper_inci")

cohort <- cprd$tables$patient %>%
  select(patid, gender, pracid, regstartdate, regenddate) %>%
  left_join((dob %>% select(patid, dob)), by="patid") %>%
  left_join((cprd$tables$practice %>% select(pracid, prac_region=region, lcd)), by="pracid") %>%
  left_join((ethnicity %>% select(patid, ethnicity_5cat, ethnicity_16cat)), by="patid") %>%
  left_join((cprd$tables$onsDeath %>% select(patid, death_date=reg_date_of_death)), by="patid") %>%
  left_join((cprd$tables$patidsWithLinkage %>% mutate(with_hes=1L)), by="patid") %>%
  left_join((cprd$tables$patientImd %>% select(-pracid)), by="patid") %>%
  
  mutate(gp_record_end=pmin(if_else(is.na(lcd), as.Date("2024-06-30"), lcd),
                            if_else(is.na(regenddate), as.Date("2024-06-30"), regenddate),
                            if_else(is.na(death_date), as.Date("2024-06-30"), death_date),
                            as.Date("2024-06-30"), na.rm=TRUE)) %>%
  
  select(patid, gender, dob, pracid, prac_region, ethnicity_5cat, ethnicity_16cat, regstartdate, gp_record_end, death_date, with_hes, imd_decile) %>%
  
  filter(with_hes==1 & regstartdate<=second_index_date & !(!is.na(death_date) & death_date<first_index_date) & !(!is.na(gp_record_end) & gp_record_end<first_index_date)) %>%
  
  analysis$cached("cohort_interim_1", unique_indexes="patid")

cohort %>% count()
# 1,807,596

# Remove patients from practices which may have duplicated data as per CPRD's guidance

analysis = cprd$analysis("diabetes_cohort")

practice_exclusion_ids <- cprd$tables$patient %>% 
  filter(pracid == "20024" | pracid == "20036" |pracid == "20091" |pracid == "20171" | pracid == "20178" |pracid == "20202" | pracid == "20254" | pracid == "20389" |pracid == "20430" |pracid == "20452" |
           pracid == "20469" | pracid == "20487" | pracid == "20552" | pracid == "20554" | pracid == "20640" | pracid == "20717" | pracid == "20734" | pracid == "20737" | pracid == "20740" | pracid == "20790" |
           pracid == "20803" | pracid == "20822" | pracid == "20868" | pracid == "20912" | pracid == "20996" | pracid == "21001" | pracid == "21015" | pracid == "21078" | pracid == "21112" | pracid == "21118" |
           pracid == "21172" | pracid == "21173" | pracid == "21277" | pracid == "21281" | pracid == "21331" | pracid == "21334" | pracid == "21390" | pracid == "21430" | pracid == "21444" | pracid == "21451" |
           pracid == "21529" | pracid == "21553" | pracid == "21558" | pracid == "21585") %>%
  distinct(patid) %>%
  analysis$cached("practice_exclusion_ids")


analysis = cprd$analysis("dpctn_paper_inci")

cohort <- cohort %>%
  anti_join(practice_exclusion_ids, by="patid") %>%
  analysis$cached("cohort_interim_2", unique_indexes="patid")

cohort %>% count()
# 1,784,963



# Just keep men and women

cohort <- cohort %>%
  filter(gender==1 | gender==2) %>%
  analysis$cached("cohort_interim_3", unique_indexes="patid")

cohort %>% count()
# 1,784,941


############################################################################################

# Add diagnosis dates and exclude if not within registration OR not between 01/03/2013 and 01/03/2023

cohort <- cohort %>%
  inner_join(diagnosis_dates, by="patid") %>%
  mutate(dm_diag_age=round((datediff(diagnosis_date, dob))/365.25, 1),
         follow_up_time=datediff(second_index_date, diagnosis_date)/365.25) %>%
  filter(diagnosis_date>regstartdate & diagnosis_date>=first_index_date & diagnosis_date<=second_index_date & diagnosis_date<=gp_record_end & (is.na(death_date) | diagnosis_date<=death_date)) %>%
  analysis$cached("cohort_interim_4", unique_indexes="patid")

cohort %>% count()
#719857


############################################################################################

# Add earliest type code and keep if type 1 or type 2

cohort <- cohort %>%
  inner_join(earliest_type_code, by="patid") %>%
  left_join(diabetes_type, by="patid") %>%
  rename(latest_diabetes_type=diabetes_type,
         diabetes_type=earliest_type_code) %>%
  mutate(t1_code_ever=ifelse(type1_code_count>0, 1L, 0L),
         t2_code_ever=ifelse(type2_code_count>0, 1L, 0L)) %>%
  analysis$cached("cohort_interim_5", unique_indexes="patid")

cohort %>% count()
#459735

719857-459735

#lose 260122 people with unspecified codes only


# Look at breakdown of earliest type code

clipr::write_clip(cohort %>%
  filter(!(diabetes_type=="type 1" | diabetes_type=="type 2")) %>%
  group_by(diabetes_type) %>%
  summarise(count=n()) %>%
  collect())

#n=18,468 gestational
#n=2,657 secondary
#n=1,410 other


cohort <- cohort %>%
  filter(diabetes_type=="type 1" | diabetes_type=="type 2") %>%
  analysis$cached("cohort_interim_6", unique_indexes="patid")

cohort %>% count()
#437200

459735-437200

cohort %>% group_by(diabetes_type) %>% count()
# 1 type 2         421215
# 2 type 1          15985


############################################################################################

# Remove if diagnosed <18 or >50

cohort <- cohort %>%
  filter(dm_diag_age>=18 & dm_diag_age<=50) %>%
  analysis$cached("cohort_interim_7", unique_indexes="patid")

cohort %>% count()
#108217

cohort %>% group_by(diabetes_type) %>% count()
# 1 type 2         100971
# 2 type 1          7246

437200-108217
421215-100971
15985-7246


############################################################################################

# Add in treatment

cohort <- cohort %>%
  inner_join((diabetes_meds %>% select(-c(regstartdate, diagnosis_date))), by="patid") %>%
  analysis$cached("cohort_interim_8", unique_indexes="patid")

cohort %>% count()
#108217

cohort %>% group_by(diabetes_type) %>% count()
# 1 type 2         100971
# 2 type 1          7246


############################################################################################

# Data quality checking: remove if close to registration

cohort <- cohort %>%
  filter((datediff(diagnosis_date, regstartdate)< -30 | datediff(diagnosis_date, regstartdate)>90)) %>%
  analysis$cached("cohort_interim_9", unique_indexes="patid")

cohort %>% count()
#96699

cohort %>% group_by(diabetes_type) %>% count()
# 1 type 2         91377
# 2 type 1          5322

108217-96699
7246-5322
100971-91377


############################################################################################

# Add in remaining features

cohort <- cohort %>%
  left_join(biomarkers, by="patid") %>%
  left_join(gad, by="patid") %>%
  left_join(ia2, by="patid") %>%
  left_join(c_pep, by="patid") %>%
  left_join(fh_positive_negative, by="patid") %>%
  left_join(at_diag_features, by="patid") %>%
  left_join(autoimmune, by="patid") %>%
  left_join(dk_hypo_contacts, by="patid") %>%
  analysis$cached("cohort_interim_10", unique_indexes="patid")


## Check missing variables for model

cohort %>% mutate(age_at_bmi=datediff(diag_bmi_date, dob)/365.25) %>% filter(is.na(diag_bmi) | age_at_bmi<18) %>% count()
#27384

cohort %>% filter(is.na(diag_totalcholesterol)) %>% count()
#20643

cohort %>% filter(is.na(diag_hdl)) %>% count()
#23524

cohort %>% filter(is.na(diag_triglyceride)) %>% count()
#36379


cohort %>% mutate(age_at_bmi=datediff(diag_bmi_date, dob)/365.25) %>% filter(is.na(diag_bmi) | age_at_bmi<18 | is.na(diag_totalcholesterol) | is.na(diag_hdl) | is.na(diag_triglyceride)) %>% count()
#51055

51055/108217 #47.2%


cohort <- cohort %>%
  mutate(age_at_bmi=datediff(diag_bmi_date, dob)/365.25) %>%
  filter(!is.na(diag_bmi) & age_at_bmi>=18 & !is.na(diag_totalcholesterol) & !is.na(diag_hdl) & !is.na(diag_triglyceride)) %>%
  analysis$cached("cohort_interim_11", unique_indexes="patid")

cohort %>% count()
#45644

cohort %>% group_by(diabetes_type) %>% count()

# 1 type 2          44823
# 2 type 1          821


############################################################################################

# Code up lipid model

cohort <- cohort %>%
  
  mutate(femalesex=ifelse(gender==2, 1, ifelse(gender==1, 0, NA)),
         
         lipid_pred_score=9.0034272-(0.1915482*diag_bmi)-(0.1686227*dm_diag_age)+(0.3026012*femalesex)-(0.2269216*diag_totalcholesterol)+(1.540850*diag_hdl)-(0.2784059*diag_triglyceride),
         lipid_pred_prob=exp(lipid_pred_score)/(1+exp(lipid_pred_score))) %>%
  
  analysis$cached("cohort", unique_index="patid")


cohort %>% count()
#45644

cohort %>% filter(!is.na(lipid_pred_prob)) %>% count()
#45644




