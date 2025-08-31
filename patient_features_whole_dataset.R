
# Pull out patient features for all in download - then can define prevalent (01/03/2023) and incident (at diagnosis) cohorts from this (ended up not using at diagnosis cohort for analysis)

## Define diagnosis date based on earliest (cleaned) diabetes code and add in age at diagnosis

## Pull out all diabetes type codes and define a) earliest, and b) diabetes type based on latest (pre-01/03/2023)

## Pull out and clean diabetes medication codes and define a) earliest insulin/OHA, b) current medication based on 01/03/2022-01/03/2023 period, c) treatment at 10-years post-diagnosis (based on data in previous year)

## Pull out and clean BMI, HDL, total cholesterol, triglyceride, HbA1c codes and define a) at diagnosis (up to 2 years previous), and b) latest post-diagnosis (pre-01/03/2023)

## Pull out and clean earliest positive/negative GAD antibody, IA2 antibody and C-peptide codes

## Define features at diagnosis x3 and autoimmune conditions

## Add hospitalisation for hypoglycaemia and DKA outcomes at diagnosis and earliest post-diagnosis


############################################################################################

# Setup
library(tidyverse)
library(aurum)
library(EHRBiomarkr)
rm(list=ls())

cprd = CPRDData$new(cprdEnv = "diabetes-jun2024", cprdConf = "~/.aurum.yaml")
codesets = cprd$codesets()
codes = codesets$getAllCodeSetVersion(v = "01/06/2024")

analysis = cprd$analysis("dpctn_paper")


############################################################################################

# Set index dates

latest_date <- as.Date("2023-03-01")


############################################################################################

# Collect all diabetes codes and use to define diagnosis date (earliest code - excluding those before DOB)

analysis = cprd$analysis("all_patid")

raw_diabetes_medcodes <- cprd$tables$observation %>%
  inner_join(codes$all_diabetes, by="medcodeid") %>%
  analysis$cached("raw_diabetes_medcodes", indexes=c("patid", "obsdate", "all_diabetes_cat"))


analysis = cprd$analysis("dpctn_paper")


#Look at codes before DOB

# test <- raw_diabetes_medcodes %>%
#   inner_join(cprd$tables$validDateLookup, by="patid") %>%
#   filter(obsdate<min_dob) %>%
#   collect()
#6,628 instance for 5,080 people
#2,821 1800, 2,324 1900, 814 1899, 180 1901, rest of years <26 instances


## Remove those before DOB or after latest index date - will reuse later for working out diabetes type
clean_diabetes_medcodes <- raw_diabetes_medcodes %>%
  inner_join(cprd$tables$validDateLookup, by="patid") %>%
  filter(obsdate>=min_dob & obsdate<=latest_date) %>%
  select(patid, date=obsdate, enterdate, category=all_diabetes_cat) %>%
  analysis$cached("clean_diabetes_medcodes", indexes=c("patid", "date", "category"))

## Find diagnosis dates - exclude gestational, gestational history and ambiguous pregnancy codes
diagnosis_dates <- clean_diabetes_medcodes %>%
  filter(category!="gestational" & category!="gestational history" & category!="pregnancy") %>%
  group_by(patid) %>%
  summarise(diagnosis_date=min(date, na.rm=TRUE)) %>%
  ungroup() %>%
  analysis$cached("diagnosis_dates", unique_indexes="patid")


############################################################################################

# Define diabetes type at diagnosis and at 01/03/2023

# Find earliest type code (type at diagnosis)

earliest_type_code <- clean_diabetes_medcodes %>%
  filter(category!="unspecified") %>%
  distinct(patid, date, category) %>%
  group_by(patid) %>%
  mutate(earliest_type_code_date=min(date, na.rm=TRUE)) %>%
  filter(date==earliest_type_code_date) %>%
  summarise(earliest_type_code=sql("group_concat(distinct category order by category separator ' & ')")) %>%
  ungroup() %>%
  analysis$cached("earliest_type_code", unique_indexes="patid")


## Find code counts at 01/03/2023 for each diabetes type

dm_code_counts <- clean_diabetes_medcodes %>%
  group_by(patid, category) %>%
  summarise(count=n()) %>%
  ungroup() %>%
  pivot_wider(id_cols=patid, names_from=category, values_from=count, values_fill=list(count=0L)) %>%
  analysis$cached("dm_code_counts", indexes="patid")

## Are a few people with insulin receptor abs 
## Have ambiguous pregnancy codes (not clear if gestational - fine, don't use for classification as not clear which category - will be unclassified if only have these codes
## Don't use gestational history codes either

## Find latest type code, excluding unspecified and pregnancy and history of gestational diabetes, and keep date

latest_type_code <- clean_diabetes_medcodes %>%
  filter(category!="unspecified" & category!="gestational - history" & category!="pregnancy") %>%
  group_by(patid) %>%
  mutate(most_recent_date=max(date, na.rm=TRUE),
         days_since_type_code=datediff(index_date, most_recent_date)) %>%
  filter(date==most_recent_date) %>%
  ungroup() %>%
  group_by(patid, days_since_type_code) %>%
  summarise(provisional_diabetes_type=sql("group_concat(distinct category order by category separator ' & ')")) %>%
  ungroup() %>%
  analysis$cached("latest_type_code", indexes="patid")


## Determine diabetes type at 01/03/2023

diabetes_type <- dm_code_counts %>%
  left_join(latest_type_code, by="patid") %>%
  mutate(diabetes_type=case_when(
    `type 1`==0 & `type 2`==0 & gestational==0 & `gestational history`==0 & malnutrition==0 & mody==0 & `other type not specified`==0 & `other/unspec genetic inc syndromic`==0 & secondary==0 & `insulin receptor abs`==0 ~ "unspecified",
    
    `type 1`>0 & `type 2`==0 & gestational==0 & `gestational history`==0 & malnutrition==0 & mody==0 & `other type not specified`==0 & `other/unspec genetic inc syndromic`==0 & secondary==0 & `insulin receptor abs`==0 ~ "type 1",
    
    `type 1`==0 & `type 2`>0 & gestational==0 & `gestational history`==0 & malnutrition==0 & mody==0 & `other type not specified`==0 & `other/unspec genetic inc syndromic`==0 & secondary==0 & `insulin receptor abs`==0 ~ "type 2",
    
    `type 1`==0 & `type 2`==0 & (gestational>0 | `gestational history`>0) & malnutrition==0 & mody==0 & `other type not specified`==0 & `other/unspec genetic inc syndromic`==0 & secondary==0 & `insulin receptor abs`==0 ~ "gestational",
    
    `type 1`==0 & `type 2`==0 & gestational==0 & `gestational history`==0 & malnutrition==0 & mody==0 & `other type not specified`==0 & `other/unspec genetic inc syndromic`==0 & secondary==0 & `insulin receptor abs`>0 ~ "insulin receptor abs",
    
    `type 1`==0 & `type 2`==0 & gestational==0 & `gestational history`==0 & malnutrition>0 & mody==0 & `other type not specified`==0 & `other/unspec genetic inc syndromic`==0 & secondary==0 & `insulin receptor abs`==0 ~ "malnutrition",
    
    `type 1`==0 & `type 2`==0 & gestational==0 & `gestational history`==0 & malnutrition==0 & mody>0 & `other type not specified`==0 & `other/unspec genetic inc syndromic`==0 & secondary==0 & `insulin receptor abs`==0 ~ "mody",
    
    `type 1`==0 & `type 2`==0 & gestational==0 & `gestational history`==0 & malnutrition==0 & mody==0 & `other type not specified`>0 & `other/unspec genetic inc syndromic`>0 & secondary==0 & `insulin receptor abs`==0 ~ "other type not specified",
    
    `type 1`==0 & `type 2`==0 & gestational==0 & `gestational history`==0 & malnutrition==0 & mody==0 & `other type not specified`==0 & `other/unspec genetic inc syndromic`>0 & secondary==0 & `insulin receptor abs`==0 ~ "other/unspec genetic inc syndromic",
    
    `type 1`==0 & `type 2`==0 & gestational==0 & `gestational history`==0 & malnutrition==0 & mody==0 & `other type not specified`==0 & `other/unspec genetic inc syndromic`==0 & secondary>0 & `insulin receptor abs`==0 ~ "secondary"),
    
    diabetes_type=ifelse(!is.na(diabetes_type), diabetes_type, paste("mixed;", provisional_diabetes_type)),
    
    type1_code_count= `type 1`,
    type2_code_count= `type 2`,
    mody_code_count=mody
    
  ) %>%
  
  mutate(gestational_dm_code_ever=ifelse(gestational>0 | `gestational history`>0, 1L, 0L),
         secondary_dm_code_ever=ifelse(secondary>0, 1L, 0L),
         other_dm_code_ever=ifelse(`other type not specified`>0 | `other/unspec genetic inc syndromic`>0 | `insulin receptor abs`>0 | malnutrition>0, 1L, 0L)) %>%
         
  select(patid, diabetes_type, type1_code_count, type2_code_count, mody_code_count, days_since_type_code, gestational_dm_code_ever, secondary_dm_code_ever, other_dm_code_ever) %>%
  
  analysis$cached("diabetes_type", unique_indexes="patid")


############################################################################################

# Add in earliest insulin and non-insulin diabetes meds
# And define current treatment (01/03/2022-01/03/2023)
# And define medication 10-year post-diagnosis

# Get clean OHA and insulin scripts

analysis = cprd$analysis("all_patid")

## All OHA scripts
raw_oha_prodcodes <- cprd$tables$drugIssue %>%
  inner_join(cprd$tables$ohaLookup, by="prodcodeid") %>%
  analysis$cached("raw_oha_prodcodes", indexes=c("patid", "issuedate"))

clean_oha_prodcodes <- raw_oha %>%
  inner_join(cprd$tables$validDateLookup, by="patid") %>%
  filter(issuedate>=min_dob & issuedate<=gp_end_date) %>%
  select(patid, date=issuedate, dosageid, quantity, quantunitid, duration, drug_class_1, drug_class_2, drug_substance_1, drug_substance_2) %>%
  analysis$cached("clean_oha_prodcodes", indexes=c("patid", "date", "drug_class_1", "drug_class_2", "drug_substance_1", "drug_substance_2"))

## All insulin scripts
raw_insulin_prodcodes <- cprd$tables$drugIssue %>%
  inner_join(codes$insulin, by="prodcodeid") %>%
  analysis$cached("raw_insulin_prodcodes", indexes=c("patid", "issuedate"))

clean_insulin_prodcodes <- raw_insulin_prodcodes %>%
  inner_join(cprd$tables$validDateLookup, by="patid") %>%
  filter(issuedate>=min_dob & issuedate<=gp_end_date) %>%
  select(patid, date=issuedate, dosageid, quantity, quantunitid, duration, insulin_cat) %>%
  analysis$cached("clean_insulin_prodcodes", indexes=c("patid", "date", "insulin_cat"))


analysis = cprd$analysis("diabetes_cohort")

## Earliest (clean) OHA script
first_oha <- clean_oha_prodcodes %>%
  group_by(patid) %>%
  summarise(dm_diag_ohadate=min(date, na.rm=TRUE)) %>%
  analysis$cached("first_oha", unique_index="patid")


## Earliest (clean) insulin script
first_insulin <- clean_insulin_prodcodes %>%
  group_by(patid) %>%
  summarise(dm_diag_insdate=min(date, na.rm=TRUE)) %>%
  analysis$cached("first_insulin", unique_index="patid")


analysis = cprd$analysis("dpctn_paper")

## Current medication
### Do drug classes separately and then combine (ignore Acarbose)

current_oha <- clean_oha_prodcodes %>%
  filter(date<=latest_date & datediff(latest_date, date)<=366) %>%
  distinct(patid, drug_class_1, drug_class_2) %>%
  pivot_longer(cols=c("drug_class_1", "drug_class_2")) %>%
  group_by(patid) %>%
  summarise(current_dpp4=any(value=="DPP4"),
            current_gipglp1=any(value=="GIPGLP1"),
            current_glinide=any(value=="Glinide"),
            current_glp1=any(value=="GLP1"),
            current_mfn=any(value=="MFN"),
            current_sglt2=any(value=="SGLT2"),
            current_su=any(value=="SU"),
            current_tzd=any(value=="TZD")) %>%
  ungroup() %>%
  analysis$cached("current_oha", unique_indexes="patid")

 
### Add insulin from insulin scripts (all prodcodes with OHA-insulin mixes are also in insulin)
current_insulin <- clean_insulin_prodcodes %>%
  filter(date<=latest_date & datediff(latest_date, date)<=366) %>%
  distinct(patid) %>%
  mutate(current_insulin=1L) %>%
  analysis$cached("current_insulin", unique_indexes="patid")


### Bolus insulin
current_bolus_insulin <- clean_insulin_prodcodes %>%
  filter(date<=latest_date & datediff(latest_date, date)<=366 & insulin_cat=="Bolus insulin") %>%
  distinct(patid) %>%
  mutate(current_bolus_insulin=1L) %>%
  analysis$cached("current_bolus_insulin", unique_indexes="patid")


### Mix insulin
current_mix_insulin <- clean_insulin_prodcodes %>%
  filter(date<=latest_date & datediff(latest_date, date)<=366 & insulin_cat=="Mix insulin") %>%
  distinct(patid) %>%
  mutate(current_mix_insulin=1L) %>%
  analysis$cached("current_mix_insulin", unique_indexes="patid")


## Add in treatment at 10 years-post-diagnosis for those registered years 9-10 and diagnosed >=1997 (2006 = year 9)
### Do drug classes separately and then combine (ignore Acarbose)

ten_yrs_post_diag <- diagnosis_dates %>%
  inner_join((cprd$tables$patient %>% select(patid, pracid, regstartdate, regenddate)), by="patid") %>%
  left_join((cprd$tables$practice %>% select(pracid, lcd)), by="pracid") %>%
  left_join((cprd$tables$onsDeath %>% select(patid, death_date=reg_date_of_death)), by="patid") %>%
  
  mutate(gp_record_end=pmin(if_else(is.na(lcd), as.Date("2024-06-30"), lcd),
                            if_else(is.na(regenddate), as.Date("2024-06-30"), regenddate),
                            if_else(is.na(death_date), as.Date("2024-06-30"), death_date),
                            as.Date("2024-06-30"), na.rm=TRUE)) %>%
  
  mutate(ten_yrs_post_diag=sql("date_add(diagnosis_date, interval 10 year)"),
         nine_yrs_post_diag=sql("date_add(diagnosis_date, interval 9 year)")) %>%
  mutate(ten_yrs_post_diag=if_else(year(diagnosis_date)>=1997 & regstartdate<=nine_yrs_post_diag & ten_yrs_post_diag<=gp_record_end & ten_yrs_post_diag<=latest_date, ten_yrs_post_diag, as.Date(NA))) %>%
  select(patid, gp_record_end, ten_yrs_post_diag) %>%
  analysis$cached("ten_yrs_post_diag", unique_indexes="patid")


ten_yrs_oha <- ten_yrs_post_diag %>%
  filter(!is.na(ten_yrs_post_diag)) %>%
  inner_join(clean_oha_prodcodes, by="patid") %>%
  filter(date<=ten_yrs_post_diag & datediff(ten_yrs_post_diag, date)<=366) %>%
  distinct(patid, drug_class_1, drug_class_2) %>%
  pivot_longer(cols=c("drug_class_1", "drug_class_2")) %>%
  group_by(patid) %>%
  summarise(ten_yrs_dpp4=any(value=="DPP4"),
            ten_yrs_gipglp1=any(value=="GIPGLP1"),
            ten_yrs_glinide=any(value=="Glinide"),
            ten_yrs_glp1=any(value=="GLP1"),
            ten_yrs_mfn=any(value=="MFN"),
            ten_yrs_sglt2=any(value=="SGLT2"),
            ten_yrs_su=any(value=="SU"),
            ten_yrs_tzd=any(value=="TZD")) %>%
  ungroup() %>%
  analysis$cached("ten_yrs_oha", unique_indexes="patid")


### Add insulin from insulin scripts (all prodcodes with OHA-insulin mixes are also in insulin)
ten_yrs_insulin <- ten_yrs_post_diag %>%
  filter(!is.na(ten_yrs_post_diag)) %>%
  inner_join(clean_insulin_prodcodes, by="patid") %>%
  filter(date<=ten_yrs_post_diag & datediff(ten_yrs_post_diag, date)<=366) %>%
  distinct(patid) %>%
  mutate(ten_yrs_insulin=1L) %>%
  analysis$cached("ten_yrs_insulin", unique_indexes="patid")


### Bolus insulin
ten_yrs_bolus_insulin <- ten_yrs_post_diag %>%
  filter(!is.na(ten_yrs_post_diag)) %>%
  inner_join(clean_insulin_prodcodes, by="patid") %>%
  filter(date<=ten_yrs_post_diag & datediff(ten_yrs_post_diag, date)<=366 & insulin_cat=="Bolus insulin") %>%
  distinct(patid) %>%
  mutate(ten_yrs_bolus_insulin=1L) %>%
  analysis$cached("ten_yrs_bolus_insulin", unique_indexes="patid")


# Combine and make current_oha and ten_yrs_oha for any OHA
# And whether on insulin at 1/2/3 years post-diagnosis

diabetes_meds <- cprd$tables$patient %>%
  select(patid, regstartdate) %>%
  inner_join(diagnosis_dates, by="patid") %>%
  left_join(first_oha, by="patid") %>%
  left_join(first_insulin, by="patid") %>%
  mutate(earliest_oha=if_else(dm_diag_ohadate>latest_date, as.Date(NA), dm_diag_ohadate),
         earliest_ins=if_else(dm_diag_insdate>latest_date, as.Date(NA), dm_diag_insdate)) %>%
  select(-c(dm_diag_ohadate, dm_diag_insdate)) %>%
  left_join(current_oha, by="patid") %>%
  left_join(current_insulin, by="patid") %>%
  left_join(current_bolus_insulin, by="patid") %>%
  left_join(current_mix_insulin, by="patid") %>%
  mutate(across(c("current_dpp4",
                  "current_glinide",
                  "current_gipglp1",
                  "current_glp1",
                  "current_mfn",
                  "current_sglt2",
                  "current_su",
                  "current_tzd",
                  "current_insulin",
                  "current_bolus_insulin",
                  "current_bolus_insulin"), coalesce, 0L)) %>%
  mutate(current_oha=ifelse(current_dpp4==1 | current_gipglp1==1 | current_glinide==1 | current_glp1==1 | current_mfn==1 | current_sglt2==1 | current_su==1 | current_tzd==1, 1L, 0L)) %>%
  left_join(ten_yrs_post_diag, by="patid") %>%
  left_join(ten_yrs_oha, by="patid") %>%
  left_join(ten_yrs_insulin, by="patid") %>%
  left_join(ten_yrs_bolus_insulin, by="patid") %>%
  mutate(across(c("ten_yrs_dpp4",
                  "ten_yrs_glinide",
                  "ten_yrs_gipglp1",
                  "ten_yrs_glp1",
                  "ten_yrs_mfn",
                  "ten_yrs_sglt2",
                  "ten_yrs_su",
                  "ten_yrs_tzd",
                  "ten_yrs_insulin",
                  "ten_yrs_bolus_insulin"), coalesce, 0L)) %>%
  mutate(ten_yrs_oha=ifelse(ten_yrs_dpp4==1 | ten_yrs_gipglp1==1 | ten_yrs_glinide==1 | ten_yrs_glp1==1 | ten_yrs_mfn==1 | ten_yrs_sglt2==1 | ten_yrs_su==1 | ten_yrs_tzd==1, 1L, 0L),
         ins_1_year=ifelse(year(diagnosis_date)<2006 | datediff(latest_date, diagnosis_date)<=366 | regstartdate>diagnosis_date, NA,
                           ifelse(!is.na(earliest_ins) & datediff(earliest_ins, diagnosis_date)<=366, 1L, 0L)),
         ins_2_years=ifelse(year(diagnosis_date)<2005 | datediff(latest_date, diagnosis_date)<=731 | datediff(regstartdate, diagnosis_date)>=366, NA,
                            ifelse(!is.na(earliest_ins) & datediff(earliest_ins, diagnosis_date)<=731, 1L, 0L)),
         ins_3_years=ifelse(year(diagnosis_date)<2004 | datediff(latest_date, diagnosis_date)<=1096 | datediff(regstartdate, diagnosis_date)>=731, NA,
                            ifelse(!is.na(earliest_ins) & datediff(earliest_ins, diagnosis_date)<=1096, 1L, 0L))) %>%
  select(-gp_record_end) %>%
  analysis$cached("diabetes_meds", unique_indexes="patid")


############################################################################################

# Add in biomarkers

# Clean biomarkers:
## Only keep those within acceptable value limits
## Only keep those with valid unit codes (numunitid)
## If multiple values on the same day, take mean
## Remove those with invalid dates (before DOB or after LCD/death/deregistration)

# Find values at diagnosis: up to 2 years prior and 7 days after
# And latest value at 01/03/2023 (as long as after diagnosis)


biomarkers <- c("bmi", "hba1c", "hdl", "totalcholesterol", "triglyceride")


# Pull out all raw biomarker values and cache

analysis = cprd$analysis("all_patid")

for (i in biomarkers) {
  
  print(i)
  
  raw_tablename <- paste0("raw_", i, "_medcodes")
  
  data <- cprd$tables$observation %>%
    inner_join(codes[[i]], by="medcodeid") %>%
    analysis$cached(raw_tablename, indexes=c("patid", "obsdate", "testvalue", "numunitid"))
  
  assign(raw_tablename, data)
  
}


# Clean biomarkers:
## Only keep those within acceptable value limits
## Only keep those with valid unit codes (numunitid)
## If multiple values on the same day, take mean
## Remove those with invalid dates (before min DOB or after LCD/death/deregistration)

## HbA1c only: remove before 1990, and convert all values to mmol/mol

analysis = cprd$analysis("all_patid")

for (i in biomarkers) {
  
  print(i)
  
  raw_tablename <- paste0("raw_", i, "_medcodes")
  clean_tablename <- paste0("clean_", i, "_medcodes")
  
  data <- get(raw_tablename)
  
  if (i=="hba1c") {
    
    data <- data %>%
      filter(year(obsdate)>=1990) %>%
      mutate(testvalue=ifelse(testvalue<=20, ((testvalue-2.152)/0.09148), testvalue))
    
  }
  
  data <- data %>%
    
    clean_biomarker_values(testvalue, i) %>%
    clean_biomarker_units(numunitid, i) %>%
    
    group_by(patid,obsdate) %>%
    summarise(testvalue=mean(testvalue, na.rm=TRUE)) %>%
    ungroup() %>%
    
    inner_join(cprd$tables$validDateLookup, by="patid") %>%
    filter(obsdate>=min_dob & obsdate<=gp_ons_end_date) %>%
    
    select(patid, date=obsdate, testvalue) %>%
    
    analysis$cached(clean_tablename, indexes=c("patid", "date", "testvalue"))
  
  assign(clean_tablename, data)
  
}


# For each biomarker, find values at diagnosis

analysis = cprd$analysis("dpctn_paper")

for (i in biomarkers) {
  
  print(i)
  
  clean_tablename <- paste0("clean_", i, "_medcodes")
  biomarker_diag_variable <- paste0("diag_",i)
  biomarker_diag_date_variable <- paste0("diag_",i, "_date")
  biomarker_diag_datediff_variable <- paste0("diag_",i, "datediff")
  interim_diag_biomarker_table <- paste0("diag_biomarkers_interim_", i)
  
  diag_value <- diagnosis_dates %>%
    inner_join(get(clean_tablename), by="patid") %>%
    mutate(diagdatediff=datediff(date, diagnosis_date)) %>%
    filter(diagdatediff<=7 & diagdatediff>=-730) %>%
    group_by(patid) %>%
    mutate(min_timediff=min(abs(diagdatediff), na.rm=TRUE)) %>%
    filter(abs(diagdatediff)==min_timediff) %>%
    mutate(pre_biomarker=min(testvalue, na.rm=TRUE)) %>%
    filter(pre_biomarker==testvalue) %>%
    dbplyr::window_order(diagdatediff) %>%
    filter(row_number()==1) %>%
    ungroup() %>%
    dbplyr::window_order() %>%
    select(patid, testvalue, date, diagdatediff) %>%
    rename({{biomarker_diag_variable}}:=testvalue,
           {{biomarker_diag_date_variable}}:=date,
           {{biomarker_diag_datediff_variable}}:=diagdatediff) %>%
    analysis$cached(interim_diag_biomarker_table, unique_indexes="patid")
  
    assign(interim_diag_biomarker_table, diag_value)
    
}


# For each biomarker, find baseline value at index date
## Use closest date to index as long as prior to this

analysis = cprd$analysis("dpctn_paper")

for (i in biomarkers) {
  
  print(i)
  
  clean_tablename <- paste0("clean_", i, "_medcodes")
  biomarker_index_variable <- paste0("index_",i)
  biomarker_index_date_variable <- paste0("index_",i, "_date")
  biomarker_index_datediff_variable <- paste0("index_",i, "datediff")
  interim_index_biomarker_table <- paste0("index_biomarkers_interim_", i)
  
  index_value <- diagnosis_dates %>%
    inner_join(get(clean_tablename), by="patid") %>%
    mutate(indexdatediff=datediff(date, latest_date)) %>%
    filter(indexdatediff<=0 & date>diagnosis_date) %>%
    group_by(patid) %>%
    mutate(min_timediff=min(abs(indexdatediff), na.rm=TRUE)) %>%
    filter(abs(indexdatediff)==min_timediff) %>%
    ungroup() %>%
    select(patid, testvalue, date, indexdatediff) %>%
    rename({{biomarker_index_variable}}:=testvalue,
           {{biomarker_index_date_variable}}:=date,
           {{biomarker_index_datediff_variable}}:=indexdatediff) %>%
    analysis$cached(interim_index_biomarker_table, unique_indexes="patid")
  
    assign(interim_index_biomarker_table, index_value)
  
}


biomarkers <- cprd$tables$patient %>%
  select(patid) %>%
  left_join(diag_biomarkers_interim_bmi, by="patid") %>%
  left_join(diag_biomarkers_interim_hba1c, by="patid") %>%
  left_join(diag_biomarkers_interim_hdl, by="patid") %>%
  left_join(diag_biomarkers_interim_totalcholesterol, by="patid") %>%
  left_join(diag_biomarkers_interim_triglyceride, by="patid") %>%
  left_join(index_biomarkers_interim_bmi, by="patid") %>%
  left_join(index_biomarkers_interim_hba1c, by="patid") %>%
  left_join(index_biomarkers_interim_hdl, by="patid") %>%
  left_join(index_biomarkers_interim_totalcholesterol, by="patid") %>%
  left_join(index_biomarkers_interim_triglyceride, by="patid") %>%
  analysis$cached("biomarkers", unique_indexes="patid")
  

############################################################################################

# Add antibodies
## Have checked and none are before DOB

analysis = cprd$analysis("all_patid")

# GAD
raw_gad <- cprd$tables$observation %>%
  inner_join(codes$gad, by="medcodeid") %>%
  analysis$cached("raw_gad_medcodes", indexes=c("patid", "obsdate", "testvalue", "numunitid"))
#n=11,416

# IA2
raw_ia2 <- cprd$tables$observation %>%
  inner_join(codes$ia2, by="medcodeid") %>%
  analysis$cached("raw_ia2_medcodes", indexes=c("patid", "obsdate", "testvalue", "numunitid"))
#n=9,042


# Assume those with missing units are in U/mL

print(raw_gad %>% group_by(numunitid) %>% summarise(count=n()) %>% left_join(cprd$tables$numUnit), n=10000)
## GAD: ~half are missing (5,972), then rest are U/mL / kU/L (equivalent) / 'units'; 2 n/ml, 1 mmol/L, 1 mu/L, 1 u/L - exclude

print(raw_ia2 %>% group_by(numunitid) %>% summarise(count=n()) %>% left_join(cprd$tables$numUnit), n=10000)
## IA2: majority missing (8,317), then rest are U/mL / kU/L (equivalent) / 'units'; 2 u/L - exclude

# Assume 0 values are 0, not missing

analysis = cprd$analysis("dpctn_paper")

### GAD

gad <- raw_gad %>%
  filter(obsdate<=latest_date & (is.na(numunitid) | (numunitid!=229 & numunitid!=11589 & numunitid!=218 & numunitid!=276))) %>%
  mutate(result=ifelse(!is.na(testvalue) & testvalue<11, "negative",
                       ifelse(!is.na(testvalue) & testvalue>=11, "positive", NA))) %>%
  filter(!is.na(result)) %>%
  distinct(patid, obsdate, result)

earliest_gad <- gad %>%
  group_by(patid, result) %>%
  mutate(earliest_gad=min(obsdate, na.rm=TRUE)) %>%
  filter(obsdate==earliest_gad) %>%
  ungroup() %>%
  pivot_wider(id_cols="patid", names_from=result, values_from=obsdate)

latest_gad <- gad %>%
  group_by(patid, result) %>%
  mutate(latest_gad=max(obsdate, na.rm=TRUE)) %>%
  filter(obsdate==latest_gad) %>%
  ungroup() %>%
  pivot_wider(id_cols="patid", names_from=result, values_from=obsdate)

gad <- earliest_gad %>%
  rename(earliest_negative_gad=negative, earliest_positive_gad=positive) %>%
  left_join((latest_gad %>% rename(latest_negative_gad=negative, latest_positive_gad=positive)), by="patid") %>%
  select(patid, earliest_negative_gad, latest_negative_gad, earliest_positive_gad, latest_positive_gad) %>%
  analysis$cached("gad", unique_indexes="patid")


### IA2

ia2 <- raw_ia2 %>%
  filter(obsdate<=latest_date & (is.na(numunitid) | numunitid!=276)) %>%
  mutate(result=ifelse(!is.na(testvalue) & testvalue<7.5, "negative",
                       ifelse(!is.na(testvalue) & testvalue>=7.5, "positive", NA))) %>%
  filter(!is.na(result)) %>%
  distinct(patid, obsdate, result)

earliest_ia2 <- ia2 %>%
  group_by(patid, result) %>%
  mutate(earliest_ia2=min(obsdate, na.rm=TRUE)) %>%
  filter(obsdate==earliest_ia2) %>%
  ungroup() %>%
  pivot_wider(id_cols="patid", names_from=result, values_from=obsdate)

latest_ia2 <- ia2 %>%
  group_by(patid, result) %>%
  mutate(latest_ia2=max(obsdate, na.rm=TRUE)) %>%
  filter(obsdate==latest_ia2) %>%
  ungroup() %>%
  pivot_wider(id_cols="patid", names_from=result, values_from=obsdate)

ia2 <- earliest_ia2 %>%
  rename(earliest_negative_ia2=negative, earliest_positive_ia2=positive) %>%
  left_join((latest_ia2 %>% rename(latest_negative_ia2=negative, latest_positive_ia2=positive)), by="patid") %>%
  select(patid, earliest_negative_ia2, latest_negative_ia2, earliest_positive_ia2, latest_positive_ia2) %>%
  analysis$cached("ia2", unique_indexes="patid")


############################################################################################

# Add in C-peptide
## Have checked and none are before DOB

analysis = cprd$analysis("all_patid")

raw_c_peptide <- cprd$tables$observation %>%
  inner_join(codes$c_peptide, by="medcodeid") %>%
  analysis$cached("raw_c_peptide_medcodes", indexes=c("patid", "obsdate", "testvalue", "numunitid"))

raw_c_peptide %>% count()
#15,050


# Look at units

print(raw_c_peptide %>% filter(c_peptide_cat=="ucpcr") %>% group_by(numunitid) %>% summarise(count=n()) %>% left_join(cprd$tables$numUnit), n=10000)
## UCPCR: most (~2000/2200) are nmol/mmol creatinine (codes 2434, 899, 1744), 185 missing, 43 umol/mol (code 959; same units)
## Also have 8 in different units (nmol/L, mmol/L, ng/mL) - exclude (codes: 235, 218, 233)


print(raw_c_peptide %>% filter(c_peptide_cat=="blood") %>% group_by(numunitid) %>% summarise(count=n()) %>% left_join(cprd$tables$numUnit), n=10000)
## Blood C-peptide: majority (~8000/13000) are pmol/L (codes 256, 23399), ~4000 missing
## 577 nmol/L (code 235) - convert to pmol/L (divide by 1000)
## 2 mmol/L (code 218) - convert to pmol/L (divide by 10^6)
## 2 pmol/ml (code 9244) - convert to pmol/L (multiply by 1000)
## Exclude ug/L (code 283), nmol/mmol (code 899), U/mL (code 277), mu/L (code 229)


# Clean and define insulin status
## Use thresholds of <0.2 and >=0.6 (UCPCR) / <200 and >=600 (blood) to define insulin status (https://www.exeterlaboratory.com/test/c-peptide-urine/)
### No indication as to whether blood samples are fasted or not so assume not
## Assume 0 values are 0, not missing

analysis = cprd$analysis("dpctn_paper")

processed_c_peptide <- raw_c_peptide %>%
  filter(obsdate<=latest_date & !is.na(testvalue) & !(c_peptide_cat=="ucpcr" & (numunitid==218 | numunitid==233 | numunitid==235)) & !(c_peptide_cat=="blood" & (numunitid==283 | numunitid==899 | numunitid==277 | numunitid==229))) %>%
  mutate(new_testvalue=ifelse(numunitid==235, testvalue/1000,
                              ifelse(numunitid==218, testvalue/1000000,
                                     ifelse(numunitid==9244, testvalue*1000, testvalue)))) %>%
  mutate(c_pep_insulin=ifelse(c_peptide_cat=="ucpcr" & new_testvalue<0.2, "c_pep_ins_deficient",
                              ifelse(c_peptide_cat=="ucpcr" & new_testvalue>=0.2 & testvalue<0.6, "c_pep_ins_intermediate",
                                     ifelse(c_peptide_cat=="ucpcr" & new_testvalue>=0.6, "c_pep_ins_normal",
                                            ifelse(c_peptide_cat=="blood" & new_testvalue<200, "c_pep_ins_deficient",
                                                   ifelse(c_peptide_cat=="blood" & new_testvalue>=200 & testvalue<600, "c_pep_ins_intermediate",
                                                          ifelse(c_peptide_cat=="blood" & new_testvalue>=600, "c_pep_ins_normal", NA))))))) %>%
  select(patid, date=obsdate, testvalue=new_testvalue, c_pep_cat=c_peptide_cat, c_pep_insulin) %>%
  distinct() %>%
  analysis$cached("processed_c_peptide", indexes=c("patid", "date", "c_pep_insulin"))


# Add earliest and latest result for each category to table

earliest_c_pep <- processed_c_peptide %>%
  group_by(patid, c_pep_insulin) %>%
  mutate(earliest=min(date, na.rm=TRUE)) %>%
  filter(date==earliest) %>%
  ungroup() %>%
  pivot_wider(id_cols="patid", names_from=c_pep_insulin, values_from=date)

latest_c_pep <- processed_c_peptide %>%
  group_by(patid, c_pep_insulin) %>%
  mutate(latest=max(date, na.rm=TRUE)) %>%
  filter(date==latest) %>%
  ungroup() %>%
  pivot_wider(id_cols="patid", names_from=c_pep_insulin, values_from=date)

c_pep <- earliest_c_pep %>%
  rename(earliest_c_pep_ins_deficient=c_pep_ins_deficient, earliest_c_pep_ins_intermediate=c_pep_ins_intermediate, earliest_c_pep_ins_normal=c_pep_ins_normal) %>%
  left_join((latest_c_pep %>% rename(latest_c_pep_ins_deficient=c_pep_ins_deficient, latest_c_pep_ins_intermediate=c_pep_ins_intermediate, latest_c_pep_ins_normal=c_pep_ins_normal)), by="patid") %>%
  select(patid, earliest_c_pep_ins_deficient, latest_c_pep_ins_deficient, earliest_c_pep_ins_intermediate, latest_c_pep_ins_intermediate, earliest_c_pep_ins_normal, latest_c_pep_ins_normal) %>%
  analysis$cached("c_pep", unique_indexes="patid")


############################################################################################

# Add in family history of diabetes (cleaned: includes 99% of raw occurrences so no difference)
# Codelist has changed since running main analysis - categories updated
## Want to include family member count and diabetes type if available
## Now including positive siblings and children

analysis = cprd$analysis("all_patid")

## Raw FH of diabetes
raw_fh_diabetes_medcodes <- cprd$tables$observation %>%
  inner_join(codes$fh_diabetes, by="medcodeid") %>%
  analysis$cached("raw_fh_diabetes_medcodes", indexes=c("patid", "obsdate"))

## Clean FH of diabetes
clean_fh_diabetes_medcodes <- raw_fh_diabetes_medcodes %>%
  inner_join(cprd$tables$validDateLookup, by="patid") %>%
  filter(obsdate>=min_dob & obsdate<=gp_end_date) %>%
  select(patid, date=obsdate, fh_diabetes_cat) %>%
  analysis$cached("clean_fh_diabetes_medcodes", indexes=c("patid", "date", "fh_diabetes_cat"))


analysis = cprd$analysis("dpctn_paper")


# First - determine whether positive or negative
# For people with positive and negative codes:
## If all negative codes are earlier than positive codes, fine - use positive
## Otherwise, treat as missing

fh_positive_negative_interim <- clean_fh_diabetes_medcodes %>%
  filter(fh_diabetes_cat!="positive - gestational" & date<=latest_date) %>%
  mutate(fh_diabetes_cat=ifelse(fh_diabetes_cat=="negative", "negative", "positive")) %>%
  group_by(patid, fh_diabetes_cat) %>%
  summarise(earliest_date=min(date, na.rm=TRUE),
            latest_date=max(date, na.rm=TRUE)) %>%
  ungroup() %>%
  group_by(patid) %>%
  pivot_wider(id_cols=patid, names_from = c(fh_diabetes_cat), names_glue = "{fh_diabetes_cat}_{.value}", values_from=c(earliest_date, latest_date)) %>%
  ungroup() %>%
  analysis$cached("fh_positive_negative_interim", unique_indexes="patid")

fh_positive_negative <- fh_positive_negative_interim %>%
  mutate(fh_diabetes=ifelse(is.na(positive_earliest_date), 0L,
                            ifelse(is.na(negative_earliest_date), 1L,
                                   ifelse(!is.na(positive_earliest_date) & !is.na(negative_earliest_date) & negative_latest_date<positive_earliest_date, 1L, NA)))) %>%
  analysis$cached("fh_positive_negative", unique_indexes="patid")


############################################################################################

# Features at diagnosis (weight loss, polydipsia, polyuria)

analysis = cprd$analysis("all_patid")

## Raw weight loss medcodes
raw_weight_loss_medcodes <- cprd$tables$observation %>%
  inner_join(codes$weight_loss, by="medcodeid") %>%
  analysis$cached("raw_weight_loss_medcodes", indexes=c("patid", "obsdate"))

## Clean weight loss medcodes
clean_weight_loss_medcodes <- raw_weight_loss_medcodes %>%
  inner_join(cprd$tables$validDateLookup, by="patid") %>%
  filter(obsdate>=min_dob & obsdate<=gp_end_date) %>%
  select(patid, date=obsdate) %>%
  analysis$cached("clean_weight_loss_medcodes", indexes=c("patid", "date"))


## Raw polydipsia medcodes
raw_polydipsia_medcodes <- cprd$tables$observation %>%
  inner_join(codes$polydipsia, by="medcodeid") %>%
  analysis$cached("raw_polydipsia_medcodes", indexes=c("patid", "obsdate"))

## Clean polydipsia medcodes
clean_polydipsia_medcodes <- raw_polydipsia_medcodes %>%
  inner_join(cprd$tables$validDateLookup, by="patid") %>%
  filter(obsdate>=min_dob & obsdate<=gp_end_date) %>%
  select(patid, date=obsdate) %>%
  analysis$cached("clean_polydipsia_medcodes", indexes=c("patid", "date"))


## Raw urinary frequency medcodes
raw_urinary_frequency_medcodes <- cprd$tables$observation %>%
  inner_join(codes$urinary_frequency, by="medcodeid") %>%
  analysis$cached("raw_urinary_frequency_medcodes", indexes=c("patid", "obsdate"))

## Clean urinary frequency medcodes
clean_urinary_frequency_medcodes <- raw_urinary_frequency_medcodes %>%
  inner_join(cprd$tables$validDateLookup, by="patid") %>%
  filter(obsdate>=min_dob & obsdate<=gp_end_date) %>%
  select(patid, date=obsdate) %>%
  analysis$cached("clean_urinary_frequency_medcodes", indexes=c("patid", "date"))


# See if any have features at diagnosis
## Use year of diagnosis

analysis = cprd$analysis("dpctn_paper")

weight_loss_at_diag <- diagnosis_dates %>%
  inner_join(clean_weight_loss_medcodes, by="patid") %>%
  filter(year(date)==year(diagnosis_date)) %>%
  distinct(patid) %>%
  mutate(weight_loss_at_diag=1L) %>%
  analysis$cached("weight_loss_at_diag", unique_indexes="patid")

polydipsia_at_diag <- diagnosis_dates %>%
  inner_join(clean_polydipsia_medcodes, by="patid") %>%
  filter(year(date)==year(diagnosis_date)) %>%
  distinct(patid) %>%
  mutate(polydipsia_at_diag=1L) %>%
  analysis$cached("polydipsia_at_diag", unique_indexes="patid")

urinary_freq_at_diag <- diagnosis_dates %>%
  inner_join(clean_urinary_frequency_medcodes, by="patid") %>%
  filter(year(date)==year(diagnosis_date)) %>%
  distinct(patid)  %>%
  mutate(urinary_freq_at_diag=1L) %>%
  analysis$cached("urinary_freq_at_diag", unique_indexes="patid")


## Combine

at_diag_features <- cprd$tables$patient %>%
  select(patid) %>%
  left_join(weight_loss_at_diag, by="patid") %>%
  left_join(polydipsia_at_diag, by="patid") %>%
  left_join(urinary_freq_at_diag, by="patid") %>%
  mutate(across(c("weight_loss_at_diag",
                  "polydipsia_at_diag",
                  "urinary_freq_at_diag"), coalesce, 0L)) %>%
  analysis$cached("at_diag_features", unique_indexes="patid")


############################################################################################

# Add autoimmune conditions

coeliac_codes <- read_delim("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2025/Julieanne progression/Final/autoimmune codelists/coeliac.txt")

thyroid_codes <- read_delim("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2025/Julieanne progression/Final/autoimmune codelists/hashimotos.txt")


raw_coeliac <- cprd$tables$observation %>%
  inner_join(coeliac_codes, by=c("medcodeid"="MedCodeId"), copy=TRUE) %>%
  analysis$cached("raw_coeliac", indexes=c("patid", "obsdate"))

raw_thyroid <- cprd$tables$observation %>%
  inner_join(thyroid_codes, by=c("medcodeid"="MedCodeId"), copy=TRUE) %>%
  analysis$cached("raw_thyroid", indexes=c("patid", "obsdate"))

clean_coeliac <- raw_coeliac %>%
  select(patid, date=obsdate) %>%
  inner_join(cprd$tables$validDateLookup, by="patid") %>%
  filter(date>=min_dob & date<=latest_date) %>%
  distinct(patid) %>%
  mutate(coeliac:=1L) %>%
  analysis$cached("clean_coeliac", unique_index="patid")

clean_thyroid <- raw_thyroid %>%
  select(patid, date=obsdate) %>%
  inner_join(cprd$tables$validDateLookup, by="patid") %>%
  filter(date>=min_dob & date<=latest_date) %>%
  distinct(patid) %>%
  mutate(autoimmune_thyroid:=1L) %>%
  analysis$cached("clean_thyroid", unique_index="patid")


autoimmune <- cprd$tables$patient %>%
  select(patid) %>%
  left_join(clean_coeliac, by="patid") %>%
  left_join(clean_thyroid, by="patid") %>%
  mutate(coeliac=ifelse(is.na(coeliac), 0L, 1L),
         autoimmune_thyroid=ifelse(is.na(autoimmune_thyroid), 0L, 1L),
         coeliac_thyroid=ifelse(coeliac==1 | autoimmune_thyroid==1, 1L, 0L)) %>%
  analysis$cached("autoimmune", unique_indexes="patid")


############################################################################################

# Add hospitalisation outcomes

# Add in DKA at diagnosis, earliest post-diagnosis DKA and hypo, and no. of DKA and hypo postdiag visits

dka_at_diagnosis <- cprd$tables$hesDiagnosisEpi %>%
  inner_join(codes$icd10_dka, by=c("ICD"="icd10")) %>%
  inner_join(diagnosis_dates, by="patid") %>%
  filter(d_order==1 & epistart<=latest_date & abs(datediff(epistart, diagnosis_date))<=30) %>%
  distinct(patid) %>%
  mutate(dka_at_diagnosis=1L) %>%
  analysis$cached("dka_at_diagnosis", unique_indexes="patid")

earliest_dka <- cprd$tables$hesDiagnosisEpi %>%
  inner_join(codes$icd10_dka, by=c("ICD"="icd10")) %>%
  inner_join(diagnosis_dates, by="patid") %>%
  filter(d_order==1 & epistart<=latest_date & datediff(epistart, diagnosis_date)>30) %>%
  group_by(patid) %>%
  summarise(earliest_dka=min(epistart, na.rm=TRUE),
            dka_count=n()) %>%
  ungroup() %>%
  analysis$cached("earliest_dka", unique_indexes="patid")

earliest_hypo <- cprd$tables$hesDiagnosisEpi %>%
  inner_join(codes$icd10_hypoglycaemia, by=c("ICD"="icd10")) %>%
  inner_join(diagnosis_dates, by="patid") %>%
  filter(d_order==1 & epistart<=latest_date & epistart>diagnosis_date) %>%
  group_by(patid) %>%
  summarise(earliest_hypo=min(epistart, na.rm=TRUE),
            hypo_count=n()) %>%
  ungroup() %>%
  analysis$cached("earliest_hypo", unique_indexes="patid")


## No. of contacts in last year (just using unique patid-date interactions as per T1 paper I reviewed)

contact_count <- cprd$tables$observation %>%
  filter(obsdate<=latest_date & datediff(latest_date, obsdate)<=366) %>%
  distinct(patid, obsdate) %>%
  group_by(patid) %>%
  summarise(contact_count=n()) %>%
  ungroup() %>%
  analysis$cached("contact_count", unique_indexes="patid")



# Combine all

dka_hypo_contacts <- cprd$tables$patient %>%
  select(patid) %>%
  left_join(dka_at_diagnosis, by="patid") %>%
  mutate(dka_at_diagnosis=ifelse(is.na(dka_at_diagnosis), 0L, 1L)) %>%
  left_join(earliest_dka, by="patid") %>%
  left_join(earliest_hypo, by="patid") %>%
  left_join(contact_count, by="patid") %>%
  mutate(dka_post_diagnosis=ifelse(!is.na(earliest_dka), 1L, 0L),
         hypo_post_diagnosis=ifelse(!is.na(earliest_hypo), 1L, 0L)) %>%
  analysis$cached("dka_hypo_contacts", unique_indexes="patid")

