
# Define prevalent cohort of insulin-treated T1 and T2 diagnosed aged 18-50 years:

## All IDs in download who are registered on 01/01/2024

## Define diabetes type based on (latest) type codes, and only keep if current diagnosis is Type 1 or 2

## Define diagnosis date based on earliest diabetes code: exclude those diagnosed between -30 and +90 days of registration and those diagnosed <18 or >50 years

## Define current treatment - at this stage don't exclude non-insulin treated as will want to flag as part of algorithm

## Pull in variables for T1T2 calculator: current BMI, HDL, total cholesterol and triglyceride

## Pull in other variables of interest: current HbA1c, antibodies, C-peptide, earliest insulin and OHA, features at diagnosis, autoimmune conditions

## (When have linkage data): add hospitalisation for hypoglycaemic and DKA outcomes

## Run lipid model


############################################################################################

# Setup
library(tidyverse)
library(aurum)
library(EHRBiomarkr)
library(lubridate)
rm(list=ls())

cprd = CPRDData$new(cprdEnv = "diabetes-jun2024", cprdConf = "~/.aurum.yaml")
codesets = cprd$codesets()
codes = codesets$getAllCodeSetVersion(v = "01/06/2024")

analysis = cprd$analysis("dpctn_paper")


############################################################################################

# Set index date

index_date <- as.Date("2023-03-01")


############################################################################################

# Find IDs of those registered on index date

# Include general patient variables as per all_diabetes_cohort:
## gender
## DOB (derived as per: https://github.com/Exeter-Diabetes/CPRD-Codelists/blob/main/readme.md#general-notes-on-implementation, cached for all IDs in 'diagnosis_date_dob' table from all_patid_dob script)
## pracid
## prac_region
## ethnicity 5-category and 16-category (derived as per: https://github.com/Exeter-Diabetes/CPRD-Codelists#ethnicity, from GP data only at the moment)
## regstartdate	
## gp_record_end (earliest of last collection date from practice, deregistration and 30/06/2024 (latest date in records))	
## death_date	('cprddeathdate')
## Calculate age at index date


# Just want those registered on the index date

analysis = cprd$analysis("diabetes_cohort")

dob <- dob %>% analysis$cached("dob")

analysis = cprd$analysis("all_patid")
ethnicity <- ethnicity %>% analysis$cached("ethnicity")


analysis = cprd$analysis("dpctn_paper")

cohort <- cprd$tables$patient %>%
  select(patid, gender, pracid, regstartdate, regenddate, cprd_ddate) %>%
  left_join((dob %>% select(patid, dob)), by="patid") %>%
  left_join((cprd$tables$practice %>% select(pracid, prac_region=region, lcd)), by="pracid") %>%
  left_join((ethnicity %>% select(patid, ethnicity_5cat, ethnicity_16cat)), by="patid") %>%
  
  mutate(gp_record_end=pmin(if_else(is.na(lcd), as.Date("2024-06-30"), lcd),
                            if_else(is.na(regenddate), as.Date("2024-06-30"), regenddate),
                            as.Date("2024-06-30"), na.rm=TRUE),
         
         death_date=cprd_ddate,
         
         age_at_index=round(datediff(index_date, dob)/365.25, 1)) %>%
  
  select(patid, gender, dob, age_at_index, pracid, prac_region, ethnicity_5cat, ethnicity_16cat, regstartdate, gp_record_end, death_date) %>%
  
  filter(regstartdate<=index_date & !(!is.na(death_date) & death_date<index_date) & !(!is.na(gp_record_end) & gp_record_end<index_date)) %>%
  
  analysis$cached("cohort_interim_1", unique_indexes="patid")

cohort %>% count()
# 1,373,214


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


analysis = cprd$analysis("dpctn_paper")

cohort <- cohort %>%
  anti_join(practice_exclusion_ids, by="patid") %>%
  analysis$cached("cohort_interim_2", unique_indexes="patid")

cohort %>% count()
# 1,373,214



# Just keep men and women

cohort <- cohort %>%
  filter(gender==1 | gender==2) %>%
  analysis$cached("cohort_interim_3", unique_indexes="patid")

cohort %>% count()
# 1,373,188


############################################################################################

# Define diabetes type (NB: not removing codes before DOB for this)

## Find all diabetes codes prior to/at index date

analysis = cprd$analysis("all_patid")

raw_diabetes_medcodes <- cprd$tables$observation %>%
  inner_join(new_codes$all_diabetes, by="medcodeid") %>%
  analysis$cached("raw_diabetes_medcodes", indexes=c("patid", "obsdate", "all_diabetes_cat"))


analysis = cprd$analysis("dpctn_paper")

dm_codes <- cohort %>%
  left_join(raw_diabetes_medcodes, by="patid") %>%
  select(patid, date=obsdate, enterdate, category=all_diabetes_cat) %>%
  filter(date<=index_date) %>%
  analysis$cached("dm_codes", indexes=c("patid", "date", "category"))

## Check how have diabetes codes before index date in download - won't be everyone
#dm_codes %>% distinct(patid) %>% count()
#1,252,239 - have lost 178,649

  
## Find code counts for each diabetes type

dm_code_counts <- dm_codes %>%
  group_by(patid, category) %>%
  summarise(count=n()) %>%
  ungroup() %>%
  pivot_wider(id_cols=patid, names_from=category, values_from=count, values_fill=list(count=0L)) %>%
  analysis$cached("dm_code_counts", indexes="patid")

## Are a few people with insulin receptor abs now 
## Now have ambiguous pregnancy codes - fine, don't use for classification as not clear which category - will be unclassified if only have these codes
## Don't use gestational history codes either

## Find latest type code, excluding unspecified and pregnancy, and keep date

latest_type_code <- dm_codes %>%
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


## Determine diabetes type and add to cohort table

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
  
  select(patid, diabetes_type, type1_code_count, type2_code_count, mody_code_count, days_since_type_code) %>%
  
  analysis$cached("diabetes_type", unique_indexes="patid")


# Keep those with Type 1 or 2 diagnosis from most code (even if mixed previously)
## Remove those where most recent code is type 1 or 2 but they have another code on the same day
## Have a look at people who are excluded

# clipr::write_clip(cohort %>%
#   inner_join(diabetes_type, by="patid") %>%
#   filter(diabetes_type!="type 1" & diabetes_type!="type 2" & diabetes_type!="mixed; type 1" & diabetes_type!="mixed; type 2") %>%
#   group_by(diabetes_type) %>%
#   count() %>%
#   collect())

cohort <- cohort %>%
  inner_join(diabetes_type, by="patid") %>%
  filter(diabetes_type=="type 1" | diabetes_type=="type 2" | diabetes_type=="mixed; type 1" | diabetes_type=="mixed; type 2") %>%
  analysis$cached("cohort_interim_4", unique_indexes="patid")

cohort %>% count()
#915,836


############################################################################################

# Define diagnosis dates
## Remove codes if before DOB
## Exclude codes if Type 2 and in year in birth
## Remove gestational, gestational history and pregnancy codes

## Check how many have diabetes codes before DOB only
#dm_codes %>% inner_join(cohort, by="patid") %>% inner_join(cprd$tables$validDateLookup, by="patid") %>% filter(date>=min_dob) %>% distinct(patid) %>% count()
#915,836
#no-one with codes before DOB only


## Check how many of T2s have codes in YOB only
#dm_codes %>% inner_join(cohort, by="patid") %>% inner_join(cprd$tables$validDateLookup, by="patid") %>% filter(!((diabetes_type=="type 2" | diabetes_type=="mixed; type 2") & year(date)==year(dob))) %>% distinct(patid) %>% count()
#915,836
#no-one with Type 2 and codes in YOB only

diagnosis_dates <- dm_codes %>%
  inner_join(cprd$tables$validDateLookup, by="patid") %>%
  filter(date>=min_dob) %>%
  inner_join(cohort, by="patid") %>%
  filter(!(diabetes_type=="type 2" & year(date)==year(dob)) & diabetes_type!="gestational" & diabetes_type!="gestational history" & diabetes_type!="pregnancy") %>%
  group_by(patid) %>%
  summarise(diagnosis_date=min(date, na.rm=TRUE)) %>%
  ungroup() %>%
  analysis$cached("diagnosis_dates", unique_indexes="patid")

diagnosis_dates %>% count()
#915,836


# Add to cohort table and calculate follow up time and age at diagnosis
## exclude if not diagnosed aged 18-50

cohort <- cohort %>%
  inner_join(diagnosis_dates, by="patid") %>%
  mutate(dm_diag_age=round((datediff(diagnosis_date, dob))/365.25, 1),
         follow_up_time=datediff(index_date, diagnosis_date)/365.25) %>%
  filter(dm_diag_age>=18 & dm_diag_age<51) %>%
  analysis$cached("cohort_interim_5", unique_indexes="patid")

cohort %>% count()
#334,304


# Exclude if diagnosis date within -30 to +90 days of registration start

cohort <- cohort %>%
  filter(datediff(diagnosis_date, regstartdate)< -30 | datediff(diagnosis_date, regstartdate)>90) %>%
  analysis$cached("cohort_interim_6", unique_indexes="patid")

cohort %>% count()
#316,076


############################################################################################

# Add in current treatment (last year)
## Not resticting to those on insulin at the moment in case need characteristics

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
  inner_join(new_codes$insulin, by="prodcodeid") %>%
  analysis$cached("raw_insulin_prodcodes", indexes=c("patid", "issuedate"))

clean_insulin_prodcodes <- raw_insulin_prodcodes %>%
  inner_join(cprd$tables$validDateLookup, by="patid") %>%
  filter(issuedate>=min_dob & issuedate<=gp_end_date) %>%
  select(patid, date=issuedate, dosageid, quantity, quantunitid, duration, insulin_cat) %>%
  analysis$cached("clean_insulin_prodcodes", indexes=c("patid", "date", "insulin_cat"))

analysis = cprd$analysis("dpctn_paper")


# Do drug classes separately and then combine (ignore Acarbose)

current_meds <- clean_oha_prodcodes %>%
  filter(date<=index_date & datediff(index_date, date)<=366) %>%
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
  analysis$cached("current_meds", unique_indexes="patid")

 
## Add insulin from insulin scripts (all prodcodes with OHA-insulin mixes are also in insulin)

current_insulin <- clean_insulin_prodcodes %>%
  filter(date<=index_date & datediff(index_date, date)<=366) %>%
  distinct(patid) %>%
  mutate(current_insulin=1L) %>%
  analysis$cached("current_insulin", unique_indexes="patid")


## Bolus insulin

current_bolus_insulin <- clean_insulin_prodcodes %>%
  filter(date<=index_date & datediff(index_date, date)<=366 & insulin_cat=="Bolus insulin") %>%
  distinct(patid) %>%
  mutate(current_bolus_insulin=1L) %>%
  analysis$cached("current_bolus_insulin", unique_indexes="patid")


# Combine and make current_oha for any OHA

cohort <- cohort %>%
  left_join(current_meds, by="patid") %>%
  left_join(current_insulin, by="patid") %>%
  left_join(current_bolus_insulin, by="patid") %>%
  mutate(across(c("current_dpp4",
                  "current_glinide",
                  "current_gipglp1",
                  "current_glp1",
                  "current_mfn",
                  "current_sglt2",
                  "current_su",
                  "current_tzd",
                  "current_insulin",
                  "current_bolus_insulin"), coalesce, 0L)) %>%
  mutate(current_oha=ifelse(current_dpp4==1 | current_gipglp1==1 | current_glinide==1 | current_glp1==1 | current_mfn==1 | current_sglt2==1 | current_su==1 | current_tzd==1, 1L, 0L)) %>%
  analysis$cached("cohort_interim_7", unique_indexes="patid")

cohort %>% count()
#316,076


# with current insulin for model
cohort %>% group_by(current_insulin) %>% count()

cohort %>% mutate(diabetes_type_new=ifelse(diabetes_type=="type 2" | diabetes_type=="mixed; type 2", "type 2", "type 1")) %>% group_by(diabetes_type_new, current_insulin) %>% count()

# 1 type 2                          0  228466
# 2 type 1                          1   27809
# 3 type 2                          1   58839
# 4 type 1                          0     962


############################################################################################

# Add in treatment at 10 years-post-diagnosis for those registered then

## Do drug classes combined (ignore Acarbose)

ten_yrs_oha <- cohort %>%
  inner_join(clean_oha_prodcodes, by="patid") %>%
  mutate(ten_yrs_post_diag=sql("date_add(diagnosis_date, interval 10 year)")) %>%
  filter(date<=ten_yrs_post_diag & datediff(ten_yrs_post_diag, date)<=366 & drug_class_1!="Acarbose") %>%
  distinct(patid) %>%
  mutate(ten_yrs_post_diag_oha=1L) %>%
  analysis$cached("ten_yrs_oha", unique_indexes="patid")


## Add insulin from insulin scripts (all prodcodes with OHA-insulin mixes are also in insulin)

ten_yrs_insulin <- cohort %>%
  inner_join(clean_insulin_prodcodes, by="patid") %>%
  mutate(ten_yrs_post_diag=sql("date_add(diagnosis_date, interval 10 year)")) %>%
  filter(date<=ten_yrs_post_diag & datediff(ten_yrs_post_diag, date)<=366) %>%
  distinct(patid) %>%
  mutate(ten_yrs_post_diag_ins=1L) %>%
  analysis$cached("ten_yrs_insulin", unique_indexes="patid")


## Bolus insulin

ten_yrs_bolus_insulin <- cohort %>%
  inner_join(clean_insulin_prodcodes, by="patid") %>%
  mutate(ten_yrs_post_diag=sql("date_add(diagnosis_date, interval 10 year)")) %>%
  filter(date<=ten_yrs_post_diag & datediff(ten_yrs_post_diag, date)<=366 & insulin_cat=="Bolus insulin") %>%
  distinct(patid) %>%
  mutate(ten_yrs_post_diag_bolus_ins=1L) %>%
  analysis$cached("ten_yrs_bolus_insulin", unique_indexes="patid")


# Combine

cohort <- cohort %>%
  mutate(ten_yrs_post_diag=sql("date_add(diagnosis_date, interval 10 year)")) %>%
  left_join(ten_yrs_oha, by="patid") %>%
  left_join(ten_yrs_insulin, by="patid") %>%
  left_join(ten_yrs_bolus_insulin, by="patid") %>%
  mutate(across(c("ten_yrs_post_diag_oha",
                  "ten_yrs_post_diag_ins",
                  "ten_yrs_post_diag_bolus_ins"), coalesce, 0L)) %>%
  analysis$cached("cohort_interim_8", unique_indexes="patid")


############################################################################################

# Add in biomarkers

# Clean biomarkers:
## Only keep those within acceptable value limits
## Only keep those with valid unit codes (numunitid)
## If multiple values on the same day, take mean
## Remove those with invalid dates (before DOB or after LCD/death/deregistration)

# Find baseline values
## Use closest date to index date as long as prior to this


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


# For each biomarker, find baseline value at index date
## Use closest date to index as long as prior to this

analysis = cprd$analysis("dpctn_paper")

for (i in biomarkers) {
  
  print(i)
  
  clean_tablename <- paste0("clean_", i, "_medcodes")
  biomarker_index_variable <- paste0("index_",i)
  biomarker_index_date_variable <- paste0("index_",i, "_date")
  biomarker_index_datediff_variable <- paste0("index_",i, "datediff")
  
  index_value <- get(clean_tablename) %>%
    mutate(indexdatediff=datediff(date, index_date)) %>%
    filter(indexdatediff<=0) %>%
    group_by(patid) %>%
    mutate(min_timediff=min(abs(indexdatediff), na.rm=TRUE)) %>%
    filter(abs(indexdatediff)==min_timediff) %>%
    ungroup() %>%
    select(patid, testvalue, date, indexdatediff) %>%
    rename({{biomarker_index_variable}}:=testvalue,
           {{biomarker_index_date_variable}}:=date,
           {{biomarker_index_datediff_variable}}:=indexdatediff)
  
  cohort <- cohort %>%
    left_join(index_value, by="patid")
  
}

cohort <- cohort %>% analysis$cached("cohort_interim_9", unique_indexes="patid")

cohort %>% count()
#316,076


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
  filter(obsdate<=index_date & (is.na(numunitid) | (numunitid!=229 & numunitid!=11589 & numunitid!=218 & numunitid!=276))) %>%
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
  filter(obsdate<=index_date & (is.na(numunitid) | numunitid!=276)) %>%
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


cohort <- cohort %>%
  left_join(gad, by="patid") %>%
  left_join(ia2, by="patid") %>%
  analysis$cached("cohort_interim_10", unique_indexes="patid")


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
  filter(obsdate<=index_date & !is.na(testvalue) & !(c_peptide_cat=="ucpcr" & (numunitid==218 | numunitid==233 | numunitid==235)) & !(c_peptide_cat=="blood" & (numunitid==283 | numunitid==899 | numunitid==277 | numunitid==229))) %>%
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


cohort <- cohort %>%
  left_join(c_pep, by="patid") %>%
  analysis$cached("cohort_interim_11", unique_indexes="patid")


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
  filter(fh_diabetes_cat!="positive - gestational" & date<=index_date) %>%
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

cohort <- cohort %>%
  left_join((fh_positive_negative %>% select(patid, fh_diabetes)), by="patid") %>%
  analysis$cached("cohort_interim_12", unique_indexes="patid")


############################################################################################

# Add earliest insulin and OHA and whether on insulin within 1, 2 or 3 years

## Earliest insulin per patient

earliest_insulin <- clean_insulin_prodcodes %>%
  filter(date<=index_date) %>%
  group_by(patid) %>%
  summarise(earliest_ins=min(date, na.rm=TRUE)) %>%
  ungroup() %>%
  analysis$cached("earliest_insulin", unique_indexes="patid")


## Earliest OHA per patient

earliest_oha <- clean_oha_prodcodes %>%
  filter(date<=index_date) %>%
  group_by(patid) %>%
  summarise(earliest_oha=min(date, na.rm=TRUE)) %>%
  ungroup() %>%
  analysis$cached("earliest_oha", unique_indexes="patid")


## Combine

cohort <- cohort %>%
  left_join(earliest_insulin, by="patid") %>%
  left_join(earliest_oha, by="patid") %>%
  mutate(ins_1_year=ifelse(year(diagnosis_date)<1995, NA,
                           ifelse(!is.na(earliest_ins) & datediff(earliest_ins, diagnosis_date)<=366, 1L,
                                  ifelse(!is.na(earliest_ins) & datediff(earliest_ins, regstartdate)<=366 & regstartdate>diagnosis_date, NA, 0))),
         ins_2_years=ifelse(year(diagnosis_date)<1995, NA,
                            ifelse(!is.na(earliest_ins) & datediff(earliest_ins, diagnosis_date)<=731, 1L,
                                   ifelse(!is.na(earliest_ins) & datediff(earliest_ins, regstartdate)<=366 & datediff(earliest_ins, regstartdate)>366, NA, 0))),
         ins_3_years=ifelse(year(diagnosis_date)<1995, NA,
                   ifelse(!is.na(earliest_ins) & datediff(earliest_ins, diagnosis_date)<=1096, 1L,
                          ifelse(!is.na(earliest_ins) & datediff(earliest_ins, regstartdate)<=366 & datediff(earliest_ins, regstartdate)>731, NA, 0)))) %>%
  analysis$cached("cohort_interim_13", unique_indexes="patid")


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

weight_loss_at_diag <- cohort %>%
  inner_join(clean_weight_loss_medcodes, by="patid") %>%
  filter(year(date)==year(diagnosis_date)) %>%
  distinct(patid) %>%
  mutate(weight_loss_at_diag=1L) %>%
  analysis$cached("weight_loss_at_diag", unique_indexes="patid")

polydipsia_at_diag <- cohort %>%
  inner_join(clean_polydipsia_medcodes, by="patid") %>%
  filter(year(date)==year(diagnosis_date)) %>%
  distinct(patid) %>%
  mutate(polydipsia_at_diag=1L) %>%
  analysis$cached("polydipsia_at_diag", unique_indexes="patid")

urinary_freq_at_diag <- cohort %>%
  inner_join(clean_urinary_frequency_medcodes, by="patid") %>%
  filter(year(date)==year(diagnosis_date)) %>%
  distinct(patid)  %>%
  mutate(urinary_freq_at_diag=1L) %>%
  analysis$cached("urinary_freq_at_diag", unique_indexes="patid")


## Combine

cohort <- cohort %>%
  left_join(weight_loss_at_diag, by="patid") %>%
  left_join(polydipsia_at_diag, by="patid") %>%
  left_join(urinary_freq_at_diag, by="patid") %>%
  mutate(across(c("weight_loss_at_diag",
                  "polydipsia_at_diag",
                  "urinary_freq_at_diag"), coalesce, 0L)) %>%
  analysis$cached("cohort_interim_14", unique_indexes="patid")


#test <- cohort %>% select(diabetes_type, weight_loss_at_diag, polydipsia_at_diag, urinary_freq_at_diag) %>% collect()

#prop.table(table(test$diabetes_type, test$weight_loss_at_diag), margin=1)
#prop.table(table(test$diabetes_type, test$polydipsia_at_diag), margin=1)
#prop.table(table(test$diabetes_type, test$urinary_freq_at_diag), margin=1)
#slightly higher in T1s for all - can have a look once model coded up


############################################################################################

# Add autoimmune conditions

# Autoimmune coding

raw_addisons <- cprd$tables$observation %>%
  filter(medcodeid==12482011000006119 | medcodeid==223491000000112 | medcodeid==3051551000006111 | medcodeid==3051561000006113 | medcodeid==3563081000006116 | medcodeid==3862341000006115 | medcodeid==3862381000006114 | medcodeid==3862401000006114 | medcodeid==3862421000006116 | medcodeid==3866961000006116 | medcodeid==3866981000006114 | medcodeid==3867001000006116 | medcodeid==3978861000006111 | medcodeid==41677012 | medcodeid==417461000006112 | medcodeid==460981000006112 | medcodeid==4741291000006117 | medcodeid==4741311000006118 | medcodeid==4741321000006114 | medcodeid==485624014 | medcodeid==6250371000006117 | medcodeid==682621000006119) %>%
  analysis$cached("raw_addisons", indexes=c("patid", "obsdate"))

raw_ankylosing <- cprd$tables$observation %>%    
    filter(medcodeid==1138031000000113 | medcodeid==11886451000006116 | medcodeid==11886461000006119 | medcodeid==12363861000006114 | medcodeid==12454581000006118 | medcodeid==1574911000006110 | medcodeid==16833013 | medcodeid==1898101000006110 | medcodeid==1909381000006115 | medcodeid==2360741000000115 | medcodeid==253940015 | medcodeid==2653421000006118 | medcodeid==2653431000006115 | medcodeid==3531707015 | medcodeid==3533616015 | medcodeid==3533627018 | medcodeid==359312010 | medcodeid==4554801000006113 | medcodeid==5139501000006110 | medcodeid==7693921000006117 | medcodeid==7709471000006115 | medcodeid==7709481000006117 | medcodeid==7709491000006119) %>%
  analysis$cached("raw_ankylosing", indexes=c("patid", "obsdate"))
    
raw_celiac <- cprd$tables$observation %>%        
    filter(medcodeid==1119251000000112 | medcodeid==131661000006114 | medcodeid==1772803017 | medcodeid==1859331000006118 | medcodeid==303686010 | medcodeid==303687018 | medcodeid==303690012 | medcodeid==3230001000006118 | medcodeid==3230011000006115 | medcodeid==3230031000006114 | medcodeid==3230051000006119 | medcodeid==3230061000006117 | medcodeid==3230071000006112 | medcodeid==3230081000006110 | medcodeid==3230091000006113 | medcodeid==3230101000006119 | medcodeid==3230111000006116 | medcodeid==3230131000006110 | medcodeid==3230141000006117 | medcodeid==411383013 | medcodeid==4787971000006112 | medcodeid==4787991000006113 | medcodeid==5572931000006112 | medcodeid==6567671000006110 | medcodeid==6567681000006113 | medcodeid==6567691000006111 | medcodeid==6567701000006111 | medcodeid==6567721000006118 | medcodeid==6567731000006115 | medcodeid==6567751000006110 | medcodeid==6567761000006112 | medcodeid==6567771000006117 | medcodeid==6567831000006111 | medcodeid==6567851000006116 | medcodeid==7281971000006119) %>%
  analysis$cached("raw_celiac", indexes=c("patid", "obsdate"))
    
raw_graves <- cprd$tables$observation %>%   
    filter(medcodeid==1842051000006119 | medcodeid==2579191000006115 | medcodeid==2579201000006117 | medcodeid==2694980014 | medcodeid==3918261000006117 | medcodeid==472452010 | medcodeid==4757781000006116 | medcodeid==490987015 | medcodeid==5112881000006118 | medcodeid==5581801000006112 | medcodeid==5581811000006110 | medcodeid==5581841000006114 | medcodeid==5581851000006111 | medcodeid==5581861000006113 | medcodeid==57561015 | medcodeid==6196451000006112 | medcodeid==6196461000006114) %>%
  analysis$cached("raw_graves", indexes=c("patid", "obsdate"))
    
raw_hashimoto <- cprd$tables$observation %>%
    filter(medcodeid==2695192017 | medcodeid==2849751000006113 | medcodeid==2849761000006110 | medcodeid==2849831000006114 | medcodeid==355964016 | medcodeid==36888019 | medcodeid==3862341000006115 | medcodeid==3862381000006114 | medcodeid==3862421000006116 | medcodeid==3872531000006114 | medcodeid==3872561000006117 | medcodeid==4197341000006114 | medcodeid==481188016 | medcodeid==5108351000006116 | medcodeid==5108361000006119) %>%
  analysis$cached("raw_hashimoto", indexes=c("patid", "obsdate"))
    
raw_ibd <- cprd$tables$observation %>%
  filter(medcodeid==107644019 | medcodeid==1222351011 | medcodeid==179501000006113 | medcodeid==2269891000000116 | medcodeid==2269901000000115 | medcodeid==2532950019 | medcodeid==2532953017 | medcodeid==2532958014 | medcodeid==2559781000006115 | medcodeid==2559801000006116 | medcodeid==2559821000006114 | medcodeid==2621151000006116 | medcodeid==2621161000006119 | medcodeid==2659581000000113 | medcodeid==2872721013 | medcodeid==2891431000006118 | medcodeid==302322010 | medcodeid==302939018 | medcodeid==302940016 | medcodeid==302941017 | medcodeid==302953018 | medcodeid==302956014 | medcodeid==303761010 | medcodeid==303762015 | medcodeid==3047391000006119 | medcodeid==3047411000006119 | medcodeid==3047421000006110 | medcodeid==309743013 | medcodeid==309744019 | medcodeid==309833017 | medcodeid==309836013 | medcodeid==3113541000006111 | medcodeid==3113551000006113 | medcodeid==3316751000006117 | medcodeid==3316801000006112 | medcodeid==3316811000006110 | medcodeid==3346681000006111 | medcodeid==3346691000006114 | medcodeid==3351341000006118 | medcodeid==3414681000006119 | medcodeid==3414701000006116 | medcodeid==3414711000006118 | medcodeid==3420881000006117 | medcodeid==3420891000006119 | medcodeid==3553391000006113 | medcodeid==396357012 | medcodeid==41137017 | medcodeid==435370011 | medcodeid==4784091000006115 | medcodeid==4784111000006112 | medcodeid==4785581000006112 | medcodeid==4808981000006112 | medcodeid==4809351000006111 | medcodeid==56765016 | medcodeid==601031000006119 | medcodeid==601091000006115 | medcodeid==655231000033115 | medcodeid==6853111000006114 | medcodeid==6853131000006115 | medcodeid==76750016 | medcodeid==85891000006115 | medcodeid==85901000006116 | medcodeid==85931000006112 | medcodeid==886291000006112) %>%
  analysis$cached("raw_ibd", indexes=c("patid", "obsdate"))
    
raw_ms <- cprd$tables$observation %>%
  filter(medcodeid==1682241000006111 | medcodeid==2674605012 | medcodeid==2692565012 | medcodeid==2894411000006110 | medcodeid==297177019 | medcodeid==297179016 | medcodeid==297180018 | medcodeid==297181019 | medcodeid==41398015 | medcodeid==4769281000006110 | medcodeid==4769341000006111 | medcodeid==641211000000118 | medcodeid==641301000000114 | medcodeid==695191000006119 | medcodeid==699671000000115 | medcodeid==699731000000110 | medcodeid==699791000000111 | medcodeid==699911000000117 | medcodeid==7045281000006112 | medcodeid==7058051000006116 | medcodeid==7092351000006112 | medcodeid==983261000006119 | medcodeid==983291000006110 | medcodeid==983301000006111 | medcodeid==983311000006114 | medcodeid==4769321000006116) %>%
  analysis$cached("raw_ms", indexes=c("patid", "obsdate"))

raw_mg <- cprd$tables$observation %>%
  filter(medcodeid==136302010 | medcodeid==151829017 | medcodeid==297580012 | medcodeid==297581011 | medcodeid==297582016 | medcodeid==3988991000006115) %>%
  analysis$cached("raw_mg", indexes=c("patid", "obsdate"))
    
raw_pernicious <- cprd$tables$observation %>%
  filter(medcodeid==11919451000006118 | medcodeid==12717091000006118 | medcodeid==293924010 | medcodeid==293948018 | medcodeid==294568012 | medcodeid==294571016 | medcodeid==294579019 | medcodeid==297018015 | medcodeid==297589013 | medcodeid==302619015 | medcodeid==3300611000006113 | medcodeid==3484001000006116 | medcodeid==3484061000006115 | medcodeid==351106013 | medcodeid==351107016 | medcodeid==3866951000006118 | medcodeid==3866961000006116 | medcodeid==3867001000006116 | medcodeid==3867031000006112 | medcodeid==3875671000006118 | medcodeid==3875681000006115 | medcodeid==399198011 | medcodeid==399201018 | medcodeid==4759041000006118 | medcodeid==4761061000006112 | medcodeid==4761161000006111 | medcodeid==4770761000006113 | medcodeid==505888017 | medcodeid==5063531000006115 | medcodeid==5063541000006113 | medcodeid==5063581000006119 | medcodeid==5063641000006114 | medcodeid==5063651000006111 | medcodeid==62261000006118 | medcodeid==62281000006111 | medcodeid==886021000006118) %>%
  analysis$cached("raw_pernicious", indexes=c("patid", "obsdate"))

raw_pr <- cprd$tables$observation %>%
  filter(medcodeid==108529013 | medcodeid==158050019 | medcodeid==3562021000006111 | medcodeid==359513016) %>%
  analysis$cached("raw_pr", indexes=c("patid", "obsdate"))

raw_pbc <- cprd$tables$observation %>%
  filter(medcodeid==21289010 | medcodeid==303452014 | medcodeid==303453016 | medcodeid==303614019 | medcodeid==303617014 | medcodeid==4044018 | medcodeid==4787861000006113 | medcodeid==52980015) %>%
  analysis$cached("raw_pbc", indexes=c("patid", "obsdate"))

raw_psoriasis <- cprd$tables$observation %>%
  filter(medcodeid==10111000033114 | medcodeid==111418010 | medcodeid==146283019 | medcodeid==1583831000006118 | medcodeid==1629311000006117 | medcodeid==17555019 | medcodeid==1887571000006111 | medcodeid==198791000006112 | medcodeid==198841000006110 | medcodeid==198881000006116 | medcodeid==2162831000033118 | medcodeid==2239001000000113 | medcodeid==243461000006111 | medcodeid==243501000006111 | medcodeid==2555281000006114 | medcodeid==2660221000006116 | medcodeid==2660241000006111 | medcodeid==2660261000006110 | medcodeid==2759071000006116 | medcodeid==2759081000006118 | medcodeid==2759091000006115 | medcodeid==2941391000006114 | medcodeid==2941451000006110 | medcodeid==2941481000006119 | medcodeid==308725015 | medcodeid==308730016 | medcodeid==308731017 | medcodeid==308732012 | medcodeid==308733019 | medcodeid==308734013 | medcodeid==308735014 | medcodeid==308739015 | medcodeid==308740018 | medcodeid==308741019 | medcodeid==308744010 | medcodeid==308745011 | medcodeid==308746012 | medcodeid==308747015 | medcodeid==308749017 | medcodeid==308750017 | medcodeid==308752013 | medcodeid==308753015 | medcodeid==308755010 | medcodeid==308756011 | medcodeid==308757019 | medcodeid==308761013 | medcodeid==308762018 | medcodeid==308769010 | medcodeid==308770011 | medcodeid==308775018 | medcodeid==308776017 | medcodeid==309331014 | medcodeid==309332019 | medcodeid==3096771000006113 | medcodeid==312521010 | medcodeid==3459061000006112 | medcodeid==357612014 | medcodeid==357623013 | medcodeid==359310019 | medcodeid==359321011 | medcodeid==3906401000006119 | medcodeid==3906501000006118 | medcodeid==40067019 | medcodeid==4805761000006111 | medcodeid==4805821000006112 | medcodeid==4805831000006110 | medcodeid==4805871000006113 | medcodeid==4805881000006111 | medcodeid==4805951000006112 | medcodeid==4805981000006116 | medcodeid==4805991000006118 | medcodeid==4806041000006115 | medcodeid==4806051000006118 | medcodeid==483572016 | medcodeid==5123861000006113 | medcodeid==5139461000006110 | medcodeid==5139581000006118 | medcodeid==55627011 | medcodeid==55628018 | medcodeid==61797015) %>%
  analysis$cached("raw_psoriasis", indexes=c("patid", "obsdate"))

raw_ra <- cprd$tables$observation %>%
  filter(medcodeid==116082011 | medcodeid==11903911000006111 | medcodeid==123542016 | medcodeid==125937016 | medcodeid==149911000006117 | medcodeid==162041000006117 | medcodeid==162051000006115 | medcodeid==162081000006111 | medcodeid==162101000006115 | medcodeid==162111000006117 | medcodeid==162131000006111 | medcodeid==162141000006118 | medcodeid==162191000006110 | medcodeid==162231000006117 | medcodeid==162301000006117 | medcodeid==162341000006115 | medcodeid==168751000006115 | medcodeid==168761000006118 | medcodeid==16911000006114 | medcodeid==1778239014 | medcodeid==1779323014 | medcodeid==1786505011 | medcodeid==1786545019 | medcodeid==2472447010 | medcodeid==255925015 | medcodeid==2653371000006113 | medcodeid==2670341000006112 | medcodeid==2839291014 | medcodeid==2842168015 | medcodeid==297544012 | medcodeid==297641010 | medcodeid==300221014 | medcodeid==3042571000006112 | medcodeid==3042581000006110 | medcodeid==309787016 | medcodeid==309788014 | medcodeid==309789018 | medcodeid==309790010 | medcodeid==309791014 | medcodeid==309792019 | medcodeid==309794018 | medcodeid==309798015 | medcodeid==309800010 | medcodeid==309802019 | medcodeid==309803012 | medcodeid==309804018 | medcodeid==309805017 | medcodeid==309812014 | medcodeid==309816012 | medcodeid==309824019 | medcodeid==309827014 | medcodeid==309828016 | medcodeid==311496011 | medcodeid==312518013 | medcodeid==312520011 | medcodeid==312534011 | medcodeid==3428891000006116 | medcodeid==3428901000006117 | medcodeid==3428911000006119 | medcodeid==359292012 | medcodeid==3636541000006112 | medcodeid==3636551000006114 | medcodeid==3636561000006111 | medcodeid==3709251000006117 | medcodeid==371261000000111 | medcodeid==3732891000006112 | medcodeid==424981000006114 | medcodeid==426510015 | medcodeid==451461014 | medcodeid==4580651000006110 | medcodeid==4809101000006114 | medcodeid==4809181000006117 | medcodeid==4809221000006114 | medcodeid==4809301000006112 | medcodeid==4809341000006114 | medcodeid==48361011 | medcodeid==485614018 | medcodeid==5721431000006114 | medcodeid==6608771000006110 | medcodeid==6608781000006113 | medcodeid==6608791000006111 | medcodeid==6802581000006114 | medcodeid==6802661000006115 | medcodeid==7099641000006113 | medcodeid==755441000006119 | medcodeid==889731000006111 | medcodeid==95067016) %>%
  analysis$cached("raw_ra", indexes=c("patid", "obsdate"))

raw_sjogren <- cprd$tables$observation %>%
  filter(medcodeid==133918015 | medcodeid==141091000006112 | medcodeid==222101000000116 | medcodeid==297645018 | medcodeid==301720019 | medcodeid==3813051000006112 | medcodeid==3865031000006117 | medcodeid==3865051000006112 | medcodeid==3865061000006114 | medcodeid==3865071000006119 | medcodeid==3865081000006116 | medcodeid==4193701000006119 | medcodeid==4193711000006116 | medcodeid==444889012 | medcodeid==4770971000006116 | medcodeid==4770991000006115 | medcodeid==4771001000006119 | medcodeid==4782461000006113 | medcodeid==4782481000006115 | medcodeid==5912751000006118) %>%
  analysis$cached("raw_sjogren", indexes=c("patid", "obsdate"))

raw_lupus <- cprd$tables$observation %>%
  filter(medcodeid==114251000006118 | medcodeid==114310015 | medcodeid==158372014 | medcodeid==158442019 | medcodeid==1728151000006118 | medcodeid==177301000006114 | medcodeid==25612013 | medcodeid==2929831000006119 | medcodeid==2929841000006112 | medcodeid==2929851000006114 | medcodeid==297542011 | medcodeid==297639014 | medcodeid==301721015 | medcodeid==308697013 | medcodeid==308700012 | medcodeid==308701011 | medcodeid==308704015 | medcodeid==308705019 | medcodeid==308706018 | medcodeid==308707010 | medcodeid==308709013 | medcodeid==309383019 | medcodeid==309405011 | medcodeid==309408013 | medcodeid==3111201000006111 | medcodeid==3111211000006114 | medcodeid==312579010 | medcodeid==3174611000006115 | medcodeid==3400601000006117 | medcodeid==359447015 | medcodeid==3619221000006118 | medcodeid==3619241000006113 | medcodeid==3789491000006111 | medcodeid==3789501000006115 | medcodeid==3862331000006113 | medcodeid==4060291000006118 | medcodeid==44918014 | medcodeid==453243010 | medcodeid==4597351000006115 | medcodeid==4597361000006118 | medcodeid==4597381000006111 | medcodeid==4597391000006114 | medcodeid==4609091000006113 | medcodeid==4805591000006111 | medcodeid==4805631000006111 | medcodeid==4805641000006118 | medcodeid==5140921000006119 | medcodeid==5140931000006116 | medcodeid==660871000006118 | medcodeid==677561000006119 | medcodeid==6909361000006110 | medcodeid==92208011 | medcodeid==92209015) %>%
  analysis$cached("raw_lupus", indexes=c("patid", "obsdate"))

raw_ss <- cprd$tables$observation %>%
  filter(medcodeid==12719701000006115 | medcodeid==147833014 | medcodeid==147834015 | medcodeid==297644019 | medcodeid==3012931000006114 | medcodeid==3012941000006116 | medcodeid==3012951000006119 | medcodeid==3012961000006117 | medcodeid==3012971000006112 | medcodeid==3012981000006110 | medcodeid==3012991000006113 | medcodeid==301716018 | medcodeid==308855016 | medcodeid==308859010 | medcodeid==308860017 | medcodeid==308871016 | medcodeid==309413012 | medcodeid==312581012 | medcodeid==316770014 | medcodeid==354512012 | medcodeid==354513019 | medcodeid==359440018 | medcodeid==38247018 | medcodeid==3949101000006111 | medcodeid==3949121000006118 | medcodeid==3949131000006115 | medcodeid==4770951000006114 | medcodeid==4782401000006112 | medcodeid==4782411000006110 | medcodeid==4806451000006111 | medcodeid==4806471000006118 | medcodeid==4806491000006117 | medcodeid==4806501000006113 | medcodeid==4806521000006115 | medcodeid==53211014 | medcodeid==5850801000006113 | medcodeid==5850811000006111 | medcodeid==709501000006112 | medcodeid==7288571000006116) %>%
  analysis$cached("raw_ss", indexes=c("patid", "obsdate"))


## Query is too long so do in two parts

raw_vasculitis_a <- cprd$tables$observation %>%
  filter(medcodeid==100606014 | medcodeid==101049016 | medcodeid==104572011 | medcodeid==11878901000006118 | medcodeid==11878911000006115 | medcodeid==124653012 | medcodeid==12704531000006115 | medcodeid==12726911000006118 | medcodeid==165061000006114 | medcodeid==1779323014 | medcodeid==1779388019 | medcodeid==1783725016 | medcodeid==1787109017 | medcodeid==1787110010 | medcodeid==189131000006116 | medcodeid==189141000006114 | medcodeid==189181000006115 | medcodeid==191841000006118 | medcodeid==191861000006119 | medcodeid==191881000006112 | medcodeid==2157271000000115 | medcodeid==2157584018 | medcodeid==2214081000000112 | medcodeid==242065019 | medcodeid==2534171012 | medcodeid==2557891000006117 | medcodeid==2658281000006111 | medcodeid==266831000006117 | medcodeid==2743591000006110 | medcodeid==2753291000006114 | medcodeid==2807131000006118 | medcodeid==2832821000006115 | medcodeid==2835591000006113 | medcodeid==286395017 | medcodeid==293373010 | medcodeid==294337012 | medcodeid==294341011 | medcodeid==2955624019 | medcodeid==297543018 | medcodeid==297640011 | medcodeid==297777014 | medcodeid==2998271000006119 | medcodeid==300503012 | medcodeid==300565019 | medcodeid==300577012 | medcodeid==300579010 | medcodeid==300584016 | medcodeid==300593015 | medcodeid==300596011 | medcodeid==3008151000006116 | medcodeid==300945017 | medcodeid==3015291000006114 | medcodeid==302199015 | medcodeid==303812017 | medcodeid==303815015 | medcodeid==303816019 | medcodeid==303820015 | medcodeid==30382015 | medcodeid==303824012 | medcodeid==303825013 | medcodeid==303826014 | medcodeid==303834015 | medcodeid==303839013 | medcodeid==303849011 | medcodeid==303863017 | medcodeid==303864011 | medcodeid==303866013 | medcodeid==303867016 | medcodeid==303868014 | medcodeid==303869018 | medcodeid==303870017 | medcodeid==303877019 | medcodeid==303881019 | medcodeid==303883016 | medcodeid==303907019 | medcodeid==304769016 | medcodeid==3082891000006115 | medcodeid==309386010 | medcodeid==309393014 | medcodeid==309483018 | medcodeid==309484012 | medcodeid==309485013 | medcodeid==309486014 | medcodeid==309488010 | medcodeid==309489019 | medcodeid==309490011 | medcodeid==309491010 | medcodeid==309492015 | medcodeid==312575016 | medcodeid==312576015 | medcodeid==3222361000006116 | medcodeid==3222401000006114 | medcodeid==3222411000006112 | medcodeid==32619017 | medcodeid==3306421000006118 | medcodeid==3319041000006114 | medcodeid==3319071000006118 | medcodeid==3319081000006115 | medcodeid==3319091000006117 | medcodeid==3319101000006111 | medcodeid==3319111000006114 | medcodeid==3319121000006118 | medcodeid==3319131000006115 | medcodeid==3319141000006113 | medcodeid==3377931000006119 | medcodeid==3480841000006118 | medcodeid==3483561000006111 | medcodeid==3483571000006116 | medcodeid==3483581000006118 | medcodeid==3483611000006114 | medcodeid==3483621000006118 | medcodeid==3483631000006115 | medcodeid==35171019 | medcodeid==352258013 | medcodeid==3522791000006117 | medcodeid==3522801000006116 | medcodeid==354346014 | medcodeid==354347017 | medcodeid==354348010) %>%
  analysis$cached("raw_vasculitis_a", indexes=c("patid", "obsdate"))

raw_vasculitis_b <- cprd$tables$observation %>%
  filter(medcodeid==354349019 | medcodeid==354352010 | medcodeid==354363013 | medcodeid==354364019 | medcodeid==354380012 | medcodeid==3544341000006115 | medcodeid==354514013 | medcodeid==354516010 | medcodeid==354520014 | medcodeid==354632016 | medcodeid==354636018 | medcodeid==357845012 | medcodeid==359513016 | medcodeid==3619221000006118 | medcodeid==3619241000006113 | medcodeid==3620071000006112 | medcodeid==3620091000006113 | medcodeid==370681000006116 | medcodeid==3720161000006112 | medcodeid==3720191000006116 | medcodeid==3720201000006118 | medcodeid==3720211000006115 | medcodeid==3755261000006110 | medcodeid==3806541000006112 | medcodeid==3806551000006114 | medcodeid==3806561000006111 | medcodeid==3838311000006114 | medcodeid==3838321000006118 | medcodeid==3838341000006113 | medcodeid==3838351000006110 | medcodeid==3838361000006112 | medcodeid==395799018 | medcodeid==396706011 | medcodeid==400165012 | medcodeid==430091000006117 | medcodeid==430101000006111 | medcodeid==430591000006111 | medcodeid==440031000006117 | medcodeid==4514891000006110 | medcodeid==4514901000006114 | medcodeid==454287012 | medcodeid==458441000006117 | medcodeid==458451000006115 | medcodeid==458521000006111 | medcodeid==47484015 | medcodeid==475138012 | medcodeid==475140019 | medcodeid==4759971000006113 | medcodeid==4762481000006118 | medcodeid==4762491000006115 | medcodeid==4762511000006114 | medcodeid==4769861000006118 | medcodeid==4769911000006111 | medcodeid==4770431000006110 | medcodeid==4772681000006110 | medcodeid==4779011000006115 | medcodeid==4779031000006114 | medcodeid==4779041000006116 | medcodeid==4788161000006115 | medcodeid==4788191000006111 | medcodeid==4788491000006115 | medcodeid==4788511000006114 | medcodeid==4791861000006115 | medcodeid==4791871000006110 | medcodeid==479396010 | medcodeid==4808051000006117 | medcodeid==4808071000006110 | medcodeid==4808091000006111 | medcodeid==4808101000006117 | medcodeid==4808121000006110 | medcodeid==4808141000006115 | medcodeid==492911000006119 | medcodeid==499415017 | medcodeid==5074381000006118 | medcodeid==508654012 | medcodeid==5093471000006118 | medcodeid==5093551000006112 | medcodeid==5093581000006116 | medcodeid==5093731000006110 | medcodeid==5095121000006112 | medcodeid==5095151000006115 | medcodeid==5095181000006111 | medcodeid==5095191000006114 | medcodeid==5096311000006113 | medcodeid==5141541000006114 | medcodeid==5141551000006111 | medcodeid==5141561000006113 | medcodeid==5141571000006118 | medcodeid==5141661000006114 | medcodeid==53457013 | medcodeid==553181000006118 | medcodeid==553191000006115 | medcodeid==553201000006117 | medcodeid==554481000006111 | medcodeid==555191000006119 | medcodeid==555201000006116 | medcodeid==555371000006118 | medcodeid==555381000006115 | medcodeid==556211000006112 | medcodeid==557151000006115 | medcodeid==6002791000006115 | medcodeid==6002811000006116 | medcodeid==6002821000006112 | medcodeid==601821000033114 | medcodeid==6206401000006113 | medcodeid==6206421000006115) %>%
  analysis$cached("raw_vasculitis_b", indexes=c("patid", "obsdate"))

raw_vasculitis_c <- cprd$tables$observation %>%
  filter(medcodeid==6206431000006117 | medcodeid==6206451000006112 | medcodeid==6206471000006119 | medcodeid==6206521000006116 | medcodeid==6206541000006111 | medcodeid==6206551000006113 | medcodeid==6206561000006110 | medcodeid==6206581000006117 | medcodeid==6206591000006119 | medcodeid==6381141000006118 | medcodeid==65971000006112 | medcodeid==6621301000006118 | medcodeid==6621311000006115 | medcodeid==6635251000006112 | medcodeid==6635261000006114 | medcodeid==666291000006116 | medcodeid==6718781000006115 | medcodeid==677461000006115 | medcodeid==677471000006110 | medcodeid==677481000006113 | medcodeid==677621000006112 | medcodeid==677971000006118 | medcodeid==677981000006115 | medcodeid==7082541000006119 | medcodeid==763111000006117 | medcodeid==7724871000006110 | medcodeid==7724911000006113 | medcodeid==7724951000006114 | medcodeid==7724991000006115 | medcodeid==7725071000006116 | medcodeid==789831000006117 | medcodeid==81251000006113 | medcodeid==81261000006110 | medcodeid==81581000006119 | medcodeid==828201000006113 | medcodeid==8412881000006110 | medcodeid==84262013 | medcodeid==87208017 | medcodeid==884681000006117 | medcodeid==91895012 | medcodeid==990621000006115) %>%
  analysis$cached("raw_vasculitis_c", indexes=c("patid", "obsdate"))

raw_vasculitis <- raw_vasculitis_a %>%
  union_all(raw_vasculitis_b) %>%
  union_all(raw_vasculitis_c) %>%
  analysis$cached("raw_vasculitis", indexes=c("patid", "obsdate"))

raw_vitiligo <- cprd$tables$observation %>%
  filter(medcodeid==145353016 | medcodeid==94346015) %>%
  analysis$cached("raw_vitiligo", indexes=c("patid", "obsdate"))


autoimmune <- c("addisons", "ankylosing", "celiac", "graves", "hashimoto", "ibd", "ms", "mg", "pernicious", "pr", "pbc", "psoriasis", "ra", "sjogren", "lupus", "ss", "vasculitis", "vitiligo")


for (i in autoimmune) {
  
  raw_tablename <- paste0("raw_", i)
  clean_tablename <- paste0("clean_", i)
  
  clean_table <- get(raw_tablename) %>%
    select(patid, date=obsdate) %>%
    inner_join(cprd$tables$validDateLookup, by="patid") %>%
    filter(date>=min_dob & date<=index_date) %>%
    distinct(patid) %>%
    mutate({{i}}:=1L) %>%
    analysis$cached(clean_tablename, unique_index="patid")
  
  assign(clean_tablename, clean_table)
  
}

# Add individual autoimmune conditions

for (i in autoimmune) {

  clean_tablename <- paste0("clean_", i)
  
  cohort <- cohort %>%
    left_join(get(clean_tablename), by="patid") %>%
    mutate({{i}}:=ifelse(is.na(!!sym(i)), 0L, 1L))
  
}


# Add overall, plus coeliac/thyroid

cohort <- cohort %>%
  mutate(any_autoimmune=ifelse(addisons==1 | ankylosing==1 | celiac==1 | graves==1 | hashimoto==1 | ibd==1 | ms==1 | mg==1 | pernicious==1 | pr==1 | pbc==1 | psoriasis==1 | ra==1 | sjogren==1 | lupus==1 | ss==1 | vasculitis==1 | vitiligo==1, 1L, 0L),
         coeliac_thyroid=ifelse(celiac==1 | graves==1 | hashimoto==1, 1L, 0L)) %>%
  analysis$cached("cohort_interim_15", unique_indexes="patid")


############################################################################################

# Add hospitalisation outcomes once have linkage data

# # Add in DKA at diagnosis, earliest post-diagnosis DKA and hypo, and no. of DKA and hypo postdiag visits
# 
# dka_at_diagnosis <- cprd$tables$hesDiagnosisEpi %>%
#   inner_join(codes$icd10_dka, by=c("ICD"="icd10")) %>%
#   inner_join((cohort %>% select(patid, diagnosis_date)), by="patid") %>%
#   filter(d_order==1 & epistart<=index_date & abs(datediff(epistart, diagnosis_date)<=30)) %>%
#   distinct(patid) %>%
#   mutate(dka_at_diagnosis=1L) %>%
#   analysis$cached("dka_at_diagnosis", unique_indexes="patid")
# 
# earliest_dka <- cprd$tables$hesDiagnosisEpi %>%
#   inner_join(codes$icd10_dka, by=c("ICD"="icd10")) %>%
#   inner_join((cohort %>% select(patid, diagnosis_date)), by="patid") %>%
#   filter(d_order==1 & epistart<=index_date & datediff(epistart, diagnosis_date)>30) %>%
#   group_by(patid) %>%
#   summarise(earliest_dka=min(epistart, na.rm=TRUE),
#             dka_count=n()) %>%
#   ungroup() %>%
#   analysis$cached("earliest_dka", unique_indexes="patid")
# 
# earliest_hypo <- cprd$tables$hesDiagnosisEpi %>%
#   inner_join(codes$icd10_hypoglycaemia, by=c("ICD"="icd10")) %>%
#   inner_join((cohort %>% select(patid, diagnosis_date)), by="patid") %>%
#   filter(d_order==1 & epistart<=index_date & datediff(epistart, diagnosis_date)>30) %>%
#   group_by(patid) %>%
#   summarise(earliest_hypo=min(epistart, na.rm=TRUE),
#             hypo_count=n()) %>%
#   ungroup() %>%
#   analysis$cached("earliest_hypo", unique_indexes="patid")
# 
# 
# ## No. of contacts in last year (just using unique patid-date interactions as per T1 paper I reviewed)
# 
# 
one_year_prior <- index_date-years(1)

contact_count <- cprd$tables$observation %>%
  filter(obsdate<=index_date & obsdate>one_year_prior) %>%
  distinct(patid, obsdate) %>%
  group_by(patid) %>%
  summarise(contact_count=n()) %>%
  ungroup() %>%
  analysis$cached("contact_count", unique_indexes="patid")
#
# 
# 
# # Add all to cohort
# 
cohort <- cohort %>%
  #left_join(dka_at_diagnosis, by="patid") %>%
  #mutate(dka_at_diagnosis=ifelse(is.na(dka_at_diagnosis), 0L, 1L)) %>%
  #left_join(earliest_dka, by="patid") %>%
  #left_join(earliest_hypo, by="patid") %>%
  left_join(contact_count, by="patid") %>%
  #mutate(dka_post_diagnosis=ifelse(!is.na(earliest_dka), 1L, 0L),
  #       hypo_post_diagnosis=ifelse(!is.na(earliest_hypo), 1L, 0L)) %>%
  analysis$cached("cohort_interim_16", unique_indexes="patid")
# 
# 
# 
# ## Check no. with DKA/hypos before diagnosis - exclude?
# ## Previously when including pre-diagnosis
# #test <- cohort %>% filter(earliest_hosp_hypo<diagnosis_date | earliest_hosp_dka<diagnosis_date) %>% select(patid, diagnosis_date, earliest_hosp_hypo, earliest_hosp_dka) %>% collect()
# #8 people - all hypos - fine


############################################################################################

# Code up lipid model

cohort <- cohort %>%
  
  mutate(age_at_bmi=datediff(index_bmi_date, dob)/365.25,
         bmi_post_diag=ifelse(index_bmi_date>=diagnosis_date & age_at_bmi>=18, index_bmi, NA),
         bmi_post_diag_datediff=ifelse(!is.na(bmi_post_diag), index_bmidatediff, NA),
         totalchol_post_diag=ifelse(index_totalcholesterol_date>=diagnosis_date, index_totalcholesterol, NA),
         totalchol_post_diag_datediff=ifelse(!is.na(totalchol_post_diag), index_totalcholesteroldatediff, NA),
         hdl_post_diag=ifelse(index_hdl_date>=diagnosis_date, index_hdl, NA),
         hdl_post_diag_datediff=ifelse(!is.na(hdl_post_diag), index_hdldatediff, NA),
         triglyceride_post_diag=ifelse(index_triglyceride_date>=diagnosis_date, index_triglyceride, NA),
         triglyceride_post_diag_datediff=ifelse(!is.na(triglyceride_post_diag), index_triglyceridedatediff, NA),
         diabetes_type_new=ifelse(diabetes_type=="type 1" | diabetes_type=="mixed; type 1", "type 1", "type 2"),
         femalesex=ifelse(gender==2, 1, ifelse(gender==1, 0, NA)),
         
         lipid_pred_score=9.0034272-(0.1915482*bmi_post_diag)-(0.1686227*dm_diag_age)+(0.3026012*femalesex)-(0.2269216*totalchol_post_diag)+(1.540850*hdl_post_diag)-(0.2784059*triglyceride_post_diag),
         lipid_pred_prob=exp(lipid_pred_score)/(1+exp(lipid_pred_score)),
         
         mixed_codes=ifelse(diabetes_type!=diabetes_type_new, 1, 0)) %>%
  
  analysis$cached("cohort", unique_index="patid")


cohort %>% count()
#316,076

