
# Define prevalent cohort of insulin-treated T1 and T2 diagnosed aged 18+ years:

## All IDs in download who are registered on 01/03/2023 (so have HES data)
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
library(table1)
rm(list=ls())

cprd=CPRDData$new(cprdEnv="diabetes-jun2024", cprdConf="~/.aurum.yaml")

analysis=cprd$analysis("dpctn_paper_prev")


############################################################################################

# Set index date

index_date <- as.Date("2023-03-01")


# Point to patient variables coded up in previous script

analysis=cprd$analysis("dpctn_paper")
diagnosis_dates <- diagnosis_dates %>% analysis$cached("diagnosis_dates")
diabetes_type <- diabetes_type %>% analysis$cached("diabetes_type")
diabetes_meds <- diabetes_meds %>% analysis$cached("diabetes_meds")
biomarkers <- biomarkers %>% analysis$cached("biomarkers")
gad <- gad %>% analysis$cached("gad")
ia2 <- ia2 %>% analysis$cached("ia2")
c_pep <- c_pep %>% analysis$cached("c_pep")
fh_positive_negative <- fh_positive_negative %>% analysis$cached("fh_positive_negative")
at_diag_features <- at_diag_features %>% analysis$cached("at_diag_features")
autoimmune <- autoimmune %>% analysis$cached("autoimmune")
dka_hypo_contacts <- dka_hypo_contacts %>% analysis$cached("dka_hypo_contacts")


############################################################################################

# Find IDs of those registered on index date and pull in static patient features

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
## Calculate age at index date

analysis=cprd$analysis("diabetes_cohort")
dob <- dob %>% analysis$cached("dob")

analysis=cprd$analysis("all_patid")
ethnicity <- ethnicity %>% analysis$cached("ethnicity")


analysis=cprd$analysis("dpctn_paper_prev")

cohort <- cprd$tables$patient %>%
  select(patid, gender, pracid, regstartdate, regenddate) %>%
  left_join((dob %>% select(patid, dob)), by="patid") %>%
  left_join((cprd$tables$practice %>% select(pracid, prac_region=region, lcd)), by="pracid") %>%
  left_join((ethnicity %>% select(patid, ethnicity_5cat, ethnicity_16cat)), by="patid") %>%
  left_join((cprd$tables$onsDeath %>% select(patid, death_date=reg_date_of_death)), by="patid") %>%
  left_join((cprd$tables$patidsWithLinkage %>% mutate(with_hes=1L)), by="patid") %>%
  left_join((cprd$tables$patientImd %>% select(-pracid)), by="patid") %>%
  
  mutate(gp_record_end=pmin(if_elsedk(is.na(lcd), as.Date("2024-06-30"), lcd),
                            if_else(is.na(regenddate), as.Date("2024-06-30"), regenddate),
                            if_else(is.na(death_date), as.Date("2024-06-30"), death_date),
                            as.Date("2024-06-30"), na.rm=TRUE),
         
         age_at_index=round(datediff(index_date, dob)/365.25, 1)) %>%
  
  select(patid, gender, dob, age_at_index, pracid, prac_region, ethnicity_5cat, ethnicity_16cat, regstartdate, gp_record_end, death_date, with_hes, imd_decile) %>%
  
  filter(with_hes==1 & regstartdate<=index_date & !(!is.na(death_date) & death_date<index_date) & !(!is.na(gp_record_end) & gp_record_end<index_date)) %>%
  
  analysis$cached("cohort_interim_1", unique_indexes="patid")

cohort %>% count()
# 1,076,165

# Remove patients from practices which may have duplicated data as per CPRD's guidance

analysis=cprd$analysis("diabetes_cohort")

practice_exclusion_ids <- cprd$tables$patient %>% 
  filter(pracid=="20024" | pracid=="20036" |pracid=="20091" |pracid=="20171" | pracid=="20178" |pracid=="20202" | pracid=="20254" | pracid=="20389" |pracid=="20430" |pracid=="20452" |
           pracid=="20469" | pracid=="20487" | pracid=="20552" | pracid=="20554" | pracid=="20640" | pracid=="20717" | pracid=="20734" | pracid=="20737" | pracid=="20740" | pracid=="20790" |
           pracid=="20803" | pracid=="20822" | pracid=="20868" | pracid=="20912" | pracid=="20996" | pracid=="21001" | pracid=="21015" | pracid=="21078" | pracid=="21112" | pracid=="21118" |
           pracid=="21172" | pracid=="21173" | pracid=="21277" | pracid=="21281" | pracid=="21331" | pracid=="21334" | pracid=="21390" | pracid=="21430" | pracid=="21444" | pracid=="21451" |
           pracid=="21529" | pracid=="21553" | pracid=="21558" | pracid=="21585") %>%
  distinct(patid) %>%
  analysis$cached("practice_exclusion_ids")


analysis=cprd$analysis("dpctn_paper_prev")

cohort <- cohort %>%
  anti_join(practice_exclusion_ids, by="patid") %>%
  analysis$cached("cohort_interim_2", unique_indexes="patid")

cohort %>% count()
# 1,076,165



# Just keep men and women

cohort <- cohort %>%
  filter(gender==1 | gender==2) %>%
  analysis$cached("cohort_interim_3", unique_indexes="patid")

cohort %>% count()
# 1,076,147


############################################################################################

# Add 'current' diabetes type and exclude if not type 1 or type 2

cohort <- cohort %>%
  inner_join(diabetes_type, by="patid") %>%
  mutate(t1_code_ever=ifelse(type1_code_count>0, 1L, 0L),
         t2_code_ever=ifelse(type2_code_count>0, 1L, 0L),
         diabetes_type=str_replace_all(diabetes_type, "mixed; ", "")) %>%
  analysis$cached("cohort_interim_4", unique_indexes="patid")

cohort %>% count()
# 984,166
## lose people diagnosed after index date


# Look at breakdown of diabetes type

# clipr::write_clip(cohort %>%
#   filter(!(diabetes_type=="type 1" | diabetes_type=="type 2")) %>%
#   group_by(diabetes_type) %>%
#   summarise(count=n()) %>%
#   collect())

#n=215,134 unspecified
#n=54,619 gestational or mixed; gestational
#n=1499 secondary or mixed; secondary
#n=1579 other


cohort <- cohort %>%
  filter(diabetes_type=="type 1" | diabetes_type=="type 2") %>%
  analysis$cached("cohort_interim_5", unique_indexes="patid")

cohort %>% count()
#711,335

cohort %>% group_by(diabetes_type) %>% count()
# 1 type 2         661273
# 2 type 1          50062


############################################################################################

# Add in treatment and exclude those not currently on insulin

cohort <- cohort %>%
  inner_join((diabetes_meds %>% select(-c(regstartdate, diagnosis_date))), by="patid") %>%
  filter(current_insulin==1) %>%
  analysis$cached("cohort_interim_6", unique_indexes="patid")

cohort %>% count()
#136,127

cohort %>% group_by(diabetes_type) %>% count()
# 1 type 2          87762
# 2 type 1          48365


# Diagnosis code in YOB
cohort %>% inner_join(diagnosis_dates, by="patid") %>% filter(year(diagnosis_date)==year(dob)) %>% count()
#643

test <- cohort %>% inner_join(diagnosis_dates, by="patid") %>% filter(year(diagnosis_date)==year(dob)) %>% select(patid, diagnosis_date, earliest_oha, earliest_ins, dob) %>% collect()
#643



############################################################################################

# Add diagnosis dates and exclude if aged<18 or >50 years

cohort <- cohort %>%
  inner_join(diagnosis_dates, by="patid") %>%
  mutate(dm_diag_age=round((datediff(diagnosis_date, dob))/365.25, 1),
         follow_up_time=datediff(index_date, diagnosis_date)/365.25) %>%
  filter(dm_diag_age>=18 & dm_diag_age<51) %>%
  analysis$cached("cohort_interim_7", unique_indexes="patid")

cohort %>% count()
#70,348

cohort %>% group_by(diabetes_type) %>% count()
# 1 type 2          47612
# 2 type 1          22736

87762-47612
48365-22736


############################################################################################

# Data quality checking: remove if close to registration or type 1 and not on insulin within 5 years

cohort %>% filter(datediff(diagnosis_date, regstartdate)>= -30 & datediff(diagnosis_date, regstartdate)<=90) %>% count()
#3859
cohort %>% filter(datediff(diagnosis_date, regstartdate)>= -30 & datediff(diagnosis_date, regstartdate)<=90) %>% group_by(diabetes_type) %>% count()
# 1 type 1           1540
# 2 type 2           2319


cohort %>% mutate(ins_5_years=ifelse(year(diagnosis_date)<2002 | datediff(index_date, diagnosis_date)<=1827 | datediff(regstartdate, diagnosis_date)>=1461, NA,
                                     ifelse(!is.na(earliest_ins) & datediff(earliest_ins, diagnosis_date)<=1827, 1L, 0L))) %>%
  filter(diabetes_type=="type 1" & !is.na(ins_5_years) & ins_5_years==0) %>%
  count()
#239  
  
  
cohort <- cohort %>%
  filter((datediff(diagnosis_date, regstartdate)< -30 | datediff(diagnosis_date, regstartdate)>90)) %>%
  mutate(ins_5_years=ifelse(year(diagnosis_date)<2002 | datediff(index_date, diagnosis_date)<=1827 | datediff(regstartdate, diagnosis_date)>=1461, NA,
                            ifelse(!is.na(earliest_ins) & datediff(earliest_ins, diagnosis_date)<=1827, 1L, 0L))) %>%
  filter(!(diabetes_type=="type 1" & !is.na(ins_5_years) & ins_5_years==0)) %>%
  analysis$cached("cohort_interim_8", unique_indexes="patid")

cohort %>% count()
#66269

cohort %>% group_by(diabetes_type) %>% count()
# 1 type 2         45293
# 2 type 1          20976

70348-66269


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
  left_join(dka_hypo_contacts, by="patid") %>%
  analysis$cached("cohort_interim_9", unique_indexes="patid")


## Check missing variables for model

cohort %>% filter(is.na(index_bmi)) %>% count()
#798

cohort %>% filter(is.na(index_totalcholesterol)) %>% count()
#896

cohort %>% filter(is.na(index_hdl)) %>% count()
#1083

cohort %>% filter(is.na(index_triglyceride)) %>% count()
#3938


cohort %>% filter(is.na(index_bmi) | is.na(index_totalcholesterol) | is.na(index_hdl) | is.na(index_triglyceride)) %>% count()
#4,256

4256/66269


# Look at missingness by sex, ethnicity, IMD

local_data <- cohort %>%
  select(index_bmi, index_totalcholesterol, index_hdl, index_triglyceride, gender, imd_decile, ethnicity_5cat, diabetes_type, pracid) %>%
  collect() %>%
  mutate(missing=as.factor(ifelse(!(is.na(index_bmi) | is.na(index_totalcholesterol) | is.na(index_hdl) | is.na(index_triglyceride)), 0, 1)),
         imd_quintiles=factor(case_when(imd_decile==1 | imd_decile==2 ~1,
                                        imd_decile==3 | imd_decile==4 ~2,
                                        imd_decile==5 | imd_decile==6 ~3,
                                        imd_decile==7 | imd_decile==8 ~4,
                                        imd_decile==9 | imd_decile==10 ~5)),
         malesex=as.factor(ifelse(gender==1, 1, 0)),
         ethnicity_decoded=factor(case_when(ethnicity_5cat==0 ~ "White",
                                            ethnicity_5cat==1 ~ "South Asian",
                                            ethnicity_5cat==2 ~ "Black",
                                            ethnicity_5cat==3 ~"Other",
                                            ethnicity_5cat==4 ~"Mixed",
                                            is.na(ethnicity_5cat) ~"Missing"), levels=c("White", "South Asian", "Black", "Mixed", "Other", "Missing")),
         diabetes_type=factor(diabetes_type))

prop_with_ci_wilson <- function(x, digits=3) {
  tab <- table(x, useNA="ifany")
  n <- sum(tab)
  if (n==0) stop("No non-missing values in x")
  
  z  <- qnorm(1 - (1 - 0.95) / 2)
  z2 <- z^2
  
  res <- vapply(seq_along(tab), function(i) {
    k <- as.integer(tab[i])
    p <- k / n
    
    # Wilson (score) interval
    denom  <- 1 + z2 / n
    center <- (p + z2 / (2 * n)) / denom
    half   <- (z * sqrt((p * (1 - p) / n) + (z2 / (4 * n^2)))) / denom
    lo <- max(0, center - half)
    hi <- min(1, center + half)
    
    sprintf(
      "%s (%.*f%% [%.*f-%.*f%%])",
      prettyNum(k, big.mark=","),
      digits, 100 * p,
      digits, 100 * lo,
      digits, 100 * hi
    )
  }, FUN.VALUE=character(1))
  
  names(res) <- names(tab)
  res
}

render.categorical.ci <- function(x, ...) {
  prop_with_ci_wilson(x, ...)
}

render_cat_ci <- function(digits=3) {
  force(digits)
  function(x, ...) render.categorical.ci(x, digits=digits)
}

new_cohort <- rbind((local_data %>% mutate(new_group="overall")),
                    (local_data %>% filter(malesex==1) %>% mutate(new_group="male")),
                    (local_data %>% filter(malesex==0) %>% mutate(new_group="female")),
                    (local_data %>% filter(ethnicity_decoded=="White") %>% mutate(new_group="White")),
                    (local_data %>% filter(ethnicity_decoded=="South Asian") %>% mutate(new_group="South Asian")),
                    (local_data %>% filter(ethnicity_decoded=="Black") %>% mutate(new_group="Black")),
                    (local_data %>% filter(ethnicity_decoded=="Mixed") %>% mutate(new_group="Mixed")),
                    (local_data %>% filter(ethnicity_decoded=="Other") %>% mutate(new_group="Other")),
                    (local_data %>% filter(ethnicity_decoded=="Missing") %>% mutate(new_group="Missing")),
                    (local_data %>% filter(imd_quintiles==1) %>% mutate(new_group="imd_1")),
                    (local_data %>% filter(imd_quintiles==2) %>% mutate(new_group="imd_2")),
                    (local_data %>% filter(imd_quintiles==3) %>% mutate(new_group="imd_3")),
                    (local_data %>% filter(imd_quintiles==4) %>% mutate(new_group="imd_4")),
                    (local_data %>% filter(imd_quintiles==5) %>% mutate(new_group="imd_5"))) %>%
  mutate(new_group=factor(new_group, levels=c("overall", "male", "female", "White", "South Asian", "Black", "Mixed", "Other", "Missing", "imd_1", "imd_2", "imd_3", "imd_4", "imd_5")))


new_cohort %>% filter(new_group=="overall") %>% count()
#66269
strat <- function (label, n, ...) {
  perc <-  paste0(round_pad(as.numeric(n)*100/66269,1), "%")
  sprintf("<span class='stratlabel'>%s<br><span class='stratn'>%s (%s)</span></span>", 
          label, prettyNum(n, big.mark=","), perc)
}
table1(~ missing | new_group, data=new_cohort, overall=F, render.categorical=render_cat_ci(digits=1), render.strat=strat)


new_cohort_t1 <- new_cohort %>% filter(diabetes_type=="type 1")
new_cohort_t1 %>% filter(new_group=="overall") %>% count()
#20,976
strat <- function (label, n, ...) {
  perc <-  paste0(round_pad(as.numeric(n)*100/20976,1), "%")
  sprintf("<span class='stratlabel'>%s<br><span class='stratn'>%s (%s)</span></span>", 
          label, prettyNum(n, big.mark=","), perc)
}
table1(~ missing | new_group, data=new_cohort_t1, overall=F, render.categorical=render_cat_ci(digits=1), render.strat=strat)

new_cohort_t2 <- new_cohort %>% filter(diabetes_type=="type 2")
new_cohort_t2 %>% filter(new_group=="overall") %>% count()
#45293
strat <- function (label, n, ...) {
  perc <-  paste0(round_pad(as.numeric(n)*100/45293,1), "%")
  sprintf("<span class='stratlabel'>%s<br><span class='stratn'>%s (%s)</span></span>", 
          label, prettyNum(n, big.mark=","), perc)
}
table1(~ missing | new_group, data=new_cohort_t2, overall=F, render.categorical=render_cat_ci(digits=1), render.strat=strat)


## Look at missingness by practice

test <- local_data %>% group_by(pracid) %>% summarise(count=n(), missing=sum(missing==1), missing_perc=missing*100/count)

summary(test$missing_perc)


## Also look at time to biomarker measurements

local_data <- cohort %>% select(index_bmidatediff, index_hdldatediff, index_totalcholesteroldatediff, index_triglyceridedatediff) %>% collect()

datediffs <- local_data %>%
  pivot_longer(
    cols = starts_with("index_"),    
    names_to = "measure",
    values_to = "datediff"
  )

datediffs %>% filter(!is.na(datediff)) %>% count()
#258361
datediffs %>% filter(!is.na(datediff) & datediff>=-(2*365.25)) %>% count()
#228054

228054/258361 #88%

datediffs <- local_data %>%
  filter(!is.na(index_bmidatediff) & !is.na(index_hdldatediff) & !is.na(index_totalcholesteroldatediff) & !is.na(index_triglyceridedatediff)) %>%
  pivot_longer(
    cols = starts_with("index_"),    
    names_to = "measure",
    values_to = "datediff"
  )

datediffs %>% filter(!is.na(datediff)) %>% count()
#248052
datediffs %>% filter(!is.na(datediff) & datediff>=-(2*365.25)) %>% count()
#219046

219046/248052 #88%
summary(datediffs)



cohort <- cohort %>%
  filter(!is.na(index_bmi) & !is.na(index_totalcholesterol) & !is.na(index_hdl) & !is.na(index_triglyceride)) %>%
  analysis$cached("cohort_interim_10", unique_indexes="patid")

cohort %>% count()
#62013

cohort %>% group_by(diabetes_type) %>% count()

# 1 type 2          43177
# 2 type 1          18836


############################################################################################

# Code up lipid model

cohort <- cohort %>%
  
  mutate(femalesex=ifelse(gender==2, 1, ifelse(gender==1, 0, NA)),
         
         lipid_pred_score=9.0034272-(0.1915482*index_bmi)-(0.1686227*dm_diag_age)+(0.3026012*femalesex)-(0.2269216*index_totalcholesterol)+(1.540850*index_hdl)-(0.2784059*index_triglyceride),
         lipid_pred_prob=exp(lipid_pred_score)/(1+exp(lipid_pred_score))) %>%
  
  analysis$cached("cohort", unique_index="patid")


cohort %>% count()
#62013

cohort %>% filter(!is.na(lipid_pred_prob)) %>% count()
#62013
