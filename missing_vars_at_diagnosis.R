
# Setup
library(tidyverse)
library(aurum)
rm(list=ls())

cprd = CPRDData$new(cprdEnv = "diabetes-jun2024",cprdConf = "~/.aurum.yaml")

analysis = cprd$analysis("dpctn_paper_inci")

cohort <- cohort %>% analysis$cached("cohort_interim_10")


############################################################################################

# Add closest biomarkers and look at dates

biomarkers <- c("bmi", "hba1c", "hdl", "totalcholesterol", "triglyceride")


# Pull out all raw biomarker values and cache

analysis = cprd$analysis("all_patid")
clean_bmi_medcodes <- clean_bmi_medcodes %>% analysis$cached("clean_bmi_medcodes")
clean_hba1c_medcodes <- clean_hba1c_medcodes %>% analysis$cached("clean_hba1c_medcodes")
clean_hdl_medcodes <- clean_hdl_medcodes %>% analysis$cached("clean_hdl_medcodes")
clean_totalcholesterol_medcodes <- clean_totalcholesterol_medcodes %>% analysis$cached("clean_totalcholesterol_medcodes")
clean_triglyceride_medcodes <- clean_triglyceride_medcodes %>% analysis$cached("clean_triglyceride_medcodes")


# For each biomarker, find values at diagnosis

analysis = cprd$analysis("dpctn_paper_inci")

for (i in biomarkers) {
  
  print(i)
  
  clean_tablename <- paste0("clean_", i, "_medcodes")
  biomarker_prediag_variable <- paste0("prediag_",i)
  biomarker_postdiag_variable <- paste0("postdiag_",i)
  biomarker_prediag_datediff_variable <- paste0("prediag_",i, "datediff")
  biomarker_postdiag_datediff_variable <- paste0("postdiag_",i, "datediff")
  interim_diag_biomarker_table <- paste0("diag_nowin_biomarkers_", i)
  
  
  prediag_value <- cohort %>%
    select(patid, diagnosis_date) %>%
    inner_join(get(clean_tablename), by="patid") %>%
    filter(date<=diagnosis_date) %>%
    mutate(diagdatediff=datediff(date, diagnosis_date)) %>%
    group_by(patid) %>%
    mutate(min_timediff=min(abs(diagdatediff), na.rm=TRUE)) %>%
    filter(abs(diagdatediff)==min_timediff) %>%
    mutate(pre_biomarker=min(testvalue, na.rm=TRUE)) %>%
    filter(pre_biomarker==testvalue) %>%
    dbplyr::window_order(diagdatediff) %>%
    filter(row_number()==1) %>%
    ungroup() %>%
    dbplyr::window_order() %>%
    select(patid, testvalue, diagdatediff) %>%
    rename({{biomarker_prediag_variable}}:=testvalue,
           {{biomarker_prediag_datediff_variable}}:=diagdatediff)
  
  postdiag_value <- cohort %>%
    select(patid, diagnosis_date) %>%
    inner_join(get(clean_tablename), by="patid") %>%
    filter(date>=diagnosis_date) %>%
    mutate(diagdatediff=datediff(date, diagnosis_date)) %>%
    group_by(patid) %>%
    mutate(min_timediff=min(abs(diagdatediff), na.rm=TRUE)) %>%
    filter(abs(diagdatediff)==min_timediff) %>%
    mutate(pre_biomarker=min(testvalue, na.rm=TRUE)) %>%
    filter(pre_biomarker==testvalue) %>%
    dbplyr::window_order(diagdatediff) %>%
    filter(row_number()==1) %>%
    ungroup() %>%
    dbplyr::window_order() %>%
    select(patid, testvalue, diagdatediff) %>%
    rename({{biomarker_postdiag_variable}}:=testvalue,
           {{biomarker_postdiag_datediff_variable}}:=diagdatediff)
  
  diag_value <- cohort %>%
    select(patid, diagnosis_date) %>%
    left_join(prediag_value, by="patid") %>%
    left_join(postdiag_value, by="patid") %>%
    analysis$cached(interim_diag_biomarker_table, unique_indexes="patid")
  
  assign(interim_diag_biomarker_table, diag_value)
  
}


############################################################################################

cohort <- cohort %>%
  left_join(diag_nowin_biomarkers_bmi, by="patid") %>%
  left_join(diag_nowin_biomarkers_hba1c, by="patid") %>%
  left_join(diag_nowin_biomarkers_hdl, by="patid") %>%
  left_join(diag_nowin_biomarkers_totalcholesterol, by="patid") %>%
  left_join(diag_nowin_biomarkers_triglyceride, by="patid") %>%
  analysis$cached("alternative_cohort_interim_11", unique_indexes="patid")


cohort %>% count()
#75,814


cohort %>% filter((prediag_bmidatediff>=-730 | postdiag_bmidatediff<=7) &
                    (prediag_hdldatediff>=-730 | postdiag_hdldatediff<=7) &
                    (prediag_totalcholesteroldatediff>=-730 | postdiag_totalcholesteroldatediff<=7) &
                    (prediag_triglyceridedatediff>=-730 | postdiag_triglyceridedatediff<=7)) %>% count()

#35,856


cohort %>% filter((prediag_bmidatediff>=-730 | postdiag_bmidatediff<=190) &
                    (prediag_hdldatediff>=-730 | postdiag_hdldatediff<=190) &
                    (prediag_totalcholesteroldatediff>=-730 | postdiag_totalcholesteroldatediff<=190) &
                    (prediag_triglyceridedatediff>=-730 | postdiag_triglyceridedatediff<=190)) %>% count()

#53,524


cohort %>% filter((prediag_bmidatediff>=-730 | postdiag_bmidatediff<=366) &
                    (prediag_hdldatediff>=-730 | postdiag_hdldatediff<=366) &
                    (prediag_totalcholesteroldatediff>=-730 | postdiag_totalcholesteroldatediff<=366) &
                    (prediag_triglyceridedatediff>=-730 | postdiag_triglyceridedatediff<=366)) %>% count()

#58,154


cohort %>% filter((prediag_bmidatediff>=-730 | postdiag_bmidatediff<=366) &
                    (prediag_hdldatediff>=-730 | postdiag_hdldatediff<=366) &
                    (prediag_totalcholesteroldatediff>=-730 | postdiag_totalcholesteroldatediff<=366) &
                    (prediag_triglyceridedatediff>=-730 | postdiag_triglyceridedatediff<=366)) %>% group_by(diabetes_type) %>% count()

# 1 type 2          55775
# 2 type 1           2379



cohort %>% filter((prediag_bmidatediff>=-366 | postdiag_bmidatediff<=366) &
                    (prediag_hdldatediff>=-366 | postdiag_hdldatediff<=366) &
                    (prediag_totalcholesteroldatediff>=-366 | postdiag_totalcholesteroldatediff<=366) &
                    (prediag_triglyceridedatediff>=-366 | postdiag_triglyceridedatediff<=366)) %>% count()

#56,841




