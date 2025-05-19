
# Analysis

############################################################################################

# Setup
library(tidyverse)
library(aurum)
library(table1)
rm(list=ls())

cprd = CPRDData$new(cprdEnv = "diabetes-jun2024",cprdConf = "~/.aurum.yaml")


############################################################################################

# Import cohort data and reformat variables 

analysis = cprd$analysis("dpctn_paper_prev")
prev_cohort <- prev_cohort %>% analysis$cached("cohort")

#analysis = cprd$analysis("dpctn_paper_inci")
#inci_cohort <- inci_cohort %>% analysis$cached("cohort")

prev_cohort_local <- prev_cohort %>%
  select(diabetes_type,
         pracid,
         lipid_pred_prob,
         age_at_index,
         gender,
         ethnicity_5cat,
         imd_decile,
         dm_diag_age,
         bmi=index_bmi,
         hdl=index_hdl,
         totalcholesterol=index_totalcholesterol,
         triglyceride=index_triglyceride,
         hba1c=index_hba1c,
         coeliac,
         autoimmune_thyroid,
         coeliac_thyroid,
         t1_code_ever,
         t2_code_ever,
         weight_loss_at_diag,
         polydipsia_at_diag,
         urinary_freq_at_diag,
         with_hes,
         dka_at_diagnosis,
         earliest_dka,
         dka_post_diagnosis,
         dka_count,
         earliest_hypo,
         hypo_post_diagnosis,
         hypo_count,
         contact_count,
         follow_up_time,
         ins_1_year,
         ins_2_years,
         ins_3_years,
         current_bolus_insulin,
         current_oha,
         regstartdate,
         diagnosis_date,
         gp_record_end,
         ten_yrs_post_diag,
         ten_yrs_oha,
         ten_yrs_insulin,
         ten_yrs_bolus_insulin) %>%
  collect() %>%
  mutate(lipid_prob_perc=lipid_pred_prob*100,
         diabetes_type=as.factor(diabetes_type),
         malesex=as.factor(ifelse(gender==1, 1, 0)),
         dka_at_diagnosis=factor(dka_at_diagnosis),
         dka_post_diagnosis=factor(dka_post_diagnosis),
         hypo_post_diagnosis=factor(hypo_post_diagnosis),
         dka_count=as.integer(dka_count),
         hypo_count=as.integer(hypo_count),
         ins_1_year=factor(ins_1_year),
         ins_2_years=factor(ins_2_years),
         ins_3_years=factor(ins_3_years),
         current_bolus_insulin=as.factor(current_bolus_insulin),
         current_oha=as.factor(current_oha),
         ten_yrs_oha=as.factor(ten_yrs_oha),
         ten_yrs_insulin=as.factor(ten_yrs_insulin),
         ten_yrs_bolus_insulin=as.factor(ten_yrs_bolus_insulin),
         ethnicity_decoded=factor(case_when(ethnicity_5cat==0 ~ "White",
                                            ethnicity_5cat==1 ~ "South Asian",
                                            ethnicity_5cat==2 ~ "Black",
                                            ethnicity_5cat==3 ~"Other",
                                            ethnicity_5cat==4 ~"Mixed",
                                            is.na(ethnicity_5cat) ~"Missing"), levels=c("White", "South Asian", "Black", "Mixed", "Other", "Missing")),
         imd_quintiles=factor(case_when(imd_decile==1 | imd_decile==2 ~1,
                                        imd_decile==3 | imd_decile==4 ~2,
                                        imd_decile==5 | imd_decile==6 ~3,
                                        imd_decile==7 | imd_decile==8 ~4,
                                        imd_decile==9 | imd_decile==10 ~5)),
         t1_code_ever=as.factor(t1_code_ever),
         t2_code_ever=as.factor(t2_code_ever),
         weight_loss_at_diag=as.factor(weight_loss_at_diag),
         polydipsia_at_diag=as.factor(polydipsia_at_diag),
         urinary_freq_at_diag=as.factor(urinary_freq_at_diag),
         coeliac=as.factor(coeliac),
         autoimmune_thyroid=as.factor(autoimmune_thyroid),
         coeliac_thyroid=as.factor(coeliac_thyroid),
         diagnosis_date=as.Date(diagnosis_date, format="%Y-%m-%d"),
         regstartdate=as.Date(regstartdate, format="%Y-%m-%d"),
         ten_yrs_post_diag=as.Date(ten_yrs_post_diag, format="%Y-%m-%d"),
         cohort="prev")


# inci_cohort <- inci_cohort %>%
#   select(diabetes_type,
#          pracid,
#          lipid_pred_prob,
#          #age_at_index,
#          gender,
#          ethnicity_5cat,
#          imd_decile,
#          dm_diag_age,
#          bmi=diag_bmi,
#          hdl=diag_hdl,
#          totalcholesterol=diag_totalcholesterol,
#          triglyceride=diag_triglyceride,
#          hba1c=diag_hba1c,
#          coeliac,
#          autoimmune_thyroid,
#          coeliac_thyroid,
#          t1_code_ever,
#          t2_code_ever,
#          weight_loss_at_diag,
#          polydipsia_at_diag,
#          urinary_freq_at_diag,
#          with_hes,
#          dka_at_diagnosis,
#          earliest_dka,
#          dka_post_diagnosis,
#          dka_count,
#          earliest_hypo,
#          hypo_post_diagnosis,
#          hypo_count,
#          contact_count,
#          follow_up_time,
#          ins_1_year,
#          ins_2_years,
#          ins_3_years,
#          current_bolus_insulin,
#          current_oha,
#          regstartdate,
#          diagnosis_date,
#          gp_record_end,
#          ten_yrs_post_diag,
#          ten_yrs_oha,
#          ten_yrs_insulin,
#          ten_yrs_bolus_insulin) %>%
#   collect() %>%
#   mutate(lipid_prob_perc=lipid_pred_prob*100,
#          diabetes_type=as.factor(diabetes_type),
#          malesex=as.factor(ifelse(gender==1, 1, 0)),
#          dka_at_diagnosis=factor(dka_at_diagnosis),
#          dka_post_diagnosis=factor(dka_post_diagnosis),
#          hypo_post_diagnosis=factor(hypo_post_diagnosis),
#          dka_count=as.integer(dka_count),
#          hypo_count=as.integer(hypo_count),
#          ins_1_year=factor(ins_1_year),
#          ins_2_years=factor(ins_2_years),
#          ins_3_years=factor(ins_3_years),
#          current_bolus_insulin=as.factor(current_bolus_insulin),
#          current_oha=as.factor(current_oha),
#          ten_yrs_oha=as.factor(ten_yrs_oha),
#          ten_yrs_insulin=as.factor(ten_yrs_insulin),
#          ten_yrs_bolus_insulin=as.factor(ten_yrs_bolus_insulin),
#          ethnicity_decoded=factor(case_when(ethnicity_5cat==0 ~ "White",
#                                             ethnicity_5cat==1 ~ "South Asian",
#                                             ethnicity_5cat==2 ~ "Black",
#                                             ethnicity_5cat==3 ~"Other",
#                                             ethnicity_5cat==4 ~"Mixed",
#                                             is.na(ethnicity_5cat) ~"Missing"), levels=c("White", "South Asian", "Black", "Mixed", "Other", "Missing")),
#          imd_quintiles=factor(case_when(imd_decile==1 | imd_decile==2 ~1,
#                                         imd_decile==3 | imd_decile==4 ~2,
#                                         imd_decile==5 | imd_decile==6 ~3,
#                                         imd_decile==7 | imd_decile==8 ~4,
#                                         imd_decile==9 | imd_decile==10 ~5)),
#          t1_code_ever=as.factor(t1_code_ever),
#          t2_code_ever=as.factor(t2_code_ever),
#          weight_loss_at_diag=as.factor(weight_loss_at_diag),
#          polydipsia_at_diag=as.factor(polydipsia_at_diag),
#          urinary_freq_at_diag=as.factor(urinary_freq_at_diag),
#          coeliac=as.factor(coeliac),
#          autoimmune_thyroid=as.factor(autoimmune_thyroid),
#          coeliac_thyroid=as.factor(coeliac_thyroid),
#          diagnosis_date=as.Date(diagnosis_date, format="%Y-%m-%d"),
#          regstartdate=as.Date(regstartdate, format="%Y-%m-%d"),
#          ten_yrs_post_diag=as.Date(ten_yrs_post_diag, format="%Y-%m-%d"),
#          age_at_index=dm_diag_age,
#          cohort="inci")


prev_cohort_no_mixed_local <- prev_cohort %>%
  filter((diabetes_type=="type 2" & t1_code_ever==0) | (diabetes_type=="type 1" & t2_code_ever==0)) %>%
  select(diabetes_type,
         pracid,
         lipid_pred_prob,
         gender,
         age_at_index,
         ethnicity_5cat,
         imd_decile,
         dm_diag_age,
         bmi=index_bmi,
         hdl=index_hdl,
         totalcholesterol=index_totalcholesterol,
         triglyceride=index_triglyceride,
         hba1c=index_hba1c,
         coeliac,
         autoimmune_thyroid,
         coeliac_thyroid,
         t1_code_ever,
         t2_code_ever,
         weight_loss_at_diag,
         polydipsia_at_diag,
         urinary_freq_at_diag,
         with_hes,
         dka_at_diagnosis,
         earliest_dka,
         dka_post_diagnosis,
         dka_count,
         earliest_hypo,
         hypo_post_diagnosis,
         hypo_count,
         contact_count,
         follow_up_time,
         ins_1_year,
         ins_2_years,
         ins_3_years,
         current_bolus_insulin,
         current_oha,
         regstartdate,
         diagnosis_date,
         ten_yrs_post_diag,
         ten_yrs_oha,
         ten_yrs_insulin,
         ten_yrs_bolus_insulin,
         gp_record_end) %>%
  collect() %>%
  mutate(lipid_prob_perc=lipid_pred_prob*100,
         diabetes_type=as.factor(diabetes_type),
         malesex=as.factor(ifelse(gender==1, 1, 0)),
         dka_at_diagnosis=factor(dka_at_diagnosis),
         dka_post_diagnosis=factor(dka_post_diagnosis),
         hypo_post_diagnosis=factor(hypo_post_diagnosis),
         dka_count=as.integer(dka_count),
         hypo_count=as.integer(hypo_count),
         ins_1_year=factor(ins_1_year),
         ins_2_years=factor(ins_2_years),
         ins_3_years=factor(ins_3_years),
         current_bolus_insulin=as.factor(current_bolus_insulin),
         current_oha=as.factor(current_oha),
         ten_yrs_oha=as.factor(ten_yrs_oha),
         ten_yrs_insulin=as.factor(ten_yrs_insulin),
         ten_yrs_bolus_insulin=as.factor(ten_yrs_bolus_insulin),
         ethnicity_decoded=factor(case_when(ethnicity_5cat==0 ~ "White",
                                            ethnicity_5cat==1 ~ "South Asian",
                                            ethnicity_5cat==2 ~ "Black",
                                            ethnicity_5cat==3 ~"Other",
                                            ethnicity_5cat==4 ~"Mixed",
                                            is.na(ethnicity_5cat) ~"Missing"), levels=c("White", "South Asian", "Black", "Mixed", "Other", "Missing")),
         imd_quintiles=factor(case_when(imd_decile==1 | imd_decile==2 ~1,
                                        imd_decile==3 | imd_decile==4 ~2,
                                        imd_decile==5 | imd_decile==6 ~3,
                                        imd_decile==7 | imd_decile==8 ~4,
                                        imd_decile==9 | imd_decile==10 ~5)),
         t1_code_ever=as.factor(t1_code_ever),
         t2_code_ever=as.factor(t2_code_ever),
         weight_loss_at_diag=as.factor(weight_loss_at_diag),
         polydipsia_at_diag=as.factor(polydipsia_at_diag),
         urinary_freq_at_diag=as.factor(urinary_freq_at_diag),
         coeliac=as.factor(coeliac),
         autoimmune_thyroid=as.factor(autoimmune_thyroid),
         coeliac_thyroid=as.factor(coeliac_thyroid),
         diagnosis_date=as.Date(diagnosis_date, format="%Y-%m-%d"),
         regstartdate=as.Date(regstartdate, format="%Y-%m-%d"),
         ten_yrs_post_diag=as.Date(ten_yrs_post_diag, format="%Y-%m-%d"),
         cohort="prev_no_mixed")


prev_cohort_nda_local <- prev_cohort %>%
  filter(diabetes_type=="type 1" & (current_bolus_insulin==1 | current_mix_insulin==1) & current_su==0 & current_dpp4==0) %>%
  select(diabetes_type,
         pracid,
         lipid_pred_prob,
         gender,
         age_at_index,
         ethnicity_5cat,
         imd_decile,
         dm_diag_age,
         bmi=index_bmi,
         hdl=index_hdl,
         totalcholesterol=index_totalcholesterol,
         triglyceride=index_triglyceride,
         hba1c=index_hba1c,
         coeliac,
         autoimmune_thyroid,
         coeliac_thyroid,
         t1_code_ever,
         t2_code_ever,
         weight_loss_at_diag,
         polydipsia_at_diag,
         urinary_freq_at_diag,
         with_hes,
         dka_at_diagnosis,
         earliest_dka,
         dka_post_diagnosis,
         dka_count,
         earliest_hypo,
         hypo_post_diagnosis,
         hypo_count,
         contact_count,
         follow_up_time,
         ins_1_year,
         ins_2_years,
         ins_3_years,
         current_bolus_insulin,
         current_oha,
         regstartdate,
         diagnosis_date,
         ten_yrs_post_diag,
         ten_yrs_oha,
         ten_yrs_insulin,
         ten_yrs_bolus_insulin,
         gp_record_end) %>%
  collect() %>%
  mutate(lipid_prob_perc=lipid_pred_prob*100,
         diabetes_type=as.factor(diabetes_type),
         malesex=as.factor(ifelse(gender==1, 1, 0)),
         dka_at_diagnosis=factor(dka_at_diagnosis),
         dka_post_diagnosis=factor(dka_post_diagnosis),
         hypo_post_diagnosis=factor(hypo_post_diagnosis),
         dka_count=as.integer(dka_count),
         hypo_count=as.integer(hypo_count),
         ins_1_year=factor(ins_1_year),
         ins_2_years=factor(ins_2_years),
         ins_3_years=factor(ins_3_years),
         current_bolus_insulin=as.factor(current_bolus_insulin),
         current_oha=as.factor(current_oha),
         ten_yrs_oha=as.factor(ten_yrs_oha),
         ten_yrs_insulin=as.factor(ten_yrs_insulin),
         ten_yrs_bolus_insulin=as.factor(ten_yrs_bolus_insulin),
         ethnicity_decoded=factor(case_when(ethnicity_5cat==0 ~ "White",
                                            ethnicity_5cat==1 ~ "South Asian",
                                            ethnicity_5cat==2 ~ "Black",
                                            ethnicity_5cat==3 ~"Other",
                                            ethnicity_5cat==4 ~"Mixed",
                                            is.na(ethnicity_5cat) ~"Missing"), levels=c("White", "South Asian", "Black", "Mixed", "Other", "Missing")),
         imd_quintiles=factor(case_when(imd_decile==1 | imd_decile==2 ~1,
                                        imd_decile==3 | imd_decile==4 ~2,
                                        imd_decile==5 | imd_decile==6 ~3,
                                        imd_decile==7 | imd_decile==8 ~4,
                                        imd_decile==9 | imd_decile==10 ~5)),
         t1_code_ever=as.factor(t1_code_ever),
         t2_code_ever=as.factor(t2_code_ever),
         weight_loss_at_diag=as.factor(weight_loss_at_diag),
         polydipsia_at_diag=as.factor(polydipsia_at_diag),
         urinary_freq_at_diag=as.factor(urinary_freq_at_diag),
         coeliac=as.factor(coeliac),
         autoimmune_thyroid=as.factor(autoimmune_thyroid),
         coeliac_thyroid=as.factor(coeliac_thyroid),
         diagnosis_date=as.Date(diagnosis_date, format="%Y-%m-%d"),
         regstartdate=as.Date(regstartdate, format="%Y-%m-%d"),
         ten_yrs_post_diag=as.Date(ten_yrs_post_diag, format="%Y-%m-%d"),
         cohort="prev_nda")


data <- rbind(prev_cohort_local, prev_cohort_no_mixed_local, prev_cohort_nda_local) %>% mutate(group=paste0(cohort, "_", diabetes_type))

data <- data %>% filter(dm_diag_age<=50)


############################################################################################

# Table 1

cat <- function(x, ...) {
  c("", sapply(stats.apply.rounding(stats.default(x)), function(y) with(y, sprintf("%s (%s%%)", prettyNum(FREQ, big.mark=","), PCT))))
}

cont <- function(x) {
  with(stats.apply.rounding(stats.default(x)), c("Median (IQR)"=sprintf("%s (%s-%s)", round_pad(as.numeric(MEDIAN),1), round_pad(as.numeric(Q1),1), round_pad(as.numeric(Q3),1))))
}

missing <- function(x, ...) {
  with(stats.apply.rounding(stats.default(x)), c("Missing"=sprintf("%s", prettyNum(NMISS, big.mark=","))))
}

strat <- function (label, n, ...) {
  sprintf("<span class='stratlabel'>%s<br><span class='stratn'>(N=%s)</span></span>", 
          label, prettyNum(n, big.mark=","))
}

rndr <- function(x, name, ...) {
  y <- render.default(x, name, ...)
  if (is.logical(x)) {
    y[2]
  } else {
    y
  }
}

# without missing category
table1(~ lipid_prob_perc + malesex + age_at_index + ethnicity_decoded + imd_quintiles + dm_diag_age + bmi + hdl + totalcholesterol + triglyceride + hba1c + coeliac_thyroid + t1_code_ever + t2_code_ever | group, data=data, overall=F, render=rndr, render.categorical=cat, render.continuous=cont, render.strat=strat, render.missing=NULL)

