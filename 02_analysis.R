
# Analysis

## Reduce cohort to currently insulin treated
## Then do quality checking step
## Then remove if missing model variables

############################################################################################

# Setup
library(tidyverse)
library(aurum)
library(table1)
library(epiR)
rm(list=ls())

cprd = CPRDData$new(cprdEnv = "diabetes-jun2024",cprdConf = "~/.aurum.yaml")

analysis = cprd$analysis("dpctn_paper")


############################################################################################

# Set index date

index_date <- as.Date("2023-03-01")


############################################################################################

# Cohort refinement

cohort <- cohort %>% analysis$cached("cohort")
cohort %>% count()
#316076  #everyone with current diagnosis of T1 or T2 diagnosed aged 18-50

ins_cohort <- cohort %>% filter(current_insulin==1) %>% analysis$cached("ins_cohort", unique_indexes="patid")
ins_cohort %>% count() #86648


# Data quality checking
## Excluded as on SGLT2/DPP4/GLP1/TZD/SU
ins_cohort %>% filter(diabetes_type_new=="type 1" & (current_sglt2==1 | current_dpp4==1 | current_glp1==1 | current_gipglp1==1 | current_tzd==1 | current_su==1)) %>% count()
#1,182
## Not on insulin within 5 years of diagnosis
ins_cohort %>% filter(diabetes_type_new=="type 1" & year(diagnosis_date)>=1995 & (datediff(earliest_ins, diagnosis_date)/365.25)>5 & (datediff(regstartdate, diagnosis_date)/365.25)<=4) %>% count()
#494

## Either
ins_cohort %>% filter((diabetes_type_new=="type 1" & (current_sglt2==1 | current_dpp4==1 | current_glp1==1 | current_gipglp1==1 | current_tzd==1 | current_su==1)) | (diabetes_type_new=="type 1" & year(diagnosis_date)>=1995 & (datediff(earliest_ins, diagnosis_date)/365.25)>5 & (datediff(regstartdate, diagnosis_date)/365.25)<=4)) %>% count()
#1,588  


ins_cohort_quality_checked <- ins_cohort %>% filter(!(diabetes_type_new=="type 1" & (current_sglt2==1 | current_dpp4==1 | current_glp1==1 | current_gipglp1==1 | current_tzd==1 | current_su==1)) & !(diabetes_type_new=="type 1" & year(diagnosis_date)>=1995 & (datediff(earliest_ins, diagnosis_date)/365.25)>5 & (datediff(regstartdate, diagnosis_date)/365.25)<=4)) %>% analysis$cached("ins_cohort_quality_checked", unique_indexes="patid")
ins_cohort_quality_checked %>% count()
#85060


# Excluding if missing model variables

ins_cohort_quality_checked %>% filter(is.na(bmi_post_diag)) %>% count()
#1080
1080*100/85060 #1.3%

ins_cohort_quality_checked %>% filter(is.na(totalchol_post_diag)) %>% count()
#1276
1276*100/85060 #1.5%

ins_cohort_quality_checked %>% filter(is.na(hdl_post_diag)) %>% count()
#1554
1554*100/85060 #1.8%

ins_cohort_quality_checked %>% filter(is.na(triglyceride_post_diag)) %>% count()
#5258
5258*100/85060 #6.2%

ins_cohort_quality_checked %>% filter(is.na(lipid_pred_prob)) %>% count()
#5751
5751*100/85060 #6.8%

                                      
model_cohort <- ins_cohort_quality_checked %>% filter(!is.na(lipid_pred_prob)) %>% analysis$cached("model_cohort", unique_indexes="patid")
model_cohort %>% count()
#79309

model_cohort %>% group_by(diabetes_type_new) %>% count()
23333/79309
55976/79309


## How recent biomarkers are

model_cohort %>% filter(bmi_post_diag_datediff>=-730 & hdl_post_diag_datediff>=-730 & totalchol_post_diag_datediff>=-730 & triglyceride_post_diag_datediff>=-730) %>% count()
#56107

56107/79309 #70.7% within 2 years


############################################################################################

# Table 1

local_vars <- model_cohort %>%
  mutate(hba1c_post_diag=ifelse(index_hba1c_date>=diagnosis_date, index_hba1c, NA),
         lipid_prob_perc=lipid_pred_prob*100,
         malesex=ifelse(gender==1, 1, 0),
         ethnicity_decoded=case_when(ethnicity_5cat==0 ~ "White",
                                     ethnicity_5cat==1 ~ "South Asian",
                                     ethnicity_5cat==2 ~ "Black",
                                     ethnicity_5cat==3 ~"Other",
                                     ethnicity_5cat==4 ~"Mixed")) %>%
         # imd_quintiles=case_when(is.na(imd2015_10) ~NULL,
         #                         imd2015_10==1 | imd2015_10==2 ~1,
         #                         imd2015_10==3 | imd2015_10==4 ~2,
         #                         imd2015_10==5 | imd2015_10==6 ~3,
         #                         imd2015_10==7 | imd2015_10==8 ~4,
         #                         imd2015_10==9 | imd2015_10==10 ~5),
         #
  select(diabetes_type_new,
         diagnosis_date,
         lipid_prob_perc,
         malesex,
         dm_diag_age,
         bmi_post_diag,
         hdl_post_diag,
         totalchol_post_diag,
         triglyceride_post_diag,
         hba1c_post_diag,
         #dka_at_diagnosis,
         #dka_post_diagnosis,
         #hypo_post_diagnosis,
         follow_up_time,
         #dka_count,
         #hypo_count,
         ins_1_year,
         ins_2_years,
         ins_3_years,
         current_bolus_insulin,
         current_oha,
         ten_yrs_post_diag_oha,
         ten_yrs_post_diag_ins,
         ten_yrs_post_diag_bolus_ins,
         ethnicity_decoded,
         #imd_quintiles,
         age_at_index,
         mixed_codes,
         regstartdate,
         ten_yrs_post_diag,
         weight_loss_at_diag,
         polydipsia_at_diag,
         urinary_freq_at_diag,
         any_autoimmune,
         coeliac_thyroid) %>%
  collect() %>%
  mutate(diabetes_type_new=as.factor(diabetes_type_new),
         malesex=as.factor(malesex),
         #dka_at_diagnosis=factor(dka_at_diagnosis),
         #dka_post_diagnosis=factor(dka_post_diagnosis),
         #hypo_post_diagnosis=factor(hypo_post_diagnosis),
         #dka_count=as.integer(dka_count),
         #hypo_count=as.integer(hypo_count),
         ins_1_year=factor(ins_1_year),
         ins_2_years=factor(ins_2_years),
         ins_3_years=factor(ins_3_years),
         current_bolus_insulin=as.factor(current_bolus_insulin),
         current_oha=as.factor(current_oha),
         ten_yrs_post_diag_oha=as.factor(ten_yrs_post_diag_oha),
         ten_yrs_post_diag_ins=as.factor(ten_yrs_post_diag_ins),
         ten_yrs_post_diag_bolus_ins=as.factor(ten_yrs_post_diag_bolus_ins),
         ethnicity_decoded=factor(ethnicity_decoded, levels=c("White", "South Asian", "Black", "Mixed", "Other")),
         #imd_quintiles=as.factor(imd_quintiles),
         mixed_codes=as.factor(mixed_codes),
         weight_loss_at_diag=as.factor(weight_loss_at_diag),
         polydipsia_at_diag=as.factor(polydipsia_at_diag),
         urinary_freq_at_diag=as.factor(urinary_freq_at_diag),
         any_autoimmune=as.factor(any_autoimmune),
         coeliac_thyroid=as.factor(coeliac_thyroid))
         

# Table 1

type_1_misclass <- local_vars %>% filter(diabetes_type_new=="type 1" & lipid_prob_perc<=5) %>% mutate(group="type_1_misclass")
type_1_correctclass <- local_vars %>% filter(diabetes_type_new=="type 1" & lipid_prob_perc>5) %>% mutate(group="type_1_correctclass")
type_2_misclass <- local_vars %>% filter(diabetes_type_new=="type 2" & lipid_prob_perc>=70) %>% mutate(group="type_2_misclass")
type_2_correctclass <- local_vars %>% filter(diabetes_type_new=="type 2" & lipid_prob_perc<70) %>% mutate(group="type_2_correctclass")
gs_type_1 <- local_vars %>% filter(diabetes_type_new=="type 1" & lipid_prob_perc>=70) %>% mutate(group="gs_type_1")
gs_type_2 <- local_vars %>% filter(diabetes_type_new=="type 2" & lipid_prob_perc<=5) %>% mutate(group="gs_type_2")
groups <- rbind(type_1_misclass, type_1_correctclass, type_2_misclass, type_2_correctclass, gs_type_1, gs_type_2)

# add empty group for gaps on graph
groups <- rbind(groups, NA) %>% mutate(group=ifelse(is.na(group), "a", group))

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

# Add labels
label(groups$lipid_prob_perc)         <- "Type 1 probability (%)"
label(groups$malesex)                 <- "Sex (% male)"
label(groups$age_at_index)            <- "Current age (years)"
label(groups$ethnicity_decoded)       <- "Ethnicity"
#label(groups$imd_quintiles)           <- "Index of Multiple Deprivation quintile"
label(groups$dm_diag_age)             <- "Diabetes diagnosis age (years)"
label(groups$bmi_post_diag)           <- "BMI (kg/m2)"
label(groups$hdl_post_diag)           <- "HDL (mmol/L)"
label(groups$totalchol_post_diag)     <- "Total cholesterol (mmol/L)"
label(groups$triglyceride_post_diag)  <- "Triglycerides (mmol/L)"
label(groups$hba1c_post_diag)         <- "HbA1c (mmol/mol)"
label(groups$mixed_codes)         <- "History of other type of diabetes codes"


table1(~ lipid_prob_perc + malesex + age_at_index + ethnicity_decoded + dm_diag_age + bmi_post_diag + hdl_post_diag + totalchol_post_diag + triglyceride_post_diag + hba1c_post_diag + mixed_codes | group, data=groups, overall=F, render=rndr, render.categorical=cat, render.continuous=cont, render.strat=strat)
#superscript kg/m2

# Add in IMD quintile


############################################################################################

# Type 2 treatment

### Treatment

#### On insulin within 1 year - cohort or those with this variable

local_vars %>% summarise(with_ins_1_year_perc=sum(!is.na(ins_1_year))/n())

groups %>% group_by(group) %>% summarise(with_ins_1_year_perc=sum(!is.na(ins_1_year))/n())
#available for 34.5% type 2 misclassed
#available for 59.6% type 2 correct class
#available for 38.5% GS type 1
#available for 64.4% GS type 2

t2_ins_1_yr_data <- groups %>% filter((group=="type_2_misclass" | group=="type_2_correctclass" | group=="a" | group=="gs_type_1" | group=="gs_type_2") & (!is.na(ins_1_year) | group=="a")) %>% mutate(group=factor(group, levels=c("type_2_misclass", "type_2_correctclass", "a", "gs_type_1", "gs_type_2")))

ins_1_year <- t2_ins_1_yr_data %>% group_by(group) %>% summarise(count=n(), outcome="Insulin within 1\nyear of diagnosis", perc=sum(ins_1_year==1)/n(), upper=perc + (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count)), lower=perc - (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count)))

tiff("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2025/DePICtion/Paper/Plots/t2_ins_1_yr.tiff", width=6, height=8, units = "in", res=800)

ggplot(ins_1_year, aes(fill=group, y=perc*100, x=outcome, ymax=35)) + 
  geom_bar(position="dodge", stat="identity", width=0.7) +
  geom_text(aes(label=paste0(round(perc*100,0), "%"), group = group, fontface = "bold"),
            position = position_dodge(width = 0.7), vjust=-2, size=7) +
  geom_errorbar(aes(ymin=lower*100, ymax=upper*100), width=.3, position=position_dodge(.7)) +
  theme_bw() +
  ylab("Percentage (%)") +
  scale_fill_manual("legend", labels=c("type_2_misclass"="type 2 potentially misclassified", "type_2_correctclass"="type 2; no query over classification", "gs_type_1"="Gold standard type 1", "gs_type_2"="Gold standard type 2"), values = c("type_2_misclass" = "darkred", "type_2_correctclass" = "seagreen3", "gs_type_1"="dodgerblue3", "gs_type_2"="darkgoldenrod2")) +
  scale_y_continuous(limits=c(0,100), breaks=c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  theme(panel.grid.major.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20, face="bold",vjust=-1.5),
        panel.grid.major.y = element_line(size=1),
        panel.grid.minor.y = element_blank(),
        axis.title.x=element_blank(),
        legend.position = "none",
        legend.text=element_text(size=18, margin = margin(t=10, b=10, unit = "pt")),
        legend.title=element_blank(),
        plot.margin = unit(c(0,0,0.7,0.7), "cm"))

dev.off()



### On OHA at 10 years post-treatment - only in cohort registered at this point

local_vars %>% count()
local_vars %>% filter(regstartdate<=ten_yrs_post_diag & ten_yrs_post_diag<=index_date) %>% count()
42757/79309 #53.9%

t2_10yr_data <- groups %>% filter((group=="type_2_misclass" | group=="type_2_correctclass" | group=="a" | group=="gs_type_1" | group=="gs_type_2") & ((regstartdate<=ten_yrs_post_diag & ten_yrs_post_diag<=index_date) | group=="a")) %>% mutate(group=factor(group, levels=c("type_2_misclass", "type_2_correctclass", "a", "gs_type_1", "gs_type_2")))

ten_yrs_oha <- t2_10yr_data %>% group_by(group) %>% summarise(count=n(), outcome="Non-insulin diabetes meds\nat 10 years post-diagnosis", perc=sum(ten_yrs_post_diag_oha==1)/n(), upper=perc + (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count)), lower=perc - (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count)))

ten_yrs_bolus_ins <- t2_10yr_data %>% group_by(group) %>% summarise(count=n(), outcome="Basal-bolus insulin regime\nat 10 years post-diagnosis", perc=sum(ten_yrs_post_diag_bolus_ins==1)/n(), upper=perc + (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count)), lower=perc - (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count)))

t2_chart_data <- rbind(ten_yrs_oha, ten_yrs_bolus_ins) %>%
  mutate(outcome=factor(outcome, levels=c("Non-insulin diabetes meds\nat 10 years post-diagnosis", "Basal-bolus insulin regime\nat 10 years post-diagnosis")))


tiff("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2025/DePICtion/Paper/Plots/t2_treatment.tiff", width=15, height=8, units = "in", res=800)

ggplot(t2_chart_data, aes(fill=group, y=perc*100, x=outcome, ymax=35)) + 
  geom_bar(position="dodge", stat="identity", width=0.7) +
  geom_text(aes(label=paste0(round(perc*100,0), "%"), group = group, fontface = "bold"),
            position = position_dodge(width = 0.7), vjust=-2, size=7) +
  geom_errorbar(aes(ymin=lower*100, ymax=upper*100), width=.3, position=position_dodge(.7)) +
  theme_bw() +
  ylab("Percentage (%)") +
  scale_fill_manual("legend", labels=c("type_2_misclass"="type 2 potentially misclassified", "type_2_correctclass"="type 2; no query over classification", "gs_type_1"="Gold standard type 1", "gs_type_2"="Gold standard type 2"), values = c("type_2_misclass" = "darkred", "type_2_correctclass" = "seagreen3", "gs_type_1"="dodgerblue3", "gs_type_2"="darkgoldenrod2")) +
  scale_y_continuous(limits=c(0,100), breaks=c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  theme(panel.grid.major.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20, face="bold",vjust=-1.5),
        panel.grid.major.y = element_line(size=1),
        panel.grid.minor.y = element_blank(),
        axis.title.x=element_blank(),
        legend.position = c(1.22, 0.85),
        legend.text=element_text(size=17, margin = margin(t=8, b=8, unit = "pt")),
        legend.title=element_blank(),
        plot.margin = unit(c(0,11,0.7,0.7), "cm"))

dev.off()


############################################################################################

# Get old data for hosp

# For now - use old data for bar charts

cprd = CPRDData$new(cprdEnv = "diabetes-2020",cprdConf = "~/.aurum.yaml")
analysis = cprd$analysis("katie_duk24")
old_cohort <- old_cohort %>% analysis$cached("cohort")
old_local_vars <- old_cohort %>%
  filter(!is.na(lipid_pred_prob)) %>%
  mutate(lipid_prob_perc=lipid_pred_prob*100,
         dka_post_diagnosis=ifelse(!is.na(earliest_dka), 1L, 0L),
         hypo_post_diagnosis=ifelse(!is.na(earliest_hypo), 1L, 0L),
         follow_up_time=datediff(index_date, diagnosis_date)/365.25) %>%
  select(patid, diagnosis_date, earliest_ins, regstartdate, diabetes_type_new, lipid_prob_perc, dka_post_diagnosis, hypo_post_diagnosis, dka_at_diagnosis, dka_count, hypo_count, follow_up_time, contact_count, age_at_index) %>%
  collect() %>%
  mutate(dka_post_diagnosis=as.factor(dka_post_diagnosis),
         hypo_post_diagnosis=as.factor(hypo_post_diagnosis),
         dka_at_diagnosis=as.factor(dka_at_diagnosis),
         dka_count=as.integer(dka_count),
         hypo_count=as.integer(hypo_count),
         contact_count=as.integer(contact_count))

type_1_misclass <- old_local_vars %>% filter(diabetes_type_new=="type 1" & lipid_prob_perc<=5) %>% mutate(group="type_1_misclass")
type_1_correctclass <- old_local_vars %>% filter(diabetes_type_new=="type 1" & lipid_prob_perc>5) %>% mutate(group="type_1_correctclass")
type_2_misclass <- old_local_vars %>% filter(diabetes_type_new=="type 2" & lipid_prob_perc>=70) %>% mutate(group="type_2_misclass")
type_2_correctclass <- old_local_vars %>% filter(diabetes_type_new=="type 2" & lipid_prob_perc<70) %>% mutate(group="type_2_correctclass")
gs_type_1 <- old_local_vars %>% filter(diabetes_type_new=="type 1" & lipid_prob_perc>=70) %>% mutate(group="gs_type_1")
gs_type_2 <- old_local_vars %>% filter(diabetes_type_new=="type 2" & lipid_prob_perc<=5) %>% mutate(group="gs_type_2")
old_groups <- rbind(type_1_misclass, type_1_correctclass, type_2_misclass, type_2_correctclass, gs_type_1, gs_type_2)

# add empty group for gaps on graph

old_groups <- rbind(old_groups, NA) %>% mutate(group=ifelse(is.na(group), "a", group))


############################################################################################

# Type 2 hosp

t2_data_old <- old_groups %>% filter(group=="type_2_misclass" | group=="type_2_correctclass" | group=="gs_type_1" | group=="gs_type_2" | group=="a") %>% mutate(group=factor(group, levels=c("type_2_misclass", "type_2_correctclass", "a", "gs_type_1", "gs_type_2")))
                                                                                                                                                                             
dka_ever <- t2_data_old %>% group_by(group) %>% summarise(count=n(), outcome="dka_ever", perc=sum(dka_post_diagnosis==1)/n(), upper=perc + (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count)), lower=perc - (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count)))

hypo_ever <- t2_data_old %>% group_by(group) %>% summarise(count=n(), outcome="hypo_ever", perc=sum(hypo_post_diagnosis==1)/n(), upper=perc + (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count)), lower=perc - (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count)))


t2_chart_data <- rbind(dka_ever, hypo_ever) %>%
  mutate(outcome=factor(ifelse(outcome=="dka_ever", "DKA hospitalisation\n(post-diagnosis)", "Hypoglycaemia\nhospitalisation"), levels=c("DKA hospitalisation\n(post-diagnosis)", "Hypoglycaemia\nhospitalisation")))


tiff("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2025/DePICtion/Paper/Plots/t2_hosp.tiff", width=10.6, height=8, units = "in", res=800)

ggplot(t2_chart_data, aes(fill=group, y=perc*100, x=outcome, ymax=35)) + 
  geom_bar(position="dodge", stat="identity", width=0.7) +
  geom_text(aes(label=paste0(round(perc*100,0), "%"), group = group, fontface = "bold"),
            position = position_dodge(width = 0.7), vjust=-3.3, size=7) +
  geom_errorbar(aes(ymin=lower*100, ymax=upper*100), width=.3, position=position_dodge(.7)) +
  theme_bw() +
  ylab("Percentage (%)") +
  scale_fill_manual("legend", labels=c("type_2_misclass"="type 2 potentially misclassified", "type_2_correctclass"="type 2; no query over classification", "gs_type_1"="Gold standard type 1", "gs_type_2"="Gold standard type 2"), values = c("type_2_misclass" = "darkred", "type_2_correctclass" = "seagreen3", "gs_type_1"="dodgerblue3", "gs_type_2"="darkgoldenrod2")) +
  scale_y_continuous(limits=c(0,30), breaks=c(0, 5, 10, 15, 20, 25, 30)) +
  theme(panel.grid.major.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20, face="bold",vjust=-1.5),
        panel.grid.major.y = element_line(size=1),
        panel.grid.minor.y = element_blank(),
        axis.title.x=element_blank(),
        legend.position = c(0.75, 0.8),
        legend.text=element_text(size=17, margin = margin(t=8, b=8, unit = "pt")),
        legend.title=element_blank(),
        plot.margin = unit(c(0,0,0.7,0.7), "cm"))

dev.off()


#### Test for differences between misclass and correctclass
data <- t2_data_old %>% filter(group=="type_2_misclass" | group=="type_2_correctclass") %>% mutate(group=factor(group, levels=c("type_2_misclass", "type_2_correctclass")))

chisq.test(table(data$group, data$dka_post_diagnosis)) #This is same as prop.test
#p<0.0001

chisq.test(table(data$group, data$hypo_post_diagnosis)) #This is same as prop.test
#p<0.0001


#### Look at incidence rates

data %>% group_by(group) %>% summarise(median=median(follow_up_time))

rate_data <- data %>% group_by(group) %>% summarise(events=sum(dka_post_diagnosis==1), person_time=sum(follow_up_time)) %>% select(-group)

epi.2by2(
  dat = as.table(as.matrix(rate_data)),
  method = "cohort.time",
  conf.level = 0.95
)
#IR: 0.42 (0.34 to 0.51)
#IR: 0.15 (0.14 to 0.16)

#IRR: 2.76 (2.21, 3.41)

#p-value:
irr <- 2.76
ci_lower <- 2.21
ci_upper <- 3.41
se_log_irr <- (log(ci_upper) - log(ci_lower)) / (2 * 1.96)
z <- (log(irr) - log(1)) / se_log_irr
2 * (1 - pnorm(abs(z)))
#<0.0001


## Actual hypo rate - normalised by follow-up time

rate_data <- data %>% group_by(group) %>% summarise(events=sum(hypo_post_diagnosis==1), person_time=sum(follow_up_time)) %>% select(-group)

epi.2by2(
  dat = as.table(as.matrix(rate_data)),
  method = "cohort.time",
  conf.level = 0.95
)
#IR: 0.29 (0.23 to 0.37)
#IR: 0.15 (0.14 to 0.16)

#IRR: 1.91 (1.46, 2.45)

#p-value:
irr <- 1.91
ci_lower <- 1.46
ci_upper <- 2.45
se_log_irr <- (log(ci_upper) - log(ci_lower)) / (2 * 1.96)
z <- (log(irr) - log(1)) / se_log_irr
2 * (1 - pnorm(abs(z)))
#p<0.0001


## DKA at diagnosis

dka_at_diag <- t2_data_old %>% group_by(group) %>% summarise(count=n(), outcome="dka_at_diag", perc=sum(dka_at_diagnosis==1)/n(), upper=perc + (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count)), lower=perc - (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count)))


dka_at_diag



############################################################################################

# Type 1 hosp

t1_data_old <- old_groups %>% filter(group=="type_1_misclass" | group=="type_1_correctclass" | group=="gs_type_1" | group=="gs_type_2" | group=="a") %>% mutate(group=factor(group, levels=c("type_1_misclass", "type_1_correctclass", "a", "gs_type_1", "gs_type_2")))

dka_ever <- t1_data_old %>% group_by(group) %>% summarise(count=n(), outcome="dka_ever", perc=sum(dka_post_diagnosis==1)/n(), upper=perc + (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count)), lower=perc - (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count)))

hypo_ever <- t1_data_old %>% group_by(group) %>% summarise(count=n(), outcome="hypo_ever", perc=sum(hypo_post_diagnosis==1)/n(), upper=perc + (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count)), lower=perc - (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count)))


t1_chart_data <- rbind(dka_ever, hypo_ever) %>%
  mutate(outcome=factor(ifelse(outcome=="dka_ever", "DKA hospitalisation\n(post-diagnosis)", "Hypoglycaemia\nhospitalisation"), levels=c("DKA hospitalisation\n(post-diagnosis)", "Hypoglycaemia\nhospitalisation")))


tiff("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2025/DePICtion/Paper/Plots/t1_hosp.tiff", width=10.6, height=8, units = "in", res=800)

ggplot(t1_chart_data, aes(fill=group, y=perc*100, x=outcome, ymax=35)) + 
  geom_bar(position="dodge", stat="identity", width=0.7) +
  geom_text(aes(label=paste0(round(perc*100,0), "%"), group = group, fontface = "bold"),
            position = position_dodge(width = 0.7), vjust=-2, size=7) +
  geom_errorbar(aes(ymin=lower*100, ymax=upper*100), width=.3, position=position_dodge(.7)) +
  theme_bw() +
  ylab("Percentage (%)") +
  scale_fill_manual("legend", labels=c("type_1_misclass"="type 1 potentially misclassified", "type_1_correctclass"="type 1; no query over classification", "gs_type_1"="Gold standard type 1", "gs_type_2"="Gold standard type 2"), values = c("type_1_misclass" = "darkred", "type_1_correctclass" = "seagreen3", "gs_type_1"="dodgerblue3", "gs_type_2"="darkgoldenrod2")) +
  scale_y_continuous(limits=c(0,30), breaks=c(0, 5, 10, 15, 20, 25, 30)) +
  theme(panel.grid.major.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20, face="bold",vjust=-1.5),
        panel.grid.major.y = element_line(size=1),
        panel.grid.minor.y = element_blank(),
        axis.title.x=element_blank(),
        legend.position = c(0.75, 0.8),
        legend.text=element_text(size=17, margin = margin(t=8, b=8, unit = "pt")),
        legend.title=element_blank(),
        plot.margin = unit(c(0,0,0.7,0.7), "cm"))

dev.off()


#### Test for differences between misclass and correctclass
data <- t1_data_old %>% filter(group=="type_1_misclass" | group=="type_1_correctclass") %>% mutate(group=factor(group, levels=c("type_1_misclass", "type_1_correctclass")))

chisq.test(table(data$group, data$dka_post_diagnosis)) #This is same as prop.test
#p<0.0001

chisq.test(table(data$group, data$hypo_post_diagnosis)) #This is same as prop.test
#p<0.0001


#### Look at incidence rates

data %>% group_by(group) %>% summarise(median=median(follow_up_time))

rate_data <- data %>% group_by(group) %>% summarise(events=sum(dka_post_diagnosis==1), person_time=sum(follow_up_time)) %>% select(-group)

epi.2by2(
  dat = as.table(as.matrix(rate_data)),
  method = "cohort.time",
  conf.level = 0.95
)
#IR: 0.59 (0.52 to 0.66)
#IR: 0.80 (0.77 to 0.83)

#IRR: 0.73 (0.64, 0.83)

#p-value:
irr <- 0.73
ci_lower <- 0.64
ci_upper <- 0.83
se_log_irr <- (log(ci_upper) - log(ci_lower)) / (2 * 1.96)
z <- (log(irr) - log(1)) / se_log_irr
2 * (1 - pnorm(abs(z)))
#<0.0001


## Actual hypo rate - normalised by follow-up time

rate_data <- data %>% group_by(group) %>% summarise(events=sum(hypo_post_diagnosis==1), person_time=sum(follow_up_time)) %>% select(-group)

epi.2by2(
  dat = as.table(as.matrix(rate_data)),
  method = "cohort.time",
  conf.level = 0.95
)
#IR: 0.23 (0.19 to 0.28)
#IR: 0.34 (0.32 to 0.36)

#IRR: 0.68 (0.55, 0.83)

#p-value:
irr <- 0.68
ci_lower <- 0.55
ci_upper <- 0.83
se_log_irr <- (log(ci_upper) - log(ci_lower)) / (2 * 1.96)
z <- (log(irr) - log(1)) / se_log_irr
2 * (1 - pnorm(abs(z)))
#p=0.0002



## DKA at diagnosis

dka_at_diag <- t1_data_old %>% group_by(group) %>% summarise(count=n(), outcome="dka_at_diag", perc=sum(dka_at_diagnosis==1)/n(), upper=perc + (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count)), lower=perc - (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count)))


dka_at_diag


############################################################################################

# Look at those on insulin within 3 years of diagnosis
## Need to use old data for outcomes

local_vars %>% filter(diabetes_type_new=="type 2") %>% count()
local_vars %>% filter(diabetes_type_new=="type 2" & !is.na(ins_3_years)) %>% count()
46065/55976 #82.3

ins_3_years <- local_vars %>% filter(diabetes_type_new=="type 2" & !is.na(ins_3_years) & ins_3_years==1) %>% mutate(group="type_2_misclass_ins")
no_ins_3_years <- local_vars %>% filter(diabetes_type_new=="type 2" & !is.na(ins_3_years) & ins_3_years==0) %>% mutate(group="type_2_correctclass_ins")
type_2_misclass <- local_vars %>% filter(diabetes_type_new=="type 2" & lipid_prob_perc>=70 & !is.na(ins_3_years)) %>% mutate(group="type_2_misclass")
type_2_correctclass <- local_vars %>% filter(diabetes_type_new=="type 2" & lipid_prob_perc<70 & !is.na(ins_3_years)) %>% mutate(group="type_2_correctclass")

new_groups <- rbind(ins_3_years, no_ins_3_years, type_2_misclass, type_2_correctclass)

table1(~ lipid_prob_perc | group, data=new_groups, overall=F, render=rndr, render.categorical=cat, render.continuous=cont, render.strat=strat)
ins_3_years %>% filter(lipid_prob_perc>=70) %>% count()

634/46065
45431/46065

9067/46065
36998/46065


# Need old data for hospitalisations

old_local_vars <- old_local_vars %>%
  mutate(ins_3_years=ifelse(year(diagnosis_date)<1995, NA,
                            ifelse(!is.na(earliest_ins) & as.numeric(difftime(earliest_ins, diagnosis_date, units="days"))<=1096, 1L,
                                   ifelse(!is.na(earliest_ins) & as.numeric(difftime(earliest_ins, regstartdate, units="days"))<=366 & as.numeric(difftime(earliest_ins, regstartdate, units="days"))>731, NA, 0))))


ins_3_years <- old_local_vars %>% filter(diabetes_type_new=="type 2" & !is.na(ins_3_years) & ins_3_years==1) %>% mutate(group="type_2_misclass_ins")
no_ins_3_years <- old_local_vars %>% filter(diabetes_type_new=="type 2" & !is.na(ins_3_years) & ins_3_years==0) %>% mutate(group="type_2_correctclass_ins")
type_2_misclass <- old_local_vars %>% filter(diabetes_type_new=="type 2" & lipid_prob_perc>=70 & !is.na(ins_3_years)) %>% mutate(group="type_2_misclass")
type_2_correctclass <- old_local_vars %>% filter(diabetes_type_new=="type 2" & lipid_prob_perc<70 & !is.na(ins_3_years)) %>% mutate(group="type_2_correctclass")

new_groups_old_data <- rbind(ins_3_years, no_ins_3_years, type_2_misclass, type_2_correctclass)



dka_ever <- new_groups_old_data %>% group_by(group) %>% summarise(count=n(), outcome="dka_ever", perc=sum(dka_post_diagnosis==1)/n(), upper=perc + (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count)), lower=perc - (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count))) %>%
  mutate(perc=round(perc*100,1),
         lower=round(lower*100,1),
         upper=round(upper*100,1))

rate_data <- new_groups_old_data %>% filter(group=="type_2_misclass" | group=="type_2_correctclass") %>% mutate(group=as.factor(group)) %>% group_by(group) %>% summarise(events=sum(dka_post_diagnosis==1), person_time=sum(follow_up_time)) %>% select(-group)

epi.2by2(
  dat = as.table(as.matrix(rate_data)),
  method = "cohort.time",
  conf.level = 0.95
)

rate_data <- new_groups_old_data %>% filter(group=="type_2_misclass_ins" | group=="type_2_correctclass_ins") %>% mutate(group=as.factor(group)) %>% group_by(group) %>% summarise(events=sum(dka_post_diagnosis==1), person_time=sum(follow_up_time)) %>% select(-group)

epi.2by2(
  dat = as.table(as.matrix(rate_data)),
  method = "cohort.time",
  conf.level = 0.95
)



hypo_ever <- new_groups_old_data %>% group_by(group) %>% summarise(count=n(), outcome="hypo_ever", perc=sum(hypo_post_diagnosis==1)/n(), upper=perc + (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count)), lower=perc - (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count))) %>%
  mutate(perc=round(perc*100,1),
         lower=round(lower*100,1),
         upper=round(upper*100,1))

rate_data <- new_groups_old_data %>% filter(group=="type_2_misclass" | group=="type_2_correctclass") %>% mutate(group=as.factor(group)) %>% group_by(group) %>% summarise(events=sum(hypo_post_diagnosis==1), person_time=sum(follow_up_time)) %>% select(-group)

epi.2by2(
  dat = as.table(as.matrix(rate_data)),
  method = "cohort.time",
  conf.level = 0.95
)

rate_data <- new_groups_old_data %>% filter(group=="type_2_misclass_ins" | group=="type_2_correctclass_ins") %>% mutate(group=as.factor(group)) %>% group_by(group) %>% summarise(events=sum(hypo_post_diagnosis==1), person_time=sum(follow_up_time)) %>% select(-group)

epi.2by2(
  dat = as.table(as.matrix(rate_data)),
  method = "cohort.time",
  conf.level = 0.95
)

local_vars %>% filter(diabetes_type_new=="type 2" & !is.na(ins_3_years)) %>% count() #46065
local_vars %>% filter(diabetes_type_new=="type 2" & !is.na(ins_3_years) & !is.na(ins_1_year)) %>% count() #33067
33067/46065 #71.8%

ins_1_year <- new_groups %>% filter(!is.na(ins_1_year)) %>% group_by(group) %>% summarise(count=n(), outcome="Insulin within 1\nyear of diagnosis", perc=sum(ins_1_year==1)/n(), upper=perc + (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count)), lower=perc - (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count))) %>%
  mutate(perc=round(perc*100,1),
         lower=round(lower*100,1),
         upper=round(upper*100,1))



local_vars %>% filter(diabetes_type_new=="type 2" & !is.na(ins_3_years) & regstartdate<=ten_yrs_post_diag & ten_yrs_post_diag<=index_date) %>% count() #27947
27947/46065 #60.7

ten_yrs_oha <- new_groups %>% filter(regstartdate<=ten_yrs_post_diag & ten_yrs_post_diag<=index_date) %>% group_by(group) %>% summarise(count=n(), perc=sum(ten_yrs_post_diag_oha==1)/n(), upper=perc + (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count)), lower=perc - (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count))) %>%
  mutate(perc=round(perc*100,1),
         lower=round(lower*100,1),
         upper=round(upper*100,1))


ten_yrs_bolus_ins <- new_groups %>% filter(regstartdate<=ten_yrs_post_diag & ten_yrs_post_diag<=index_date) %>% group_by(group) %>% summarise(count=n(), perc=sum(ten_yrs_post_diag_bolus_ins==1)/n(), upper=perc + (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count)), lower=perc - (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count))) %>%
  mutate(perc=round(perc*100,1),
         lower=round(lower*100,1),
         upper=round(upper*100,1))


############################################################################################

# Look at contact count in old data

test <- old_groups %>% group_by(group) %>% mutate(contact_count=ifelse(is.na(contact_count), 0, contact_count)) %>% summarise(mean_contact_count=mean(contact_count), median_contact_count=median(contact_count), median_age=median(age_at_index))


############################################################################################


# C-peptide and Abs
 
model_cohort %>% count() #79309
model_cohort %>% filter(!is.na(earliest_c_pep_ins_deficient) | !is.na(earliest_c_pep_ins_intermediate) | !is.na(earliest_c_pep_ins_normal)) %>% count()#1286
1286/79309

model_cohort %>% count() #79309
model_cohort %>% filter(!is.na(earliest_negative_gad) | !is.na(earliest_positive_gad)) %>% count()#887
887/79309
model_cohort %>% filter(!is.na(earliest_negative_ia2) | !is.na(earliest_positive_ia2)) %>% count()#106
106/79309



# AI conditions

any_ai <- groups %>% group_by(group) %>% summarise(count=n(), perc=sum(any_autoimmune==1)/n(), upper=perc + (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count)), lower=perc - (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count))) %>%
  mutate(perc=round(perc*100,1),
         lower=round(lower*100,1),
         upper=round(upper*100,1)) %>%
  select(group, perc, lower, upper)

any_ai

coeliac_thyroid <- groups %>% group_by(group) %>% summarise(count=n(), perc=sum(coeliac_thyroid==1)/n(), upper=perc + (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count)), lower=perc - (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count))) %>%
  mutate(perc=round(perc*100,1),
         lower=round(lower*100,1),
         upper=round(upper*100,1)) %>%
  select(group, perc, lower, upper)

coeliac_thyroid

data <- groups %>% filter(group=="type_1_misclass" | group=="type_1_correctclass") %>% mutate(group=factor(group, levels=c("type_1_misclass", "type_1_correctclass")))

chisq.test(table(data$group, data$coeliac_thyroid)) #This is same as prop.test
#p<0.0001




# DKA at diagnosis

dka_at_diagnosis <- old_groups %>% group_by(group) %>% summarise(count=n(), outcome="dka_ever", perc=sum(dka_at_diagnosis==1)/n(), upper=perc + (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count)), lower=perc - (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count))) %>%
  mutate(perc=round(perc*100,1),
         lower=round(lower*100,1),
         upper=round(upper*100,1))
dka_at_diagnosis


# Primary care features
local_vars %>% count() #79309
local_vars %>% filter(regstartdate<=diagnosis_date) %>% count() #34252
34252/79309

weight_loss <- groups %>% filter(regstartdate<=diagnosis_date) %>% group_by(group) %>% summarise(count=n(), perc=sum(weight_loss_at_diag==1)/n(), upper=perc + (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count)), lower=perc - (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count))) %>%
  mutate(perc=round(perc*100,1),
         lower=round(lower*100,1),
         upper=round(upper*100,1))
weight_loss


polydipsia <- groups %>% filter(regstartdate<=diagnosis_date) %>% group_by(group) %>% summarise(count=n(), perc=sum(polydipsia_at_diag==1)/n(), upper=perc + (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count)), lower=perc - (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count))) %>%
  mutate(perc=round(perc*100,1),
         lower=round(lower*100,1),
         upper=round(upper*100,1))
polydipsia

polyuria <- groups %>% filter(regstartdate<=diagnosis_date) %>% group_by(group) %>% summarise(count=n(), perc=sum(urinary_freq_at_diag==1)/n(), upper=perc + (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count)), lower=perc - (qnorm(1-(0.05/2))*sqrt(perc*(1-perc)/count))) %>%
  mutate(perc=round(perc*100,1),
         lower=round(lower*100,1),
         upper=round(upper*100,1))
polyuria


############################################################################################

# Patient counts

#5.6 million people in the UK with diabetes out of 68.4 million
# Our download: ~1,252,239-263,444 with diabetes (excludes those with no type-specific codes) = 988795
# So to get total out of population: divide by (990,000*(68.4/5.6))
(990000*(68.4/5.6)) #12,092,143 (previously used 9,900,000)


# High scorers
model_cohort %>% filter(diabetes_type_new=="type 2" & !is.na(lipid_pred_prob)) %>% count()
#55976
model_cohort %>% filter(diabetes_type_new=="type 2" & !is.na(lipid_pred_prob) & lipid_pred_prob>=0.8) %>% count()
#598
598*100/55976 #1.1%
(598/12092143)*10000
(598/12092143)*40000
model_cohort %>% filter(diabetes_type_new=="type 2" & !is.na(lipid_pred_prob) & lipid_pred_prob>=0.7) %>% count()
#1099
1099*100/55976 #2.0%
(1099/12092143)*10000
(1099/12092143)*40000
model_cohort %>% filter(diabetes_type_new=="type 2" & !is.na(lipid_pred_prob) & lipid_pred_prob>=0.6) %>% count()
#1097
1765*100/55976 #3.2%
(1765/12092143)*10000
(1765/12092143)*40000
model_cohort %>% filter(diabetes_type_new=="type 2" & !is.na(lipid_pred_prob) & lipid_pred_prob>=0.5) %>% count()
#1611
2594*100/55976 #4.6%
(2594/12092143)*10000
(2594/12092143)*40000



model_cohort %>% filter(diabetes_type_new=="type 1" & !is.na(lipid_pred_prob)) %>% count()
#23333
model_cohort %>% filter(diabetes_type_new=="type 1" & !is.na(lipid_pred_prob) & lipid_pred_prob<=0.025) %>% count()
#1688
1688*100/23333 #7.2%
(1688/9900000)*10000
(1688/9900000)*40000
model_cohort %>% filter(diabetes_type_new=="type 1" & !is.na(lipid_pred_prob) & lipid_pred_prob<=0.05) %>% count()
#2933
2933*100/23333 #12.6%
(2933/9900000)*10000
(2933/9900000)*40000
model_cohort %>% filter(diabetes_type_new=="type 1" & !is.na(lipid_pred_prob) & lipid_pred_prob<=0.1) %>% count()
#4758
4758*100/23333 #20.4%
(4758/9900000)*10000
(4758/9900000)*40000



# Range per practice

summary((model_cohort %>% filter(diabetes_type_new=="type 2" & !is.na(lipid_pred_prob) & lipid_pred_prob>=0.8) %>% group_by(pracid) %>% summarise(count=n()) %>% collect())$count)

summary((model_cohort %>% filter(diabetes_type_new=="type 2" & !is.na(lipid_pred_prob) & lipid_pred_prob>=0.7) %>% group_by(pracid) %>% summarise(count=n()) %>% collect())$count)

summary((model_cohort %>% filter(diabetes_type_new=="type 2" & !is.na(lipid_pred_prob) & lipid_pred_prob>=0.6) %>% group_by(pracid) %>% summarise(count=n()) %>% collect())$count)

summary((model_cohort %>% filter(diabetes_type_new=="type 2" & !is.na(lipid_pred_prob) & lipid_pred_prob>=0.5) %>% group_by(pracid) %>% summarise(count=n()) %>% collect())$count)


summary((model_cohort %>% filter(diabetes_type_new=="type 1" & !is.na(lipid_pred_prob) & lipid_pred_prob<=0.025) %>% group_by(pracid) %>% summarise(count=n()) %>% collect())$count)

summary((model_cohort %>% filter(diabetes_type_new=="type 1" & !is.na(lipid_pred_prob) & lipid_pred_prob<=0.05) %>% group_by(pracid) %>% summarise(count=n()) %>% collect())$count)

summary((model_cohort %>% filter(diabetes_type_new=="type 1" & !is.na(lipid_pred_prob) & lipid_pred_prob<=0.1) %>% group_by(pracid) %>% summarise(count=n()) %>% collect())$count)
