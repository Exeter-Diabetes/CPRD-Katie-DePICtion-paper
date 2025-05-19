
# Analysis

############################################################################################

# Setup
library(tidyverse)
library(aurum)
library(table1)
library(epiR)
rm(list=ls())

cprd = CPRDData$new(cprdEnv = "diabetes-jun2024",cprdConf = "~/.aurum.yaml")

analysis = cprd$analysis("dpctn_paper_inci")


############################################################################################

# Import cohort data and reformat variables 

cohort <- cohort %>% analysis$cached("cohort")

local_vars <- cohort %>%
  select(diabetes_type,
         pracid,
         lipid_pred_prob,
         gender,
         ethnicity_5cat,
         imd_decile,
         dm_diag_age,
         diag_bmi,
         diag_hdl,
         diag_totalcholesterol,
         diag_triglyceride,
         diag_hba1c,
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
         ten_yrs_post_diag=as.Date(ten_yrs_post_diag, format="%Y-%m-%d"))


############################################################################################

# Define groups based on thresholds

## As per Nick's paper: >=80% for defining gold standard type 1, <5% for defining gold standard type 2

type_1_misclass <- local_vars %>% filter(diabetes_type=="type 1" & lipid_prob_perc<5) %>% mutate(group="type_1_misclass")
type_1_correctclass <- local_vars %>% filter(diabetes_type=="type 1" & lipid_prob_perc>=5) %>% mutate(group="type_1_correctclass")

type_2_misclass <- local_vars %>% filter(diabetes_type=="type 2" & lipid_prob_perc>=80) %>% mutate(group="type_2_misclass")
type_2_correctclass <- local_vars %>% filter(diabetes_type=="type 2" & lipid_prob_perc<80) %>% mutate(group="type_2_correctclass")

ref_type_1 <- local_vars %>% filter(diabetes_type=="type 1" & lipid_prob_perc>=80) %>% mutate(group="ref_type_1")
ref_type_2 <- local_vars %>% filter(diabetes_type=="type 2" & lipid_prob_perc<5) %>% mutate(group="ref_type_2")

groups <- rbind(type_1_misclass, type_1_correctclass, type_2_misclass, type_2_correctclass, ref_type_1, ref_type_2)


main_cohort <- groups #already removed those diagnosed >50


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
table1(~ lipid_prob_perc + malesex + ethnicity_decoded + imd_quintiles + dm_diag_age + diag_bmi + diag_hdl + diag_totalcholesterol + diag_triglyceride + diag_hba1c + coeliac_thyroid + t1_code_ever + t2_code_ever | group, data=main_cohort, overall=F, render=rndr, render.categorical=cat, render.continuous=cont, render.strat=strat, render.missing=NULL)


############################################################################################

# Type 2 hospitalisation

t2_main_cohort <- main_cohort %>%
  filter(group!="type_1_correctclass" & group!="type_1_misclass")


### DKA at diagnosis

dka_at_diag <- t2_main_cohort %>%
  group_by(group) %>%
  summarise(total=n(),
            outcome="DKA hospitalisation\nat diagnosis",
            with_outcome=sum(dka_at_diagnosis==1)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(
    test = list(prop.test(with_outcome, total)),
    proportion = test$estimate,
    conf.low = test$conf.int[1],
    conf.high = test$conf.int[2]
  ) %>%
  union(data.frame(group="a", outcome="DKA hospitalisation\nat diagnosis", total=NA, with_outcome=NA, test=NA, proportion=NA, conf.low=NA, conf.high=NA)) %>%
  mutate(group=factor(group, levels=c("type_2_misclass", "type_2_correctclass", "a", "ref_type_1", "ref_type_2")))


t2_main_cohort %>%
  summarise(p = prop.test(
    c(sum(dka_at_diagnosis[group == "type_2_misclass"]==1), sum(dka_at_diagnosis[group == "type_2_correctclass"]==1)),
    c(sum(group == "type_2_misclass"   ), sum(group == "type_2_correctclass"    ))
  )$p.value) %>%
  pull(p)
#p=1


# tiff("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2025/DePICtion/Paper/Plots/diag_t2_dka_diag.tiff", width=6, height=8, units = "in", res=800)
# 
# ggplot(dka_at_diag, aes(fill=group, y=proportion*100, x=outcome, ymax=20)) + 
#   geom_bar(position="dodge", stat="identity", width=0.7) +
#   geom_text(aes(label=paste0(round(proportion*100,1), "%"), group = group, fontface = "bold"),
#             position = position_dodge(width = 0.7), vjust=-6, size=7) +
#   geom_errorbar(aes(ymin=conf.low*100, ymax=conf.high*100), width=.3, position=position_dodge(.7)) +
#   theme_bw() +
#   ylab("Percentage (%)") +
#   scale_fill_manual("legend", labels=c("type_2_misclass"="type 2 potentially misclassified", "type_2_correctclass"="type 2; no query over classification", "ref_type_1"="Reference type 1", "ref_type_2"="Reference type 2"), values = c("type_2_misclass" = "darkred", "type_2_correctclass" = "seagreen3", "ref_type_1"="dodgerblue3", "ref_type_2"="darkgoldenrod2")) +
#   scale_y_continuous(limits=c(0,15), breaks=c(0, 2.5, 5, 7.5, 10, 12.5, 15)) +
#   scale_x_discrete(expand = expansion(add = 0.5)) +
#   theme(panel.grid.major.x = element_blank(),
#         axis.title.y = element_text(size = 20),
#         axis.text.y = element_text(size = 20),
#         axis.text.x = element_text(size = 20, face="bold",vjust=-1.5),
#         panel.grid.major.y = element_line(size=1),
#         panel.grid.minor.y = element_blank(),
#         axis.title.x=element_blank(),
#         legend.position="none",
#         #legend.position = c(0.8, 0.8),
#         #legend.text=element_text(size=17, margin = margin(t=8, b=8, l=4, unit = "pt")),
#         #legend.title=element_blank(),
#         plot.margin = unit(c(0,0,0.7,0.7), "cm"))
# 
# dev.off()



### DKA and hypos post-hosp

dka_post_diag <- t2_main_cohort %>%
  group_by(group) %>%
  summarise(
    outcome="DKA hospitalisation\n(post-diagnosis)",
    IR    = poisson.test(sum(dka_post_diagnosis == 1), sum(follow_up_time))$estimate,
    lower = poisson.test(sum(dka_post_diagnosis == 1), sum(follow_up_time))$conf.int[1],
    upper = poisson.test(sum(dka_post_diagnosis == 1), sum(follow_up_time))$conf.int[2],
    .groups = "drop"
  ) %>%
  union(data.frame(group="a", outcome="DKA hospitalisation\n(post-diagnosis)", IR=NA, lower=NA, upper=NA)) %>%
  mutate(group=factor(group, levels=c("type_2_misclass", "type_2_correctclass", "a", "ref_type_1", "ref_type_2")))

t2_main_cohort %>%
  summarise(
    p = poisson.test(
      c(
        sum(dka_post_diagnosis[group == "type_2_misclass"] == 1),
        sum(dka_post_diagnosis[group == "type_2_correctclass"] == 1)
      ),
      c(
        sum(follow_up_time[group == "type_2_misclass"]),
        sum(follow_up_time[group == "type_2_correctclass"])
      )
    )$p.value
  ) %>%
  pull(p)
#p=1


hypo_post_diag <- t2_main_cohort %>%
  group_by(group) %>%
  summarise(
    outcome="Hypoglycaemia\nhospitalisation",
    IR    = poisson.test(sum(hypo_post_diagnosis == 1), sum(follow_up_time))$estimate,
    lower = poisson.test(sum(hypo_post_diagnosis == 1), sum(follow_up_time))$conf.int[1],
    upper = poisson.test(sum(hypo_post_diagnosis == 1), sum(follow_up_time))$conf.int[2],
    .groups = "drop"
  ) %>%
  union(data.frame(group="a", outcome="Hypoglycaemia\nhospitalisation", IR=NA, lower=NA, upper=NA)) %>%
  mutate(group=factor(group, levels=c("type_2_misclass", "type_2_correctclass", "a", "ref_type_1", "ref_type_2")))

t2_main_cohort %>%
  summarise(
    p = poisson.test(
      c(
        sum(hypo_post_diagnosis[group == "type_2_misclass"] == 1),
        sum(hypo_post_diagnosis[group == "type_2_correctclass"] == 1)
      ),
      c(
        sum(follow_up_time[group == "type_2_misclass"]),
        sum(follow_up_time[group == "type_2_correctclass"])
      )
    )$p.value
  ) %>%
  pull(p)
#p=1

t2_chart_data <- rbind(dka_post_diag, hypo_post_diag) #%>%
#mutate(outcome=factor(outcome, levels=c("DKA hospitalisation\nat diagnosis", "DKA hospitalisation\n(post-diagnosis)", "Hypoglycaemia\nhospitalisation")))


# tiff("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2025/DePICtion/Paper/Plots/diag_t2_postdiaghosp.tiff", width=11, height=8, units = "in", res=800)
# 
# ggplot(t2_chart_data, aes(fill=group, y=IR*100, x=outcome, ymax=1.5)) + 
#   geom_bar(position="dodge", stat="identity", width=0.7) +
#   geom_text(aes(label=round(IR*100,2), group = group, fontface = "bold"),
#             position = position_dodge(width = 0.7), vjust=-5, size=7) +
#   geom_errorbar(aes(ymin=lower*100, ymax=upper*100), width=.3, position=position_dodge(.7)) +
#   theme_bw() +
#   ylab("Incident rate (patients per 100 patient-years)") +
#   scale_x_discrete(expand = expansion(add = 0.5)) +
#   scale_fill_manual("legend", labels=c("type_2_misclass"="type 2 potentially misclassified", "type_2_correctclass"="type 2; no query over classification", "ref_type_1"="Reference type 1", "ref_type_2"="Reference type 2"), values = c("type_2_misclass" = "darkred", "type_2_correctclass" = "seagreen3", "ref_type_1"="dodgerblue3", "ref_type_2"="darkgoldenrod2")) +
#   scale_y_continuous(limits=c(0,1.25), breaks=c(0, 0.25, 0.5, 0.75, 1, 1.25)) +
#   theme(panel.grid.major.x = element_blank(),
#         axis.title.y = element_text(size = 20),
#         axis.text.y = element_text(size = 20),
#         axis.text.x = element_text(size = 20, face="bold",vjust=-1.5),
#         panel.grid.major.y = element_line(size=1),
#         panel.grid.minor.y = element_blank(),
#         axis.title.x=element_blank(),
#         legend.position = c(0.75, 0.8),
#         legend.text=element_text(size=17, margin = margin(t=8, b=8, l=4, unit = "pt")),
#         legend.title=element_blank(),
#         plot.margin = unit(c(0,0,0.7,0.7), "cm"))
# 
# dev.off()


############################################################################################

## Type 2 treatment

### On insulin within 1 year - just use those with this available

local_vars %>% summarise(with_ins_1_year_perc=sum(!is.na(ins_1_year))/n())
#94.4% of main study cohort

t2_ins_1_yr_data <- t2_main_cohort %>% filter(!is.na(ins_1_year))

ins_1_year <- t2_ins_1_yr_data %>%
  group_by(group) %>%
  summarise(total=n(),
            outcome="Insulin within 1\nyear of diagnosis",
            with_outcome=sum(ins_1_year==1)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(
    test = list(prop.test(with_outcome, total)),
    proportion = test$estimate,
    conf.low = test$conf.int[1],
    conf.high = test$conf.int[2]
  ) %>%
  union(data.frame(group="a", outcome="Insulin within 1\nyear of diagnosis", total=NA, with_outcome=NA, test=NA, proportion=NA, conf.low=NA, conf.high=NA)) %>%
  mutate(group=factor(group, levels=c("type_2_misclass", "type_2_correctclass", "a", "ref_type_1", "ref_type_2")))

t2_ins_1_yr_data %>%
  summarise(p = prop.test(
    c(sum(ins_1_year[group == "type_2_misclass"]==1), sum(ins_1_year[group == "type_2_correctclass"]==1)),
    c(sum(group == "type_2_misclass"   ), sum(group == "type_2_correctclass"    ))
  )$p.value) %>%
  pull(p)
#p=0.001


tiff("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2025/DePICtion/Paper/Plots/diag_t2_ins_1_yr.tiff", width=6, height=8, units = "in", res=800)

ggplot(ins_1_year, aes(fill=group, y=proportion*100, x=outcome, ymax=35)) + 
  geom_bar(position="dodge", stat="identity", width=0.7) +
  geom_text(aes(label=paste0(round(proportion*100,0), "%"), group = group, fontface = "bold"),
            position = position_dodge(width = 0.7), vjust=-5.2, size=7) +
  geom_errorbar(aes(ymin=conf.low*100, ymax=conf.high*100), width=.3, position=position_dodge(.7)) +
  theme_bw() +
  ylab("Percentage (%)") +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  scale_fill_manual("legend", labels=c("type_2_misclass"="type 2 potentially misclassified", "type_2_correctclass"="type 2; no query over classification", "ref_type_1"="Reference type 1", "ref_type_2"="Reference type 2"), values = c("type_2_misclass" = "darkred", "type_2_correctclass" = "seagreen3", "ref_type_1"="dodgerblue3", "ref_type_2"="darkgoldenrod2")) +
  scale_y_continuous(limits=c(0,115), breaks=c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
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



### On OHA at 10 years post-treatment - only in cohort registered at this point and post-1995

local_vars %>% count()
local_vars %>% filter(!is.na(ten_yrs_post_diag)) %>% count()
15/45644 

# t2_10yr_data <- t2_main_cohort %>% filter(!is.na(ten_yrs_post_diag))
# 
# ten_yrs_oha <- t2_10yr_data %>%
#   group_by(group) %>%
#   summarise(total=n(),
#             outcome="Non-insulin diabetes meds\nat 10 years post-diagnosis",
#             with_outcome=sum(ten_yrs_oha==1)) %>%
#   ungroup() %>%
#   rowwise() %>%
#   mutate(
#     test = list(prop.test(with_outcome, total)),
#     proportion = test$estimate,
#     conf.low = test$conf.int[1],
#     conf.high = test$conf.int[2]
#   ) %>%
#   union(data.frame(group="a", outcome="Non-insulin diabetes meds\nat 10 years post-diagnosis", total=NA, with_outcome=NA, test=NA, proportion=NA, conf.low=NA, conf.high=NA)) %>%
#   mutate(group=factor(group, levels=c("type_2_misclass", "type_2_correctclass", "a", "ref_type_1", "ref_type_2")))
# 
# t2_10yr_data %>%
#   summarise(p = prop.test(
#     c(sum(ten_yrs_oha[group == "type_2_misclass"]==1), sum(ten_yrs_oha[group == "type_2_correctclass"]==1)),
#     c(sum(group == "type_2_misclass"   ), sum(group == "type_2_correctclass"    ))
#   )$p.value) %>%
#   pull(p)
# #p<0.001
# 
# 
# ten_yrs_bolus_insulin <- t2_10yr_data %>%
#   group_by(group) %>%
#   summarise(total=n(),
#             outcome="Basal-bolus insulin regime\nat 10 years post-diagnosis",
#             with_outcome=sum(ten_yrs_bolus_insulin==1)) %>%
#   ungroup() %>%
#   rowwise() %>%
#   mutate(
#     test = list(prop.test(with_outcome, total)),
#     proportion = test$estimate,
#     conf.low = test$conf.int[1],
#     conf.high = test$conf.int[2]
#   ) %>%
#   union(data.frame(group="a", outcome="Basal-bolus insulin regime\nat 10 years post-diagnosis", total=NA, with_outcome=NA, test=NA, proportion=NA, conf.low=NA, conf.high=NA)) %>%
#   mutate(group=factor(group, levels=c("type_2_misclass", "type_2_correctclass", "a", "ref_type_1", "ref_type_2")))
# 
# t2_10yr_data %>%
#   summarise(p = prop.test(
#     c(sum(ten_yrs_bolus_insulin[group == "type_2_misclass"]==1), sum(ten_yrs_bolus_insulin[group == "type_2_correctclass"]==1)),
#     c(sum(group == "type_2_misclass"   ), sum(group == "type_2_correctclass"    ))
#   )$p.value) %>%
#   pull(p)
# #p<0.001
# 
# 
# t2_chart_data <- rbind(ten_yrs_oha, ten_yrs_bolus_insulin)
# 
# t2_chart_data$outcome <- factor(t2_chart_data$outcome, levels=c("Non-insulin diabetes meds\nat 10 years post-diagnosis", "Basal-bolus insulin regime\nat 10 years post-diagnosis"))
# 
# tiff("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2025/DePICtion/Paper/Plots/diag_t2_treatment.tiff", width=15, height=8, units = "in", res=800)
# 
# ggplot(t2_chart_data, aes(fill=group, y=proportion*100, x=outcome, ymax=35)) + 
#   geom_bar(position="dodge", stat="identity", width=0.7) +
#   geom_text(aes(label=paste0(round(proportion*100,0), "%"), group = group, fontface = "bold"),
#             position = position_dodge(width = 0.7), vjust=-3, size=7) +
#   geom_errorbar(aes(ymin=conf.low*100, ymax=conf.high*100), width=.3, position=position_dodge(.7)) +
#   theme_bw() +
#   ylab("Percentage (%)") +
#   scale_x_discrete(expand = expansion(add = 0.5)) +
#   scale_fill_manual("legend", labels=c("type_2_misclass"="type 2 potentially misclassified", "type_2_correctclass"="type 2; no query over classification", "ref_type_1"="Reference type 1", "ref_type_2"="Reference type 2"), values = c("type_2_misclass" = "darkred", "type_2_correctclass" = "seagreen3", "ref_type_1"="dodgerblue3", "ref_type_2"="darkgoldenrod2")) +
#   scale_y_continuous(limits=c(0,100), breaks=c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
#   theme(panel.grid.major.x = element_blank(),
#         axis.title.y = element_text(size = 20),
#         axis.text.y = element_text(size = 20),
#         axis.text.x = element_text(size = 20, face="bold",vjust=-1.5),
#         panel.grid.major.y = element_line(size=1),
#         panel.grid.minor.y = element_blank(),
#         axis.title.x=element_blank(),
#         legend.position = c(1.22, 0.85),
#         legend.text=element_text(size=17, margin = margin(t=8, b=8, l=4, unit = "pt")),
#         legend.title=element_blank(),
#         plot.margin = unit(c(0,11,0.7,0.7), "cm"))
# 
# dev.off()


############################################################################################

# # Type 1 hospitalisation
# 
# t1_main_cohort <- main_cohort %>%
#   filter(group!="type_2_correctclass" & group!="type_2_misclass")
# 
# 
# ### DKA at diagnosis
# 
# dka_at_diag <- t1_main_cohort %>%
#   group_by(group) %>%
#   summarise(total=n(),
#             outcome="DKA hospitalisation\nat diagnosis",
#             with_outcome=sum(dka_at_diagnosis==1)) %>%
#   ungroup() %>%
#   rowwise() %>%
#   mutate(
#     test = list(prop.test(with_outcome, total)),
#     proportion = test$estimate,
#     conf.low = test$conf.int[1],
#     conf.high = test$conf.int[2]
#   ) %>%
#   union(data.frame(group="a", outcome="DKA hospitalisation\nat diagnosis", total=NA, with_outcome=NA, test=NA, proportion=NA, conf.low=NA, conf.high=NA)) %>%
#   mutate(group=factor(group, levels=c("type_1_misclass", "type_1_correctclass", "a", "ref_type_1", "ref_type_2")))
# 
# 
# t1_main_cohort %>%
#   summarise(p = prop.test(
#     c(sum(dka_at_diagnosis[group == "type_1_misclass"]==1), sum(dka_at_diagnosis[group == "type_1_correctclass"]==1)),
#     c(sum(group == "type_1_misclass"   ), sum(group == "type_1_correctclass"    ))
#   )$p.value) %>%
#   pull(p)
# #p=0.016
# 
# 
# tiff("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2025/DePICtion/Paper/Plots/diag_t1_dka_diag.tiff", width=6, height=8, units = "in", res=800)
# 
# ggplot(dka_at_diag, aes(fill=group, y=proportion*100, x=outcome, ymax=20)) + 
#   geom_bar(position="dodge", stat="identity", width=0.7) +
#   geom_text(aes(label=paste0(round(proportion*100,1), "%"), group = group, fontface = "bold"),
#             position = position_dodge(width = 0.7), vjust=-6, size=7) +
#   geom_errorbar(aes(ymin=conf.low*100, ymax=conf.high*100), width=.3, position=position_dodge(.7)) +
#   theme_bw() +
#   ylab("Percentage (%)") +
#   scale_fill_manual("legend", labels=c("type_1_misclass"="type 1 potentially misclassified", "type_1_correctclass"="type 1; no query over classification", "ref_type_1"="Reference type 1", "ref_type_2"="Reference type 2"), values = c("type_1_misclass" = "darkred", "type_1_correctclass" = "seagreen3", "ref_type_1"="dodgerblue3", "ref_type_2"="darkgoldenrod2")) +
#   scale_y_continuous(limits=c(0,15), breaks=c(0, 2.5, 5, 7.5, 10, 12.5, 15)) +
#   scale_x_discrete(expand = expansion(add = 0.5)) +
#   theme(panel.grid.major.x = element_blank(),
#         axis.title.y = element_text(size = 20),
#         axis.text.y = element_text(size = 20),
#         axis.text.x = element_text(size = 20, face="bold",vjust=-1.5),
#         panel.grid.major.y = element_line(size=1),
#         panel.grid.minor.y = element_blank(),
#         axis.title.x=element_blank(),
#         legend.position="none",
#         #legend.position = c(0.8, 0.8),
#         #legend.text=element_text(size=17, margin = margin(t=8, b=8, l=4, unit = "pt")),
#         #legend.title=element_blank(),
#         plot.margin = unit(c(0,0,0.7,0.7), "cm"))
# 
# dev.off()
# 
# 
# 
# ### DKA and hypos post-hosp
# 
# dka_post_diag <- t1_main_cohort %>%
#   group_by(group) %>%
#   summarise(
#     outcome="DKA hospitalisation\n(post-diagnosis)",
#     IR    = poisson.test(sum(dka_post_diagnosis == 1), sum(follow_up_time))$estimate,
#     lower = poisson.test(sum(dka_post_diagnosis == 1), sum(follow_up_time))$conf.int[1],
#     upper = poisson.test(sum(dka_post_diagnosis == 1), sum(follow_up_time))$conf.int[2],
#     .groups = "drop"
#   ) %>%
#   union(data.frame(group="a", outcome="DKA hospitalisation\n(post-diagnosis)", IR=NA, lower=NA, upper=NA)) %>%
#   mutate(group=factor(group, levels=c("type_1_misclass", "type_1_correctclass", "a", "ref_type_1", "ref_type_2")))
# 
# t1_main_cohort %>%
#   summarise(
#     p = poisson.test(
#       c(
#         sum(dka_post_diagnosis[group == "type_1_misclass"] == 1),
#         sum(dka_post_diagnosis[group == "type_1_correctclass"] == 1)
#       ),
#       c(
#         sum(follow_up_time[group == "type_1_misclass"]),
#         sum(follow_up_time[group == "type_1_correctclass"])
#       )
#     )$p.value
#   ) %>%
#   pull(p)
# #p=0.009
# 
# 
# hypo_post_diag <- t1_main_cohort %>%
#   group_by(group) %>%
#   summarise(
#     outcome="Hypoglycaemia\nhospitalisation",
#     IR    = poisson.test(sum(hypo_post_diagnosis == 1), sum(follow_up_time))$estimate,
#     lower = poisson.test(sum(hypo_post_diagnosis == 1), sum(follow_up_time))$conf.int[1],
#     upper = poisson.test(sum(hypo_post_diagnosis == 1), sum(follow_up_time))$conf.int[2],
#     .groups = "drop"
#   ) %>%
#   union(data.frame(group="a", outcome="Hypoglycaemia\nhospitalisation", IR=NA, lower=NA, upper=NA)) %>%
#   mutate(group=factor(group, levels=c("type_1_misclass", "type_1_correctclass", "a", "ref_type_1", "ref_type_2")))
# 
# t1_main_cohort %>%
#   summarise(
#     p = poisson.test(
#       c(
#         sum(hypo_post_diagnosis[group == "type_1_misclass"] == 1),
#         sum(hypo_post_diagnosis[group == "type_1_correctclass"] == 1)
#       ),
#       c(
#         sum(follow_up_time[group == "type_1_misclass"]),
#         sum(follow_up_time[group == "type_1_correctclass"])
#       )
#     )$p.value
#   ) %>%
#   pull(p)
# #p=0.13
# 
# t1_chart_data <- rbind(dka_post_diag, hypo_post_diag) #%>%
# #mutate(outcome=factor(outcome, levels=c("DKA hospitalisation\nat diagnosis", "DKA hospitalisation\n(post-diagnosis)", "Hypoglycaemia\nhospitalisation")))
# 
# 
# tiff("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2025/DePICtion/Paper/Plots/diag_t1_postdiaghosp.tiff", width=11, height=8, units = "in", res=800)
# 
# ggplot(t1_chart_data, aes(fill=group, y=IR*100, x=outcome, ymax=1.5)) + 
#   geom_bar(position="dodge", stat="identity", width=0.7) +
#   geom_text(aes(label=round_pad(IR*100,2), group = group, fontface = "bold"),
#             position = position_dodge(width = 0.7), vjust=-4, size=7) +
#   geom_errorbar(aes(ymin=lower*100, ymax=upper*100), width=.3, position=position_dodge(.7)) +
#   theme_bw() +
#   ylab("Incident rate (patients per 100 patient-years)") +
#   scale_x_discrete(expand = expansion(add = 0.5)) +
#   scale_fill_manual("legend", labels=c("type_1_misclass"="type 1 potentially misclassified", "type_1_correctclass"="type 1; no query over classification", "ref_type_1"="Reference type 1", "ref_type_2"="Reference type 2"), values = c("type_1_misclass" = "darkred", "type_1_correctclass" = "seagreen3", "ref_type_1"="dodgerblue3", "ref_type_2"="darkgoldenrod2")) +
#   scale_y_continuous(limits=c(0,1.25), breaks=c(0, 0.25, 0.5, 0.75, 1, 1.25)) +
#   theme(panel.grid.major.x = element_blank(),
#         axis.title.y = element_text(size = 20),
#         axis.text.y = element_text(size = 20),
#         axis.text.x = element_text(size = 20, face="bold",vjust=-1.5),
#         panel.grid.major.y = element_line(size=1),
#         panel.grid.minor.y = element_blank(),
#         axis.title.x=element_blank(),
#         legend.position = c(0.75, 0.8),
#         legend.text=element_text(size=17, margin = margin(t=8, b=8, l=4, unit = "pt")),
#         legend.title=element_blank(),
#         plot.margin = unit(c(0,0,0.7,0.7), "cm"))
# 
# dev.off()



