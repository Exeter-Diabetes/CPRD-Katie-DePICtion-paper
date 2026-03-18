
# Analysis

## Display items for paper:
# 1 Inclusion flowchart (Figure 1)
# 2. Missingness of model variables overall and by sex, ethnicity, and deprivation (Supplementary)
# 3. Characteristics of study population (those with model variables; Table 1)
# 
# Current diagnosis of type 2 diabetes:
# 4. 	a) Number and characteristics of patients identified as potentially misclassified or miscoded (Supplementary)
#     b) As above by sex, ethnicity and deprivation (Supplementary)
# 5. 	a) Outcomes of potentially misclassified patients compared to reference type 2s (Figure 2)
#     b) As above by sex, ethnicity and deprivation (Supplementary)
# 6.	a) Number and characteristics of patients on insulin within 3 years of diagnosis (alternative to using model threshold; Supplementary)
#     b) Outcomes of potentially misclassified patients compared to reference type 2s (Supplementary)
# 
# Current diagnosis of type 1 diabetes:
# 7. 	a) Number and characteristics of patients identified as potentially misclassified or miscoded (Supplementary)
#     b) As above by sex, ethnicity and deprivation (Supplementary)
# 8. 	a) Outcomes of potentially misclassified patients compared to reference type 1s (Figure 3)
#     b) As above by sex, ethnicity and deprivation (Supplementary)
# 
# 9. 	For both type 1 and type 2, number identified using different model thresholds (per N patients, per average practice) and outcomes of potentially misclassified patients compared to reference group (Table 2; additional note that using insulin within 3 years to flag extra type	2s adds ~4 people per practice at any of the thresholds used)


############################################################################################

# Setup
library(tidyverse)
library(aurum)
library(table1)
library(epitools)
library(gridExtra)
library(cowplot)
library(grid)
library(binom)
library(ggpubr)
rm(list=ls())

cprd=CPRDData$new(cprdEnv="diabetes-jun2024", cprdConf="~/.aurum.yaml")

analysis=cprd$analysis("dpctn_paper_prev")


############################################################################################

# Setup for formatting of tables

cat <- function(x, ...) {
  c("", sapply(stats.apply.rounding(stats.default(x)), function(y) with(y, sprintf("%s (%s%%)", prettyNum(FREQ, big.mark=","), PCT))))
}

cont <- function(x) {
  with(stats.apply.rounding(stats.default(x)), c("Median (IQR)"=sprintf("%s (%s-%s)", round_pad(as.numeric(MEDIAN),1), round_pad(as.numeric(Q1),1), round_pad(as.numeric(Q3),1))))
}

rndr <- function(x, name, ...) {
  y <- render.default(x, name, ...)
  if (is.logical(x)) {
    y[2]
  } else {
    y
  }
}

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

############################################################################################

# Import cohort data and reformat variables 

cohort <- cohort %>% analysis$cached("cohort")

local_vars <- cohort %>%
  select(diabetes_type,
         pracid,
         lipid_pred_prob,
         gender,
         age_at_index,
         ethnicity_5cat,
         imd_decile,
         dm_diag_age,
         index_bmi,
         index_hdl,
         index_totalcholesterol,
         index_triglyceride,
         index_hba1c,
         coeliac,
         autoimmune_thyroid,
         coeliac_thyroid,
         t1_code_ever,
         t2_code_ever,
         type1_code_count,
         type2_code_count,
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
         current_mix_insulin,
         current_su,
         current_dpp4,
         current_oha,
         regstartdate,
         diagnosis_date,
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
         current_mix_insulin=as.factor(current_mix_insulin),
         current_oha=as.factor(current_oha),
         current_su=as.factor(current_su),
         current_dpp4=as.factor(current_dpp4),

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

type_1_all <- local_vars %>% filter(diabetes_type=="type 1") %>% mutate(group="type_1_all")
type_1_misclass_miscode <- local_vars %>% filter(diabetes_type=="type 1" & lipid_prob_perc<5) %>% mutate(group="type_1_misclass_miscode")
type_1_misclass <- local_vars %>% filter(diabetes_type=="type 1" & t2_code_ever==0 & lipid_prob_perc<5) %>% mutate(group="type_1_misclass")

type_2_all <- local_vars %>% filter(diabetes_type=="type 2") %>% mutate(group="type_2_all")
type_2_misclass_miscode <- local_vars %>% filter(diabetes_type=="type 2" & lipid_prob_perc>=80) %>% mutate(group="type_2_misclass_miscode")
type_2_misclass <- local_vars %>% filter(diabetes_type=="type 2" & t1_code_ever==0 & lipid_prob_perc>=80) %>% mutate(group="type_2_misclass")

ref_type_1 <- local_vars %>% filter(diabetes_type=="type 1" & t2_code_ever==0 & lipid_prob_perc>=80) %>% mutate(group="ref_type_1")
ref_type_2 <- local_vars %>% filter(diabetes_type=="type 2" & t1_code_ever==0 & lipid_prob_perc<5) %>% mutate(group="ref_type_2")


############################################################################################

# Current type 2

## Characteristics

t2_table_cohort <- rbind(type_2_all, type_2_misclass_miscode, type_2_misclass, ref_type_2)

table1(~ lipid_prob_perc + malesex + age_at_index + ethnicity_decoded + imd_quintiles + dm_diag_age + index_bmi + index_hdl + index_totalcholesterol + index_triglyceride + index_hba1c + t1_code_ever | group, data=t2_table_cohort, overall=F, render=rndr, render.categorical=cat, render.continuous=cont)


### Look at pairwise p-values

get_pvalue <- function(var, group, data) {
  x <- data[[var]]
  g <- data[[group]]
  
  if (is.numeric(x)) {
    p <- kruskal.test(x ~ g)$p.value
  } else if (is.factor(x) || is.character(x)) {
    tbl <- table(x, g)
    p <- chisq.test(tbl)$p.value
  } else {
    return(NA)
  }
  return(p)
}

vars <- c("lipid_prob_perc", "malesex", "age_at_index", "ethnicity_decoded", "imd_quintiles", "dm_diag_age", "index_bmi", "index_hdl", "index_totalcholesterol", "index_triglyceride", "index_hba1c", "t1_code_ever")

## Compare flagged to correcly classified
cohort <- t2_table_cohort %>% filter(group!="type_2_all" & group!="type_2_misclass" & ethnicity_decoded!="Missing") %>% mutate(ethnicity_decoded=factor(ethnicity_decoded))
cohort %>% distinct(group)
pvals <- data.frame(p=sapply(vars, function(v) get_pvalue(v, "group", cohort))) %>%
  mutate(sig=ifelse(p<0.05, 1, NA))
pvals
# all are significantly different except IMD

## Compare flagged to misclassified
cohort <- t2_table_cohort %>% filter(group!="type_2_all" & group!="ref_type_2" & ethnicity_decoded!="Missing") %>% mutate(ethnicity_decoded=factor(ethnicity_decoded))
cohort %>% distinct(group)
pvals <- data.frame(p=sapply(vars, function(v) get_pvalue(v, "group", cohort))) %>%
  mutate(sig=ifelse(p<0.05, 1, NA))
pvals
# none are different except type 1 code ever - obviously

## Compare reference to all
cohort <- t2_table_cohort %>% filter(group!="type_2_misclass" & group!="type_2_misclass_miscode")
cohort %>% distinct(group)
pvals <- data.frame(p=sapply(vars, function(v) get_pvalue(v, "group", cohort))) %>%
  mutate(sig=ifelse(p<0.05, 1, NA))
pvals
# all are sig


t2_cohort_by_group <- rbind((type_2_all %>% mutate(new_group="overall")),
                              (type_2_all %>% filter(malesex==1) %>% mutate(new_group="male")),
                              (type_2_all %>% filter(malesex==0) %>% mutate(new_group="female")),
                              (type_2_all %>% filter(ethnicity_decoded=="White") %>% mutate(new_group="White")),
                              (type_2_all %>% filter(ethnicity_decoded=="South Asian") %>% mutate(new_group="South Asian")),
                              (type_2_all %>% filter(ethnicity_decoded=="Black") %>% mutate(new_group="Black")),
                              (type_2_all %>% filter(ethnicity_decoded=="Mixed") %>% mutate(new_group="Mixed")),
                              (type_2_all %>% filter(ethnicity_decoded=="Other") %>% mutate(new_group="Other")),
                              (type_2_all %>% filter(ethnicity_decoded=="Missing") %>% mutate(new_group="Missing")),
                              (type_2_all %>% filter(imd_quintiles==1) %>% mutate(new_group="imd_1")),
                              (type_2_all %>% filter(imd_quintiles==2) %>% mutate(new_group="imd_2")),
                              (type_2_all %>% filter(imd_quintiles==3) %>% mutate(new_group="imd_3")),
                              (type_2_all %>% filter(imd_quintiles==4) %>% mutate(new_group="imd_4")),
                              (type_2_all %>% filter(imd_quintiles==5) %>% mutate(new_group="imd_5"))) %>%
  mutate(new_group=factor(new_group, levels=c("overall", "male", "female", "White", "South Asian", "Black", "Mixed", "Other", "Missing", "imd_1", "imd_2", "imd_3", "imd_4", "imd_5")),
         flagged=as.factor(ifelse(lipid_prob_perc>=80, 1, 0)),
         misclass=as.factor(ifelse(flagged==1 & t1_code_ever==0, 1, 0)))


# Without CIs for proportions
table1(~ lipid_prob_perc + malesex + imd_quintiles + dm_diag_age + index_bmi + index_hdl + index_totalcholesterol + index_triglyceride | new_group, data=t2_cohort_by_group, overall=F, render=rndr, render.categorical=cat, render.continuous=cont)


t2_cohort_by_group_plot <- t2_cohort_by_group %>% mutate(miscode=as.factor(ifelse(flagged==1 & misclass==0, 1, 0)))

t2_cohort_by_group_plot <- t2_cohort_by_group_plot %>%
  filter(new_group!="Mixed" & new_group!="Other" & new_group!="Missing" & new_group!="imd_2" & new_group!="imd_3" & new_group!="imd_4") %>%
  mutate(new_group=case_when(new_group=="male" ~ "Male",
                              new_group=="female" ~ "Female",
                              new_group=="overall" ~ "Overall",
                              new_group=="imd_1" ~ "IMD quintile 1",
                              new_group=="imd_5" ~ "IMD quintile 5",
                              TRUE ~ new_group))
  

t2_cohort_by_group_plot$new_group <- factor(t2_cohort_by_group_plot$new_group,
                              levels = c("IMD quintile 5", "IMD quintile 1","Black","South Asian",  "White", "Female","Male", "Overall"))


plot_df <- t2_cohort_by_group_plot %>%
  group_by(new_group) %>%
  summarise(
    flagged_n = sum(flagged == 1, na.rm = TRUE),
    misclass_n = sum(misclass == 1, na.rm = TRUE),
    miscode_n = sum(miscode == 1, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

library(binom)

flagged_ci  <- binom.confint(plot_df$flagged_n,  plot_df$n, method = "wilson")
misclass_ci <- binom.confint(plot_df$misclass_n, plot_df$n, method = "wilson")
miscode_ci  <- binom.confint(plot_df$miscode_n,  plot_df$n, method = "wilson")


plot_df <- plot_df %>%
  mutate(
    flagged_est  = flagged_ci$mean,
    flagged_low  = flagged_ci$lower,
    flagged_upp  = flagged_ci$upper,
    
    misclass_est = misclass_ci$mean,
    misclass_low = misclass_ci$lower,
    misclass_upp = misclass_ci$upper,
    
    miscode_est  = miscode_ci$mean,
    miscode_low  = miscode_ci$lower,
    miscode_upp  = miscode_ci$upper
  )

plot_long <- plot_df %>%
  select(new_group,
         flagged_est, flagged_low, flagged_upp,
         misclass_est, misclass_low, misclass_upp,
         miscode_est, miscode_low, miscode_upp) %>%
  pivot_longer(
    cols = -new_group,
    names_to = c("outcome", ".value"),
    names_pattern = "(flagged|misclass|miscode)_(est|low|upp)"
  )

plot_long <- plot_long %>%
  mutate(outcome = recode(outcome,
                          flagged  = "Potentially misclassified or miscoded (model=T1)",
                          misclass = "Potentially misclassified (model=T1, no history of T1 codes)",
                          miscode  = "Potentially miscoded (model=T1, history of T1 codes)"
  ))

plot_long$outcome <- factor(plot_long$outcome, levels=c("Potentially misclassified or miscoded (model=T1)", "Potentially misclassified (model=T1, no history of T1 codes)", "Potentially miscoded (model=T1, history of T1 codes)"))



x_max <- 1.6
arrow_length <- 0.03  # small offset inside axis

# adjust truncated upper CI so arrow fits inside
plot_long <- plot_long %>%
  mutate(
    upp_trunc2 = ifelse(upp * 100 > x_max, x_max - arrow_length, upp * 100),
    arrow_needed = (upp * 100 > x_max),
    est_pct = est * 100,
    label_text = paste0(round(est_pct, 2), "%"),
    row_id = as.numeric(factor(new_group)),                # numeric row
    stripe = factor(row_id %% 2)                            # 0/1 as factor
  )

label_x <- -0.02 * x_max

# plotting

tiff("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2026/DePICtion/Plots/t2_proportion_subgroup.tiff",
     width = 10, height = 9, units = "in", res = 800)

ggplot(plot_long, aes(x = est_pct, y = new_group)) +
  # zebra stripes
  geom_rect(
    data = plot_long %>% distinct(new_group, row_id, stripe),
    aes(
      ymin = row_id - 0.5,
      ymax = row_id + 0.5,
      xmin = -Inf,
      xmax = Inf,
      fill = stripe
    ),
    inherit.aes = FALSE,
    alpha = 0.1
  ) +
  scale_fill_manual(values = c("white", "grey70"), guide = "none") +  # now works
  # points
  geom_point() +
  # error bars
  geom_errorbarh(aes(xmin = low * 100, xmax = upp_trunc2), height = 0) +
  # arrows for truncated CIs
  geom_segment(
    data = subset(plot_long, arrow_needed),
    aes(
      x = x_max - arrow_length,
      xend = x_max,
      y = new_group,
      yend = new_group
    ),
    arrow = arrow(length = unit(0.15, "cm"), ends = "last", type = "closed"),
    inherit.aes = FALSE
  ) +
  # left-aligned % labels
  geom_text(aes(y = new_group, label = label_text),
            x = label_x,
            hjust = 0,
            size = 5) +
  facet_wrap(~ outcome, ncol = 1) +
  xlim(label_x, x_max) +
  labs(x = "Percentage (%)", y = NULL) +
  theme_classic(base_size = 16) +
  theme(axis.text.y = element_text(hjust = 1))

dev.off()






## Outcomes in misclassified compared to reference groups

setwd("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2026/DePICtion/Scripts/Functions/")
source("outcomes.R")

outcomes_t2_cohort <- rbind(type_2_misclass, ref_type_2)

outcomes_t2_cohort_by_group <- rbind((outcomes_t2_cohort %>% mutate(new_group="overall")),
                                 (outcomes_t2_cohort %>% filter(malesex==1) %>% mutate(new_group="male")),
                                 (outcomes_t2_cohort %>% filter(malesex==0) %>% mutate(new_group="female")),
                                 (outcomes_t2_cohort %>% filter(ethnicity_decoded=="White") %>% mutate(new_group="White")),
                                 (outcomes_t2_cohort %>% filter(ethnicity_decoded=="South Asian") %>% mutate(new_group="South Asian")),
                                 (outcomes_t2_cohort %>% filter(ethnicity_decoded=="Black") %>% mutate(new_group="Black")),
                                 (outcomes_t2_cohort %>% filter(ethnicity_decoded=="Mixed") %>% mutate(new_group="Mixed")),
                                 (outcomes_t2_cohort %>% filter(ethnicity_decoded=="Other") %>% mutate(new_group="Other")),
                                 (outcomes_t2_cohort %>% filter(ethnicity_decoded=="Missing") %>% mutate(new_group="Missing")),
                                 (outcomes_t2_cohort %>% filter(imd_quintiles==1) %>% mutate(new_group="imd_1")),
                                 (outcomes_t2_cohort %>% filter(imd_quintiles==2) %>% mutate(new_group="imd_2")),
                                 (outcomes_t2_cohort %>% filter(imd_quintiles==3) %>% mutate(new_group="imd_3")),
                                 (outcomes_t2_cohort %>% filter(imd_quintiles==4) %>% mutate(new_group="imd_4")),
                                 (outcomes_t2_cohort %>% filter(imd_quintiles==5) %>% mutate(new_group="imd_5")))

new_groups <- c("overall", "male", "female", "White", "South Asian", "Black", "imd_1", "imd_5")


for (i in new_groups) {
  
  data <- outcomes_t2_cohort_by_group %>% filter(new_group==i)
  
  outcomes_t2(data, i)
  
}

compare_t2_sex(outcomes_t2_cohort)

compare_t2_ethnicity(outcomes_t2_cohort)

compare_t2_deprivation(outcomes_t2_cohort)


# check autoimmune diffs
outcomes_t2_cohort %>%
  summarise(p=prop.test(
    c(sum(coeliac_thyroid[group=="type_2_misclass"]==1), sum(coeliac_thyroid[group=="ref_type_2"]==1)),
    c(sum(group=="type_2_misclass"   ), sum(group=="ref_type_2"    ))
  )$p.value) %>%
  pull(p)
#p=0.9



# Type 2 and on insulin within 3 years

local_vars %>% filter(diabetes_type=="type 2") %>% count()
local_vars %>% filter(diabetes_type=="type 2" & !is.na(ins_3_years)) %>% count()
12496/43180 #28.9%

type_2_ins_3yrs_cohort <- local_vars %>%
  filter(diabetes_type=="type 2" & !is.na(ins_3_years)) %>%
  mutate(model_prob_over_80=factor(ifelse(lipid_prob_perc>=80, 1, 0)))


type_2_ins_3yrs_cohort %>% filter(ins_3_years==1) %>% count()
#3,778
3778/12496

type_2_ins_3yrs_cohort %>% filter(ins_3_years==1 & t1_code_ever==0) %>% count()
#3,311

type_2_ins_3yrs_cohort %>% filter(ins_3_years==1 & t1_code_ever==1) %>% count()
#3,311


table1(~ lipid_prob_perc + model_prob_over_80 + malesex + age_at_index + ethnicity_decoded + imd_quintiles + dm_diag_age + index_bmi + index_hdl + index_totalcholesterol + index_triglyceride + index_hba1c + t1_code_ever | ins_3_years, data=type_2_ins_3yrs_cohort, overall=F, render=rndr, render.categorical=cat, render.continuous=cont)

type_2_ins_3yrs_cohort_misclass_miscod <- type_2_ins_3yrs_cohort %>% filter(ins_3_years==1)

table1(~ lipid_prob_perc + model_prob_over_80 + malesex + age_at_index + ethnicity_decoded + imd_quintiles + dm_diag_age + index_bmi + index_hdl + index_totalcholesterol + index_triglyceride + index_hba1c | t1_code_ever, data=type_2_ins_3yrs_cohort_misclass_miscod, overall=F, render=rndr, render.categorical=cat, render.continuous=cont)

## Compare groups
type_2_ins_3yrs_cohort %>% distinct(ins_3_years)
pvals <- data.frame(p=sapply(vars, function(v) get_pvalue(v, "ins_3_years", cohort))) %>%
  mutate(sig=ifelse(p<0.05, 1, NA))
pvals

type_2_misclass_ins <- local_vars %>% filter(diabetes_type=="type 2" & t1_code_ever==0 & !is.na(ins_3_years) & ins_3_years==1) %>% mutate(group="type_2_misclass")
outcomes_t2_cohort <- rbind(type_2_misclass_ins, ref_type_2)

outcomes_t2_ins3yrs(outcomes_t2_cohort, group_name="overall_ins3yrs")

ref_type_2 %>% count() #25615
ref_type_2 %>% filter(!is.na(ins_3_years)) %>% count() #8481
8481/25615
#33.1%
ref_type_2 %>% filter(!is.na(ins_3_years) & ins_3_years==1) %>% count() #2076
2076/8481
#24.5%




# Type 2 and on insulin within 1 year

setwd("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2026/DePICtion/Scripts/Functions/")
source("outcomes.R")



local_vars %>% filter(diabetes_type=="type 2") %>% count()
local_vars %>% filter(diabetes_type=="type 2" & !is.na(ins_1_year)) %>% count()
9835/43180 #22.8%

type_2_ins_1yr_cohort <- local_vars %>%
  filter(diabetes_type=="type 2" & !is.na(ins_1_year)) %>%
  mutate(model_prob_over_80=factor(ifelse(lipid_prob_perc>=80, 1, 0)))


type_2_ins_1yr_cohort %>% filter(ins_1_year==1) %>% count()
#2150
2150/9835 #21.9%

type_2_ins_1yr_cohort %>% filter(ins_1_year==1 & t1_code_ever==0) %>% count()
#1854/9835 #21.9%

type_2_ins_1yr_cohort %>% filter(ins_1_year==1 & t1_code_ever==1) %>% count()
#296


table1(~ lipid_prob_perc + model_prob_over_80 + malesex + age_at_index + ethnicity_decoded + imd_quintiles + dm_diag_age + index_bmi + index_hdl + index_totalcholesterol + index_triglyceride + index_hba1c + t1_code_ever | ins_1_year, data=type_2_ins_1yr_cohort, overall=F, render=rndr, render.categorical=cat, render.continuous=cont)

type_2_ins_1yr_cohort_misclass_miscod <- type_2_ins_1yr_cohort %>% filter(ins_1_year==1)

table1(~ lipid_prob_perc + model_prob_over_80 + malesex + age_at_index + ethnicity_decoded + imd_quintiles + dm_diag_age + index_bmi + index_hdl + index_totalcholesterol + index_triglyceride + index_hba1c | t1_code_ever, data=type_2_ins_1yr_cohort_misclass_miscod, overall=F, render=rndr, render.categorical=cat, render.continuous=cont)


type_2_misclass_ins <- local_vars %>% filter(diabetes_type=="type 2" & t1_code_ever==0 & !is.na(ins_1_year) & ins_1_year==1) %>% mutate(group="type_2_misclass")
outcomes_t2_cohort <- rbind(type_2_misclass_ins, ref_type_2)

outcomes_t2_ins1yr(outcomes_t2_cohort, group_name="overall_ins1yr")

ref_type_2 %>% count() #25615
ref_type_2 %>% filter(!is.na(ins_3_years)) %>% count() #8481
8481/25615
#33.1%
ref_type_2 %>% filter(!is.na(ins_3_years) & ins_3_years==1) %>% count() #2076
2076/8481
#24.5%

############################################################################################

# Current type 1

## Misclassified/miscoded counts and predictors

t1_table_cohort <- rbind(type_1_all, type_1_misclass_miscode, type_1_misclass, ref_type_1)

table1(~ lipid_prob_perc + malesex + age_at_index + ethnicity_decoded + imd_quintiles + dm_diag_age + index_bmi + index_hdl + index_totalcholesterol + index_triglyceride + index_hba1c + t2_code_ever | group, data=t1_table_cohort, overall=F, render=rndr, render.categorical=cat, render.continuous=cont)

# Look at ethnicity-deprivation relationship
table1(~ imd_quintiles | ethnicity_decoded, data=type_1_all, overall=F, render=rndr, render.categorical=cat, render.continuous=cont)
#Overall, South Asian, Black, Mixed and Other all have more deprivation that White - true for type 2 too


### Look at pairwise p-values

get_pvalue <- function(var, group, data) {
  x <- data[[var]]
  g <- data[[group]]
  
  if (is.numeric(x)) {
    p <- kruskal.test(x ~ g)$p.value
  } else if (is.factor(x) || is.character(x)) {
    tbl <- table(x, g)
    p <- chisq.test(tbl)$p.value
  } else {
    return(NA)
  }
  return(p)
}

vars <- c("lipid_prob_perc", "malesex", "age_at_index", "ethnicity_decoded", "imd_quintiles", "dm_diag_age", "index_bmi", "index_hdl", "index_totalcholesterol", "index_triglyceride", "index_hba1c", "t2_code_ever")


### Compare flagged to overall
cohort <- t1_table_cohort %>% filter(group!="type_1_misclass" & group!="ref_type_1")
cohort %>% distinct(group)
pvals <- data.frame(p=sapply(vars, function(v) get_pvalue(v, "group", cohort))) %>%
  mutate(sig=ifelse(p<0.05, 1, NA))
pvals
# all are significantly different except total cholesterol

## Compare flagged to misclassified
cohort <- t1_table_cohort %>% filter(group!="ref_type_1" & group!="type_1_all")
cohort %>% distinct(group)
pvals <- data.frame(p=sapply(vars, function(v) get_pvalue(v, "group", cohort))) %>%
  mutate(sig=ifelse(p<0.05, 1, NA))
pvals
# None are sig except current age and diagnosis age and T2 code ever

## Compare reference to all
cohort <- t1_table_cohort %>% filter(group!="type_1_misclass" & group!="type_1_misclass_miscode")
cohort %>% distinct(group)
pvals <- data.frame(p=sapply(vars, function(v) get_pvalue(v, "group", cohort))) %>%
  mutate(sig=ifelse(p<0.05, 1, NA))
pvals
# all are sig except total cholesterol


t1_cohort_by_group <- rbind((type_1_all %>% mutate(new_group="overall")),
                            (type_1_all %>% filter(malesex==1) %>% mutate(new_group="male")),
                            (type_1_all %>% filter(malesex==0) %>% mutate(new_group="female")),
                            (type_1_all %>% filter(ethnicity_decoded=="White") %>% mutate(new_group="White")),
                            (type_1_all %>% filter(ethnicity_decoded=="South Asian") %>% mutate(new_group="South Asian")),
                            (type_1_all %>% filter(ethnicity_decoded=="Black") %>% mutate(new_group="Black")),
                            (type_1_all %>% filter(ethnicity_decoded=="Mixed") %>% mutate(new_group="Mixed")),
                            (type_1_all %>% filter(ethnicity_decoded=="Other") %>% mutate(new_group="Other")),
                            (type_1_all %>% filter(ethnicity_decoded=="Missing") %>% mutate(new_group="Missing")),
                            (type_1_all %>% filter(imd_quintiles==1) %>% mutate(new_group="imd_1")),
                            (type_1_all %>% filter(imd_quintiles==2) %>% mutate(new_group="imd_2")),
                            (type_1_all %>% filter(imd_quintiles==3) %>% mutate(new_group="imd_3")),
                            (type_1_all %>% filter(imd_quintiles==4) %>% mutate(new_group="imd_4")),
                            (type_1_all %>% filter(imd_quintiles==5) %>% mutate(new_group="imd_5"))) %>%
  mutate(new_group=factor(new_group, levels=c("overall", "male", "female", "White", "South Asian", "Black", "Mixed", "Other", "Missing", "imd_1", "imd_2", "imd_3", "imd_4", "imd_5")),
         flagged=as.factor(ifelse(lipid_prob_perc<5, 1, 0)),
         misclass=as.factor(ifelse(flagged==1 & t2_code_ever==0, 1, 0)))

# With CIs for proportions
table1(~ flagged + misclass | new_group, data=t1_cohort_by_group, overall=F, render=rndr, render.categorical=render_cat_ci(digits=1), render.continuous=cont)

# Without CIs for proportions
table1(~ lipid_prob_perc + malesex + imd_quintiles + dm_diag_age + index_bmi + index_hdl + index_totalcholesterol + index_triglyceride | new_group, data=t1_cohort_by_group, overall=F, render=rndr, render.categorical=cat, render.continuous=cont)

flagged_t1_cohort_by_group <- t1_cohort_by_group %>% filter(flagged==1)

table1(~ t2_code_ever | new_group, data=flagged_t1_cohort_by_group, overall=F, render=rndr, render.categorical=render_cat_ci(digits=1), render.continuous=cont)





t1_cohort_by_group_plot <- t1_cohort_by_group %>% mutate(miscode=as.factor(ifelse(flagged==1 & misclass==0, 1, 0)))

t1_cohort_by_group_plot <- t1_cohort_by_group_plot %>%
  filter(new_group!="Mixed" & new_group!="Other" & new_group!="Missing" & new_group!="imd_2" & new_group!="imd_3" & new_group!="imd_4") %>%
  mutate(new_group=case_when(new_group=="male" ~ "Male",
                             new_group=="female" ~ "Female",
                             new_group=="overall" ~ "Overall",
                             new_group=="imd_1" ~ "IMD quintile 1",
                             new_group=="imd_5" ~ "IMD quintile 5",
                             TRUE ~ new_group))


t1_cohort_by_group_plot$new_group <- factor(t1_cohort_by_group_plot$new_group,
                                            levels = c("IMD quintile 5", "IMD quintile 1","Black","South Asian",  "White", "Female","Male", "Overall"))


plot_df <- t1_cohort_by_group_plot %>%
  group_by(new_group) %>%
  summarise(
    flagged_n = sum(flagged == 1, na.rm = TRUE),
    misclass_n = sum(misclass == 1, na.rm = TRUE),
    miscode_n = sum(miscode == 1, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

library(binom)

flagged_ci  <- binom.confint(plot_df$flagged_n,  plot_df$n, method = "wilson")
misclass_ci <- binom.confint(plot_df$misclass_n, plot_df$n, method = "wilson")
miscode_ci  <- binom.confint(plot_df$miscode_n,  plot_df$n, method = "wilson")


plot_df <- plot_df %>%
  mutate(
    flagged_est  = flagged_ci$mean,
    flagged_low  = flagged_ci$lower,
    flagged_upp  = flagged_ci$upper,
    
    misclass_est = misclass_ci$mean,
    misclass_low = misclass_ci$lower,
    misclass_upp = misclass_ci$upper,
    
    miscode_est  = miscode_ci$mean,
    miscode_low  = miscode_ci$lower,
    miscode_upp  = miscode_ci$upper
  )

plot_long <- plot_df %>%
  select(new_group,
         flagged_est, flagged_low, flagged_upp,
         misclass_est, misclass_low, misclass_upp,
         miscode_est, miscode_low, miscode_upp) %>%
  pivot_longer(
    cols = -new_group,
    names_to = c("outcome", ".value"),
    names_pattern = "(flagged|misclass|miscode)_(est|low|upp)"
  )

plot_long <- plot_long %>%
  mutate(outcome = recode(outcome,
                          flagged  = "Potentially misclassified or miscoded (model=T2)",
                          misclass = "Potentially misclassified (model=T2, no history of T2 codes)",
                          miscode  = "Potentially miscoded (model=T2, history of T2 codes)"
  ))

plot_long$outcome <- factor(plot_long$outcome, levels=c("Potentially misclassified or miscoded (model=T2)", "Potentially misclassified (model=T2, no history of T2 codes)", "Potentially miscoded (model=T2, history of T2 codes)"))



x_max <- 20
arrow_length <- 0.03  # small offset inside axis

# adjust truncated upper CI so arrow fits inside
plot_long <- plot_long %>%
  mutate(
    upp_trunc2 = ifelse(upp * 100 > x_max, x_max - arrow_length, upp * 100),
    arrow_needed = (upp * 100 > x_max),
    est_pct = est * 100,
    label_text = paste0(round(est_pct, 1), "%"),
    row_id = as.numeric(factor(new_group)),                # numeric row
    stripe = factor(row_id %% 2)                            # 0/1 as factor
  )

label_x <- -0.02 * x_max

# plotting

tiff("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2026/DePICtion/Plots/t1_proportion_subgroup.tiff",
     width = 10, height = 9, units = "in", res = 800)

ggplot(plot_long, aes(x = est_pct, y = new_group)) +
  # zebra stripes
  geom_rect(
    data = plot_long %>% distinct(new_group, row_id, stripe),
    aes(
      ymin = row_id - 0.5,
      ymax = row_id + 0.5,
      xmin = -Inf,
      xmax = Inf,
      fill = stripe
    ),
    inherit.aes = FALSE,
    alpha = 0.1
  ) +
  scale_fill_manual(values = c("white", "grey70"), guide = "none") +  # now works
  # points
  geom_point() +
  # error bars
  geom_errorbarh(aes(xmin = low * 100, xmax = upp_trunc2), height = 0) +
  # arrows for truncated CIs
  geom_segment(
    data = subset(plot_long, arrow_needed),
    aes(
      x = x_max - arrow_length,
      xend = x_max,
      y = new_group,
      yend = new_group
    ),
    arrow = arrow(length = unit(0.15, "cm"), ends = "last", type = "closed"),
    inherit.aes = FALSE
  ) +
  # left-aligned % labels
  geom_text(aes(y = new_group, label = label_text),
            x = label_x,
            hjust = 0,
            size = 5) +
  facet_wrap(~ outcome, ncol = 1) +
  xlim(label_x, x_max) +
  labs(x = "Percentage (%)", y = NULL) +
  theme_classic(base_size = 16) +
  theme(axis.text.y = element_text(hjust = 1))

dev.off()






## Outcomes

setwd("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2026/DePICtion/Scripts/Functions/")
source("outcomes.R")

outcomes_t1_cohort <- rbind(type_1_misclass, ref_type_1)

outcomes_t1_cohort_by_group <- rbind((outcomes_t1_cohort %>% mutate(new_group="overall")),
                                     (outcomes_t1_cohort %>% filter(malesex==1) %>% mutate(new_group="male")),
                                     (outcomes_t1_cohort %>% filter(malesex==0) %>% mutate(new_group="female")),
                                     (outcomes_t1_cohort %>% filter(ethnicity_decoded=="White") %>% mutate(new_group="White")),
                                     (outcomes_t1_cohort %>% filter(ethnicity_decoded=="South Asian") %>% mutate(new_group="South Asian")),
                                     (outcomes_t1_cohort %>% filter(ethnicity_decoded=="Black") %>% mutate(new_group="Black")),
                                     (outcomes_t1_cohort %>% filter(ethnicity_decoded=="Mixed") %>% mutate(new_group="Mixed")),
                                     (outcomes_t1_cohort %>% filter(ethnicity_decoded=="Other") %>% mutate(new_group="Other")),
                                     (outcomes_t1_cohort %>% filter(ethnicity_decoded=="Missing") %>% mutate(new_group="Missing")),
                                     (outcomes_t1_cohort %>% filter(imd_quintiles==1) %>% mutate(new_group="imd_1")),
                                     (outcomes_t1_cohort %>% filter(imd_quintiles==2) %>% mutate(new_group="imd_2")),
                                     (outcomes_t1_cohort %>% filter(imd_quintiles==3) %>% mutate(new_group="imd_3")),
                                     (outcomes_t1_cohort %>% filter(imd_quintiles==4) %>% mutate(new_group="imd_4")),
                                     (outcomes_t1_cohort %>% filter(imd_quintiles==5) %>% mutate(new_group="imd_5")))

new_groups <- c("overall", "male", "female", "White", "South Asian", "Black", "imd_1", "imd_5")

for (i in new_groups) {
  
  data <- outcomes_t1_cohort_by_group %>% filter(new_group==i)
  
  outcomes_t1(data, i)
  
}

compare_t1_sex(outcomes_t1_cohort)

compare_t1_ethnicity(outcomes_t1_cohort)

compare_t1_deprivation(outcomes_t1_cohort)



# Check whether BMI is making difference for IMD groups

type_1_all %>% filter(imd_quintiles==1) %>% summarise(mean_bmi=mean(index_bmi),  #27.3
                                                      flagged=sum(lipid_prob_perc<5),
                                                      count=n())

type_1_all %>% filter(imd_quintiles==5) %>% summarise(mean_bmi=mean(index_bmi), #28.1
                                                      flagged=sum(lipid_prob_perc<5),
                                                      count=n())

type_1_all %>%
  filter(imd_quintiles==1) %>%
  mutate(index_bmi=27.3,
         femalesex=ifelse(gender==2, 1, ifelse(gender==1, 0, NA)),
         lipid_pred_score=9.0034272-(0.1915482*index_bmi)-(0.1686227*dm_diag_age)+(0.3026012*femalesex)-(0.2269216*index_totalcholesterol)+(1.540850*index_hdl)-(0.2784059*index_triglyceride),
         lipid_pred_prob=exp(lipid_pred_score)/(1+exp(lipid_pred_score)),
         lipid_prob_perc=lipid_pred_prob*100) %>%
  summarise(flagged=sum(lipid_prob_perc<5),
            count=n())
#328/3896 - many fewer flagged, so could be



# check autoimmune diffs
outcomes_t1_cohort %>%
  summarise(p=prop.test(
    c(sum(coeliac_thyroid[group=="type_1_misclass"]==1), sum(coeliac_thyroid[group=="ref_type_1"]==1)),
    c(sum(group=="type_1_misclass"   ), sum(group=="ref_type_1"    ))
  )$p.value) %>%
  pull(p)
#p=0.8


############################################################################################

# Flagged by different thresholds

#4.6 million people in the UK with diabetes out of 68.3 million
# Our download: 984,166 with diabetes
# So to get total out of population: divide by (984,166*(68.3/4.6))
(984166*(68.3/4.6)) #14612726

model_cohort <- rbind(type_1_all, type_2_all)


## Type 2

### Flagged

# High scorers
model_cohort %>% filter(diabetes_type=="type 2" & !is.na(lipid_pred_prob)) %>% count()
#43180
model_cohort %>% filter(diabetes_type=="type 2" & !is.na(lipid_pred_prob) & lipid_pred_prob>=0.8) %>% count()
#406
406*100/43180
(406/14612726)*10000
model_cohort %>% filter(diabetes_type=="type 2" & !is.na(lipid_pred_prob) & lipid_pred_prob>=0.7) %>% count()
#759
759*100/43180
(759/14612726)*10000
model_cohort %>% filter(diabetes_type=="type 2" & !is.na(lipid_pred_prob) & lipid_pred_prob>=0.6) %>% count()
#1239
1239*100/43180
(1239/14612726)*10000
model_cohort %>% filter(diabetes_type=="type 2" & !is.na(lipid_pred_prob) & lipid_pred_prob>=0.5) %>% count()
#1853
1853*100/43180
(1853/14612726)*10000
model_cohort %>% filter(diabetes_type=="type 2" & !is.na(lipid_pred_prob) & lipid_pred_prob>=0.4) %>% count()
#2724
2724*100/43180
(2724/14612726)*10000
model_cohort %>% filter(diabetes_type=="type 2" & !is.na(ins_1_year) & ins_1_year==1) %>% count()
#2150
2150*100/43180
(2150/14612726)*10000
model_cohort %>% filter(diabetes_type=="type 2" & !is.na(ins_3_years) & ins_3_years==1) %>% count()
#3778
3778*100/43180
(3778/14612726)*10000


# Range per practice
summary((model_cohort %>% group_by(pracid) %>% summarise(count=sum(diabetes_type=="type 2" & !is.na(lipid_pred_prob) & lipid_pred_prob>=0.8)) %>% collect())$count)

summary((model_cohort %>% group_by(pracid) %>% summarise(count=sum(diabetes_type=="type 2" & !is.na(lipid_pred_prob) & lipid_pred_prob>=0.7)) %>% collect())$count)

summary((model_cohort %>% group_by(pracid) %>% summarise(count=sum(diabetes_type=="type 2" & !is.na(lipid_pred_prob) & lipid_pred_prob>=0.6)) %>% collect())$count)

summary((model_cohort %>% group_by(pracid) %>% summarise(count=sum(diabetes_type=="type 2" & !is.na(lipid_pred_prob) & lipid_pred_prob>=0.5)) %>% collect())$count)

summary((model_cohort %>% group_by(pracid) %>% summarise(count=sum(diabetes_type=="type 2" & !is.na(lipid_pred_prob) & lipid_pred_prob>=0.4)) %>% collect())$count)

summary((model_cohort %>% group_by(pracid) %>% summarise(count=sum(diabetes_type=="type 2" & !is.na(ins_3_years) & ins_3_years==1)) %>% collect())$count)




### Misclassified

# High scorers
model_cohort %>% filter(diabetes_type=="type 2" & !is.na(lipid_pred_prob)) %>% count()
#43180
model_cohort %>% filter(diabetes_type=="type 2" & !is.na(lipid_pred_prob) & lipid_pred_prob>=0.8 & t1_code_ever==0) %>% count()
#248
248*100/43180
(248/14612726)*10000
model_cohort %>% filter(diabetes_type=="type 2" & !is.na(lipid_pred_prob) & lipid_pred_prob>=0.7 & t1_code_ever==0) %>% count()
#513
513*100/43180
(513/14612726)*10000
model_cohort %>% filter(diabetes_type=="type 2" & !is.na(lipid_pred_prob) & lipid_pred_prob>=0.6 & t1_code_ever==0) %>% count()
#876
876*100/43180
(876/14612726)*10000
model_cohort %>% filter(diabetes_type=="type 2" & !is.na(lipid_pred_prob) & lipid_pred_prob>=0.5 & t1_code_ever==0) %>% count()
#1372
1372*100/43180
(1372/14612726)*10000
model_cohort %>% filter(diabetes_type=="type 2" & !is.na(lipid_pred_prob) & lipid_pred_prob>=0.4 & t1_code_ever==0) %>% count()
#2111
2111*100/43180
(2111/14612726)*10000
model_cohort %>% filter(diabetes_type=="type 2" & !is.na(ins_1_year) & ins_1_year==1 & t1_code_ever==0) %>% count()
#1854
1854*100/43180
(1854/14612726)*10000
model_cohort %>% filter(diabetes_type=="type 2" & !is.na(ins_3_years) & ins_3_years==1 & t1_code_ever==0) %>% count()
#3778
3311*100/43180
(3311/14612726)*10000


# Range per practice
summary((model_cohort %>% group_by(pracid) %>% summarise(count=sum(diabetes_type=="type 2" & !is.na(lipid_pred_prob) & lipid_pred_prob>=0.8 & t1_code_ever==0)) %>% collect())$count)

summary((model_cohort %>% group_by(pracid) %>% summarise(count=sum(diabetes_type=="type 2" & !is.na(lipid_pred_prob) & lipid_pred_prob>=0.7 & t1_code_ever==0)) %>% collect())$count)

summary((model_cohort %>% group_by(pracid) %>% summarise(count=sum(diabetes_type=="type 2" & !is.na(lipid_pred_prob) & lipid_pred_prob>=0.6 & t1_code_ever==0)) %>% collect())$count)

summary((model_cohort %>% group_by(pracid) %>% summarise(count=sum(diabetes_type=="type 2" & !is.na(lipid_pred_prob) & lipid_pred_prob>=0.5 & t1_code_ever==0)) %>% collect())$count)

summary((model_cohort %>% group_by(pracid) %>% summarise(count=sum(diabetes_type=="type 2" & !is.na(lipid_pred_prob) & lipid_pred_prob>=0.4 & t1_code_ever==0)) %>% collect())$count)

summary((model_cohort %>% group_by(pracid) %>% summarise(count=sum(diabetes_type=="type 2" & !is.na(ins_3_years) & ins_3_years==1 & t1_code_ever==0)) %>% collect())$count)






### Outcomes for misclassified

setwd("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2026/DePICtion/Scripts/Functions/")
source("outcomes_thresholds.R")

outcomes_t2_thresholds_cohort <- rbind(ref_type_2,
                                       (model_cohort %>% filter(diabetes_type=="type 2" & !is.na(lipid_pred_prob) & lipid_pred_prob>=0.8 & t1_code_ever==0) %>% mutate(group="type_2_misclass_0.8")),
                                       (model_cohort %>% filter(diabetes_type=="type 2" & !is.na(lipid_pred_prob) & lipid_pred_prob>=0.7 & t1_code_ever==0) %>% mutate(group="type_2_misclass_0.7")),
                                       (model_cohort %>% filter(diabetes_type=="type 2" & !is.na(lipid_pred_prob) & lipid_pred_prob>=0.6 & t1_code_ever==0) %>% mutate(group="type_2_misclass_0.6")),
                                       (model_cohort %>% filter(diabetes_type=="type 2" & !is.na(lipid_pred_prob) & lipid_pred_prob>=0.5 & t1_code_ever==0) %>% mutate(group="type_2_misclass_0.5")),
                                       (model_cohort %>% filter(diabetes_type=="type 2" & !is.na(lipid_pred_prob) & lipid_pred_prob>=0.4 & t1_code_ever==0) %>% mutate(group="type_2_misclass_0.4")))

outcomes_thresholds_t2(outcomes_t2_thresholds_cohort)





## Type 1

### Flagged

# High scorers
model_cohort %>% filter(diabetes_type=="type 1" & !is.na(lipid_pred_prob)) %>% count()
#18836
model_cohort %>% filter(diabetes_type=="type 1" & !is.na(lipid_pred_prob) & lipid_pred_prob<0.025) %>% count()
#1492
1492*100/18836
(1492/14612726)*10000
model_cohort %>% filter(diabetes_type=="type 1" & !is.na(lipid_pred_prob) & lipid_pred_prob<0.05) %>% count()
#2570
2570*100/18836
(2570/14612726)*10000
model_cohort %>% filter(diabetes_type=="type 1" & !is.na(lipid_pred_prob) & lipid_pred_prob<0.1) %>% count()
#4082
4082*100/18836
(4082/14612726)*10000



# Range per practice
summary((model_cohort %>% group_by(pracid) %>% summarise(count=sum(diabetes_type=="type 1" & !is.na(lipid_pred_prob) & lipid_pred_prob<0.025)) %>% collect())$count)

summary((model_cohort %>% group_by(pracid) %>% summarise(count=sum(diabetes_type=="type 1" & !is.na(lipid_pred_prob) & lipid_pred_prob<0.05)) %>% collect())$count)

summary((model_cohort %>% group_by(pracid) %>% summarise(count=sum(diabetes_type=="type 1" & !is.na(lipid_pred_prob) & lipid_pred_prob<0.1)) %>% collect())$count)




### Misclassified

# High scorers
model_cohort %>% filter(diabetes_type=="type 1" & !is.na(lipid_pred_prob)) %>% count()
#18836
model_cohort %>% filter(diabetes_type=="type 1" & !is.na(lipid_pred_prob) & lipid_pred_prob<0.025 & t2_code_ever==0) %>% count()
#827
827*100/18836
(827/14612726)*10000
model_cohort %>% filter(diabetes_type=="type 1" & !is.na(lipid_pred_prob) & lipid_pred_prob<0.05 & t2_code_ever==0) %>% count()
#1460
1460*100/18836
(1460/14612726)*10000
model_cohort %>% filter(diabetes_type=="type 1" & !is.na(lipid_pred_prob) & lipid_pred_prob<0.1 & t2_code_ever==0) %>% count()
#2396
2396*100/18836
(2396/14612726)*10000

# Range per practice
summary((model_cohort %>% group_by(pracid) %>% summarise(count=sum(diabetes_type=="type 1" & !is.na(lipid_pred_prob) & lipid_pred_prob<0.025 & t2_code_ever==0)) %>% collect())$count)

summary((model_cohort %>% group_by(pracid) %>% summarise(count=sum(diabetes_type=="type 1" & !is.na(lipid_pred_prob) & lipid_pred_prob<0.05 & t2_code_ever==0)) %>% collect())$count)

summary((model_cohort %>% group_by(pracid) %>% summarise(count=sum(diabetes_type=="type 1" & !is.na(lipid_pred_prob) & lipid_pred_prob<0.1 & t2_code_ever==0)) %>% collect())$count)




### Outcomes for misclassified

setwd("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2026/DePICtion/Scripts/Functions/")
source("outcomes_thresholds.R")

outcomes_t1_thresholds_cohort <- rbind(ref_type_1,
                                       (model_cohort %>% filter(diabetes_type=="type 1" & !is.na(lipid_pred_prob) & lipid_pred_prob<0.025 & t2_code_ever==0) %>% mutate(group="type_1_misclass_0.025")),
                                       (model_cohort %>% filter(diabetes_type=="type 1" & !is.na(lipid_pred_prob) & lipid_pred_prob<0.05 & t2_code_ever==0) %>% mutate(group="type_1_misclass_0.05")),
                                       (model_cohort %>% filter(diabetes_type=="type 1" & !is.na(lipid_pred_prob) & lipid_pred_prob<0.1 & t2_code_ever==0) %>% mutate(group="type_1_misclass_0.1")))

outcomes_thresholds_t1(outcomes_t1_thresholds_cohort)



###############################################

# Data quality checking

# prescription or Hba1c earlier

cohort <- cohort %>% analysis$cached("cohort")

analysis = cprd$analysis("all_patid")
hba1c <- hba1c %>% analysis$cached("clean_hba1c_medcodes")
earliest_high_hba1c <- hba1c %>% filter(testvalue>=48) %>% group_by(patid) %>% summarise(earliest_high_hba1c=min(date, na.rm=TRUE)) %>% ungroup()

analysis = cprd$analysis("dpctn_paper_prev")

test <- cohort %>%
  left_join(earliest_high_hba1c, by="patid") %>%
  mutate(earliest_ins_oha_hba1c=pmin(ifelse(is.na(earliest_high_hba1c), as.Date("2050-01-01"), earliest_high_hba1c),
                                     ifelse(is.na(earliest_ins), as.Date("2050-01-01"), earliest_ins),
                                     ifelse(is.na(earliest_oha), as.Date("2050-01-01"), earliest_oha), na.rm=TRUE),
         earliest_ins_oha_hba1c=ifelse(earliest_ins_oha_hba1c==as.Date("2050-01-01"), as.Date(NA), earliest_ins_oha_hba1c)) %>%
  select(patid, earliest_high_hba1c, earliest_ins, earliest_oha, earliest_ins_oha_hba1c, diagnosis_date) %>%
  collect()

test %>% filter(!is.na(earliest_ins_oha_hba1c) & earliest_ins_oha_hba1c<diagnosis_date) %>% count()
10837/62013

summary((test %>% filter(!is.na(earliest_ins_oha_hba1c) & earliest_ins_oha_hba1c<diagnosis_date) %>% mutate(datediff=difftime(diagnosis_date, earliest_ins_oha_hba1c, units="days")))$datediff)

test %>% filter(!is.na(earliest_ins_oha_hba1c) & earliest_ins_oha_hba1c<diagnosis_date) %>% mutate(datediff=as.numeric(difftime(diagnosis_date, earliest_ins_oha_hba1c, units="days"))) %>% filter(datediff>365.25) %>% count()
                                                                                                   
                                                                                                   

# earliest in HES records

earliest_hes_diabetes <- cprd$tables$hesDiagnosisEpi %>%
  filter(ICD %like% "E10%" | ICD %like% "E11%" | ICD %like% "E12%" | ICD %like% "E13%" | ICD %like% "E14%") %>%
  group_by(patid) %>%
  summarise(earliest_hes_diabetes=min(epistart, na.rm=TRUE)) %>%
  collect()

cohort <- cohort %>% select(patid, diagnosis_date) %>% collect() %>% inner_join(earliest_hes_diabetes, by="patid")
#54,112

cohort %>% filter(earliest_hes_diabetes<diagnosis_date) %>% count()
#4353

summary((cohort %>% filter(earliest_hes_diabetes<diagnosis_date) %>% mutate(time=as.numeric(difftime(diagnosis_date, earliest_hes_diabetes, units="days"))) %>% select(patid, time) %>% collect)$time)
#Median 13 days(max 8576)


test2 <- test %>% mutate(earliest_ins_oha_code_hba1c=pmin(earliest_ins_oha_hba1c, diagnosis_date, na.rm=TRUE)) %>% inner_join(earliest_hes_diabetes, by="patid")

test2 %>% filter(earliest_hes_diabetes<earliest_ins_oha_code_hba1c) %>% count()
#3748
3748/62013

summary((test2 %>% filter(earliest_hes_diabetes<earliest_ins_oha_code_hba1c) %>% mutate(time=as.numeric(difftime(earliest_ins_oha_code_hba1c, earliest_hes_diabetes, units="days"))) %>% select(patid, time))$time)
#Median 13 days(max 8482)
8482/365.25



###############################################

#C-peptide and antibodies

cohort <- cohort %>% analysis$cached("cohort")

cohort %>% filter(!is.na(earliest_c_pep_ins_deficient) | !is.na(earliest_c_pep_ins_intermediate) | !is.na(earliest_c_pep_ins_normal)) %>% count()
#1,058

cohort %>% filter(!is.na(earliest_positive_gad) | !is.na(earliest_negative_gad) | !is.na(earliest_positive_ia2) | !is.na(earliest_negative_ia2)) %>% count()
#762

cohort %>% filter(!is.na(earliest_c_pep_ins_deficient) | !is.na(earliest_c_pep_ins_intermediate) | !is.na(earliest_c_pep_ins_normal) | !is.na(earliest_positive_gad) | !is.na(earliest_negative_gad) | !is.na(earliest_positive_ia2) | !is.na(earliest_negative_ia2)) %>% count()
#1602





