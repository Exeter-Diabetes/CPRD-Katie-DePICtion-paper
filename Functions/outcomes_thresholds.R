
# Output incidence rates for hospitalisation for DKA and hypoglycaemia, and %s on insulin at 1 year post diagnosis and on basal-bolus insulin at 10 years post diagnosis
## And test for differences between potentially misclassified and reference groups

outcomes_thresholds_t2 <- function(cohort_dataset) {
  
  p_val_dka_0.8 <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(dka_post_diagnosis[group == "type_2_misclass_0.8"] == 1), sum(dka_post_diagnosis[group == "ref_type_2"] == 1)),
      c(sum(follow_up_time[group == "type_2_misclass_0.8"]), sum(follow_up_time[group == "ref_type_2"]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_dka_0.7 <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(dka_post_diagnosis[group == "type_2_misclass_0.7"] == 1), sum(dka_post_diagnosis[group == "ref_type_2"] == 1)),
      c(sum(follow_up_time[group == "type_2_misclass_0.7"]), sum(follow_up_time[group == "ref_type_2"]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_dka_0.6 <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(dka_post_diagnosis[group == "type_2_misclass_0.6"] == 1), sum(dka_post_diagnosis[group == "ref_type_2"] == 1)),
      c(sum(follow_up_time[group == "type_2_misclass_0.6"]), sum(follow_up_time[group == "ref_type_2"]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_dka_0.5 <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(dka_post_diagnosis[group == "type_2_misclass_0.5"] == 1), sum(dka_post_diagnosis[group == "ref_type_2"] == 1)),
      c(sum(follow_up_time[group == "type_2_misclass_0.5"]), sum(follow_up_time[group == "ref_type_2"]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_dka_0.4 <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(dka_post_diagnosis[group == "type_2_misclass_0.4"] == 1), sum(dka_post_diagnosis[group == "ref_type_2"] == 1)),
      c(sum(follow_up_time[group == "type_2_misclass_0.4"]), sum(follow_up_time[group == "ref_type_2"]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  
  p_val_hypo_0.8 <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(hypo_post_diagnosis[group == "type_2_misclass_0.8"] == 1), sum(hypo_post_diagnosis[group == "ref_type_2"] == 1)),
      c(sum(follow_up_time[group == "type_2_misclass_0.8"]), sum(follow_up_time[group == "ref_type_2"]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_hypo_0.7 <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(hypo_post_diagnosis[group == "type_2_misclass_0.7"] == 1), sum(hypo_post_diagnosis[group == "ref_type_2"] == 1)),
      c(sum(follow_up_time[group == "type_2_misclass_0.7"]), sum(follow_up_time[group == "ref_type_2"]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_hypo_0.6 <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(hypo_post_diagnosis[group == "type_2_misclass_0.6"] == 1), sum(hypo_post_diagnosis[group == "ref_type_2"] == 1)),
      c(sum(follow_up_time[group == "type_2_misclass_0.5"]), sum(follow_up_time[group == "ref_type_2"]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_hypo_0.5 <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(hypo_post_diagnosis[group == "type_2_misclass_0.5"] == 1), sum(hypo_post_diagnosis[group == "ref_type_2"] == 1)),
      c(sum(follow_up_time[group == "type_2_misclass_0.5"]), sum(follow_up_time[group == "ref_type_2"]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_hypo_0.4 <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(hypo_post_diagnosis[group == "type_2_misclass_0.4"] == 1), sum(hypo_post_diagnosis[group == "ref_type_2"] == 1)),
      c(sum(follow_up_time[group == "type_2_misclass_0.4"]), sum(follow_up_time[group == "ref_type_2"]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  
  p_val_ins_1_year_0.8 <- cohort_dataset %>%
    filter(!is.na(ins_1_year)) %>%
    summarise(p = ifelse(
      sum(group == "type_2_misclass_0.8") == 0 | sum(group == "ref_type_2") == 0,
      NA,
      prop.test(
        c(sum(ins_1_year[group == "type_2_misclass_0.8"] == 1), sum(ins_1_year[group == "ref_type_2"] == 1)),
        c(sum(group == "type_2_misclass_0.8"), sum(group == "ref_type_2"))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_ins_1_year_0.7 <- cohort_dataset %>%
    filter(!is.na(ins_1_year)) %>%
    summarise(p = ifelse(
      sum(group == "type_2_misclass_0.7") == 0 | sum(group == "ref_type_2") == 0,
      NA,
      prop.test(
        c(sum(ins_1_year[group == "type_2_misclass_0.7"] == 1), sum(ins_1_year[group == "ref_type_2"] == 1)),
        c(sum(group == "type_2_misclass_0.7"), sum(group == "ref_type_2"))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_ins_1_year_0.6 <- cohort_dataset %>%
    filter(!is.na(ins_1_year)) %>%
    summarise(p = ifelse(
      sum(group == "type_2_misclass_0.6") == 0 | sum(group == "ref_type_2") == 0,
      NA,
      prop.test(
        c(sum(ins_1_year[group == "type_2_misclass_0.6"] == 1), sum(ins_1_year[group == "ref_type_2"] == 1)),
        c(sum(group == "type_2_misclass_0.6"), sum(group == "ref_type_2"))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_ins_1_year_0.5 <- cohort_dataset %>%
    filter(!is.na(ins_1_year)) %>%
    summarise(p = ifelse(
      sum(group == "type_2_misclass_0.5") == 0 | sum(group == "ref_type_2") == 0,
      NA,
      prop.test(
        c(sum(ins_1_year[group == "type_2_misclass_0.5"] == 1), sum(ins_1_year[group == "ref_type_2"] == 1)),
        c(sum(group == "type_2_misclass_0.5"), sum(group == "ref_type_2"))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_ins_1_year_0.4 <- cohort_dataset %>%
    filter(!is.na(ins_1_year)) %>%
    summarise(p = ifelse(
      sum(group == "type_2_misclass_0.4") == 0 | sum(group == "ref_type_2") == 0,
      NA,
      prop.test(
        c(sum(ins_1_year[group == "type_2_misclass_0.4"] == 1), sum(ins_1_year[group == "ref_type_2"] == 1)),
        c(sum(group == "type_2_misclass_0.4"), sum(group == "ref_type_2"))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  
  p_val_bolus_10_year_0.8 <- cohort_dataset %>%
    filter(!is.na(ten_yrs_post_diag)) %>%
    summarise(p = ifelse(
      sum(group == "type_2_misclass_0.8") == 0 | sum(group == "ref_type_2") == 0,
      NA,
      prop.test(
        c(sum(ten_yrs_bolus_insulin[group == "type_2_misclass_0.8"] == 1), sum(ten_yrs_bolus_insulin[group == "ref_type_2"] == 1)),
        c(sum(group == "type_2_misclass_0.8"), sum(group == "ref_type_2"))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_bolus_10_year_0.7 <- cohort_dataset %>%
    filter(!is.na(ten_yrs_post_diag)) %>%
    summarise(p = ifelse(
      sum(group == "type_2_misclass_0.7") == 0 | sum(group == "ref_type_2") == 0,
      NA,
      prop.test(
        c(sum(ten_yrs_bolus_insulin[group == "type_2_misclass_0.7"] == 1), sum(ten_yrs_bolus_insulin[group == "ref_type_2"] == 1)),
        c(sum(group == "type_2_misclass_0.7"), sum(group == "ref_type_2"))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_bolus_10_year_0.6 <- cohort_dataset %>%
    filter(!is.na(ten_yrs_post_diag)) %>%
    summarise(p = ifelse(
      sum(group == "type_2_misclass_0.6") == 0 | sum(group == "ref_type_2") == 0,
      NA,
      prop.test(
        c(sum(ten_yrs_bolus_insulin[group == "type_2_misclass_0.6"] == 1), sum(ten_yrs_bolus_insulin[group == "ref_type_2"] == 1)),
        c(sum(group == "type_2_misclass_0.6"), sum(group == "ref_type_2"))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_bolus_10_year_0.5 <- cohort_dataset %>%
    filter(!is.na(ten_yrs_post_diag)) %>%
    summarise(p = ifelse(
      sum(group == "type_2_misclass_0.5") == 0 | sum(group == "ref_type_2") == 0,
      NA,
      prop.test(
        c(sum(ten_yrs_bolus_insulin[group == "type_2_misclass_0.5"] == 1), sum(ten_yrs_bolus_insulin[group == "ref_type_2"] == 1)),
        c(sum(group == "type_2_misclass_0.5"), sum(group == "ref_type_2"))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_bolus_10_year_0.4 <- cohort_dataset %>%
    filter(!is.na(ten_yrs_post_diag)) %>%
    summarise(p = ifelse(
      sum(group == "type_2_misclass_0.4") == 0 | sum(group == "ref_type_2") == 0,
      NA,
      prop.test(
        c(sum(ten_yrs_bolus_insulin[group == "type_2_misclass_0.4"] == 1), sum(ten_yrs_bolus_insulin[group == "ref_type_2"] == 1)),
        c(sum(group == "type_2_misclass_0.4"), sum(group == "ref_type_2"))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  print(paste0("p-value DKA for 80/70/60/50/40% thresholds: ", p_val_dka_0.8, ", ", p_val_dka_0.7, ", ", p_val_dka_0.6, ", ", p_val_dka_0.5, ", ", p_val_dka_0.4))
  print(paste0("p-value hypoglycaemia for 80/70/60/50/40% thresholds: ", p_val_hypo_0.8, ", ", p_val_hypo_0.7, ", ", p_val_hypo_0.6, ", ", p_val_hypo_0.5, ", ", p_val_hypo_0.4))
  print(paste0("p-value insulin 1 year for 80/70/60/50/40% thresholds: ", p_val_ins_1_year_0.8, ", ", p_val_ins_1_year_0.7, ", ", p_val_ins_1_year_0.6, ", ", p_val_ins_1_year_0.5, ", ", p_val_ins_1_year_0.4))
  print(paste0("p-value bolus insulin 10 years for 80/70/60/50/40% thresholds: ", p_val_bolus_10_year_0.8, ", ", p_val_bolus_10_year_0.7, ", ", p_val_bolus_10_year_0.6, ", ", p_val_bolus_10_year_0.5, ", ", p_val_bolus_10_year_0.4))


  missing_summary <- cohort_dataset %>%
    group_by(group) %>%
    summarise(
      n = n(),
      ins_1_year_non_missing = sum(!is.na(ins_1_year)),
      ten_yrs_post_diag_non_missing = sum(!is.na(ten_yrs_post_diag)),
      .groups = "drop"
    ) %>%
    mutate(
      ins_1_year_avail = paste0(round_pad(ins_1_year_non_missing / n * 100, 0), "%"),
      ten_yrs_post_diag_avail = paste0(round_pad(ten_yrs_post_diag_non_missing / n * 100, 0), "%")
    )
  
  missing_ins_1_year <- paste0(missing_summary$ins_1_year_avail[missing_summary$group == "type_2_misclass_0.8"], "/",
                               missing_summary$ins_1_year_avail[missing_summary$group == "type_2_misclass_0.7"], "/",
                               missing_summary$ins_1_year_avail[missing_summary$group == "type_2_misclass_0.6"], "/",
                               missing_summary$ins_1_year_avail[missing_summary$group == "type_2_misclass_0.5"], "/",
                               missing_summary$ins_1_year_avail[missing_summary$group == "type_2_misclass_0.4"], "/",
                               missing_summary$ins_1_year_avail[missing_summary$group == "ref_type_1"], "/",
                               missing_summary$ins_1_year_avail[missing_summary$group == "ref_type_2"])
  missing_bolus_10_year <- paste0(missing_summary$ten_yrs_post_diag_avail[missing_summary$group == "type_2_misclass_0.8"], "/",
                                  missing_summary$ten_yrs_post_diag_avail[missing_summary$group == "type_2_misclass_0.7"], "/",
                                  missing_summary$ten_yrs_post_diag_avail[missing_summary$group == "type_2_misclass_0.6"], "/",
                                  missing_summary$ten_yrs_post_diag_avail[missing_summary$group == "type_2_misclass_0.5"], "/",
                                  missing_summary$ten_yrs_post_diag_avail[missing_summary$group == "type_2_misclass_0.4"], "/",
                               missing_summary$ten_yrs_post_diag_avail[missing_summary$group == "ref_type_1"], "/",
                               missing_summary$ten_yrs_post_diag_avail[missing_summary$group == "ref_type_2"])
  
  
  print(paste0("Insulin within 1 year of diagnosis available for ", missing_ins_1_year))
  print(paste0("Basal-bolus insulin regime at 10 years post-diagnosis available for ", missing_bolus_10_year))
  
  
  ## Plots
  
  dka_post_diag <- cohort_dataset %>%
    group_by(group) %>%
    summarise(
      outcome="DKA\nhospitalisation",
      estimate    = poisson.test(sum(dka_post_diagnosis == 1), sum(follow_up_time))$estimate*100,
      lower = poisson.test(sum(dka_post_diagnosis == 1), sum(follow_up_time))$conf.int[1]*100,
      upper = poisson.test(sum(dka_post_diagnosis == 1), sum(follow_up_time))$conf.int[2]*100,
      .groups = "drop"
    ) %>%
    union(data.frame(group="a", outcome="DKA\nhospitalisation", estimate=NA, lower=NA, upper=NA)) %>%
    mutate(group=factor(group, levels=c("type_2_misclass_0.8", "type_2_misclass_0.7", "type_2_misclass_0.6", "type_2_misclass_0.5", "type_2_misclass_0.4", "a", "ref_type_1", "ref_type_2")))
  
  hypo_post_diag <- cohort_dataset %>%
    group_by(group) %>%
    summarise(
      outcome="Hypoglycaemia\nhospitalisation",
      estimate    = poisson.test(sum(hypo_post_diagnosis == 1), sum(follow_up_time))$estimate*100,
      lower = poisson.test(sum(hypo_post_diagnosis == 1), sum(follow_up_time))$conf.int[1]*100,
      upper = poisson.test(sum(hypo_post_diagnosis == 1), sum(follow_up_time))$conf.int[2]*100,
      .groups = "drop"
    ) %>%
    union(data.frame(group="a", outcome="Hypoglycaemia\nhospitalisation", estimate=NA, lower=NA, upper=NA)) %>%
    mutate(group=factor(group, levels=c("type_2_misclass_0.8", "type_2_misclass_0.7", "type_2_misclass_0.6", "type_2_misclass_0.5", "type_2_misclass_0.4", "a", "ref_type_1", "ref_type_2")))
  
  ins_1_year <- cohort_dataset %>%
    filter(!is.na(ins_1_year)) %>%
    group_by(group) %>%
    summarise(total=n(),
              outcome="Insulin within 1\nyear of diagnosis",
              with_outcome=sum(ins_1_year==1)) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(
      test = list(binom::binom.confint(with_outcome, total, methods = "wilson")),
      estimate = test$mean*100,
      lower = test$lower*100,
      upper = test$upper*100
    ) %>%
    select(-c(test, total, with_outcome)) %>%
    union(data.frame(group="a", outcome="Insulin within 1\nyear of diagnosis", estimate=NA, lower=NA, upper=NA)) %>%
    mutate(group=factor(group, levels=c("type_2_misclass_0.8", "type_2_misclass_0.7", "type_2_misclass_0.6", "type_2_misclass_0.5", "type_2_misclass_0.4", "a", "ref_type_1", "ref_type_2")))
  
  ten_yrs_bolus_insulin <- cohort_dataset %>%
    filter(!is.na(ten_yrs_post_diag)) %>%
    group_by(group) %>%
    summarise(total=n(),
              outcome="Basal-bolus insulin regime\nat 10 years post-diagnosis",
              with_outcome=sum(ten_yrs_bolus_insulin==1)) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(
      test = list(binom::binom.confint(with_outcome, total, methods = "wilson")),
      estimate = test$mean*100,
      lower = test$lower*100,
      upper = test$upper*100
    ) %>%
    select(-c(test, total, with_outcome)) %>%
    union(data.frame(group="a", outcome="Basal-bolus insulin regime\nat 10 years post-diagnosis", estimate=NA, lower=NA, upper=NA)) %>%
    mutate(group=factor(group, levels=c("type_2_misclass_0.8", "type_2_misclass_0.7", "type_2_misclass_0.6", "type_2_misclass_0.5", "type_2_misclass_0.4", "a", "ref_type_1", "ref_type_2")))
  
  
  hosp_t2_chart_data <- rbind(dka_post_diag, hypo_post_diag) %>% rowwise() %>% mutate(label_y = ifelse(group=="type_2_misclass_0.8" | group=="type_2_misclass_0.6" | group=="type_2_misclass_0.4", upper + 0.09, upper + 0.03))
  
  treatment_t2_chart_data <- rbind(ins_1_year, ten_yrs_bolus_insulin) %>% rowwise() %>% mutate(label_y = ifelse(!is.na(upper) & (group=="type_2_misclass_0.7" | group=="type_2_misclass_0.5" | group=="type_2_misclass_0.4"), upper + 3, ifelse(!is.na(upper), upper + 4.2, NA_real_))) %>% ungroup()
  treatment_t2_chart_data$outcome <- factor(treatment_t2_chart_data$outcome, levels=c("Insulin within 1\nyear of diagnosis", "Basal-bolus insulin regime\nat 10 years post-diagnosis"))
 
  
  group_counts <- cohort_dataset %>% count(group, name = "n")
  
  legend_labels <- c(
    "type_2_misclass_0.8" = "\nPotentially misclassified\nT2D using ≥80%\nmodel threshold\n",
    "type_2_misclass_0.7" = "\nPotentially misclassified\nT2D using ≥70%\nmodel threshold\n",
    "type_2_misclass_0.6" = "\nPotentially misclassified\nT2D using ≥60%\nmodel threshold\n",
    "type_2_misclass_0.5" = "\nPotentially misclassified\nT2D using ≥50%\nmodel threshold\n",
    "type_2_misclass_0.4" = "\nPotentially misclassified\nT2D using ≥40%\nmodel threshold\n",
    "ref_type_1" = "\nCorrectly classified\nT1D",
    "ref_type_2" = "\nCorrectly classified\nT2D"
  )
  
  legend_labels <- map2_chr(
    names(legend_labels), legend_labels,
    ~ paste0(.y, " (n = ", format(group_counts$n[group_counts$group == .x], big.mark = ","), ")")
  )
  
  max_y <- 1.5

  hosp_plot <- ggplot(hosp_t2_chart_data, aes(fill=group, y=estimate, x=outcome)) +
    geom_bar(position="dodge", stat="identity", width=0.7) +
    geom_text(aes(y=label_y, label=round_pad(estimate,2), group = group),
              position = position_dodge(width = 0.7), size=6, fontface="bold") +
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.3, position=position_dodge(.7)) +
    theme_bw() +
    ylab("Incidence rate (patients per 100 patient-years)") +
    scale_x_discrete(expand = expansion(add = 0.5)) +
    scale_fill_manual("legend", labels=legend_labels, values = c("type_2_misclass_0.8" = "darkred", "type_2_misclass_0.7" = "red3", "type_2_misclass_0.6" = "red1", "type_2_misclass_0.5" = "salmon", "type_2_misclass_0.4" = "lightpink", "ref_type_1"="dodgerblue3", "ref_type_2"="darkgoldenrod2")) +
    scale_y_continuous(limits = c(0, max_y)) +
    theme(panel.grid.major.x = element_blank(),
          axis.title.y = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 20, face="bold", vjust=-0.2),
          panel.grid.major.y = element_line(size=1),
          panel.grid.minor.y = element_blank(),
          axis.title.x=element_blank(),
          legend.text=element_text(size=20),
          legend.title = element_blank(),
          plot.margin = unit(c(0,0.7,0.7,0.7), "cm"))
  
  treatment_plot <- ggplot(treatment_t2_chart_data, aes(fill=group, y=estimate, x=outcome, ymax=35)) +
    geom_bar(position="dodge", stat="identity", width=0.7) +
    geom_text(aes(y=label_y, label=paste0(round_pad(estimate,0), "%"), group = group),
              position = position_dodge(width = 0.7), size=6, fontface="bold") +
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.3, position=position_dodge(.7)) +
    theme_bw() +
    ylab("Percentage (%)") +
    scale_x_discrete(expand = expansion(add = 0.5)) +
    scale_fill_manual("legend", labels=legend_labels, values = c("type_2_misclass_0.8" = "darkred", "type_2_misclass_0.7" = "red3", "type_2_misclass_0.6" = "red1", "type_2_misclass_0.5" = "salmon", "type_2_misclass_0.4" = "lightpink", "ref_type_1"="dodgerblue3", "ref_type_2"="darkgoldenrod2")) +
    scale_y_continuous(limits=c(0,105), breaks=c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
    theme(panel.grid.major.x = element_blank(),
          axis.title.y = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 20, face="bold", vjust=-0.2),
          panel.grid.major.y = element_line(size=1),
          panel.grid.minor.y = element_blank(),
          axis.title.x=element_blank(),
          plot.margin = unit(c(0,0.7,0.7,0.7), "cm"))
  
  legend <- get_legend(hosp_plot)
  hosp_plot <- hosp_plot + theme(legend.position = "none")
  treatment_plot <- treatment_plot + theme(legend.position = "none")
  
  plot <- grid.arrange(hosp_plot, treatment_plot,  legend,
                       ncol=3, nrow = 1,
                       widths = c(1.8, 1.8, 0.7))
  
  
  plot_file_name <- paste0("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2025/DePICtion/Paper/Plots/t2_outcomes_by_threshold.tiff")
  
  plot_width <- 22
  plot_height <- 10
  
  tiff(plot_file_name, width=plot_width, height=plot_height, units = "in", res=800)
  
  print(ggpubr::as_ggplot(plot)) #+
    #draw_plot_label(label = c("a)", "b)"), size = 20,  hjust=0, x = c(0, 0.4), y = c(1, 1)))
  
  dev.off()
  
  
}



outcomes_thresholds_t1 <- function(cohort_dataset) {
  
  p_val_dka_0.025 <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(dka_post_diagnosis[group == "type_1_misclass_0.025"] == 1), sum(dka_post_diagnosis[group == "ref_type_1"] == 1)),
      c(sum(follow_up_time[group == "type_1_misclass_0.025"]), sum(follow_up_time[group == "ref_type_1"]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_dka_0.05 <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(dka_post_diagnosis[group == "type_1_misclass_0.05"] == 1), sum(dka_post_diagnosis[group == "ref_type_1"] == 1)),
      c(sum(follow_up_time[group == "type_1_misclass_0.05"]), sum(follow_up_time[group == "ref_type_1"]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_dka_0.1 <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(dka_post_diagnosis[group == "type_1_misclass_0.1"] == 1), sum(dka_post_diagnosis[group == "ref_type_1"] == 1)),
      c(sum(follow_up_time[group == "type_1_misclass_0.1"]), sum(follow_up_time[group == "ref_type_1"]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  
  p_val_hypo_0.025 <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(hypo_post_diagnosis[group == "type_1_misclass_0.025"] == 1), sum(hypo_post_diagnosis[group == "ref_type_1"] == 1)),
      c(sum(follow_up_time[group == "type_1_misclass_0.025"]), sum(follow_up_time[group == "ref_type_1"]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_hypo_0.05 <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(hypo_post_diagnosis[group == "type_1_misclass_0.05"] == 1), sum(hypo_post_diagnosis[group == "ref_type_1"] == 1)),
      c(sum(follow_up_time[group == "type_1_misclass_0.05"]), sum(follow_up_time[group == "ref_type_1"]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_hypo_0.1 <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(hypo_post_diagnosis[group == "type_1_misclass_0.1"] == 1), sum(hypo_post_diagnosis[group == "ref_type_1"] == 1)),
      c(sum(follow_up_time[group == "type_1_misclass_0.1"]), sum(follow_up_time[group == "ref_type_1"]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  
  p_val_ins_1_year_0.025 <- cohort_dataset %>%
    filter(!is.na(ins_1_year)) %>%
    summarise(p = ifelse(
      sum(group == "type_1_misclass_0.025") == 0 | sum(group == "ref_type_1") == 0,
      NA,
      prop.test(
        c(sum(ins_1_year[group == "type_1_misclass_0.025"] == 1), sum(ins_1_year[group == "ref_type_1"] == 1)),
        c(sum(group == "type_1_misclass_0.025"), sum(group == "ref_type_1"))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_ins_1_year_0.05 <- cohort_dataset %>%
    filter(!is.na(ins_1_year)) %>%
    summarise(p = ifelse(
      sum(group == "type_1_misclass_0.05") == 0 | sum(group == "ref_type_1") == 0,
      NA,
      prop.test(
        c(sum(ins_1_year[group == "type_1_misclass_0.05"] == 1), sum(ins_1_year[group == "ref_type_1"] == 1)),
        c(sum(group == "type_1_misclass_0.05"), sum(group == "ref_type_1"))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_ins_1_year_0.1 <- cohort_dataset %>%
    filter(!is.na(ins_1_year)) %>%
    summarise(p = ifelse(
      sum(group == "type_1_misclass_0.1") == 0 | sum(group == "ref_type_1") == 0,
      NA,
      prop.test(
        c(sum(ins_1_year[group == "type_1_misclass_0.1"] == 1), sum(ins_1_year[group == "ref_type_1"] == 1)),
        c(sum(group == "type_1_misclass_0.1"), sum(group == "ref_type_1"))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  
  p_val_bolus_10_year_0.025 <- cohort_dataset %>%
    filter(!is.na(ten_yrs_post_diag)) %>%
    summarise(p = ifelse(
      sum(group == "type_1_misclass_0.025") == 0 | sum(group == "ref_type_1") == 0,
      NA,
      prop.test(
        c(sum(ten_yrs_bolus_insulin[group == "type_1_misclass_0.025"] == 1), sum(ten_yrs_bolus_insulin[group == "ref_type_1"] == 1)),
        c(sum(group == "type_1_misclass_0.025"), sum(group == "ref_type_1"))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_bolus_10_year_0.05 <- cohort_dataset %>%
    filter(!is.na(ten_yrs_post_diag)) %>%
    summarise(p = ifelse(
      sum(group == "type_1_misclass_0.05") == 0 | sum(group == "ref_type_1") == 0,
      NA,
      prop.test(
        c(sum(ten_yrs_bolus_insulin[group == "type_1_misclass_0.05"] == 1), sum(ten_yrs_bolus_insulin[group == "ref_type_1"] == 1)),
        c(sum(group == "type_1_misclass_0.05"), sum(group == "ref_type_1"))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_bolus_10_year_0.1 <- cohort_dataset %>%
    filter(!is.na(ten_yrs_post_diag)) %>%
    summarise(p = ifelse(
      sum(group == "type_1_misclass_0.1") == 0 | sum(group == "ref_type_1") == 0,
      NA,
      prop.test(
        c(sum(ten_yrs_bolus_insulin[group == "type_1_misclass_0.1"] == 1), sum(ten_yrs_bolus_insulin[group == "ref_type_1"] == 1)),
        c(sum(group == "type_1_misclass_0.1"), sum(group == "ref_type_1"))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  
  print(paste0("p-value DKA for 2.5/5/10% thresholds: ", p_val_dka_0.025, ", ", p_val_dka_0.05, ", ", p_val_dka_0.1))
  print(paste0("p-value hypoglycaemia for 2.5/5/10% thresholds: ", p_val_hypo_0.025, ", ", p_val_hypo_0.05, ", ", p_val_hypo_0.1))
  print(paste0("p-value insulin 1 year for 2.5/5/10% thresholds: ", p_val_ins_1_year_0.025, ", ", p_val_ins_1_year_0.05, ", ", p_val_ins_1_year_0.1))
  print(paste0("p-value bolus insulin 10 years for 2.5/5/10% thresholds: ", p_val_bolus_10_year_0.025, ", ", p_val_bolus_10_year_0.05, ", ", p_val_bolus_10_year_0.1))
  
  
  missing_summary <- cohort_dataset %>%
    group_by(group) %>%
    summarise(
      n = n(),
      ins_1_year_non_missing = sum(!is.na(ins_1_year)),
      ten_yrs_post_diag_non_missing = sum(!is.na(ten_yrs_post_diag)),
      .groups = "drop"
    ) %>%
    mutate(
      ins_1_year_avail = paste0(round_pad(ins_1_year_non_missing / n * 100, 0), "%"),
      ten_yrs_post_diag_avail = paste0(round_pad(ten_yrs_post_diag_non_missing / n * 100, 0), "%")
    )
  
  missing_ins_1_year <- paste0(missing_summary$ins_1_year_avail[missing_summary$group == "type_1_misclass_0.025"], "/",
                               missing_summary$ins_1_year_avail[missing_summary$group == "type_1_misclass_0.05"], "/",
                               missing_summary$ins_1_year_avail[missing_summary$group == "type_1_misclass_0.1"], "/",
                               missing_summary$ins_1_year_avail[missing_summary$group == "ref_type_1"], "/",
                               missing_summary$ins_1_year_avail[missing_summary$group == "ref_type_2"])
  missing_bolus_10_year <- paste0(missing_summary$ten_yrs_post_diag_avail[missing_summary$group == "type_1_misclass_0.025"], "/",
                                  missing_summary$ten_yrs_post_diag_avail[missing_summary$group == "type_1_misclass_0.05"], "/",
                                  missing_summary$ten_yrs_post_diag_avail[missing_summary$group == "type_1_misclass_0.1"], "/",
                                  missing_summary$ten_yrs_post_diag_avail[missing_summary$group == "ref_type_1"], "/",
                                  missing_summary$ten_yrs_post_diag_avail[missing_summary$group == "ref_type_2"])
  
  
  print(paste0("Insulin within 1 year of diagnosis available for ", missing_ins_1_year))
  print(paste0("Basal-bolus insulin regime at 10 years post-diagnosis available for ", missing_bolus_10_year))
  
  
  ## Plots
  
  dka_post_diag <- cohort_dataset %>%
    group_by(group) %>%
    summarise(
      outcome="DKA\nhospitalisation",
      estimate    = poisson.test(sum(dka_post_diagnosis == 1), sum(follow_up_time))$estimate*100,
      lower = poisson.test(sum(dka_post_diagnosis == 1), sum(follow_up_time))$conf.int[1]*100,
      upper = poisson.test(sum(dka_post_diagnosis == 1), sum(follow_up_time))$conf.int[2]*100,
      .groups = "drop"
    ) %>%
    union(data.frame(group="a", outcome="DKA\nhospitalisation", estimate=NA, lower=NA, upper=NA)) %>%
    mutate(group=factor(group, levels=c("type_1_misclass_0.025", "type_1_misclass_0.05", "type_1_misclass_0.1", "a", "ref_type_1", "ref_type_2")))
  
  hypo_post_diag <- cohort_dataset %>%
    group_by(group) %>%
    summarise(
      outcome="Hypoglycaemia\nhospitalisation",
      estimate    = poisson.test(sum(hypo_post_diagnosis == 1), sum(follow_up_time))$estimate*100,
      lower = poisson.test(sum(hypo_post_diagnosis == 1), sum(follow_up_time))$conf.int[1]*100,
      upper = poisson.test(sum(hypo_post_diagnosis == 1), sum(follow_up_time))$conf.int[2]*100,
      .groups = "drop"
    ) %>%
    union(data.frame(group="a", outcome="Hypoglycaemia\nhospitalisation", estimate=NA, lower=NA, upper=NA)) %>%
    mutate(group=factor(group, levels=c("type_1_misclass_0.025", "type_1_misclass_0.05", "type_1_misclass_0.1", "a", "ref_type_1", "ref_type_2")))
  
  ins_1_year <- cohort_dataset %>%
    filter(!is.na(ins_1_year)) %>%
    group_by(group) %>%
    summarise(total=n(),
              outcome="Insulin within 1\nyear of diagnosis",
              with_outcome=sum(ins_1_year==1)) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(
      test = list(binom::binom.confint(with_outcome, total, methods = "wilson")),
      estimate = test$mean*100,
      lower = test$lower*100,
      upper = test$upper*100
    ) %>%
    select(-c(test, total, with_outcome)) %>%
    union(data.frame(group="a", outcome="Insulin within 1\nyear of diagnosis", estimate=NA, lower=NA, upper=NA)) %>%
    mutate(group=factor(group, levels=c("type_1_misclass_0.025", "type_1_misclass_0.05", "type_1_misclass_0.1", "a", "ref_type_1", "ref_type_2")))
  
  ten_yrs_bolus_insulin <- cohort_dataset %>%
    filter(!is.na(ten_yrs_post_diag)) %>%
    group_by(group) %>%
    summarise(total=n(),
              outcome="Basal-bolus insulin regime\nat 10 years post-diagnosis",
              with_outcome=sum(ten_yrs_bolus_insulin==1)) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(
      test = list(binom::binom.confint(with_outcome, total, methods = "wilson")),
      estimate = test$mean*100,
      lower = test$lower*100,
      upper = test$upper*100
    ) %>%
    select(-c(test, total, with_outcome)) %>%
    union(data.frame(group="a", outcome="Basal-bolus insulin regime\nat 10 years post-diagnosis", estimate=NA, lower=NA, upper=NA)) %>%
    mutate(group=factor(group, levels=c("type_1_misclass_0.025", "type_1_misclass_0.05", "type_1_misclass_0.1", "a", "ref_type_1", "ref_type_2")))
  
  
  hosp_t2_chart_data <- rbind(dka_post_diag, hypo_post_diag) %>% rowwise() %>% mutate(label_y = upper + 0.07)
  
  treatment_t2_chart_data <- rbind(ins_1_year, ten_yrs_bolus_insulin) %>% rowwise() %>% mutate(label_y = ifelse(!is.na(upper), upper + 5, NA_real_)) %>% ungroup()
  treatment_t2_chart_data$outcome <- factor(treatment_t2_chart_data$outcome, levels=c("Insulin within 1\nyear of diagnosis", "Basal-bolus insulin regime\nat 10 years post-diagnosis"))
  
  
  group_counts <- cohort_dataset %>% count(group, name = "n")
  
  legend_labels <- c(
    "type_1_misclass_0.025" = "\nPotentially misclassified\nT1D using <2.5%\nmodel threshold\n",
    "type_1_misclass_0.05" = "\nPotentially misclassified\nT1D using <5%\nmodel threshold\n",
    "type_1_misclass_0.1" = "\nPotentially misclassified\nT1D using <10%\nmodel threshold\n",
    "ref_type_1" = "\nCorrectly classified\nT1D",
    "ref_type_2" = "\nCorrectly classified\nT2D"
  )
  
  legend_labels <- map2_chr(
    names(legend_labels), legend_labels,
    ~ paste0(.y, " (n = ", format(group_counts$n[group_counts$group == .x], big.mark = ","), ")\n")
  )
  
  max_y <- 1.5
  
  hosp_plot <- ggplot(hosp_t2_chart_data, aes(fill=group, y=estimate, x=outcome)) +
    geom_bar(position="dodge", stat="identity", width=0.7) +
    geom_text(aes(y=label_y, label=round_pad(estimate,2), group = group),
              position = position_dodge(width = 0.7), size=5, fontface="bold") +
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.3, position=position_dodge(.7)) +
    theme_bw() +
    ylab("Incidence rate (patients per 100 patient-years)") +
    scale_x_discrete(expand = expansion(add = 0.5)) +
    scale_fill_manual("legend", labels=legend_labels, values = c("type_1_misclass_0.025" = "darkred", "type_1_misclass_0.05" = "red1", "type_1_misclass_0.1" = "lightpink", "ref_type_1"="dodgerblue3", "ref_type_2"="darkgoldenrod2")) +
    scale_y_continuous(limits = c(0, max_y)) +
    theme(panel.grid.major.x = element_blank(),
          axis.title.y = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 20, face="bold", vjust=-0.2),
          panel.grid.major.y = element_line(size=1),
          panel.grid.minor.y = element_blank(),
          axis.title.x=element_blank(),
          legend.text=element_text(size=20),
          legend.title = element_blank(),
          plot.margin = unit(c(0,0.7,0.7,0.7), "cm"))
  
  treatment_plot <- ggplot(treatment_t2_chart_data, aes(fill=group, y=estimate, x=outcome, ymax=35)) +
    geom_bar(position="dodge", stat="identity", width=0.7) +
    geom_text(aes(y=label_y, label=paste0(round_pad(estimate,0), "%"), group = group),
              position = position_dodge(width = 0.7), size=5, fontface="bold") +
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.3, position=position_dodge(.7)) +
    theme_bw() +
    ylab("Percentage (%)") +
    scale_x_discrete(expand = expansion(add = 0.5)) +
    scale_fill_manual("legend", labels=legend_labels, values = c("type_1_misclass_0.025" = "darkred", "type_1_misclass_0.05" = "red1", "type_1_misclass_0.1" = "lightpink", "ref_type_1"="dodgerblue3", "ref_type_2"="darkgoldenrod2")) +
    scale_y_continuous(limits=c(0,105), breaks=c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
    theme(panel.grid.major.x = element_blank(),
          axis.title.y = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 20, face="bold", vjust=-0.2),
          panel.grid.major.y = element_line(size=1),
          panel.grid.minor.y = element_blank(),
          axis.title.x=element_blank(),
          plot.margin = unit(c(0,0.7,0.7,0.7), "cm"))
  
  legend <- get_legend(hosp_plot)
  hosp_plot <- hosp_plot + theme(legend.position = "none")
  treatment_plot <- treatment_plot + theme(legend.position = "none")
  
  plot <- grid.arrange(hosp_plot, treatment_plot,  legend,
                       ncol=3, nrow = 1,
                       widths = c(1.8, 1.8, 0.7))
  
  
  plot_file_name <- paste0("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2025/DePICtion/Paper/Plots/t1_outcomes_by_threshold.tiff")
  
  plot_width <- 22
  plot_height <- 10
  
  tiff(plot_file_name, width=plot_width, height=plot_height, units = "in", res=800)
  
  print(ggpubr::as_ggplot(plot)) #+
          #draw_plot_label(label = c("a)", "b)"), size = 20,  hjust=0, x = c(0, 0.4), y = c(1, 1)))
  
  dev.off()
  
  
}