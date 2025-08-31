
# Output incidence rates for hospitalisation for DKA and hypoglycaemia, and %s on insulin at 1 year post diagnosis and on basal-bolus insulin at 10 years post diagnosis
## And test for differences between potentially misclassified and reference groups

outcomes_t2 <- function(cohort_dataset, group_name) {
  
  print(paste0("group = ", group_name))
  
  p_val_dka <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(dka_post_diagnosis[group == "type_2_misclass"] == 1), sum(dka_post_diagnosis[group == "ref_type_2"] == 1)),
      c(sum(follow_up_time[group == "type_2_misclass"]), sum(follow_up_time[group == "ref_type_2"]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_hypo <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(hypo_post_diagnosis[group == "type_2_misclass"] == 1), sum(hypo_post_diagnosis[group == "ref_type_2"] == 1)),
      c(sum(follow_up_time[group == "type_2_misclass"]), sum(follow_up_time[group == "ref_type_2"]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_ins_1_year <- cohort_dataset %>%
    filter(!is.na(ins_1_year)) %>%
    summarise(p = ifelse(
      sum(group == "type_2_misclass") == 0 | sum(group == "ref_type_2") == 0,
      NA,
      prop.test(
        c(sum(ins_1_year[group == "type_2_misclass"] == 1), sum(ins_1_year[group == "ref_type_2"] == 1)),
        c(sum(group == "type_2_misclass"), sum(group == "ref_type_2"))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_bolus_10_year <- cohort_dataset %>%
    filter(!is.na(ten_yrs_post_diag)) %>%
    summarise(p = ifelse(
      sum(group == "type_2_misclass") == 0 | sum(group == "ref_type_2") == 0,
      NA,
      prop.test(
        c(sum(ten_yrs_bolus_insulin[group == "type_2_misclass"] == 1), sum(ten_yrs_bolus_insulin[group == "ref_type_2"] == 1)),
        c(sum(group == "type_2_misclass"), sum(group == "ref_type_2"))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  print(paste0("p-value DKA: ", p_val_dka))
  print(paste0("p-value hypoglycaemia: ", p_val_hypo))
  print(paste0("p-value insulin 1 year: ", p_val_ins_1_year))
  print(paste0("p-value bolus insulin 10 years: ", p_val_bolus_10_year))

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
  
  missing_ins_1_year <- paste0(missing_summary$ins_1_year_avail[missing_summary$group == "type_2_misclass"], "/",
                               missing_summary$ins_1_year_avail[missing_summary$group == "ref_type_1"], "/",
                               missing_summary$ins_1_year_avail[missing_summary$group == "ref_type_2"])
  missing_bolus_10_year <- paste0(missing_summary$ten_yrs_post_diag_avail[missing_summary$group == "type_2_misclass"], "/",
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
    mutate(group=factor(group, levels=c("type_2_misclass", "a", "ref_type_1", "ref_type_2")))
  
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
    mutate(group=factor(group, levels=c("type_2_misclass", "a", "ref_type_1", "ref_type_2")))
  
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
    mutate(group=factor(group, levels=c("type_2_misclass", "a", "ref_type_1", "ref_type_2")))
  
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
    mutate(group=factor(group, levels=c("type_2_misclass", "a", "ref_type_1", "ref_type_2")))
  
  
  hosp_t2_chart_data <- rbind(dka_post_diag, hypo_post_diag) %>% rowwise() %>% mutate(label_y = ifelse(group_name=="overall" | group_name=="overall_ins3yrs" | group_name=="male" | group_name=="female", upper + 0.07, upper + 0.2))
  
  treatment_t2_chart_data <- rbind(ins_1_year, ten_yrs_bolus_insulin) %>% rowwise() %>% mutate(label_y = if (!is.na(upper)) upper + 5 else NA_real_) %>% ungroup()
  treatment_t2_chart_data$outcome <- factor(treatment_t2_chart_data$outcome, levels=c("Insulin within 1\nyear of diagnosis", "Basal-bolus insulin regime\nat 10 years post-diagnosis"))
  
  group_counts <- cohort_dataset %>% count(group, name = "n")
  
  legend_labels <- c(
    "type_2_misclass" = "\nPotentially misclassified\nT2D",
    "ref_type_1" = "\nCorrectly classified\nT1D",
    "ref_type_2" = "\nCorrectly classified\nT2D"
  )
  
  legend_labels <- map2_chr(
    names(legend_labels), legend_labels,
    ~ paste0(.y, " (n = ", format(group_counts$n[group_counts$group == .x], big.mark = ","), ")\n")
  )
  
  max_y <- ifelse(group_name=="overall" | group_name=="overall_ins3yrs" | group_name=="male" | group_name=="female", 1.5, ifelse(group_name=="imd_1" | group_name=="imd_2" | group_name=="imd_3" | group_name=="imd_4" | group_name=="imd_5", 3, 4))

  hosp_plot <- ggplot(hosp_t2_chart_data, aes(fill=group, y=estimate, x=outcome)) +
    geom_bar(position="dodge", stat="identity", width=0.7) +
    geom_text(aes(y=label_y, label=round_pad(estimate,2), group = group),
              position = position_dodge(width = 0.7), size=7, fontface="bold") +
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.3, position=position_dodge(.7)) +
    theme_bw() +
    ylab("Incidence rate (patients per 100 patient-years)") +
    scale_x_discrete(expand = expansion(add = 0.5)) +
    scale_fill_manual("legend", labels=legend_labels, values = c("type_2_misclass" = "darkred", "ref_type_1"="dodgerblue3", "ref_type_2"="darkgoldenrod2")) +
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
              position = position_dodge(width = 0.7), size=7, fontface="bold") +
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.3, position=position_dodge(.7)) +
    theme_bw() +
    ylab("Percentage (%)") +
    scale_x_discrete(expand = expansion(add = 0.5)) +
    scale_fill_manual("legend", labels=c("type_2_misclass"="type 2 potentially misclassified", "ref_type_1"="Reference type 1", "ref_type_2"="Reference type 2"), values = c("type_2_misclass" = "darkred", "ref_type_1"="dodgerblue3", "ref_type_2"="darkgoldenrod2")) +
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
  
  
  plot_file_name <- paste0("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2025/DePICtion/Paper/Plots/t2_outcomes_", group_name, ".tiff")
  
  plot_width <- 22
  plot_height <- 10
  
  tiff(plot_file_name, width=plot_width, height=plot_height, units = "in", res=800)
  
  print(ggpubr::as_ggplot(plot)) #+
    #draw_plot_label(label = c("a)", "b)"), size = 20,  hjust=0, x = c(0, 0.43), y = c(1, 1)))
  
  dev.off()
  
  
}



compare_t2_sex <- function(cohort_dataset) {
  
  print("Misclassified")
  
  p_val_dka <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(dka_post_diagnosis[group == "type_2_misclass" & malesex==1] == 1), sum(dka_post_diagnosis[group == "type_2_misclass" & malesex==0] == 1)),
      c(sum(follow_up_time[group == "type_2_misclass" & malesex==1]), sum(follow_up_time[group == "type_2_misclass" & malesex==0]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_hypo <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(hypo_post_diagnosis[group == "type_2_misclass" & malesex==1] == 1), sum(hypo_post_diagnosis[group == "type_2_misclass" & malesex==0] == 1)),
      c(sum(follow_up_time[group == "type_2_misclass" & malesex==1]), sum(follow_up_time[group == "type_2_misclass" & malesex==0]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_ins_1_year <- cohort_dataset %>%
    filter(!is.na(ins_1_year)) %>%
    summarise(p = ifelse(
      sum(group == "type_2_misclass" & malesex==1) == 0 | sum(group == "type_2_misclass"  & malesex==0) == 0,
      NA,
      prop.test(
        c(sum(ins_1_year[group == "type_2_misclass" & malesex==1] == 1), sum(ins_1_year[group == "type_2_misclass" & malesex==0] == 1)),
        c(sum(group == "type_2_misclass" & malesex==1), sum(group == "type_2_misclass" & malesex==0))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_bolus_10_year <- cohort_dataset %>%
    filter(!is.na(ten_yrs_post_diag)) %>%
    summarise(p = ifelse(
      sum(group == "type_2_misclass" & malesex==1) == 0 | sum(group == "type_2_misclass"  & malesex==0) == 0,
      NA,
      prop.test(
        c(sum(ten_yrs_bolus_insulin[group == "type_2_misclass" & malesex==1] == 1), sum(ten_yrs_bolus_insulin[group == "type_2_misclass" & malesex==0] == 1)),
        c(sum(group == "type_2_misclass" & malesex==1), sum(group == "type_2_misclass" & malesex==0))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  print(paste0("p-value male vs female DKA: ", p_val_dka))
  print(paste0("p-value male vs female hypoglycaemia: ", p_val_hypo))
  print(paste0("p-value male vs female insulin 1 year: ", p_val_ins_1_year))
  print(paste0("p-value male vs female bolus insulin 10 years: ", p_val_bolus_10_year))
  
  
  print("Reference type 2")
  
  p_val_dka <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(dka_post_diagnosis[group == "ref_type_2" & malesex==1] == 1), sum(dka_post_diagnosis[group == "ref_type_2" & malesex==0] == 1)),
      c(sum(follow_up_time[group == "ref_type_2" & malesex==1]), sum(follow_up_time[group == "ref_type_2" & malesex==0]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_hypo <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(hypo_post_diagnosis[group == "ref_type_2" & malesex==1] == 1), sum(hypo_post_diagnosis[group == "ref_type_2" & malesex==0] == 1)),
      c(sum(follow_up_time[group == "ref_type_2" & malesex==1]), sum(follow_up_time[group == "ref_type_2" & malesex==0]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_ins_1_year <- cohort_dataset %>%
    filter(!is.na(ins_1_year)) %>%
    summarise(p = ifelse(
      sum(group == "ref_type_2" & malesex==1) == 0 | sum(group == "ref_type_2"  & malesex==0) == 0,
      NA,
      prop.test(
        c(sum(ins_1_year[group == "ref_type_2" & malesex==1] == 1), sum(ins_1_year[group == "ref_type_2" & malesex==0] == 1)),
        c(sum(group == "ref_type_2" & malesex==1), sum(group == "ref_type_2" & malesex==0))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_bolus_10_year <- cohort_dataset %>%
    filter(!is.na(ten_yrs_post_diag)) %>%
    summarise(p = ifelse(
      sum(group == "ref_type_2" & malesex==1) == 0 | sum(group == "ref_type_2"  & malesex==0) == 0,
      NA,
      prop.test(
        c(sum(ten_yrs_bolus_insulin[group == "ref_type_2" & malesex==1] == 1), sum(ten_yrs_bolus_insulin[group == "ref_type_2" & malesex==0] == 1)),
        c(sum(group == "ref_type_2" & malesex==1), sum(group == "ref_type_2" & malesex==0))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  print(paste0("p-value male vs female DKA: ", p_val_dka))
  print(paste0("p-value male vs female hypoglycaemia: ", p_val_hypo))
  print(paste0("p-value male vs female insulin 1 year: ", p_val_ins_1_year))
  print(paste0("p-value male vs female bolus insulin 10 years: ", p_val_bolus_10_year))
  
}



compare_t2_ethnicity <- function(cohort_dataset) {
  
  print("Misclassified")
  
  
  ## South Asian vs White
  p_val_dka <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(dka_post_diagnosis[group == "type_2_misclass" & ethnicity_decoded=="White"] == 1), sum(dka_post_diagnosis[group == "type_2_misclass" & ethnicity_decoded=="South Asian"] == 1)),
      c(sum(follow_up_time[group == "type_2_misclass" & ethnicity_decoded=="White"]), sum(follow_up_time[group == "type_2_misclass" & ethnicity_decoded=="South Asian"]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_hypo <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(hypo_post_diagnosis[group == "type_2_misclass" & ethnicity_decoded=="White"] == 1), sum(hypo_post_diagnosis[group == "type_2_misclass" & ethnicity_decoded=="South Asian"] == 1)),
      c(sum(follow_up_time[group == "type_2_misclass" & ethnicity_decoded=="White"]), sum(follow_up_time[group == "type_2_misclass" & ethnicity_decoded=="South Asian"]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_ins_1_year <- cohort_dataset %>%
    filter(!is.na(ins_1_year)) %>%
    summarise(p = ifelse(
      sum(group == "type_2_misclass" & ethnicity_decoded=="White") == 0 | sum(group == "type_2_misclass"  & ethnicity_decoded=="South Asian") == 0,
      NA,
      prop.test(
        c(sum(ins_1_year[group == "type_2_misclass" & ethnicity_decoded=="White"] == 1), sum(ins_1_year[group == "type_2_misclass" & ethnicity_decoded=="South Asian"] == 1)),
        c(sum(group == "type_2_misclass" & ethnicity_decoded=="White"), sum(group == "type_2_misclass" & ethnicity_decoded=="South Asian"))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_bolus_10_year <- cohort_dataset %>%
    filter(!is.na(ten_yrs_post_diag)) %>%
    summarise(p = ifelse(
      sum(group == "type_2_misclass" & ethnicity_decoded=="White") == 0 | sum(group == "type_2_misclass"  & ethnicity_decoded=="South Asian") == 0,
      NA,
      prop.test(
        c(sum(ten_yrs_bolus_insulin[group == "type_2_misclass" & ethnicity_decoded=="White"] == 1), sum(ten_yrs_bolus_insulin[group == "type_2_misclass" & ethnicity_decoded=="South Asian"] == 1)),
        c(sum(group == "type_2_misclass" & ethnicity_decoded=="White"), sum(group == "type_2_misclass" & ethnicity_decoded=="South Asian"))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  print(paste0("p-value South Asian vs White DKA: ", p_val_dka))
  print(paste0("p-value South Asian vs White hypoglycaemia: ", p_val_hypo))
  print(paste0("p-value South Asian vs White insulin 1 year: ", p_val_ins_1_year))
  print(paste0("p-value South Asian vs White bolus insulin 10 years: ", p_val_bolus_10_year))
  
  
  ## Black vs White
  p_val_dka <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(dka_post_diagnosis[group == "type_2_misclass" & ethnicity_decoded=="White"] == 1), sum(dka_post_diagnosis[group == "type_2_misclass" & ethnicity_decoded=="Black"] == 1)),
      c(sum(follow_up_time[group == "type_2_misclass" & ethnicity_decoded=="White"]), sum(follow_up_time[group == "type_2_misclass" & ethnicity_decoded=="Black"]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_hypo <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(hypo_post_diagnosis[group == "type_2_misclass" & ethnicity_decoded=="White"] == 1), sum(hypo_post_diagnosis[group == "type_2_misclass" & ethnicity_decoded=="Black"] == 1)),
      c(sum(follow_up_time[group == "type_2_misclass" & ethnicity_decoded=="White"]), sum(follow_up_time[group == "type_2_misclass" & ethnicity_decoded=="Black"]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_ins_1_year <- cohort_dataset %>%
    filter(!is.na(ins_1_year)) %>%
    summarise(p = ifelse(
      sum(group == "type_2_misclass" & ethnicity_decoded=="White") == 0 | sum(group == "type_2_misclass"  & ethnicity_decoded=="Black") == 0,
      NA,
      prop.test(
        c(sum(ins_1_year[group == "type_2_misclass" & ethnicity_decoded=="White"] == 1), sum(ins_1_year[group == "type_2_misclass" & ethnicity_decoded=="Black"] == 1)),
        c(sum(group == "type_2_misclass" & ethnicity_decoded=="White"), sum(group == "type_2_misclass" & ethnicity_decoded=="Black"))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_bolus_10_year <- cohort_dataset %>%
    filter(!is.na(ten_yrs_post_diag)) %>%
    summarise(p = ifelse(
      sum(group == "type_2_misclass" & ethnicity_decoded=="White") == 0 | sum(group == "type_2_misclass"  & ethnicity_decoded=="Black") == 0,
      NA,
      prop.test(
        c(sum(ten_yrs_bolus_insulin[group == "type_2_misclass" & ethnicity_decoded=="White"] == 1), sum(ten_yrs_bolus_insulin[group == "type_2_misclass" & ethnicity_decoded=="Black"] == 1)),
        c(sum(group == "type_2_misclass" & ethnicity_decoded=="White"), sum(group == "type_2_misclass" & ethnicity_decoded=="Black"))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  print(paste0("p-value Black vs White DKA: ", p_val_dka))
  print(paste0("p-value Black vs White hypoglycaemia: ", p_val_hypo))
  print(paste0("p-value Black vs White insulin 1 year: ", p_val_ins_1_year))
  print(paste0("p-value Black vs White bolus insulin 10 years: ", p_val_bolus_10_year))

  
  print("Reference type 2")
  
  
  ## South Asian vs White
  p_val_dka <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(dka_post_diagnosis[group == "ref_type_2" & ethnicity_decoded=="White"] == 1), sum(dka_post_diagnosis[group == "ref_type_2" & ethnicity_decoded=="South Asian"] == 1)),
      c(sum(follow_up_time[group == "ref_type_2" & ethnicity_decoded=="White"]), sum(follow_up_time[group == "ref_type_2" & ethnicity_decoded=="South Asian"]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_hypo <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(hypo_post_diagnosis[group == "ref_type_2" & ethnicity_decoded=="White"] == 1), sum(hypo_post_diagnosis[group == "ref_type_2" & ethnicity_decoded=="South Asian"] == 1)),
      c(sum(follow_up_time[group == "ref_type_2" & ethnicity_decoded=="White"]), sum(follow_up_time[group == "ref_type_2" & ethnicity_decoded=="South Asian"]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_ins_1_year <- cohort_dataset %>%
    filter(!is.na(ins_1_year)) %>%
    summarise(p = ifelse(
      sum(group == "ref_type_2" & ethnicity_decoded=="White") == 0 | sum(group == "ref_type_2"  & ethnicity_decoded=="South Asian") == 0,
      NA,
      prop.test(
        c(sum(ins_1_year[group == "ref_type_2" & ethnicity_decoded=="White"] == 1), sum(ins_1_year[group == "ref_type_2" & ethnicity_decoded=="South Asian"] == 1)),
        c(sum(group == "ref_type_2" & ethnicity_decoded=="White"), sum(group == "ref_type_2" & ethnicity_decoded=="South Asian"))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_bolus_10_year <- cohort_dataset %>%
    filter(!is.na(ten_yrs_post_diag)) %>%
    summarise(p = ifelse(
      sum(group == "ref_type_2" & ethnicity_decoded=="White") == 0 | sum(group == "ref_type_2"  & ethnicity_decoded=="South Asian") == 0,
      NA,
      prop.test(
        c(sum(ten_yrs_bolus_insulin[group == "ref_type_2" & ethnicity_decoded=="White"] == 1), sum(ten_yrs_bolus_insulin[group == "ref_type_2" & ethnicity_decoded=="South Asian"] == 1)),
        c(sum(group == "ref_type_2" & ethnicity_decoded=="White"), sum(group == "ref_type_2" & ethnicity_decoded=="South Asian"))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  print(paste0("p-value South Asian vs White DKA: ", p_val_dka))
  print(paste0("p-value South Asian vs White hypoglycaemia: ", p_val_hypo))
  print(paste0("p-value South Asian vs White insulin 1 year: ", p_val_ins_1_year))
  print(paste0("p-value South Asian vs White bolus insulin 10 years: ", p_val_bolus_10_year))
  
  
  ## Black vs White
  p_val_dka <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(dka_post_diagnosis[group == "ref_type_2" & ethnicity_decoded=="White"] == 1), sum(dka_post_diagnosis[group == "ref_type_2" & ethnicity_decoded=="Black"] == 1)),
      c(sum(follow_up_time[group == "ref_type_2" & ethnicity_decoded=="White"]), sum(follow_up_time[group == "ref_type_2" & ethnicity_decoded=="Black"]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_hypo <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(hypo_post_diagnosis[group == "ref_type_2" & ethnicity_decoded=="White"] == 1), sum(hypo_post_diagnosis[group == "ref_type_2" & ethnicity_decoded=="Black"] == 1)),
      c(sum(follow_up_time[group == "ref_type_2" & ethnicity_decoded=="White"]), sum(follow_up_time[group == "ref_type_2" & ethnicity_decoded=="Black"]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_ins_1_year <- cohort_dataset %>%
    filter(!is.na(ins_1_year)) %>%
    summarise(p = ifelse(
      sum(group == "ref_type_2" & ethnicity_decoded=="White") == 0 | sum(group == "ref_type_2"  & ethnicity_decoded=="Black") == 0,
      NA,
      prop.test(
        c(sum(ins_1_year[group == "ref_type_2" & ethnicity_decoded=="White"] == 1), sum(ins_1_year[group == "ref_type_2" & ethnicity_decoded=="Black"] == 1)),
        c(sum(group == "ref_type_2" & ethnicity_decoded=="White"), sum(group == "ref_type_2" & ethnicity_decoded=="Black"))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_bolus_10_year <- cohort_dataset %>%
    filter(!is.na(ten_yrs_post_diag)) %>%
    summarise(p = ifelse(
      sum(group == "ref_type_2" & ethnicity_decoded=="White") == 0 | sum(group == "ref_type_2"  & ethnicity_decoded=="Black") == 0,
      NA,
      prop.test(
        c(sum(ten_yrs_bolus_insulin[group == "ref_type_2" & ethnicity_decoded=="White"] == 1), sum(ten_yrs_bolus_insulin[group == "ref_type_2" & ethnicity_decoded=="Black"] == 1)),
        c(sum(group == "ref_type_2" & ethnicity_decoded=="White"), sum(group == "ref_type_2" & ethnicity_decoded=="Black"))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  print(paste0("p-value Black vs White DKA: ", p_val_dka))
  print(paste0("p-value Black vs White hypoglycaemia: ", p_val_hypo))
  print(paste0("p-value Black vs White insulin 1 year: ", p_val_ins_1_year))
  print(paste0("p-value Black vs White bolus insulin 10 years: ", p_val_bolus_10_year))
}



compare_t2_deprivation <- function(cohort_dataset) {
  
  cohort_dataset <- cohort_dataset %>% filter(!is.na(imd_quintiles))
  
  print("Misclassified")
  
  p_val_dka <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(dka_post_diagnosis[group == "type_2_misclass" & imd_quintiles==1] == 1), sum(dka_post_diagnosis[group == "type_2_misclass" & imd_quintiles==5] == 1)),
      c(sum(follow_up_time[group == "type_2_misclass" & imd_quintiles==1]), sum(follow_up_time[group == "type_2_misclass" & imd_quintiles==5]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_hypo <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(hypo_post_diagnosis[group == "type_2_misclass" & imd_quintiles==1] == 1), sum(hypo_post_diagnosis[group == "type_2_misclass" & imd_quintiles==5] == 1)),
      c(sum(follow_up_time[group == "type_2_misclass" & imd_quintiles==1]), sum(follow_up_time[group == "type_2_misclass" & imd_quintiles==5]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_ins_1_year <- cohort_dataset %>%
    filter(!is.na(ins_1_year)) %>%
    summarise(p = ifelse(
      sum(group == "type_2_misclass" & imd_quintiles==1) == 0 | sum(group == "type_2_misclass"  & imd_quintiles==5) == 0,
      NA,
      prop.test(
        c(sum(ins_1_year[group == "type_2_misclass" & imd_quintiles==1] == 1), sum(ins_1_year[group == "type_2_misclass" & imd_quintiles==5] == 1)),
        c(sum(group == "type_2_misclass" & imd_quintiles==1), sum(group == "type_2_misclass" & imd_quintiles==5))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_bolus_10_year <- cohort_dataset %>%
    filter(!is.na(ten_yrs_post_diag)) %>%
    summarise(p = ifelse(
      sum(group == "type_2_misclass" & imd_quintiles==1) == 0 | sum(group == "type_2_misclass"  & imd_quintiles==5) == 0,
      NA,
      prop.test(
        c(sum(ten_yrs_bolus_insulin[group == "type_2_misclass" & imd_quintiles==1] == 1), sum(ten_yrs_bolus_insulin[group == "type_2_misclass" & imd_quintiles==5] == 1)),
        c(sum(group == "type_2_misclass" & imd_quintiles==1), sum(group == "type_2_misclass" & imd_quintiles==5))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  print(paste0("p-value IMD quintile 1 vs 5 DKA: ", p_val_dka))
  print(paste0("p-value IMD quintile 1 vs 5 hypoglycaemia: ", p_val_hypo))
  print(paste0("p-value IMD quintile 1 vs 5 insulin 1 year: ", p_val_ins_1_year))
  print(paste0("p-value IMD quintile 1 vs 5 bolus insulin 10 years: ", p_val_bolus_10_year))
  
  
  print("Reference type 2")
  
  p_val_dka <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(dka_post_diagnosis[group == "ref_type_2" & imd_quintiles==1] == 1), sum(dka_post_diagnosis[group == "ref_type_2" & imd_quintiles==5] == 1)),
      c(sum(follow_up_time[group == "ref_type_2" & imd_quintiles==1]), sum(follow_up_time[group == "ref_type_2" & imd_quintiles==5]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_hypo <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(hypo_post_diagnosis[group == "ref_type_2" & imd_quintiles==1] == 1), sum(hypo_post_diagnosis[group == "ref_type_2" & imd_quintiles==5] == 1)),
      c(sum(follow_up_time[group == "ref_type_2" & imd_quintiles==1]), sum(follow_up_time[group == "ref_type_2" & imd_quintiles==5]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_ins_1_year <- cohort_dataset %>%
    filter(!is.na(ins_1_year)) %>%
    summarise(p = ifelse(
      sum(group == "ref_type_2" & imd_quintiles==1) == 0 | sum(group == "ref_type_2"  & imd_quintiles==5) == 0,
      NA,
      prop.test(
        c(sum(ins_1_year[group == "ref_type_2" & imd_quintiles==1] == 1), sum(ins_1_year[group == "ref_type_2" & imd_quintiles==5] == 1)),
        c(sum(group == "ref_type_2" & imd_quintiles==1), sum(group == "ref_type_2" & imd_quintiles==5))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_bolus_10_year <- cohort_dataset %>%
    filter(!is.na(ten_yrs_post_diag)) %>%
    summarise(p = ifelse(
      sum(group == "ref_type_2" & imd_quintiles==1) == 0 | sum(group == "ref_type_2"  & imd_quintiles==5) == 0,
      NA,
      prop.test(
        c(sum(ten_yrs_bolus_insulin[group == "ref_type_2" & imd_quintiles==1] == 1), sum(ten_yrs_bolus_insulin[group == "ref_type_2" & imd_quintiles==5] == 1)),
        c(sum(group == "ref_type_2" & imd_quintiles==1), sum(group == "ref_type_2" & imd_quintiles==5))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  print(paste0("p-value IMD quintile 1 vs 5 DKA: ", p_val_dka))
  print(paste0("p-value IMD quintile 1 vs 5 hypoglycaemia: ", p_val_hypo))
  print(paste0("p-value IMD quintile 1 vs 5 insulin 1 year: ", p_val_ins_1_year))
  print(paste0("p-value IMD quintile 1 vs 5 bolus insulin 10 years: ", p_val_bolus_10_year))
  
}



outcomes_t1 <- function(cohort_dataset, group_name) {
  
  print(paste0("group = ", group_name))
  
  p_val_dka <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(dka_post_diagnosis[group == "type_1_misclass"] == 1), sum(dka_post_diagnosis[group == "ref_type_1"] == 1)),
      c(sum(follow_up_time[group == "type_1_misclass"]), sum(follow_up_time[group == "ref_type_1"]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_hypo <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(hypo_post_diagnosis[group == "type_1_misclass"] == 1), sum(hypo_post_diagnosis[group == "ref_type_1"] == 1)),
      c(sum(follow_up_time[group == "type_1_misclass"]), sum(follow_up_time[group == "ref_type_1"]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_ins_1_year <- cohort_dataset %>%
    filter(!is.na(ins_1_year)) %>%
    summarise(p = ifelse(
      sum(group == "type_1_misclass") == 0 | sum(group == "ref_type_1") == 0,
      NA,
      prop.test(
        c(sum(ins_1_year[group == "type_1_misclass"] == 1), sum(ins_1_year[group == "ref_type_1"] == 1)),
        c(sum(group == "type_1_misclass"), sum(group == "ref_type_1"))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_bolus_10_year <- cohort_dataset %>%
    filter(!is.na(ten_yrs_post_diag)) %>%
    summarise(p = ifelse(
      sum(group == "type_1_misclass") == 0 | sum(group == "ref_type_1") == 0,
      NA,
      prop.test(
        c(sum(ten_yrs_bolus_insulin[group == "type_1_misclass"] == 1), sum(ten_yrs_bolus_insulin[group == "ref_type_1"] == 1)),
        c(sum(group == "type_1_misclass"), sum(group == "ref_type_1"))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  print(paste0("p-value DKA: ", p_val_dka))
  print(paste0("p-value hypoglycaemia: ", p_val_hypo))
  print(paste0("p-value insulin 1 year: ", p_val_ins_1_year))
  print(paste0("p-value bolus insulin 10 years: ", p_val_bolus_10_year))
  
  
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
  
  missing_ins_1_year <- paste0(missing_summary$ins_1_year_avail[missing_summary$group == "type_1_misclass"], "/",
                               missing_summary$ins_1_year_avail[missing_summary$group == "ref_type_1"], "/",
                               missing_summary$ins_1_year_avail[missing_summary$group == "ref_type_2"])
  missing_bolus_10_year <- paste0(missing_summary$ten_yrs_post_diag_avail[missing_summary$group == "type_1_misclass"], "/",
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
    mutate(group=factor(group, levels=c("type_1_misclass", "a", "ref_type_1", "ref_type_2")))
  
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
    mutate(group=factor(group, levels=c("type_1_misclass", "a", "ref_type_1", "ref_type_2")))
  
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
    mutate(group=factor(group, levels=c("type_1_misclass", "a", "ref_type_1", "ref_type_2")))
  
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
    mutate(group=factor(group, levels=c("type_1_misclass", "a", "ref_type_1", "ref_type_2")))
  
  
  hosp_t1_chart_data <- rbind(dka_post_diag, hypo_post_diag) %>% rowwise() %>% mutate(label_y = ifelse(group_name=="overall" | group_name=="overall_ins3yrs" | group_name=="male" | group_name=="female", upper + 0.07, upper + 0.2))
  
  treatment_t1_chart_data <- rbind(ins_1_year, ten_yrs_bolus_insulin) %>% rowwise() %>% mutate(label_y = if (!is.na(upper)) upper + 5 else NA_real_) %>% ungroup()
  treatment_t1_chart_data$outcome <- factor(treatment_t1_chart_data$outcome, levels=c("Insulin within 1\nyear of diagnosis", "Basal-bolus insulin regime\nat 10 years post-diagnosis"))
  
  
  
  group_counts <- cohort_dataset %>% count(group, name = "n")
  
  legend_labels <- c(
    "type_1_misclass" = "\nPotentially misclassified\nT1D",
    "ref_type_1" = "\nCorrectly classified\nT1D",
    "ref_type_2" = "\nCorrectly classified\nT2D"
  )
  
  legend_labels <- map2_chr(
    names(legend_labels), legend_labels,
    ~ paste0(.y, " (n = ", format(group_counts$n[group_counts$group == .x], big.mark = ","), ")\n")
  )
  
  max_y <- ifelse(group_name=="overall" | group_name=="overall_ins3yrs" | group_name=="male" | group_name=="female", 1.5, ifelse(group_name=="imd_1" | group_name=="imd_2" | group_name=="imd_3" | group_name=="imd_4" | group_name=="imd_5", 3, 4))
  
  hosp_plot <- ggplot(hosp_t1_chart_data, aes(fill=group, y=estimate, x=outcome)) +
    geom_bar(position="dodge", stat="identity", width=0.7) +
    geom_text(aes(y=label_y, label=round_pad(estimate,2), group = group),
              position = position_dodge(width = 0.7), size=7, fontface="bold") +
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.3, position=position_dodge(.7)) +
    theme_bw() +
    ylab("Incidence rate (patients per 100 patient-years)") +
    scale_x_discrete(expand = expansion(add = 0.5)) +
    scale_fill_manual("legend", labels=legend_labels, values = c("type_1_misclass" = "darkred", "ref_type_1"="dodgerblue3", "ref_type_2"="darkgoldenrod2")) +
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
  
  treatment_plot <- ggplot(treatment_t1_chart_data, aes(fill=group, y=estimate, x=outcome, ymax=35)) +
    geom_bar(position="dodge", stat="identity", width=0.7) +
    geom_text(aes(y=label_y, label=paste0(round_pad(estimate,0), "%"), group = group),
              position = position_dodge(width = 0.7), size=7, fontface="bold") +
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.3, position=position_dodge(.7)) +
    theme_bw() +
    ylab("Percentage (%)") +
    scale_x_discrete(expand = expansion(add = 0.5)) +
    scale_fill_manual("legend", labels=c("type_1_misclass"="type 1 potentially misclassified", "ref_type_1"="Reference type 1", "ref_type_2"="Reference type 2"), values = c("type_1_misclass" = "darkred", "ref_type_1"="dodgerblue3", "ref_type_2"="darkgoldenrod2")) +
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
  
  
  plot_file_name <- paste0("C:/Users/ky279/OneDrive - University of Exeter/CPRD/2025/DePICtion/Paper/Plots/t1_outcomes_", group_name, ".tiff")
  
  plot_width <- 22
  plot_height <- 10
  
  tiff(plot_file_name, width=plot_width, height=plot_height, units = "in", res=800)
  
  print(ggpubr::as_ggplot(plot)) #+
          #draw_plot_label(label = c("a)", "b)"), size = 20,  hjust=0, x = c(0, 0.43), y = c(1, 1)))
  
  dev.off()
  
  
}



compare_t1_sex <- function(cohort_dataset) {
  
  print("Misclassified")
  
  p_val_dka <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(dka_post_diagnosis[group == "type_1_misclass" & malesex==1] == 1), sum(dka_post_diagnosis[group == "type_1_misclass" & malesex==0] == 1)),
      c(sum(follow_up_time[group == "type_1_misclass" & malesex==1]), sum(follow_up_time[group == "type_1_misclass" & malesex==0]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_hypo <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(hypo_post_diagnosis[group == "type_1_misclass" & malesex==1] == 1), sum(hypo_post_diagnosis[group == "type_1_misclass" & malesex==0] == 1)),
      c(sum(follow_up_time[group == "type_1_misclass" & malesex==1]), sum(follow_up_time[group == "type_1_misclass" & malesex==0]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_ins_1_year <- cohort_dataset %>%
    filter(!is.na(ins_1_year)) %>%
    summarise(p = ifelse(
      sum(group == "type_1_misclass" & malesex==1) == 0 | sum(group == "type_1_misclass"  & malesex==0) == 0,
      NA,
      prop.test(
        c(sum(ins_1_year[group == "type_1_misclass" & malesex==1] == 1), sum(ins_1_year[group == "type_1_misclass" & malesex==0] == 1)),
        c(sum(group == "type_1_misclass" & malesex==1), sum(group == "type_1_misclass" & malesex==0))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_bolus_10_year <- cohort_dataset %>%
    filter(!is.na(ten_yrs_post_diag)) %>%
    summarise(p = ifelse(
      sum(group == "type_1_misclass" & malesex==1) == 0 | sum(group == "type_1_misclass"  & malesex==0) == 0,
      NA,
      prop.test(
        c(sum(ten_yrs_bolus_insulin[group == "type_1_misclass" & malesex==1] == 1), sum(ten_yrs_bolus_insulin[group == "type_1_misclass" & malesex==0] == 1)),
        c(sum(group == "type_1_misclass" & malesex==1), sum(group == "type_1_misclass" & malesex==0))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  print(paste0("p-value male vs female DKA: ", p_val_dka))
  print(paste0("p-value male vs female hypoglycaemia: ", p_val_hypo))
  print(paste0("p-value male vs female insulin 1 year: ", p_val_ins_1_year))
  print(paste0("p-value male vs female bolus insulin 10 years: ", p_val_bolus_10_year))
  
  
  print("Reference type 1")
  
  p_val_dka <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(dka_post_diagnosis[group == "ref_type_1" & malesex==1] == 1), sum(dka_post_diagnosis[group == "ref_type_1" & malesex==0] == 1)),
      c(sum(follow_up_time[group == "ref_type_1" & malesex==1]), sum(follow_up_time[group == "ref_type_1" & malesex==0]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_hypo <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(hypo_post_diagnosis[group == "ref_type_1" & malesex==1] == 1), sum(hypo_post_diagnosis[group == "ref_type_1" & malesex==0] == 1)),
      c(sum(follow_up_time[group == "ref_type_1" & malesex==1]), sum(follow_up_time[group == "ref_type_1" & malesex==0]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_ins_1_year <- cohort_dataset %>%
    filter(!is.na(ins_1_year)) %>%
    summarise(p = ifelse(
      sum(group == "ref_type_1" & malesex==1) == 0 | sum(group == "ref_type_1"  & malesex==0) == 0,
      NA,
      prop.test(
        c(sum(ins_1_year[group == "ref_type_1" & malesex==1] == 1), sum(ins_1_year[group == "ref_type_1" & malesex==0] == 1)),
        c(sum(group == "ref_type_1" & malesex==1), sum(group == "ref_type_1" & malesex==0))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_bolus_10_year <- cohort_dataset %>%
    filter(!is.na(ten_yrs_post_diag)) %>%
    summarise(p = ifelse(
      sum(group == "ref_type_1" & malesex==1) == 0 | sum(group == "ref_type_1"  & malesex==0) == 0,
      NA,
      prop.test(
        c(sum(ten_yrs_bolus_insulin[group == "ref_type_1" & malesex==1] == 1), sum(ten_yrs_bolus_insulin[group == "ref_type_1" & malesex==0] == 1)),
        c(sum(group == "ref_type_1" & malesex==1), sum(group == "ref_type_1" & malesex==0))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  print(paste0("p-value male vs female DKA: ", p_val_dka))
  print(paste0("p-value male vs female hypoglycaemia: ", p_val_hypo))
  print(paste0("p-value male vs female insulin 1 year: ", p_val_ins_1_year))
  print(paste0("p-value male vs female bolus insulin 10 years: ", p_val_bolus_10_year))
  
}



compare_t1_ethnicity <- function(cohort_dataset) {
  
  print("Misclassified")
  
  
  ## South Asian vs White
  p_val_dka <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(dka_post_diagnosis[group == "type_1_misclass" & ethnicity_decoded=="White"] == 1), sum(dka_post_diagnosis[group == "type_1_misclass" & ethnicity_decoded=="South Asian"] == 1)),
      c(sum(follow_up_time[group == "type_1_misclass" & ethnicity_decoded=="White"]), sum(follow_up_time[group == "type_1_misclass" & ethnicity_decoded=="South Asian"]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_hypo <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(hypo_post_diagnosis[group == "type_1_misclass" & ethnicity_decoded=="White"] == 1), sum(hypo_post_diagnosis[group == "type_1_misclass" & ethnicity_decoded=="South Asian"] == 1)),
      c(sum(follow_up_time[group == "type_1_misclass" & ethnicity_decoded=="White"]), sum(follow_up_time[group == "type_1_misclass" & ethnicity_decoded=="South Asian"]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_ins_1_year <- cohort_dataset %>%
    filter(!is.na(ins_1_year)) %>%
    summarise(p = ifelse(
      sum(group == "type_1_misclass" & ethnicity_decoded=="White") == 0 | sum(group == "type_1_misclass"  & ethnicity_decoded=="South Asian") == 0,
      NA,
      prop.test(
        c(sum(ins_1_year[group == "type_1_misclass" & ethnicity_decoded=="White"] == 1), sum(ins_1_year[group == "type_1_misclass" & ethnicity_decoded=="South Asian"] == 1)),
        c(sum(group == "type_1_misclass" & ethnicity_decoded=="White"), sum(group == "type_1_misclass" & ethnicity_decoded=="South Asian"))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_bolus_10_year <- cohort_dataset %>%
    filter(!is.na(ten_yrs_post_diag)) %>%
    summarise(p = ifelse(
      sum(group == "type_1_misclass" & ethnicity_decoded=="White") == 0 | sum(group == "type_1_misclass"  & ethnicity_decoded=="South Asian") == 0,
      NA,
      prop.test(
        c(sum(ten_yrs_bolus_insulin[group == "type_1_misclass" & ethnicity_decoded=="White"] == 1), sum(ten_yrs_bolus_insulin[group == "type_1_misclass" & ethnicity_decoded=="South Asian"] == 1)),
        c(sum(group == "type_1_misclass" & ethnicity_decoded=="White"), sum(group == "type_1_misclass" & ethnicity_decoded=="South Asian"))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  print(paste0("p-value South Asian vs White DKA: ", p_val_dka))
  print(paste0("p-value South Asian vs White hypoglycaemia: ", p_val_hypo))
  print(paste0("p-value South Asian vs White insulin 1 year: ", p_val_ins_1_year))
  print(paste0("p-value South Asian vs White bolus insulin 10 years: ", p_val_bolus_10_year))
  
  
  ## Black vs White
  p_val_dka <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(dka_post_diagnosis[group == "type_1_misclass" & ethnicity_decoded=="White"] == 1), sum(dka_post_diagnosis[group == "type_1_misclass" & ethnicity_decoded=="Black"] == 1)),
      c(sum(follow_up_time[group == "type_1_misclass" & ethnicity_decoded=="White"]), sum(follow_up_time[group == "type_1_misclass" & ethnicity_decoded=="Black"]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_hypo <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(hypo_post_diagnosis[group == "type_1_misclass" & ethnicity_decoded=="White"] == 1), sum(hypo_post_diagnosis[group == "type_1_misclass" & ethnicity_decoded=="Black"] == 1)),
      c(sum(follow_up_time[group == "type_1_misclass" & ethnicity_decoded=="White"]), sum(follow_up_time[group == "type_1_misclass" & ethnicity_decoded=="Black"]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_ins_1_year <- cohort_dataset %>%
    filter(!is.na(ins_1_year)) %>%
    summarise(p = ifelse(
      sum(group == "type_1_misclass" & ethnicity_decoded=="White") == 0 | sum(group == "type_1_misclass"  & ethnicity_decoded=="Black") == 0,
      NA,
      prop.test(
        c(sum(ins_1_year[group == "type_1_misclass" & ethnicity_decoded=="White"] == 1), sum(ins_1_year[group == "type_1_misclass" & ethnicity_decoded=="Black"] == 1)),
        c(sum(group == "type_1_misclass" & ethnicity_decoded=="White"), sum(group == "type_1_misclass" & ethnicity_decoded=="Black"))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_bolus_10_year <- cohort_dataset %>%
    filter(!is.na(ten_yrs_post_diag)) %>%
    summarise(p = ifelse(
      sum(group == "type_1_misclass" & ethnicity_decoded=="White") == 0 | sum(group == "type_1_misclass"  & ethnicity_decoded=="Black") == 0,
      NA,
      prop.test(
        c(sum(ten_yrs_bolus_insulin[group == "type_1_misclass" & ethnicity_decoded=="White"] == 1), sum(ten_yrs_bolus_insulin[group == "type_1_misclass" & ethnicity_decoded=="Black"] == 1)),
        c(sum(group == "type_1_misclass" & ethnicity_decoded=="White"), sum(group == "type_1_misclass" & ethnicity_decoded=="Black"))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  print(paste0("p-value Black vs White DKA: ", p_val_dka))
  print(paste0("p-value Black vs White hypoglycaemia: ", p_val_hypo))
  print(paste0("p-value Black vs White insulin 1 year: ", p_val_ins_1_year))
  print(paste0("p-value Black vs White bolus insulin 10 years: ", p_val_bolus_10_year))
  
  
  print("Reference type 1")
  
  
  ## South Asian vs White
  p_val_dka <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(dka_post_diagnosis[group == "ref_type_1" & ethnicity_decoded=="White"] == 1), sum(dka_post_diagnosis[group == "ref_type_1" & ethnicity_decoded=="South Asian"] == 1)),
      c(sum(follow_up_time[group == "ref_type_1" & ethnicity_decoded=="White"]), sum(follow_up_time[group == "ref_type_1" & ethnicity_decoded=="South Asian"]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_hypo <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(hypo_post_diagnosis[group == "ref_type_1" & ethnicity_decoded=="White"] == 1), sum(hypo_post_diagnosis[group == "ref_type_1" & ethnicity_decoded=="South Asian"] == 1)),
      c(sum(follow_up_time[group == "ref_type_1" & ethnicity_decoded=="White"]), sum(follow_up_time[group == "ref_type_1" & ethnicity_decoded=="South Asian"]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_ins_1_year <- cohort_dataset %>%
    filter(!is.na(ins_1_year)) %>%
    summarise(p = ifelse(
      sum(group == "ref_type_1" & ethnicity_decoded=="White") == 0 | sum(group == "ref_type_1"  & ethnicity_decoded=="South Asian") == 0,
      NA,
      prop.test(
        c(sum(ins_1_year[group == "ref_type_1" & ethnicity_decoded=="White"] == 1), sum(ins_1_year[group == "ref_type_1" & ethnicity_decoded=="South Asian"] == 1)),
        c(sum(group == "ref_type_1" & ethnicity_decoded=="White"), sum(group == "ref_type_1" & ethnicity_decoded=="South Asian"))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_bolus_10_year <- cohort_dataset %>%
    filter(!is.na(ten_yrs_post_diag)) %>%
    summarise(p = ifelse(
      sum(group == "ref_type_1" & ethnicity_decoded=="White") == 0 | sum(group == "ref_type_1"  & ethnicity_decoded=="South Asian") == 0,
      NA,
      prop.test(
        c(sum(ten_yrs_bolus_insulin[group == "ref_type_1" & ethnicity_decoded=="White"] == 1), sum(ten_yrs_bolus_insulin[group == "ref_type_1" & ethnicity_decoded=="South Asian"] == 1)),
        c(sum(group == "ref_type_1" & ethnicity_decoded=="White"), sum(group == "ref_type_1" & ethnicity_decoded=="South Asian"))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  print(paste0("p-value South Asian vs White DKA: ", p_val_dka))
  print(paste0("p-value South Asian vs White hypoglycaemia: ", p_val_hypo))
  print(paste0("p-value South Asian vs White insulin 1 year: ", p_val_ins_1_year))
  print(paste0("p-value South Asian vs White bolus insulin 10 years: ", p_val_bolus_10_year))
  
  
  ## Black vs White
  p_val_dka <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(dka_post_diagnosis[group == "ref_type_1" & ethnicity_decoded=="White"] == 1), sum(dka_post_diagnosis[group == "ref_type_1" & ethnicity_decoded=="Black"] == 1)),
      c(sum(follow_up_time[group == "ref_type_1" & ethnicity_decoded=="White"]), sum(follow_up_time[group == "ref_type_1" & ethnicity_decoded=="Black"]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_hypo <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(hypo_post_diagnosis[group == "ref_type_1" & ethnicity_decoded=="White"] == 1), sum(hypo_post_diagnosis[group == "ref_type_1" & ethnicity_decoded=="Black"] == 1)),
      c(sum(follow_up_time[group == "ref_type_1" & ethnicity_decoded=="White"]), sum(follow_up_time[group == "ref_type_1" & ethnicity_decoded=="Black"]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_ins_1_year <- cohort_dataset %>%
    filter(!is.na(ins_1_year)) %>%
    summarise(p = ifelse(
      sum(group == "ref_type_1" & ethnicity_decoded=="White") == 0 | sum(group == "ref_type_1"  & ethnicity_decoded=="Black") == 0,
      NA,
      prop.test(
        c(sum(ins_1_year[group == "ref_type_1" & ethnicity_decoded=="White"] == 1), sum(ins_1_year[group == "ref_type_1" & ethnicity_decoded=="Black"] == 1)),
        c(sum(group == "ref_type_1" & ethnicity_decoded=="White"), sum(group == "ref_type_1" & ethnicity_decoded=="Black"))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_bolus_10_year <- cohort_dataset %>%
    filter(!is.na(ten_yrs_post_diag)) %>%
    summarise(p = ifelse(
      sum(group == "ref_type_1" & ethnicity_decoded=="White") == 0 | sum(group == "ref_type_1"  & ethnicity_decoded=="Black") == 0,
      NA,
      prop.test(
        c(sum(ten_yrs_bolus_insulin[group == "ref_type_1" & ethnicity_decoded=="White"] == 1), sum(ten_yrs_bolus_insulin[group == "ref_type_1" & ethnicity_decoded=="Black"] == 1)),
        c(sum(group == "ref_type_1" & ethnicity_decoded=="White"), sum(group == "ref_type_1" & ethnicity_decoded=="Black"))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  print(paste0("p-value Black vs White DKA: ", p_val_dka))
  print(paste0("p-value Black vs White hypoglycaemia: ", p_val_hypo))
  print(paste0("p-value Black vs White insulin 1 year: ", p_val_ins_1_year))
  print(paste0("p-value Black vs White bolus insulin 10 years: ", p_val_bolus_10_year))
}



compare_t1_deprivation <- function(cohort_dataset) {
  
  cohort_dataset <- cohort_dataset %>% filter(!is.na(imd_quintiles))
  
  print("Misclassified")
  
  p_val_dka <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(dka_post_diagnosis[group == "type_1_misclass" & imd_quintiles==1] == 1), sum(dka_post_diagnosis[group == "type_1_misclass" & imd_quintiles==5] == 1)),
      c(sum(follow_up_time[group == "type_1_misclass" & imd_quintiles==1]), sum(follow_up_time[group == "type_1_misclass" & imd_quintiles==5]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_hypo <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(hypo_post_diagnosis[group == "type_1_misclass" & imd_quintiles==1] == 1), sum(hypo_post_diagnosis[group == "type_1_misclass" & imd_quintiles==5] == 1)),
      c(sum(follow_up_time[group == "type_1_misclass" & imd_quintiles==1]), sum(follow_up_time[group == "type_1_misclass" & imd_quintiles==5]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_ins_1_year <- cohort_dataset %>%
    filter(!is.na(ins_1_year)) %>%
    summarise(p = ifelse(
      sum(group == "type_1_misclass" & imd_quintiles==1) == 0 | sum(group == "type_1_misclass"  & imd_quintiles==5) == 0,
      NA,
      prop.test(
        c(sum(ins_1_year[group == "type_1_misclass" & imd_quintiles==1] == 1), sum(ins_1_year[group == "type_1_misclass" & imd_quintiles==5] == 1)),
        c(sum(group == "type_1_misclass" & imd_quintiles==1), sum(group == "type_1_misclass" & imd_quintiles==5))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_bolus_10_year <- cohort_dataset %>%
    filter(!is.na(ten_yrs_post_diag)) %>%
    summarise(p = ifelse(
      sum(group == "type_1_misclass" & imd_quintiles==1) == 0 | sum(group == "type_1_misclass"  & imd_quintiles==5) == 0,
      NA,
      prop.test(
        c(sum(ten_yrs_bolus_insulin[group == "type_1_misclass" & imd_quintiles==1] == 1), sum(ten_yrs_bolus_insulin[group == "type_1_misclass" & imd_quintiles==5] == 1)),
        c(sum(group == "type_1_misclass" & imd_quintiles==1), sum(group == "type_1_misclass" & imd_quintiles==5))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  print(paste0("p-value IMD quintile 1 vs 5 DKA: ", p_val_dka))
  print(paste0("p-value IMD quintile 1 vs 5 hypoglycaemia: ", p_val_hypo))
  print(paste0("p-value IMD quintile 1 vs 5 insulin 1 year: ", p_val_ins_1_year))
  print(paste0("p-value IMD quintile 1 vs 5 bolus insulin 10 years: ", p_val_bolus_10_year))
  
  
  print("Reference type 1")
  
  p_val_dka <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(dka_post_diagnosis[group == "ref_type_1" & imd_quintiles==1] == 1), sum(dka_post_diagnosis[group == "ref_type_1" & imd_quintiles==5] == 1)),
      c(sum(follow_up_time[group == "ref_type_1" & imd_quintiles==1]), sum(follow_up_time[group == "ref_type_1" & imd_quintiles==5]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_hypo <- cohort_dataset %>%
    summarise(p = poisson.test(
      c(sum(hypo_post_diagnosis[group == "ref_type_1" & imd_quintiles==1] == 1), sum(hypo_post_diagnosis[group == "ref_type_1" & imd_quintiles==5] == 1)),
      c(sum(follow_up_time[group == "ref_type_1" & imd_quintiles==1]), sum(follow_up_time[group == "ref_type_1" & imd_quintiles==5]))
    )$p.value) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_ins_1_year <- cohort_dataset %>%
    filter(!is.na(ins_1_year)) %>%
    summarise(p = ifelse(
      sum(group == "ref_type_1" & imd_quintiles==1) == 0 | sum(group == "ref_type_1"  & imd_quintiles==5) == 0,
      NA,
      prop.test(
        c(sum(ins_1_year[group == "ref_type_1" & imd_quintiles==1] == 1), sum(ins_1_year[group == "ref_type_1" & imd_quintiles==5] == 1)),
        c(sum(group == "ref_type_1" & imd_quintiles==1), sum(group == "ref_type_1" & imd_quintiles==5))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  p_val_bolus_10_year <- cohort_dataset %>%
    filter(!is.na(ten_yrs_post_diag)) %>%
    summarise(p = ifelse(
      sum(group == "ref_type_1" & imd_quintiles==1) == 0 | sum(group == "ref_type_1"  & imd_quintiles==5) == 0,
      NA,
      prop.test(
        c(sum(ten_yrs_bolus_insulin[group == "ref_type_1" & imd_quintiles==1] == 1), sum(ten_yrs_bolus_insulin[group == "ref_type_1" & imd_quintiles==5] == 1)),
        c(sum(group == "ref_type_1" & imd_quintiles==1), sum(group == "ref_type_1" & imd_quintiles==5))
      )$p.value)) %>%
    mutate(p=ifelse(p<0.001, "<0.001", round_pad(p, 3))) %>%
    pull(p)
  
  print(paste0("p-value IMD quintile 1 vs 5 DKA: ", p_val_dka))
  print(paste0("p-value IMD quintile 1 vs 5 hypoglycaemia: ", p_val_hypo))
  print(paste0("p-value IMD quintile 1 vs 5 insulin 1 year: ", p_val_ins_1_year))
  print(paste0("p-value IMD quintile 1 vs 5 bolus insulin 10 years: ", p_val_bolus_10_year))
  
}