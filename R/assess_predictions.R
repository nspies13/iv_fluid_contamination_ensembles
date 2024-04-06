aggregateTargetObjects <- function(train_cohort = "BJH", panel = "BMP", fluid_list = c("NS", "LR", "D5W", "Tech"), step_name = "wf_regression_tuned", recipe_name = "Base", feature_set = "results_with_priors", models_string = "XGB"){
  
  meta <- tar_meta()
  
  tar_names <- 
    meta %>% 
      filter(
        grepl(train_cohort, meta$name) & 
          grepl(panel, meta$name) &
          grepl(paste0(paste0("_", fluid_list, "_"), collapse = "|"), meta$name) &
          grepl(step_name, meta$name) &
          grepl(recipe_name, meta$name) &
          grepl(feature_set, name) &
          grepl(models_string, name)) %>%
      pluck("name")
  
  map(tar_names, ~qs::qread(paste0("_targets/objects/", .x))) |> set_names(tar_names)
  
}

getManualReviewLabels <- function(preds = tar_read(wf_prospective_BJH_Tech_BMP_Base_NNet_results_with_priors), consensus = readxl::read_xlsx("../../Data/results_to_review_formatted_preconsensus.xlsx")){
  
  first <- consensus %>% arrange(drawn_dt_tm) %>% filter(time_point == "Current") %>% slice_head(n = 1) %>% pluck("drawn_dt_tm")
  last <-  consensus %>% arrange(drawn_dt_tm) %>% filter(time_point == "Current") %>% slice_tail(n = 1) %>% pluck("drawn_dt_tm")
  
  reviewed <- 
    preds %>% 
      filter(between(drawn_dt_tm, first, last)) |>
      left_join(consensus %>% filter(time_point == "Current") %>% select(epic_id, drawn_dt_tm, final_prediction), by = c("epic_id", "drawn_dt_tm")) |> 
      mutate(
        final_prediction = ifelse(
          calcium > 0.95 * calcium_prior & co2_totl > 0.95 * co2_totl_prior, 0, final_prediction    
        )) |>
      drop_na() |>
      mutate(final_prediction = factor(final_prediction, levels = c(0, 1)))

  reviewed
  
}

makePipelinePredictions <- function(sim_preds_raw = aggregateTargetObjects(train_cohort = "BJH", panel = "BMP", fluid_list = c("NS", "LR", "D5W"), step_name = "wf_yale_reviewed", recipe_name = "Base", feature_set = "results_with_priors", models_string = "XGB|NNet"), tech_preds_raw = aggregateTargetObjects(train_cohort = "BJH", panel = "BMP", fluid_list = c("Tech"), step_name = "wf_yale_reviewed", recipe_name = "Base", feature_set = "results_with_priors", models_string = "XGB|NNet")){
  
  sim_preds <- map2(sim_preds_raw, names(sim_preds_raw), ~select(.x, .pred_class) |> set_names(paste0(.y, "_pred_class"))) |> bind_cols() |> 
    mutate(across(everything(), ~as.numeric(as.character(.))))
  sim_pred_count <- rowSums(sim_preds, na.rm = T)
  sim_pred_final <- factor(as.numeric(sim_pred_count > 1))
  
  tech_preds <- map2(tech_preds_raw, names(tech_preds_raw), ~select(.x, .pred_class) |> set_names(paste0(.y, "_pred_class"))) |> bind_cols() |> 
    mutate(across(everything(), ~as.numeric(as.character(.))))
  tech_pred_count <- rowSums(tech_preds, na.rm = T)
  tech_pred_final <- factor(as.numeric(tech_pred_count > 1))
  
  
  base <- sim_preds_raw[[1]] |> select(-.pred_class, -.pred_0, -.pred_1)
  
  preds <-
    bind_cols(base, sim_pred = sim_pred_final, tech_pred = tech_pred_final) |> 
      mutate(.pred_class = factor(as.numeric(sim_pred == 1 | tech_pred == 1)))
  
  preds
  
}

getMetricsAndConfusionMatrix <- function(input = pipeline_preds){
  
  list(
    data = input, 
    preds = input$.pred_class,
    metrics = metrics_binary(input, truth = final_prediction, estimate = .pred_class, event_level = "second") |> transmute(Metric = .metric, Value = round(.estimate, digits = 3)),
    conf_mat = conf_mat(input, truth = final_prediction, estimate = .pred_class))
  
}

makeMixtureRatioPredictions <- function(input = tar_read(wf_reviewed_BJH_Sim_BMP_NS_Base_XGB_all), fluid_list = fluid_names){
  
  models <- map(fluid_list, ~aggregateTargetObjects(fluid_list = .x, step_name = "wf_regression_fit", feature_set = "all", models_string = "XGB") |> set_names(.x)) |> flatten() |> map(bundle::unbundle)
  
  preds <- map2(models, names(models), ~predict(.x, input) |> set_names(.y)) |> bind_cols()
  
  max_mixtures <- apply(preds, 1, max)
  max_fluids <- apply(preds, 1, function(x) fluid_list[which.max(x)])
  
  final_calls <- bind_cols(max_mixture = max_mixtures, max_fluid = max_fluids) |> transmute(mix_ratio_pred = ifelse(max_mixture < 0, NA, max_mixture), contaminating_fluid_pred = max_fluid)
  
  output <- bind_cols(input, final_calls)
  
  output
  
}
  
joinResultMetadata <- function(input = read_csv("../../Results/Predictions/Supervised/final_pipeline_internal_preds.csv"), result_metadata = read_csv("../../Data/prospective_bmp_result_metadata.csv")){
  
    preds_meta <- 
      left_join(
        input |> select(matches("_id|draw|comment"), lab_strings_bmp_no_gap, matches("pred")) |> pivot_longer(cols = lab_strings_bmp_no_gap, names_to = "task_assay", values_to = "result_val"), 
        result_metadata |> distinct()
      )
 
    preds_meta
     
}

getResultFlags <- function(input = read_csv("../../Results/Predictions/Supervised/final_pipeline_prospective_preds.csv"), result_metadata = read_csv("../../Data/prospective_bmp_result_metadata.csv"), out_path = "../../Results/Predictions/Supervised/full_pipeline_prospective_with_result_metadata_flags.csv"){
  
  preds_with_result_metadata <- joinResultMetadata(input, result_metadata)
  
  pred_flags <- 
    preds_with_result_metadata %>% 
      mutate(flagged_by_technologist = factor(as.numeric(unlikely_comment|contam_comment)), 
             flag = case_when(criticality == "Critical" ~ "Critical", normality %in% c("Low", "High") ~ "Abnormal", T ~ "Normal"))
  
  result_preds_meta_raw <- 
    pred_flags %>%
      group_by(patient_id, task_assay) %>%
      mutate(time_to_next = dplyr::lead(drawn_dt_tm) - drawn_dt_tm, 
             next_flag = dplyr::lead(flag), 
             time_from_prior = drawn_dt_tm - dplyr::lag(drawn_dt_tm),
             prior_flag = dplyr::lag(flag)) %>%
      ungroup()
  
  result_preds_meta_flags <- 
    result_preds_meta_raw %>% 
      mutate(time_to_next = time_to_next / 60 / 60) %>% 
      mutate(flag = factor(flag, levels = c("Normal", "Abnormal", "Critical")), 
             next_flag = factor(next_flag, levels = c("Normal", "Abnormal", "Critical")), 
             flag_bin = factor(case_when(.pred_class == 1 & flagged_by_technologist == 1 ~ "Both", 
                                         .pred_class == 1 & flagged_by_technologist == 0 ~ "Only ML", 
                                         .pred_class == 0 & flagged_by_technologist == 1 ~ "Only Tech", 
                                         .pred_class == 0 & flagged_by_technologist == 0 ~ "Neither"), 
                               levels = c("Neither", "Only ML", "Only Tech", "Both")))
  
  write_csv(result_preds_meta_flags, out_path)
  
  result_preds_meta_flags

}

tmp2 <- function(){  
  
  floors <- preds_meta$draw_nursing
  preds_meta$floor_type <- case_when(floors %in% c("I044", "NI104", "I084", "94ICU", "I056", "I104", "7800S", "7800M", "83CTI", "I155") ~ "Intensive Care", 
                                     floors %in% c("58LD", "58AP", "68") ~ "Obstetrics",
                                     floors %in% c("8800", "9800", "10800", "11800", "12800", "0164", "0102", "0092", "6900", "0112", "0082", "0122") ~ "Inpatient",
                                     grepl("POD", floors) ~ "Operating Room",
                                     grepl("ER|ED", floors) ~ "Emergency",
                                     floors %in% c("BJHOP", "BJH HIGHLANDS", "BJH CAM8", "BJH OPLAB", "BJH SCCSC", "BJH CAM") ~ "Outpatient")  
  
  floor_realtime_props <- 
    preds_meta %>% select(floor_type, pred_realtime) %>% drop_na(floor_type, pred_realtime) %>% group_by(floor_type) %>% count(pred_realtime) %>% mutate(prop = n/sum(n) * 10000) %>% filter(pred_realtime == "Contaminated") %>% arrange(desc(prop)) %>% mutate(floor_type = factor(floor_type), type = "Real-Time") %>% select(floor_type, type, prop)
  
  floor_retro_props <- 
    preds_meta %>% select(floor_type, pred_retrospective) %>% drop_na(floor_type, pred_retrospective) %>% group_by(floor_type) %>% count(pred_retrospective) %>% mutate(prop = n/sum(n) * 10000) %>% filter(pred_retrospective == "Contaminated") %>% arrange(desc(prop)) %>% mutate(floor_type = factor(floor_type), type = "Retrospective") %>% select(floor_type, type, prop)
  
  floor_prop_input <- bind_rows(floor_realtime_props, floor_retro_props, as_tibble_row(list(floor_type = "Outpatient", type = "Retrospective", prop = 0)))
  prop_segments <- pivot_wider(floor_prop_input, names_from = "type", values_from = "prop") %>% janitor::clean_names()
  
  
  gg_floor_props <- 
    ggplot(floor_prop_input, aes(x = prop, y = fct_reorder(floor_type, prop, .desc = F), color = type)) + 
    geom_segment(data = prop_segments, aes(x = retrospective, xend = real_time, y = fct_reorder(floor_type, retrospective, .desc = F), yend = fct_reorder(floor_type, retrospective, .desc = F)), color = scico::scico(n = 1, palette = "bilbao", begin = 0.8), linewidth = 2) +
    geom_point(size = 4) +
    geom_text(aes(label = round(prop), x = prop + 1, group = type), position = position_dodge(width = 1.25), fontface = "bold.italic", hjust = 0.5, size = 3) +
    scale_x_continuous(name = "Predicted Positives per 10000 Specimens", breaks = c(0, 100, 200), limits = c(0, 200)) +
    scale_y_discrete(name = "") +
    annotate("text", x = 74, y = 5, label = "Real-Time", fontface = "bold.italic", hjust = 1, color = scico::scico(n = 1, palette = "bilbao", begin = 0.8)) +
    annotate("text", x = 113, y = 5, label = "Retrospective", fontface = "bold.italic", hjust = 0, color = scico::scico(n = 1, palette = "bilbao", begin = 0.2)) +
    scico::scale_color_scico_d(palette = "bilbao", begin = 0.2, end = 0.8, direction = -1) +
    facet_wrap(~"Context") + 
    theme(legend.position = "none", axis.text.y.left = element_text(face = "bold", color = "black"))
  
  floor_realtime_counts <- 
    preds_meta %>% 
    group_by(draw_nursing) %>%
    count(pred_retrospective) %>% 
    drop_na() %>% 
    mutate(prop = n/sum(n)) %>% 
    filter(pred_retrospective == "Contaminated" & n > 5) %>%
    arrange(desc(prop))
  write_csv(floor_realtime_counts, "../../Results/Predictions/Supervised/floor_contamination_rates.csv")
  
  service_realtime_counts <- 
    preds_meta %>% 
    group_by(draw_med_service) %>%
    count(pred_retrospective) %>% 
    drop_na() %>% 
    mutate(prop = n/sum(n)) %>% 
    filter(pred_retrospective == "Contaminated" & n > 5) %>%
    arrange(desc(prop))
  write_csv(service_realtime_counts, "../../Results/Predictions/Supervised/service_contamination_rates.csv")
  
  encounter_realtime_counts <- 
    preds_meta %>% 
    group_by(encounter_type) %>%
    count(pred_retrospective) %>% 
    drop_na() %>% 
    mutate(prop = n/sum(n)) %>% 
    filter(pred_retrospective == "Contaminated") %>%
    arrange(desc(prop))
  write_csv(encounter_realtime_counts, "../../Results/Predictions/Supervised/encounter_contamination_rates.csv")
  
  service_realtime_counts <- 
    preds_meta %>% 
    group_by(admit_type) %>%
    count(pred_retrospective) %>% 
    drop_na() %>% 
    mutate(prop = n/sum(n)) %>% 
    filter(pred_retrospective == "Contaminated") %>%
    arrange(desc(prop))
  write_csv(floor_realtime_counts, "../../Results/Predictions/Supervised/service_contamination_rates.csv")
  
  floor_realtime_counts <- 
    preds_meta %>% 
    group_by(draw_nursing) %>%
    count(pred_realtime) %>% 
    drop_na() %>% 
    mutate(prop = n/sum(n)) %>% 
    filter(pred_realtime == "Contaminated" & n > 5) %>%
    arrange(desc(prop))
  write_csv(floor_realtime_counts, "../../Results/Predictions/Supervised/floor_contamination_realtime_rates.csv")
  
  service_realtime_counts <- 
    preds_meta %>% 
    group_by(draw_med_service) %>%
    count(pred_realtime) %>% 
    drop_na() %>% 
    mutate(prop = n/sum(n)) %>% 
    filter(pred_realtime == "Contaminated" & n > 5) %>%
    arrange(desc(prop))
  write_csv(service_realtime_counts, "../../Results/Predictions/Supervised/service_contamination_realtime_rates.csv")
  
  encounter_realtime_counts <- 
    preds_meta %>% 
    group_by(encounter_type) %>%
    count(pred_realtime) %>% 
    drop_na() %>% 
    mutate(prop = n/sum(n)) %>% 
    filter(pred_realtime == "Contaminated") %>%
    arrange(desc(prop))
  write_csv(encounter_realtime_counts, "../../Results/Predictions/Supervised/encounter_contamination_realtime_rates.csv")
  
  service_realtime_counts <- 
    preds_meta %>% 
    group_by(admit_type) %>%
    count(pred_realtime) %>% 
    drop_na() %>% 
    mutate(prop = n/sum(n)) %>% 
    filter(pred_realtime == "Contaminated") %>%
    arrange(desc(prop))
  write_csv(floor_realtime_counts, "../../Results/Predictions/Supervised/service_contamination_realtime_rates.csv")
  
}

tmp <- function(){
  
  manual_review_final <- joined %>% drop_na(expert_label)
  
  fit_workflows <- aggregateFitWorkflows()
  
  probs_reviewed <- map(fit_workflows, ~predict(.x, manual_review_final, type = "prob") %>% select(.pred_1) %>% bind_cols(., expert_label = manual_review_final$expert_label)) %>% setNames(paste0("prob_", fit_workflows %>% map("metadata") %>% map("pred_set") %>% unlist()))
  
  
  
  contam_sim_all <- read_csv("../../Data/contam_sim_all.csv") 
  sim_preds <- predict(fit_workflows[[1]], contam_sim_all)
  
  gg_sim_input <- bind_cols(contam_sim_all %>% select(label, mix_ratio), sim_preds) %>%
    mutate(label = ifelse(grepl("D5", label), "D5", label), label = ifelse(label == "Water", "water", label)) %>%
    group_by(label, mix_ratio) %>% 
    count(.pred_class) %>% 
    mutate(prop = n/sum(n)) %>% 
    filter(.pred_class == 1 & mix_ratio < 0.5)
  
  lines <- gg_sim_input %>% filter(mix_ratio == 0.1) %>% arrange(mix_ratio) %>% ungroup() %>% slice_head(n = 1, by = label) %>% mutate(label = fct_reorder(label, prop, .desc = T))
  
  gg_sim_input <- gg_sim_input %>% ungroup() %>%  mutate(label = factor(label, levels = c(lines$label)))
  
  library(geomtextpath)
  gg_sim_sensitivity <- 
    ggplot(gg_sim_input, aes(mix_ratio, prop, color = label, label = label)) +
      geom_vline(xintercept = 0.1, linetype = "dashed", alpha = 0.25) +  
      geom_textsegment(data = lines, aes(x = 0, xend = 0.1, y = prop, yend = prop), hjust = 0.05, size = 3, fontface = "bold.italic", linetype = "dashed") + 
      geom_smooth(linewidth = 2, size = 6, hjust = 0.4, fontface = "bold", se = F, span = 0.25) +
      scico::scale_color_scico_d(palette = "davos", end = 0.8) + 
      scale_x_continuous(name = "Simulated Mixture Ratio", expand = c(0, 0), limits = c(0, 0.3)) +
      scale_y_continuous(name = "Sensitivity", expand = c(0, 0), limits = c(0, 1.02)) +
      theme(legend.position = "none", plot.margin = margin(0,0,0,0)) + coord_equal(ratio = .3)
  
  gg_pre_post_manual <- map(lab_strings_bmp_no_gap, ~makePrePostPlots(input = manual_review_final, analytes_to_show = .x, rows_to_show = nrow(filter(manual_review_final, expert_label == "Contaminated")), group_col = "expert_label", out_path = paste0("manual_review_pre_post_", .x, "_retrospective.pdf")))
  
  library(ggpubr)
  gg_pre_post_manual_all <- ggarrange(plotlist = map(gg_pre_post_manual, ~.x + theme(axis.title.y.left = element_blank())), nrow = 1, ncol = length(gg_pre_post_manual), common.legend = T, legend = "bottom")
  ggsave("../../Figures/Supervised/pre_post_tiles_all_manual_review_retrospective.pdf", gg_pre_post_manual_all, width = 17, height = 6)

  library(ggpubr)
  gg_pre_post_manual_subset <- ggarrange(plotlist = map(gg_pre_post_manual[c(2,7,8)], ~.x + theme(axis.title.y.left = element_blank())), nrow = 1, ncol = 3, common.legend = T, legend = "bottom")
  ggsave("../../Figures/Supervised/pre_post_tiles_subset_manual_review_retrospective.pdf", gg_pre_post_manual_subset, width = 8.5, height = 4)
  
  
  ggpubr::ggarrange(gg_rocs, gg_prs, gg_sim_sensitivity, align = "hv", nrow = 1, labels = "AUTO")
  ggsave("../../Figures/Supervised/curves_final_reviewed.pdf", width = 12, height = 6)
  
  yale <- read_csv("../../Data/yale_wide_with_delta_flags.csv")
  yale_preds_list <- map(fit_workflows, ~predict(.x, yale)) %>% bind_cols() %>% setNames(paste0("pred_", fit_workflows %>% map("metadata") %>% map("pred_set") %>% unlist()))
  
  yale_preds <- 
    bind_cols(yale, yale_preds_list %>% 
    mutate(pred_realtime = ifelse(is.na(pred_current_with_delta), pred_results, pred_current_with_delta), 
              pred_retrospective = ifelse(is.na(pred_all), pred_current_with_post, pred_all)) %>% 
    mutate(across(any_of(paste0("pred_", fit_workflows %>% map("metadata") %>% map("pred_set") %>% unlist())), ~factor(ifelse(. == 0, "Real", "Contaminated"), levels = c("Real", "Contaminated"), labels = c("Real", "Contaminated"))))) %>% 
    mutate(across(matches("retro|realt"), ~factor(ifelse(. == 1, "Real", "Contaminated"), levels = c("Real", "Contaminated"), labels = c("Real", "Contaminated"))))
  
  
  yale_preds <- yale_preds %>% rename(epic_id = patient_id)
  
  yale_preds %>% 
    filter(!is.na(pred_all)) %>%
    arrange(epic_id) %>% 
    ungroup() %>% 
    reformatReviews() %>% 
    select(epic_id, drawn_dt_tm, time_point, any_of(lab_strings_bmp), final_prediction, contaminating_fluid) %>%
    left_join(yale_preds %>% select(epic_id, specimen_id, any_of(lab_strings))) %>% 
    select(epic_id, specimen_id, time_point, any_of(lab_strings_bmp), final_prediction, contaminating_fluid) %>% 
    write_csv("../../Results/Predictions/Supervised/yale_results_to_review_all.csv")
  
  yale_pre_post_realtime <- map(lab_strings_bmp_no_gap, ~makePrePostPlots(input = yale_preds, group_col = "pred_realtime", analytes_to_show = .x, rows_to_show = nrow(filter(yale_preds, pred_realtime == "Contaminated") %>% drop_na()), out_path = paste0("../../Figures/Supervised/pre_post_tiles_", .x, "_yale_realtime.pdf")))
  yale_pre_post_retrospective <- map(lab_strings_bmp_no_gap, ~makePrePostPlots(input = yale_preds, group_col = "pred_retrospective", analytes_to_show = .x, rows_to_show = nrow(filter(yale_preds, pred_retrospective == "Contaminated") %>% drop_na()), out_path = paste0("../../Figures/Supervised/pre_post_tiles_", .x, "_yale_retrospective.pdf")))

  library(ggpubr)
  gg_yale_realtime <- ggarrange(plotlist = map(yale_pre_post_realtime, ~.x + theme(axis.title.y.left = element_blank())), nrow = 1, ncol = length(yale_pre_post_realtime), common.legend = T, legend = "bottom")
  ggsave("../../Figures/Supervised/pre_post_tiles_all_yale_realtime.pdf", gg_yale_realtime, width = 17, height = 6)
  
  gg_yale_retro <- ggarrange(plotlist = map(yale_pre_post_retrospective, ~.x + theme(axis.title.y.left = element_blank())), nrow = 1, ncol = length(yale_pre_post_retrospective), common.legend = T, legend = "bottom")
  ggsave("../../Figures/Supervised/pre_post_tiles_all_yale_retrospective.pdf", gg_yale_retro, width = 17, height = 6)

  gg_yale_retro_subset <- ggarrange(plotlist = map(yale_pre_post_retrospective[c(2,7,8)], ~.x + theme(axis.title.y.left = element_blank())), nrow = 1, ncol = 3, common.legend = T, legend = "bottom")
  ggsave("../../Figures/Supervised/pre_post_tiles_subset_yale_retrospective.pdf", gg_yale_retro_subset, width = 8.5, height = 4)
  
    
  
  
  
  
  
  library(zoo)
  preds_rolling <- 
    preds_meta %>% 
      mutate(month = floor_date(drawn_dt_tm, unit = "month")) %>%
      group_by(month) %>%
      count(pred_retrospective) %>%
      drop_na() %>%
      mutate(prop = n/sum(n)) %>%
      filter(pred_retrospective == "Contaminated")
  
  ggplot(preds_rolling, aes(month, prop)) +
    geom_smooth()
  ggsave("../../Figures/Supervised/pred_rate_over_time.tiff", width = 8.5, height = 3)
  
  
  library(applicable)
  train <- tar_read(train_input)[[1]][["data"]]
  train_pca <- apd_pca(train %>% select(lab_strings_bmp_no_gap))
  
  scores <- score(train_pca, train %>% select(lab_strings_bmp_no_gap))
  
  
  
  
  
  timings <- 
    result_preds_meta %>% mutate(drawn_to_reported = perform_result_updt_dt_tm - drawn_dt_tm,
                                    performed_to_reported = perform_result_updt_dt_tm - perform_dt_tm)
  
  ed_timings_last_update <- timings %>% filter(grepl("ED", draw_nursing) & task_assay == "potassium_plas") %>% arrange(perform_result_updt_dt_tm) %>% group_by(specimen_id, task_assay) %>% summarise(drawn_to_last_reported = max(drawn_to_reported), performed_to_last_reported = max(performed_to_reported)) %>% left_join(preds_meta %>% select(specimen_id, pred_realtime, pred_retrospective))
  
  ed_timings_last_update %>% group_by(pred_realtime) %>% summarize(median = median(performed_to_last_reported), mean = mean(performed_to_last_reported))
  
}

assessDemographicFlagRates <- function(input = prospective_preds, reviewed = pipeline_preds, specimen_metadata = read_csv("../../Data/prospective_bmp_result_metadata.csv") |> select(person_id = patient_id, matches("id|draw_|encounter|admit"), drawn_dt_tm) |>  distinct(), demographics = read_csv("../../../BMP/Data/20230822/demographics.csv") %>% janitor::clean_names() %>% select(-x1)){

specimen_preds <- 
  input %>% 
    select(specimen_id, drawn_dt_tm, matches("pred")) %>%
    left_join(reviewed |> select(specimen_id, final_prediction)) |> 
    left_join(specimen_metadata) %>% 
    left_join(demographics) %>%
    left_join(input |> transmute(specimen_id, flagged_by_tech = unlikely_comment|contam_comment)) |> 
    mutate(age = as.numeric((drawn_dt_tm - birth_dt_tm) / 365.25))

demo_input <- 
  specimen_preds %>% 
    mutate(race = factor(fct_other(race, keep = c("White", "Black")), levels = c("White", "Black")), 
           sex = factor(sex, levels = c("Male", "Female")),
           age_bin = factor(case_when(age <= 65 ~ "<65", age > 65 ~ ">=65"), levels = c("<65", ">=65")),
           encounter = fct_infreq(fct_lump_n(encounter_type, 3)), 
           medical_service = fct_infreq(fct_lump_n(draw_med_service, 10)))

floors <- demo_input$draw_nursing
demo_input$floor_type <- fct_infreq(fct_other(factor(case_when(floors %in% c("I044", "NI104", "I084", "94ICU", "I056", "I104", "7800S", "7800M", "83CTI", "I155") ~ "Intensive Care", 
                                   floors %in% c("58LD", "58AP", "68") ~ "Obstetrics",
                                   floors %in% c("8800", "9800", "10800", "11800", "12800", "0164", "0102", "0092", "6900", "0112", "0082", "0122") ~ "Inpatient",
                                   grepl("POD", floors) ~ "Operating Room",
                                   grepl("ER|ED", floors) ~ "Emergency")), keep = c("Intensive Care", "Obstetrics", "Inpatient", "Operating Room", "Emergency"))) 

library(gt)
library(gtsummary)

gt_univariate <- 
  tbl_uvregression(demo_input |> select(.pred_class, age, race, sex, floor_type) |> drop_na(), y = .pred_class, method = glm, method.args = list(family = binomial), exponentiate = T, add_estimate_to_reference_rows = T, label = c("race" ~ "Race", "age" ~ "Age (years)", "sex" ~ "Sex", "floor_type" ~ "Unit Type")) |> 
  add_n(location = "level") |> 
  add_nevent(location = "level") |> 
  bold_labels() |> 
  italicize_levels() |> 
  add_q(method = "BH") |> 
  bold_p(q = T) |>
  modify_table_body(
    ~ .x %>% 
      dplyr::mutate(
        stat_nevent = 
          ifelse(
            !is.na(stat_nevent),
            paste0(style_sigfig(stat_nevent / stat_n, scale = 100), "%"),
            NA
          ),
        ci = ifelse(
          !is.na(conf.low) & !is.na(conf.high),
          paste0("[", style_sigfig(conf.low, digits = 2), " - ", style_sigfig(conf.high, digits = 2), "]"),
          NA
        ))
  ) %>%
  modify_cols_merge(
    pattern = "{estimate}, {ci}",
    rows = !is.na(conf.low)
  ) |>
  modify_header(stat_nevent = "**Flag Rate**", estimate = "**OR, [95% CI]**")
gt::gtsave(gt_univariate |> gtsummary::as_gt(), filename = "../../Figures/Supervised/final_figures/supp_table_univariate_demographics_by_prediction.docx")
gt::gtsave(gt_univariate |> gtsummary::as_gt(), filename = "../../Figures/Supervised/final_figures/supp_table_univariate_demographics_by_prediction.html")

gt_univariate_tech <- 
  tbl_uvregression(demo_input |> select(flagged_by_tech, age, race, sex, floor_type) |> drop_na(), y = flagged_by_tech, method = glm, method.args = list(family = binomial), exponentiate = T, add_estimate_to_reference_rows = T, label = c("race" ~ "Race", "age" ~ "Age (years)", "sex" ~ "Sex", "floor_type" ~ "Unit Type")) |> 
  add_n(location = "level") |> 
  add_nevent(location = "level") |> 
  bold_labels() |> 
  italicize_levels() |> 
  add_q(method = "BH") |>
  bold_p(q = T) |>
  modify_table_body(
    ~ .x %>% 
      dplyr::mutate(
        stat_nevent = 
          ifelse(
            !is.na(stat_nevent),
            paste0(style_sigfig(stat_nevent / stat_n, scale = 100), "%"),
            NA
          ),
        ci = ifelse(
          !is.na(conf.low) & !is.na(conf.high),
          paste0("[", style_sigfig(conf.low, digits = 2), " - ", style_sigfig(conf.high, digits = 2), "]"),
          NA
        ))
  ) %>%
  modify_cols_merge(
    pattern = "{estimate}, {ci}",
    rows = !is.na(conf.low)
  ) |>
  modify_header(stat_nevent = "**Flag Rate**", estimate = "**OR, [95% CI]**")
gt::gtsave(gt_univariate_tech |> gtsummary::as_gt(), filename = "../../Figures/Supervised/final_figures/supp_table_univariate_demographics_by_prediction_tech_flags.docx")
gt::gtsave(gt_univariate_tech |> gtsummary::as_gt(), filename = "../../Figures/Supervised/final_figures/supp_table_univariate_demographics_by_prediction_tech_flags.html")

gt_univariate_review <- 
  tbl_uvregression(demo_input |> select(final_prediction, age, race, sex, floor_type) |> drop_na(), y = final_prediction, method = glm, method.args = list(family = binomial), exponentiate = T, add_estimate_to_reference_rows = T, label = c("race" ~ "Race", "age" ~ "Age (years)", "sex" ~ "Sex", "floor_type" ~ "Unit Type")) |> 
  add_n(location = "level") |> 
  add_nevent(location = "level") |> 
  bold_labels() |> 
  italicize_levels() |> 
  add_q(method = "BH") |>
  bold_p(q = T) |>
  modify_table_body(
    ~ .x %>% 
      dplyr::mutate(
        stat_nevent = 
          ifelse(
            !is.na(stat_nevent),
            paste0(style_sigfig(stat_nevent / stat_n, scale = 100), "%"),
            NA
          ),
        ci = ifelse(
          !is.na(conf.low) & !is.na(conf.high),
          paste0("[", style_sigfig(conf.low, digits = 2), " - ", style_sigfig(conf.high, digits = 2), "]"),
          NA
        ))
  ) %>%
  modify_cols_merge(
    pattern = "{estimate}, {ci}",
    rows = !is.na(conf.low)
  ) |>
  modify_header(stat_nevent = "**Flag Rate**", estimate = "**OR, [95% CI]**")
gt::gtsave(gt_univariate_review |> gtsummary::as_gt(), filename = "../../Figures/Supervised/final_figures/supp_table_univariate_demographics_by_prediction_expert_review.docx")
gt::gtsave(gt_univariate_review |> gtsummary::as_gt(), filename = "../../Figures/Supervised/final_figures/supp_table_univariate_demographics_by_prediction_expert_review.html")


gt_combo_uni <- tbl_merge(list(gt_univariate_tech |> modify_column_hide(c("ci", "p.value")), gt_univariate |> modify_column_hide(c("stat_n", "ci", "p.value")), gt_univariate_review |> modify_column_hide(c("stat_n", "ci", "p.value"))), tab_spanner = c("Flagged by Current", "Predicted by Pipeline", "Reviewed by Judges")) |>
  modify_header(label = "") |> 
  as_gt() |>
  tab_header(title = "Univariate Analysis: Demographic Breakdown of Flagged Results", subtitle = "Internal Prospective Validation Set")

gt::gtsave(gt_combo_uni, filename = "../../Figures/Supervised/final_figures/table_uniivariate_demographics_combined.docx")
gt::gtsave(gt_combo_uni, filename = "../../Figures/Supervised/final_figures/table_univariate_demographics_combined.html")



glm_pred <- glm(.pred_class ~ age + race + sex + floor_type, demo_input, family = binomial(link = "logit"))
gt_multi_raw <- tbl_regression(glm_pred, exponentiate = T, add_estimate_to_reference_rows = T, label = c("race" ~ "Race", "age" ~ "Age (years)", "sex" ~ "Sex", "floor_type" ~ "Unit Type"))

gt_multivariate <- 
   gt_multi_raw |> 
    add_n(location = "level") |>
    add_nevent(location = "level") |> 
    add_q(method = "fdr") |> 
    bold_p(q = T) |> 
    bold_labels() |> 
    italicize_levels()
gt::gtsave(gt_multivariate |> gtsummary::as_gt(), filename = "../../Figures/Supervised/final_figures/table_multivariate_demographics_by_prediction.docx")
gt::gtsave(gt_multivariate |> gtsummary::as_gt(), filename = "../../Figures/Supervised/final_figures/table_multivariate_demographics_by_prediction.html")



glm_tech <- glm(flagged_by_tech ~ age + race + sex + floor_type, demo_input, family = binomial(link = "logit")) 
gt_tech_raw <- tbl_regression(glm_tech, exponentiate = T, add_estimate_to_reference_rows = T, label = c("race" ~ "Race", "age" ~ "Age (years)", "sex" ~ "Sex", "floor_type" ~ "Unit Type"))

gt_tech_multivariate <- 
  glm(flagged_by_tech ~ age + race + sex + floor_type, demo_input, family = binomial(link = "logit")) |> 
  tbl_regression(exponentiate = T, add_estimate_to_reference_rows = T, label = c("race" ~ "Race", "age" ~ "Age (years)", "sex" ~ "Sex", "floor_type" ~ "Unit Type")) |> 
  add_n(location = "level") |>
  add_nevent(location = "level") |> 
  add_q(method = "BH") |>
  bold_p(q = T) |>
  bold_labels() |> 
  italicize_levels()
gt::gtsave(gt_tech_multivariate |> gtsummary::as_gt(), filename = "../../Figures/Supervised/final_figures/table_tech_multivariate_demographics_by_prediction.docx")
gt::gtsave(gt_tech_multivariate |> gtsummary::as_gt(), filename = "../../Figures/Supervised/final_figures/table_tech_multivariate_demographics_by_prediction.html")

glm_review <- glm(final_prediction ~ age + race + sex + floor_type, demo_input |> drop_na(final_prediction), family = binomial(link = "logit")) 
gt_review_raw <- tbl_regression(glm_review, exponentiate = T, add_estimate_to_reference_rows = T, label = c("race" ~ "Race", "age" ~ "Age (years)", "sex" ~ "Sex", "floor_type" ~ "Unit Type"))

gt_review_multivariate <- 
  gt_review_raw |> 
    add_n(location = "level") |>
    add_nevent(location = "level") |> 
    add_q(method = "BH") |>
    bold_p(q = T) |>
    bold_labels() |> 
    italicize_levels()
gt::gtsave(gt_review_multivariate |> gtsummary::as_gt(), filename = "../../Figures/Supervised/final_figures/table_review_multivariate_demographics_by_prediction.docx")
gt::gtsave(gt_review_multivariate |> gtsummary::as_gt(), filename = "../../Figures/Supervised/final_figures/table_review_multivariate_demographics_by_prediction.html")





gt_combo <- 
  tbl_merge(list(
   gt_tech_raw |> 
    add_n(location = "level") |> 
    add_nevent(location = "level") |> 
    bold_labels() |> 
    italicize_levels() |> 
    add_q(method = "BH") |>
    bold_p(q = T) |>
     modify_table_body(
       ~ .x %>% 
         dplyr::mutate(
           stat_nevent = 
             ifelse(
               !is.na(stat_nevent),
               paste0(style_sigfig(stat_nevent / stat_n, scale = 100), "%"),
               NA
             ),
           ci = ifelse(
             !is.na(conf.low) & !is.na(conf.high),
             paste0("[", style_sigfig(conf.low, digits = 2), " - ", style_sigfig(conf.high, digits = 2), "]"),
             NA
           ))
     ) %>%
     modify_cols_merge(
       pattern = "{estimate}, {ci}",
       rows = !is.na(conf.low)
     ) |>
     modify_header(stat_nevent = "**Flag Rate**", estimate = "**OR, [95% CI]**") |>
     modify_column_hide(c("ci", "p.value")), 
  gt_multi_raw |> 
    add_n(location = "level") |>
    add_nevent(location = "level") |> 
    bold_labels() |> 
    italicize_levels() |> 
    add_q(method = "BH") |>
    bold_p(q = T) |>
    modify_table_body(
      ~ .x %>% 
        dplyr::mutate(
          stat_nevent = 
            ifelse(
              !is.na(stat_nevent),
              paste0(style_sigfig(stat_nevent / stat_n, scale = 100), "%"),
              NA
            ),
          ci = ifelse(
            !is.na(conf.low) & !is.na(conf.high),
            paste0("[", style_sigfig(conf.low, digits = 2), " - ", style_sigfig(conf.high, digits = 2), "]"),
            NA
          ))
    ) |>
    modify_header(stat_nevent = "**Flag Rate**", estimate = "**OR, [95% CI]**") |>
    modify_column_merge(
      pattern = "{estimate}, {ci}",
      rows = estimate != 1
    ) |> 
    modify_column_hide(c("stat_n", "ci", "p.value")),
  gt_review_raw |> 
    add_n(location = "level") |>
    add_nevent(location = "level") |> 
    bold_labels() |> 
    italicize_levels() |> 
    add_q(method = "BH") |>
    bold_p(q = T) |>
    modify_table_body(
      ~ .x %>% 
        dplyr::mutate(
          stat_nevent = 
            ifelse(
              !is.na(stat_nevent),
              paste0(style_sigfig(stat_nevent / stat_n, scale = 100), "%"),
              NA
            ),
          ci = ifelse(
            !is.na(conf.low) & !is.na(conf.high),
            paste0("[", style_sigfig(conf.low, digits = 2), " - ", style_sigfig(conf.high, digits = 2), "]"),
            NA
          ))
    ) |>
    modify_header(stat_nevent = "**Flag Rate**", estimate = "**OR, [95% CI]**") |>
    modify_column_merge(
      pattern = "{estimate}, {ci}",
      rows = estimate != 1
    ) |> 
    modify_column_hide(c("stat_n", "ci", "p.value"))), 
      tab_spanner = c("Flagged by Current", "Predicted by Pipeline", "Reviewed by Judges"))

  gt_combo_out <- 
    gt_combo |>
    modify_header(label = "") |> 
      as_gt() |>
      tab_header(title = "Multivariate Analysis: Demographic Breakdown of Flagged Results", subtitle = "Internal Prospective Validation Set")

gt::gtsave(gt_combo_out, filename = "../../Figures/Supervised/final_figures/table_multivariate_demographics_combined.docx")
gt::gtsave(gt_combo_out, filename = "../../Figures/Supervised/final_figures/table_multivariate_demographics_combined.html")

demo_input

}

getYaleReviewedFalsePositives <- function(preds = makeMixtureRatioPredictions(input = pipeline_preds_yale, fluid_list = fluid_names)){
  
  preds_with_review <- preds |> mutate(reviewed = calcium < calcium_prior * 0.95 & co2_totl < co2_totl_prior * 0.95)
  
  input <- 
    preds_with_review |> 
      filter(reviewed & final_prediction == 0 & .pred_class == 1 & sim_pred == 1) |> 
      mutate(epic_id = patient_id_hash, drawn_dt_tm = drawn_dt_time_relative, anion_gap = sodium - chloride - co2_totl)
    
  output <- 
    reformatReviews(input) |> left_join(input |> select(epic_id, matches("hash|pred"), -final_prediction))
  
  output_slim <- output[,-c(12:20)] |> select(matches("hash"), drawn_dt_tm, time_point, any_of(lab_strings), matches("pred"), -sim_pred, -tech_pred) 
  output_above_ten <- output_slim |> filter(mix_ratio_pred > 0.10)
  
  write_csv(output_above_ten, "../../Results/Predictions/Supervised/yale_false_positive_reviewed.csv")
  
  output_above_ten
  
}
  
comparePredsToFlags <- function(input = tar_read(test)[[1]], fluid_names = c("NS", "LR", "D5NS", "D5LR", "D5W", "hyperNS", "D5halfNSwK", "D5halfNS")){
  
  preds <- map(fluid_names, ~aggregateFinalFluidPrediction(input, fluid_name = .x))
  
  output <- bind_cols(input, bind_cols(preds)) %>% bind_cols(binarizePredictions(.))
  
  flags <- getAbnormalFlags(output)
  
  output_with_flags <- output %>% left_join(getAbnormalFlags(output))
  
  output_with_flags %>%
    filter(any_realtime != "Real" | any_retrospective != "Real") %>%
    select(contam_comment, any_realtime, any_retrospective) %>% 
    ggplot(aes(y = , fill = contam_comment)) +
      geom_bar() + 
      scale_fill_manual("") +
      scale_x_continuous()

}

comparePredsToCurrentWorkflow <- function(input = pipeline_preds){
  
  flagged_with_current <-
    input |> 
      mutate(flagged_by_current = factor(as.numeric(unlikely_comment | contam_comment)))
  
  multi_conf_mat <-  
    list(
      metrics = metrics_binary(flagged_with_current, truth = final_prediction, estimate = flagged_by_current, event_level = "second") |> transmute(Metric = .metric, Value = round(.estimate, digits = 3)),
      conf_mat = conf_mat(flagged_with_current, truth = final_prediction, estimate = flagged_by_current))
  
  list(data = flagged_with_current, preds = flagged_with_current$flagged_by_current,  metrics = multi_conf_mat$metrics, conf_mat = multi_conf_mat$conf_mat)
  
  
}


comparePredsToYaleRules <- function(input = pipeline_preds){
  
  flagged_with_yale <-
    input |> 
      bind_cols(getYaleDeltaCheckRules(input)) |> 
      mutate(yale_any = factor(as.numeric(yale_NS | yale_LR | yale_D5W)))
  
  multi_conf_mat <-  
    list(
      metrics = metrics_binary(flagged_with_yale, truth = final_prediction, estimate = yale_any, event_level = "second") |> transmute(Metric = .metric, Value = round(.estimate, digits = 3)),
      conf_mat = conf_mat(flagged_with_yale, truth = final_prediction, estimate = yale_any))
  
  list(data = flagged_with_yale, preds = flagged_with_yale$yale_any, metrics = multi_conf_mat$metrics, conf_mat = multi_conf_mat$conf_mat)
  
    
}

getReferenceChangeValueFlags <- function(input = pipeline_preds){
  
  output <- 
    input |> 
      mutate(across(any_of(lab_strings_bmp_no_gap), ~ (. > (1 + RCV_increase[[cur_column()]]) * get(paste0(cur_column(), "_prior"))) | (. < (1 - RCV_decrease[[cur_column()]])) * get(paste0(cur_column(), "_prior")), .names = "{col}_exceeded_RCV"))
  
  total_count_over_RCV <- rowSums(output |> select(matches("RCV")))
  
  output$total_exceeding_RCV <- total_count_over_RCV
  
  output
  
}

getAbnormalFlags <- function(input, flags = qs::qread("../Preprocessing/_targets/objects/raw_bmp_data")){
  
  flag_input <- 
    flags %>% 
      filter(task_assay %in% (critical_ranges %>% names()) & container_id %in% input$specimen_id & !is.na(result_value_numeric)) %>% 
      transmute(specimen_id = container_id, task_assay = task_assay, result_value_numeric, abnormal = (result_value_numeric > normal_high | result_value_numeric < normal_low))
                
  critical = map2(flag_input$task_assay, flag_input$result_value_numeric, ~(.y < critical_ranges[[.x]][["min"]] | .y > critical_ranges[[.x]][["max"]])) %>% unlist()
  
  flag_input$critical <- critical
  
  output <- flag_input %>% group_by(specimen_id) %>% summarise(abnormal = any(abnormal), critical = any(critical))
  
  output
  
}

assessPredictions <- function(fluids_to_include = c("NS", "SW", "D5W", "D5NS", "D5LR"), models_string = "XGB|NNet", out_tag = paste0("Predictions/Supervised/xgb_nnet_", paste0(fluids_to_include, collapse = "_"))){
  
##### Assess Internal Preds
  prospective_preds <- 
    makePipelinePredictions(
      sim_preds_raw = aggregateTargetObjects(fluid_list = fluids_to_include, step_name = "wf_prospective", feature_set = "results_with_priors", models_string = models_string), 
      tech_preds_raw = aggregateTargetObjects(fluid_list = c("Tech"), step_name = "wf_prospective", feature_set = "results_with_priors", models_string = models_string)
    )
  write_csv(prospective_preds, paste0("../../Results/", out_tag, "_final_pipeline_prospective_preds.csv"))
  
  pipeline_preds <- 
    makePipelinePredictions(
      sim_preds_raw = aggregateTargetObjects(fluid_list = fluids_to_include, step_name = "wf_reviewed", feature_set = "results_with_priors", models_string = models_string), 
      tech_preds_raw = aggregateTargetObjects(fluid_list = c("Tech"), step_name = "wf_calibrated", feature_set = "results_with_priors", models_string = models_string)
    )
  write_csv(prospective_preds, paste0("../../Results/", out_tag, "_final_pipeline_internal_preds.csv"))
  
  pipeline_preds |> 
    filter(tech_pred == 1 & sim_pred == 0) |> 
    write_csv("../../Results/Predictions/Supervised/internal_tech_pred_pos_but_sim_pred_negative.csv")
  
  pipeline_preds |> 
    filter(.pred_class == 1 & final_prediction == 0)
  
  pipeline_preds_reported_results_only <- 
    pipeline_preds |> 
      filter(!reported_as_comment|contam_comment|code_comment|unlikely_comment)
  
  pipeline_pre_post_scatter <- 
    makePrePostScatterplots(
      pipeline_preds, target = "final_prediction", path = "../../Figures/Supervised/full_pipeline_pre_post_scatter", analytes_to_include = c("sodium", "chloride", "co2_totl", "calcium", "glucose")
    )
  
  pipeline_conf_mat <- getMetricsAndConfusionMatrix(pipeline_preds)
  write_csv(pipeline_conf_mat$metrics, paste0("../../Results/", out_tag, "_full_pipeline_all_model_metrics_compared_to_expert_review.csv"))

  pipeline_reported_only_conf_mat <-
    list(
      metrics = metrics_binary(pipeline_preds_reported_results_only, truth = final_prediction, estimate = .pred_class, event_level = "second") |> transmute(Metric = .metric, Value = round(.estimate, digits = 3)),
      conf_mat = conf_mat(pipeline_preds_reported_results_only, truth = final_prediction, estimate = .pred_class))
  write_csv(pipeline_conf_mat$metrics, "../../Results/Predictions/Supervised/full_pipeline_reported_results_only_all_model_metrics_compared_to_expert_review.csv")
  
    
  pipeline_umap_label <- 
    plotPredictionsOnUMAP(
      input = pipeline_preds, umap_model = pin_read(model_board, "UMAP_ClinChem_Final") |> bundle::unbundle(), 
      target = "final_prediction", pred = ".pred_class", 
      path = "../../Figures/Supervised/full_pipeline_internal", 
      subtitle_label = "WashU Model Predicting on WashU Data")
  
  pipeline_mixture_ratios <- makeMixtureRatioPredictions(input = pipeline_preds, fluid_list = fluids_to_include)
  pipeline_false_pos <- pipeline_mixture_ratios |> filter(.pred_class == 1 & final_prediction == 0)
  reformatReviews(pipeline_false_pos) |> write_csv("../../Results/Predictions/Supervised/full_pipeline_false_positives_to_review.csv")
  
  pipeline_mixture_distributions <- 
    plotMixtureDistributions(
      input = pipeline_preds, 
      target = "final_prediction", pred = ".pred_class", 
      path = "../../Figures/Supervised/full_pipeline_internal", 
      subtitle_label = "WashU Model Predicting on WashU Data")
  
  prospective_apd <- 
    getApplicabilityScores(input = pipeline_preds, apd_model = tar_read(apd_pca))
  
  preds_exceeding_RCV <- getReferenceChangeValueFlags(pipeline_mixture_ratios)
  

###### Assess External Preds
  
  pipeline_preds_yale <- 
    makePipelinePredictions(
      sim_preds_raw = aggregateTargetObjects(fluid_list = fluids_to_include, step_name = "wf_yale_reviewed", feature_set = "results_with_priors", models_string = "XGB|NNet"), 
      tech_preds_raw = aggregateTargetObjects(fluid_list = c("Tech"), step_name = "wf_yale_calibrated", feature_set = "results_with_priors", models_string = "XGB|NNet")
    )
  write_csv(pipeline_preds_yale, "../../Results/Predictions/Supervised/final_pipeline_yale_preds.csv")
  
  pipeline_preds_yale |> 
    filter(tech_pred == 1 & sim_pred == 0) |> 
    write_csv("../../Results/Predictions/Supervised/external_tech_pred_pos_but_sim_pred_negative.csv")
  
  
  pipeline_pre_post_scatter <- 
    makePrePostScatterplots(
      pipeline_preds_yale, target = "final_prediction", path = "../../Figures/Supervised/full_pipeline_pre_post_scatter_yale", analytes_to_include = c("sodium", "chloride", "co2_totl", "calcium", "glucose")
    )
  
  pipeline_conf_mat_yale <-
    list(
      metrics = metrics_binary(pipeline_preds_yale, truth = final_prediction, estimate = .pred_class, event_level = "second") |> transmute(Metric = .metric, Value = round(.estimate, digits = 3)),
      conf_mat = conf_mat(pipeline_preds_yale, truth = final_prediction, estimate = .pred_class))
  write_csv(pipeline_conf_mat_yale$metrics, "../../Results/Predictions/Supervised/full_pipeline_yale_metrics_compared_to_expert_review.csv")
  
  pipeline_mixture_distributions_yale <- 
    plotMixtureDistributions(input = pipeline_preds_yale, 
                             target = "final_prediction", pred = ".pred_class", 
                             path = "../../Figures/Supervised/full_pipeline_yale", 
                             subtitle_label = "WashU Model Predicting on Yale Data")
  
  pipeline_yale_mixture_ratios <- makeMixtureRatioPredictions(input = pipeline_preds_yale, fluid_list = fluids_to_include)
  pipeline_yale_false_pos <- pipeline_yale_mixture_ratios |> filter(.pred_class == 1 & final_prediction == 0)
  
  pipeline_umap_label_yale <- 
    plotPredictionsOnUMAP(
      input = pipeline_yale_mixture_ratios |> mutate(anion_gap = sodium - chloride - co2_totl) |> filter(calcium < calcium_prior * 0.95 & co2_totl < co2_totl_prior * 0.95 & mix_ratio_pred > 0.10), umap_model = pin_read(model_board, "UMAP_ClinChem_Final") |> bundle::unbundle(), 
      target = "final_prediction", pred = ".pred_class", 
      path = "../../Figures/Supervised/full_pipeline_yale", 
      subtitle_label = "WashU Model Predicting on Yale Data")
  
  reformatReviews(pipeline_yale_false_pos |> mutate(epic_id = patient_id_hash, drawn_dt_tm = drawn_dt_time_relative)) |> write_csv("../../Results/Predictions/Supervised/full_pipeline_false_positives_to_review_yale.csv")
  
  #preds_yale_exceeding_RCV_table <- compareRCVbyClassLabel(input = pipeline_yale_mixture_ratios |> getReferenceChangeValueFlags(), out_path = "../../Figures/Supervised/final_figures/prospective_RCV_table_yale")
  
##### Compare WashU and Yale #####
  library(ggpubr)
  
  distributions <- compareResultDistributions(internal = pipeline_preds, external = pipeline_preds_yale)
  
  gg_false_positive_distributions <- ggpubr::ggarrange(pipeline_mixture_distributions$plot + ggtitle("", subtitle = "WashU Expert Review Set") + theme(plot.title = element_blank(), plot.subtitle = element_text(face = "italic")), pipeline_mixture_distributions_yale$plot + ggtitle("", subtitle = "Yale Expert Review Set") + theme(axis.title.y.left = element_blank(), plot.title = element_blank(), plot.subtitle = element_text(face = "italic")), nrow = 1, ncol = 2) |> ggpubr::annotate_figure(plot, top = text_grob("Comparison of Contamination Severity by Label", face = "bold", size = 14))
  ggsave("../../Figures/Supervised/final_figures/false_vs_true_positive_mixture_ratio_distribution.svg", gg_false_positive_distributions, width = 8.5, height = 4, dpi = 600)
  ggsave("../../Figures/Supervised/final_figures/false_vs_true_positive_mixture_ratio_distribution.pdf", gg_false_positive_distributions, width = 8.5, height = 4, dpi = 600)
  ggsave("../../Figures/Supervised/final_figures/false_vs_true_positive_mixture_ratio_distribution.png", gg_false_positive_distributions, width = 8.5, height = 4, dpi = 600)
  
  sim_vs_tech_comp <- compareSimToTechPreds(input = pipeline_preds)
  #preds_exceeding_RCV_table <- compareRCVbyClassLabel()
  
  fluid_list_comparison
    

}

