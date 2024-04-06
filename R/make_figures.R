makeCrossValMetricFigures <- function(metrics_to_include = c("sensitivity", "specificity", "roc_auc"), fluids_to_include = c(fluid_names, "Tech"), out_path = "../../Figures/Supervised/cv_metrics/"){
  
  cv_tune_results <- aggregateTargetObjects(fluid_list = fluids_to_include, step_name = "wf_cv_metrics", models_string = "XGB|NNet|CART|LogReg")
  
  tuned_metrics <- 
    cv_tune_results |> 
      bind_rows() |> 
      select(.metric:label) |>
      transmute(config = .config, metric = .metric, value = mean, std_err = std_err, fluid = str_match(label, "_BMP_(.*?)_Base")[,2], model = str_match(label, "Base_(.*?)_")[,2]) |> 
      mutate(fluid = ifelse(is.na(fluid), "Tech", fluid))
  
  best_metrics <- 
    tuned_metrics |> 
      group_by(metric, model, fluid) |>
      summarise(value = round(max(value), digits = 3), sd = max(std_err)) |> 
      ungroup() |> 
      filter(metric %in% metrics_to_include)
  
  ggplot(best_metrics, aes(y = model, fill = value)) +
    geom_tile(aes(x = 0)) +
    geom_text(aes(x = 0, label = value)) +
    scico::scale_fill_scico(palette = "berlin", limits = c(0, 1)) + 
    xlab("Fluid") + ylab("Model") + ggtitle("Model Performance on Simulated Cross-Validation Sets") +
    facet_grid(metric~fluid, scales = "free") + 
    coord_cartesian(clip = "off") +
    theme(axis.text.x = element_blank(), strip.clip = "off")
  ggsave(paste0(out_path, "cv_metrics_grid.pdf"), width = 12, height = 6, dpi = 600)
  ggsave(paste0(out_path, "cv_metrics_grid.png"), width = 12, height = 6, dpi = 600)
  ggsave(paste0(out_path, "cv_metrics_grid.svg"), width = 12, height = 6)
  
  tuned_metrics
  
}

makeExpertLabelMetricFigures <- function(metrics_to_include = c("ppv", "mcc", "pr_auc"), fluids_to_include = c(fluid_names, "Tech"), out_path = "../../Figures/Supervised/cv_metrics/"){
  
  test_results <- 
    aggregateTargetObjects(fluid_list = fluids_to_include, step_name = "wf_test_metrics", models_string = "XGB|NNet|CART|LogReg") |> 
    map("metrics")
  
  test_results <- map2(test_results, names(test_results), ~mutate(.x, label = .y))
  
  test_metrics <- 
    test_results |> 
      bind_rows() |> 
      transmute(metric = .metric, value = round(.estimate, digits = 3), fluid = str_match(label, "_BMP_(.*?)_Base")[,2], model = str_match(label, "Base_(.*?)_")[,2]) |> 
      mutate(fluid = ifelse(is.na(fluid), "Tech", fluid)) |> 
      filter(metric %in% c(metrics_to_include))
  
  ggplot(test_metrics, aes(y = model, fill = value)) +
    geom_tile(aes(x = 0)) +
    geom_text(aes(x = 0, label = value, color = (value > 0.05 & value < 0.95)), show.legend = F) +
    scico::scale_fill_scico(palette = "berlin", limits = c(0, 1)) + 
    scale_color_manual(values = c("black", "white")) +
    xlab("Fluid") + ylab("Model") + ggtitle("Model Performance on Simulated Cross-Validation Sets") +
    facet_grid(metric~fluid, scales = "free") + 
    coord_cartesian(clip = "off") +
    theme(axis.text.x = element_blank(), strip.clip = "off")
  ggsave(paste0(out_path, "expert_label_metrics_grid.pdf"), width = 12, height = 6, dpi = 600)
  ggsave(paste0(out_path, "expert_label_metrics_grid.pdf"), width = 12, height = 6, dpi = 600)
  ggsave(paste0(out_path, "expert_label_metrics_grid.pdf"), width = 12, height = 6)
  
  tuned_metrics
  
}

makeSimulatedSensitivityComparisonCurvesHelper <- function(fluid_name = "SW"){
  
  input = aggregateTargetObjects(fluid_list = c(fluid_name, "Tech"), step_name = "wf_sim_sensitivity", models_string = "XGB")
  
  input[[1]]$model <- "Fluid-Specific"
  input[[2]]$model <- "General"
  
  gg_input <- 
    input |> 
      bind_rows() |> 
      filter(label == fluid_name) |> select(.pred_class, mix_ratio, model) |> 
      group_by(mix_ratio, model) |> 
      summarise(n = sum(as.numeric(as.character(.pred_class)))) |> 
      mutate(prop = n / 1000)
  
  library(geomtextpath)
  
  binomial_smooth <- function(...) {
    geom_textsmooth(method = "glm", method.args = list(family = "binomial"), ...)
  }
  
  ggplot(gg_input, aes(mix_ratio, prop, color = fct_rev(model))) + 
    binomial_smooth(aes(label = model), size = 2, linewidth = 1.25, fontface = "bold", span = 0.25, hjust = 0.15, se = F) + 
    scale_x_continuous(name = "Mixture Ratio", breaks = seq(0, 0.5, by = 0.1), limits = c(0, 0.5)) +
    scale_y_continuous(name = "Sensitivity", breaks = c(0, 0.5, 1), limits = c(0, 1)) + 
    scale_color_grey() + 
    ggtitle(fluid_name) + 
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold.italic"))
  
}

makeSimulatedSensitivityComparisonCurves <- function(){
  
  curves <- map(c("NS", "D5NS", "D5LR", "SW", "D5W", "D5halfNSwK"), ~makeSimulatedSensitivityComparisonCurvesHelper(.x))
  
  half_sensitive <- map(curves, ~.x |> pluck("data") |> filter(prop > 0.5) |> ungroup() |> slice_head(n = 1, by = model)) |> bind_rows()
  
  tmp <- map(curves[c(2,3,5,6)], ~.x + theme(axis.line.y.left = element_blank(), axis.text.y.left = element_blank(), axis.title.y.left = element_blank()))
  tmp_no_labels <- map(tmp[c(1,2)], ~.x + theme(axis.line.x.bottom = element_blank(), axis.text.x.bottom = element_blank(), axis.title.x.bottom = element_blank()))
  
  gg_curves_combined <- 
    ggpubr::ggarrange(
      curves[[1]] + theme(axis.line.x.bottom = element_blank(), axis.text.x.bottom = element_blank(), axis.title.x.bottom = element_blank()), 
      tmp_no_labels[[1]],
      tmp_no_labels[[2]], 
      curves[[4]], 
      tmp[[3]],
      tmp[[4]], 
      nrow = 2, ncol = 3)
  ggsave("../../Figures/Supervised/final_figures/fig2b_sensitivity_curves.png", gg_curves_combined, width = 8, height = 4, dpi = 600)
  ggsave("../../Figures/Supervised/final_figures/fig2b_sensitivity_curves.pdf", gg_curves_combined, width = 8, height = 4)
  ggsave("../../Figures/Supervised/final_figures/fig2b_sensitivity_curves.svg", gg_curves_combined, width = 8, height = 4)
    
  gg_curves_combined
  
}

makeCurves <- function(preds = tar_read(wf_reviewed_BJH_Tech_BMP_Base_NNet_results_with_priors), target = "final_prediction", path = "../../Figures/Supervised/Tech_model_NNet_curves"){
  
  library(gridExtra)
  metric_list <- getMetrics(preds[[".pred_1"]], preds[[target]], threshold = 0.5)
  conf_table <- tableGrob(metric_list[[2]] %>% filter(!.metric %in% c("binary_classification_cost", "mn_log_loss")) %>% transmute(Metric = factor(.metric, levels = names(metric_labels), labels = metric_labels), Value = round(.estimate, digits = 2)), theme = ttheme_minimal(base_size = 8, padding = unit(c(4, 4), units = "pt")), rows = NULL)
  
  pr = preds %>% pr_curve(.pred_1, truth = !!target, event_level = "second") %>% mutate(curve = "pr") %>% bind_rows(tibble_row(.threshold = 0, recall = 1, precision = 0))
  roc = preds %>% roc_curve(.pred_1, truth = !!target, event_level = "second") %>% mutate(curve = "roc")
  
  gg_pr = 
    ggplot(pr, aes(recall, precision)) + 
    geom_path() + 
    xlim(0,1) + ylim(0,1) + 
    xlab("Recall") + ylab("Precision")
  
  gg_roc = 
    ggplot(roc, aes(1-specificity, sensitivity)) + 
    geom_line() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
    xlim(0,1) + ylim(0,1) + 
    xlab("1 - Specificity") + ylab("Sensitivity") +
    annotation_custom(conf_table, xmin = 0.7, ymin = 0, xmax = 1, ymax = 0.4)
  
  library(ggpubr)  
  gg_curves <- ggarrange(gg_roc, gg_pr, nrow = 1, ncol = 2)
  ggsave(paste0(path, ".pdf"), gg_curves, width = 8.5, height = 4, dpi = 600)  
  ggsave(paste0(path, ".png"), gg_curves, width = 8.5, height = 4, dpi = 600) 
  ggsave(paste0(path, ".svg"), gg_curves, width = 8.5, height = 4, dpi = 600) 
  
  list(gg_roc, gg_pr)
  
}

makeSimSensitivityPlots <- function(contam_sim_all = read_csv("../../Data/contam_sim_all.csv"), wf_fit = tar_read(wf_fit_BJH_Tech_BMP_Base_NNet_results_with_priors), path = "../../Figures/Supervised/Tech_model_NNet_sim_sensitivity"){
  
  contam_sim_all <- contam_sim_all |> slice_sample(prop = 0.2, by = c("label", "mix_ratio"))
  
  sim_preds <- predict(wf_fit |> bundle::unbundle(), contam_sim_all)
  
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
  ggsave(paste0(path, ".pdf"), gg_sim_sensitivity, width = 8.5, height = 4, dpi = 600)  
  ggsave(paste0(path, ".png"), gg_sim_sensitivity, width = 8.5, height = 4, dpi = 600) 
  ggsave(paste0(path, ".svg"), gg_sim_sensitivity, width = 8.5, height = 4, dpi = 600) 
  
  list(input = gg_sim_input, plot = gg_sim_sensitivity)
  
}

makePrePostTiles <- function(preds = read_csv("../../Results/Predictions/Supervised/final_pipeline_internal_preds.csv"), target = "final_prediction", cap = 5, path = "../../Figures/Supervised/full_pipeline", subtitle_label = "WashU Pipeline on Prospective WashU Data", analytes_to_include = c("chloride", "calcium", "sodium", "glucose")){


gg_input <- 
  preds |> 
  select(target = !!target, pred = .pred_class, matches("delta")) |> 
  mutate(tag = paste(pred, target, sep = "_"),
         across(matches(paste0(lab_strings, collapse = "|")), ~scale(.)), 
         across(matches(paste0(lab_strings, collapse = "|")), ~case_when(. > cap ~ cap, . < -cap ~ -cap, T ~ .))) |>
  group_by(tag) |>
  mutate(row_id = row_number()) |> 
  pivot_longer(-matches("row_id|target|pred|tag"), names_to = c("var", "time"), names_sep = "_delta_", values_to = "result") |>
  mutate(var = factor(var, levels = lab_strings)) |> 
  ungroup()


gg_pre_post_tiles <- 
  ggplot(gg_input |> filter(var %in% analytes_to_include), aes(x = var, y = row_id, fill = result)) + 
    geom_tile(color = NA) + 
    scale_x_discrete(labels = analyte_labels) + 
    scico::scale_fill_scico(palette = "vik", begin = 0.1, end = 0.9) +
    guides(fill = guide_colourbar(title = "Normalized Delta", direction = "horizontal", title.position = "top", title.theme = element_text(size = 8, face = "italic", hjust = 0.5), label.theme = element_text(size = 6, face = "bold", hjust = 0.5))) + 
    xlab("Result Deltas") + ylab("Specimen") + 
    ggtitle("Comparison of Deltas by Model and Expert Label Group", subtitle = subtitle_label) +
    facet_grid(tag ~ fct_rev(time), scales = "free", labeller = as_labeller(c(analyte_labels, "0_0" = "True Neg", "0_1" = "False Neg", "1_0" = "False Pos", "1_1" = "True Pos", "prior" = "Priors", "post" = "Re-draws"))) + 
    theme(legend.position = c(0.91, 1.07), legend.background = element_blank(), strip.text.x.top = element_text(size = 14), strip.clip = "off", plot.title = element_text(face = "bold"), legend.key.size = unit(.2, .2, units = "in"), axis.text.x.bottom = element_text(face = "bold"))
ggsave(paste0(path, "_pre_post_tiles.svg"), gg_pre_post_tiles, width = 8.5, height = 8)  
ggsave(paste0(path, "_pre_post_tiles.pdf"), gg_pre_post_tiles, width = 8.5, height = 8)  
ggsave(paste0(path, "_pre_post_tiles.png"), gg_pre_post_tiles, width = 8.5, height = 8)

list(input = gg_input, plot = gg_pre_post_tiles)

}

makePrePostScatterplots <- function(preds = tar_read(wf_reviewed_BJH_Tech_BMP_Base_NNet_results_with_priors), target = "final_prediction", path = "../../Figures/Supervised/Tech_model_NNet_pre_post_scatter_results_with_priors", subtitle_label = "WashU Pipeline on Prospective WashU Data", analytes_to_include = c("chloride", "calcium", "sodium", "glucose")){
  
  gg_input <- 
    preds |> 
      select(target = !!target, pred = .pred_class, matches("delta")) |> 
      mutate(
        across(matches("delta"), ~scale(.)),
        tag = paste(pred, target, sep = "_")) |>
      group_by(tag) |>
      mutate(row_id = row_number()) |> 
      pivot_longer(-matches("row_id|target|pred|tag"), names_to = c("var", "time"), names_sep = "_delta_", values_to = "result") |>
      mutate(var = factor(var, levels = lab_strings)) |> 
      ungroup() |>
      pivot_wider(id_cols = -matches("time"), names_from = "time", values_from = "result") |> 
      drop_na()
  
  corrs <- gg_input |> group_by(var, tag) |> summarise(corr = cor(prior, post)^2)
  gg_input <- gg_input |> filter(var %in% analytes_to_include) |> left_join(corrs)
  
  gg_pre_post_scatters <- 
    ggplot(gg_input, aes(prior, post)) + 
      geom_vline(xintercept = 0, alpha = 0.2, linetype = "dashed") + 
      geom_hline(yintercept = 0, alpha = 0.2, linetype = "dashed") + 
      ggpmisc::stat_poly_line(data = gg_input |> filter(tag != "0_0"), aes(color = corr), se = F, na.rm = T, alpha = 1, linewidth = 1.25) + 
      geom_point(data = gg_input |> filter(tag == "0_0") |> slice_sample(n = 1000), aes(prior, post), alpha = 0.25) +
      geom_point(data = gg_input |> filter(tag != "0_0"), aes(prior, post), alpha = 0.25) +
      ggpmisc::stat_poly_eq(aes(color = after_stat(r.squared), label = after_stat(rr.label)), label.x = 0.05, label.y = 0.05, rr.digits = 1, na.rm = T) +
      #    scale_x_discrete(limits = c("prior", "post"), labels = c("Prior", "Post")) + 
      scico::scale_color_scico(palette = "lajolla", limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
      guides(color = guide_colorbar(title = "Correlation", title.position = "top", title.hjust = 0.5, label.theme = element_text(size = 6), barwidth = unit(1, "in"), barheight = unit(0.15, "in"), direction = "horizontal", title.theme = element_text(size = 8, face = "bold.italic", margin = margin(0,0,0,0)))) +
      coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) + 
      xlab("Delta From Prior") + ylab("Delta To Post") + 
      ggtitle("Comparison of Deltas by Label Group", subtitle = subtitle_label) +
      ggh4x::facet_grid2(tag ~ var, strip = ggh4x::strip_vanilla(clip = "off"), labeller = as_labeller(c(analyte_labels, "0_0" = "True Neg", "0_1" = "False Neg", "1_0" = "False Pos", "1_1" = "True Pos"))) + 
      theme(legend.position = c(0.92, 1.09), legend.background = element_blank(), plot.title = element_text(face = "bold"), axis.text = element_text(size = 8), panel.background = element_rect(fill = "gray95", color = NA))
  #  ggsave(paste0(path, ".svg"), gg_pre_post_scatters, width = 8.5, height = 7)  
  ggsave(paste0(path, ".pdf"), gg_pre_post_scatters, width = 8.5, height = 6, dpi = 600)  
  ggsave(paste0(path, ".png"), gg_pre_post_scatters, width = 8.5, height = 6, dpi = 600) 
  
  list(input = gg_input, plot = gg_pre_post_scatters)
  
}

makeOldPrePostPlots <- function(input = manual_review_final, group_col = "expert_label", analytes_to_show = c("calcium"), rows_to_show = 100, out_path = "../../Figures/Supervised/pre_post_calcium_expert_review_yale.pdf"){
  
  gg_input <- 
    input %>% 
      drop_na(!!group_col) %>%
      mutate(index = factor(row_number()), target = factor(eval(parse(text = group_col)), labels = c("Real", "Contaminated"))) %>%
      mutate(across(matches("delta") & matches(paste0(lab_strings, collapse = "|")), ~scale(.) %>% as.vector())) %>% 
      select(index, matches(paste0(analytes_to_show, collapse = "|")) & matches("delta"), target) %>%
      drop_na() %>%
      pivot_longer(matches("delta"), names_to = c("analyte", "time_point"), names_sep = "_delta_", values_to = "delta") %>%
      left_join(., filter(., time_point == "prior") %>% transmute(index, prior_delta = delta), by = "index") %>% 
      group_by(index, analyte) %>% 
      arrange(prior_delta) %>%
      group_by(target, time_point) %>%
      mutate(delta = ifelse(delta > 5, 5, delta), delta = ifelse(delta < -5,  -5, delta), row = fct_reorder(index, prior_delta, .desc = T)) %>% ungroup()
  
  rows_to_keep <- gg_input %>% distinct(index, target) %>% slice_sample(n = rows_to_show, by = target) %>% pluck("index")
  gg_input_sample <- gg_input %>% filter(index %in% rows_to_keep) %>% mutate(row = row_number(), .by = c(target, time_point))
  
  gg <- 
    ggplot() + 
      geom_tile(data = gg_input_sample, aes(x = fct_rev(time_point), y = max(row) + 1 - row, fill = delta)) + 
      scico::scale_fill_scico(palette = "vik", begin = 0.05, end = 0.95, direction = 1) +
      annotate("segment", x = 1.5, xend = 1.5, y = 0, yend = rows_to_show + 1, linewidth = 1) +
      facet_wrap(~target, scales = "free") +
      guides(fill = guide_colorbar(title = "Normalized Delta", direction = "horizontal", label.position = "bottom", title.position = "left")) + 
      scale_y_discrete(name = "Specimen", expand = c(0, 0.05)) + 
      scale_x_discrete(name = "Deltas", breaks = c("prior", "post"), labels = c("Prior", "Post")) +
      ggtitle(analyte_labels[[analytes_to_show]]) +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom", axis.text.y.left = element_blank(), axis.text.x.bottom = element_text(face = "bold", color = "black"), axis.title.x.bottom = element_blank(), strip.text = element_text(face = "bold"))
  ggsave(out_path, gg, width = 4, height = 4)
  
  gg
  
}

plotRedrawFigures <- function(){
  
  gg_redraw_ecdf <- 
    ggplot(result_preds_meta %>% filter(!is.na(.pred_class)) %>% distinct(specimen_id, time_to_next, flag_bin), aes(time_to_next, color = fct_rev(flag_bin))) + 
    geom_vline(xintercept = 24, linetype = "dashed", alpha = 0.25) +
    geom_vline(xintercept = 12, linetype = "dashed", alpha = 0.25) +
    stat_ecdf(linewidth = 1.5) + 
    scale_x_continuous(name = "Hours Until Next Result", limits = c(0, 26), breaks = c(0, 6, 12, 18, 24)) + 
    scale_y_continuous(name = "Cumulative Proportion", breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) + 
    scico::scale_color_scico_d(palette = "davos", begin = 0, end = 0.8) +
    theme()
  ggsave("../../Figures/Supervised/contaminated_ecdf_redraw_time.pdf", gg_redraw_ecdf, width = 8.5, height = 4)
  
}

plotFlagFigures <- function(result_preds_meta = read_csv("../../Results/Predictions/Supervised/full_pipeline_prospective_with_result_metadata_flags.csv")){
  
  theme_set(theme_ns)
  
  flag_counts <- 
    result_preds_meta %>% 
      group_by(task_assay, .pred_class) %>% 
      select(task_assay, flag, flagged_by_technologist, .pred_class) %>% 
      count(task_assay, flag, flagged_by_technologist, .pred_class) %>% 
    filter(!task_assay %in% c("anion_gap", "bun", "creatinine")) %>% 
    mutate(task_assay = factor(task_assay, levels = c("sodium", "chloride", "bun", "calcium", "potassium_plas", "co2_totl", "creatinine", "glucose")), 
           flag = factor(flag, levels = c("Normal", "Abnormal", "Critical")),
           flagged_by_technologist = factor(flagged_by_technologist, levels = c(0, 1), labels= c("Unflagged", "Flagged")))
  
  gg_flag_counts <- 
    flag_counts %>% filter(.pred_class == 1) %>%
      ggplot(aes(x = n, y = flagged_by_technologist, fill = flag)) +
      geom_col(position = position_stack(reverse = T)) + 
      scico::scale_fill_scico_d(palette = "bilbao", begin = 0.2, end = 0.8, direction = -1) +
      facet_wrap(~task_assay, nrow = 2, ncol = 3, labeller = as_labeller(analyte_labels)) +
      scale_x_continuous(name = "Count") + 
      ylab("Technologist Flag") + 
      theme(legend.position = c(0.18, 0.9), legend.background = element_blank(), legend.key.size = unit(0.3, "cm"), legend.title = element_blank())
  ggsave("../../Figures/Supervised/final_figures/prospective_flag_counts_pred_vs_tech.pdf", gg_flag_counts, width = 8.5, height = 5)
  
  

  
  flag_props <- 
    result_preds_meta %>% 
    mutate(flagged_by_technologist = unlikely_comment|contam_comment, 
           flag_tag = paste0(.pred_class, flagged_by_technologist),
           flag_tag = factor(flag_tag, labels = c("Neither", "Tech Only", "ML Only", "Both"))) |> 
    count(task_assay, flag_tag, flag) |> 
    group_by(task_assay, flag_tag) |> 
    mutate(prop = n / sum(n), 
           task_assay = factor(task_assay, levels = c("sodium", "chloride", "bun", "calcium", "potassium_plas", "co2_totl", "creatinine", "glucose"))) |> 
    filter(!task_assay %in% c("anion_gap", "bun", "creatinine")) %>% 
    mutate(flag = factor(flag, levels = c("Normal", "Abnormal", "Critical")))
  
  gg_flag_props <- 
    flag_props %>%
    ggplot(aes(x = prop, y = fct_rev(flag_tag), fill = flag)) +
    geom_col() + 
    scico::scale_fill_scico_d(palette = "bilbao", begin = 0.05, end = 0.95, direction = -1) +
    facet_wrap(~task_assay, nrow = 2, ncol = 3, labeller = as_labeller(analyte_labels)) +
    scale_x_continuous(name = "Count") + 
    ylab("Technologist Flag") + 
    theme(legend.position = c(0.18, 0.9), legend.background = element_blank(), legend.key.size = unit(0.3, "cm"), legend.title = element_blank())
  ggsave("../../Figures/Supervised/final_figures/prospective_flag_props_pred_vs_tech.pdf", gg_flag_props, width = 8.5, height = 5)
  
  
  
  flag_input <- 
    result_preds_meta %>% 
    select(task_assay, flagged_by_technologist, .pred_class, flag, next_flag) %>%
    arrange(flag, next_flag)
  
  gg_current_next_retro <- map(c("calcium", "glucose"), ~plotCurrentPostFlags(flag_input, group_col = ".pred_class", analytes_to_show = .x, rows_to_show = 50, out_path = paste0("../../Figures/Supervised/impact_plots/current_next_", .x, "_pred_flags.pdf")))
  
  gg_cur_next <- ggpubr::ggarrange(plotlist = gg_current_next_retro, nrow = 1)
  ggsave("../../Figures/Supervised/impact_plots/current_next_all_pred_flags.pdf", gg_cur_next, width = 5, height = 4)
  
  gg_current_next_tech <- map(c("calcium", "glucose"), ~plotCurrentPostFlags(flag_input, group_col = "flagged_by_technologist", analytes_to_show = .x, rows_to_show = 50, out_path = paste0("../../Figures/Supervised/impact_plots/current_next_", .x, "_tech_flags.pdf")))
  
  ggpubr::ggarrange(plotlist = gg_current_next_tech, nrow = 1)
  ggsave("../../Figures/Supervised/impact_plots/current_next_all_tech_flags.pdf", width = 5, height = 4)
  
  
  fig5 <- ggarrange(gg_flag_counts, ggarrange(plotlist = list(gg_time_ecdf, gg_cur_next), nrow = 1, ncol = 2, labels = c("B", "C")), nrow = 2, labels = c("A", ""))
  ggsave("../../Figures/Supervised/impact_plots/counts_time_ecdfs_cur_next.pdf", fig5, width = 12, height = 6)
  
}

plotReferenceChangeValueFlags <- function(input = read_csv("../../Results/Predictions/Supervised/final_pipeline_prospective_preds.csv") |> mutate(flagged_by_technologist = unlikely_comment|contam_comment)){
    
    prior_prop <- 
      input %>% 
        reframe(across(any_of(lab_strings) & !matches("anion_gap"), ~ (. - get(paste0(cur_column(), "_prior")))/get(paste0(cur_column(), "_prior")))) |> 
        bind_cols(input |> select(matches("_id|_dt|pred|flag")) |> mutate(row_id = row_number()))
    
    post_prop <- 
      input %>% 
        reframe(across(any_of(lab_strings) & !matches("anion_gap"), ~ (get(paste0(cur_column(), "_post")) - .)/.)) |>
        bind_cols(input |> select(matches("_id|_dt|pred|flag")) |> mutate(row_id = row_number()))
    
    RCV_prior_flags <- 
      prior_prop %>% 
        mutate(across(any_of(lab_strings), ~ (. > RCV_increase[[cur_column()]]) | (. < -RCV_decrease[[cur_column()]]))) %>% 
        pivot_longer(cols = any_of(lab_strings), names_to = "task_assay", values_to = "prior_RCV_flag")

    RCV_post_flags <- 
      post_prop %>% 
        mutate(across(any_of(lab_strings), ~ (. > RCV_increase[[cur_column()]]) | (. < -RCV_decrease[[cur_column()]]))) %>% 
        pivot_longer(cols = any_of(lab_strings), names_to = "task_assay", values_to = "post_RCV_flag")
    
    
    gg_input <- 
      left_join(RCV_prior_flags, RCV_post_flags) |> 
        mutate(flag_tag = paste0(.pred_class, flagged_by_technologist),
               flag_tag = factor(flag_tag, labels = c("Neither", "Tech Only", "ML Only", "Both")),
               RCV_tag = paste0(prior_RCV_flag, post_RCV_flag),
               RCV_tag = factor(RCV_tag, labels = c("Neither", "Prior Only", "Post Only", "Both"))) |> 
        count(task_assay, flag_tag, RCV_tag) |> 
        group_by(task_assay, flag_tag) |> 
        mutate(prop = n / sum(n), 
               task_assay = factor(task_assay, levels = c("sodium", "chloride", "bun", "calcium", "potassium_plas", "co2_totl", "creatinine", "glucose"))) 
      
    ggplot(gg_input |> filter(!task_assay %in% c("creatinine", "bun")), aes(x = prop, y = fct_rev(flag_tag), fill = RCV_tag)) + 
      geom_col() + 
      scale_x_continuous(name = "Proportion", breaks = c(0, 0.5, 1)) + 
      scico::scale_fill_scico_d(palette = "bilbao", direction = -1, begin = 0.05, end = 0.95) + 
      facet_wrap(~task_assay, labeller = as_labeller(analyte_labels))
    
}

plotCurrentPostFlags <- function(input = read_csv("../../Results/Predictions/Supervised/final_pipeline_internal_preds.csv"), group_col = "final_prediction", analytes_to_show = c("sodium"), rows_to_show = 100, out_path = "../../Figures/Supervised/impact_plots/current_next_sodium_pred_flags.pdf"){
  
  gg_input <- 
    input %>% 
      filter(task_assay == analytes_to_show) %>%
      select(task_assay, !!group_col, flag, next_flag) %>% 
      drop_na(!!group_col) %>%
      mutate(target = factor(eval(parse(text = group_col)), labels = c("Real", "Contaminated"))) %>%
      group_by(target, flag) %>%
      slice_sample(n = rows_to_show) %>%
      arrange(flag, next_flag) %>%
      group_by(target) %>%
      mutate(index = row_number()) %>%
      select(index, task_assay, target, flag, next_flag) %>% 
      drop_na() %>%
      ungroup()
    
  rows_to_keep <- gg_input %>% distinct(index, target) %>% slice_sample(n = rows_to_show, by = c(target)) %>% pluck("index")
  gg_input_sample <- gg_input %>% filter(index %in% rows_to_keep) %>% mutate(row = row_number(), .by = target)
  
  gg <- 
    ggplot() + 
    geom_rect(data = gg_input, aes(ymin = index, ymax = index + 1, xmin = 0, xmax = 1, fill = flag, color = flag)) + 
    geom_rect(data = gg_input, aes(ymin = index, ymax = index + 1, xmin = 1, xmax = 2, fill = next_flag, color = next_flag)) +
    facet_wrap(~target, scales = "free") +
    scico::scale_fill_scico_d(palette = "bilbao", begin = 0.2, end = 0.8, direction = -1) +
    scico::scale_color_scico_d(palette = "bilbao", begin = 0.2, end = 0.8, direction = -1) +
    annotate("segment", x = 1, xend = 1, y = 1, yend = rows_to_show * 3 + 1, linewidth = 0.5) +
    annotate("segment", x = 0, xend = 2, y = rows_to_show + 1, yend = rows_to_show + 1, linewidth = 0.5) +
    annotate("segment", x = 0, xend = 2, y = rows_to_show * 2 + 1, yend = rows_to_show * 2 + 1, linewidth = 0.5) +
    annotate("text", x = 0.5, y = rows_to_show * 2.5, label = "CRITICAL", color = "white", fontface = "bold", size = 2) +
    annotate("text", x = 0.5, y = rows_to_show * 1.5, label = "ABNORMAL", color = "white", fontface = "bold", size = 2) +
    annotate("text", x = 0.5, y = rows_to_show * 0.5, label = "NORMAL", color = "black", fontface = "bold", size = 2) +
    scale_x_continuous(name = "", breaks = c(0.5, 1.5), labels = c("Current", "Next"), expand = c(0,0)) +
    scale_y_continuous(name = "", expand = c(0,0)) +
    ggtitle(analyte_labels[analytes_to_show]) +
    theme(plot.title = element_text(hjust = 0.5, margin = margin(0, 0, 4, 0), face = "bold.italic", size = 12), strip.text = element_text(size = 10, face = "bold", margin = margin(0, 0, 4, 0)), axis.line = element_blank(), axis.text.y.left = element_blank(), legend.position = "none", axis.text.x.bottom = element_text(face = "bold", color = "black"))
  ggsave(out_path, gg, width = 3, height = 4)
  
  list(input = gg_input, plot = gg)
  
}

plotPredictionsOnUMAP <- function(input = read_csv("../../Results/Predictions/Supervised/final_pipeline_external_preds.csv"), umap_model = pin_read(model_board, "UMAP_ClinChem_Final") |> bundle::unbundle(), target = "final_prediction", pred = ".pred_class", path = "../../Figures/Supervised/full_pipeline_internal", subtitle_label = "WashU Model Predicting on WashU Data"){
  
  library(embed)
  theme_set(theme_ns)
  
  umap_coords <- bake(umap_model, input) |> select(UMAP1, UMAP2) |> bind_cols(input)
  
  gg_umap_input <- 
    umap_coords |> 
      select(matches("UMAP"), target = !!target, pred = !!pred) |> 
      mutate(tag = factor(paste(pred, target, sep = "_"), levels = c("0_0", "1_0", "1_1", "0_1")))
  
  gg_umap_labels <- 
    ggplot() +
      geom_point(data = gg_umap_input |> filter(tag == "0_0"), aes(UMAP1, UMAP2), color = "grey75", shape = ".") + 
      geom_point(data = gg_umap_input |> filter(tag != "0_0" & tag != "0_1"), aes(UMAP1, UMAP2, color = tag, shape = tag), size = 3) +
      scale_color_manual(values = c(scico::scico(3, begin = 0.1, end = 0.9, palette = "berlin", direction = -1)), labels = c("0_0" = "True Neg", "0_1" = "False Negative", "1_0" = "False Positive", "1_1" = "True Positive")) +     
      scale_shape_manual(values = c(17, 16, 15), labels = c("0_0" = "True Neg", "0_1" = "False Negative", "1_0" = "False Positive", "1_1" = "True Positive")) +
      labs(shape = "", color = "") + 
      ggtitle("UMAP Embedding of Basic Metabolic Panel Results", subtitle = subtitle_label) +
      theme(axis.text = element_blank(), legend.position = c(0.2, 0.1), legend.background = element_blank(), legend.key = element_blank(), plot.title = element_text(face = "bold", margin = margin(0, 0, 2, 0)))
  ggsave(paste0(path, "_umap_labels.png"), gg_umap_labels, width = 8, height = 8, dpi = 600)
  ggsave(paste0(path, "_umap_labels.pdf"), gg_umap_labels, width = 8, height = 8, dpi = 600)
  
  list(input = gg_umap_input, plot = gg_umap_labels)
  
}

plotTrainingMixRatios <- function(input = tar_read(ml_inputs)){
  
  gg_input <- 
    input |> 
      filter(fluid_name %in% c("NS", "D5NS")) |> 
      pluck("train") |> 
      bind_rows() |> 
      filter(label != "Patient")
  
  medians <- gg_input |> group_by(label) |> summarise(Median = median(mix_ratio), Mean = mean(mix_ratio)) |> pivot_longer(c(Median, Mean), names_to = "metric", values_to = "value")
  
  gg_mix_ratio <- 
    ggplot(gg_input, aes(y = mix_ratio, fill = label)) + 
      geom_density(show.legend = F) + 
      geom_text(data = medians, aes(x = 1.05, y = value, label = metric), size = 4, hjust = 0, fontface = "bold") +
      geom_segment(data = medians, aes(y = value, yend = value, x = 0, xend = 1), linetype = "dashed", show.legend = F) +
      scico::scale_fill_scico_d(alpha = 0.5) + 
      scale_y_continuous(name = "Mixture Ratio", breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) + 
      xlab("Density") + 
      facet_wrap(~fct_rev(label), labeller = as_labeller(c("NS" = "Non-D5", "D5NS" = "D5")), nrow = 1, ncol = 2, scales = "free") + 
      theme(legend.position = "none", axis.text.x.bottom = element_blank())
  ggsave("../../Figures/Supervised/final_figures/sup_fig_2_distribution_of_sim_mixture_ratios.pdf", gg_mix_ratio, height = 8, width = 8)
  ggsave("../../Figures/Supervised/final_figures/sup_fig_2_distribution_of_sim_mixture_ratios.svg", gg_mix_ratio, height = 8, width = 8)
  ggsave("../../Figures/Supervised/final_figures/sup_fig_2_distribution_of_sim_mixture_ratios.png", gg_mix_ratio, height = 8, width = 8, dpi = 600)
  
}

plotMixtureRatioLoss <- function(){
  
  tuned <- aggregateTargetObjects(train_cohort = "BJH", panel = "BMP", fluid_list = c(fluid_names), step_name = "wf_regression_tuned", recipe_name = "Base", feature_set = "all", models_string = "XGB")
  
  metrics <- tuned |> map("cv_metrics") |> bind_rows()
  
  loss <- metrics |> group_by(label) |> filter(.metric == "huber_loss") |> arrange(mean) |> slice_head(n = 1) |> ungroup() |> transmute(Fluid = str_match(label, "BMP_(.*?)_Base")[,2], HuberLoss = signif(mean, digits = 2), StdErr = signif(std_err, digits = 2)) |> arrange(HuberLoss)

  gt_loss <- gt::gt(loss) 
  gt::gtsave(gt_loss, "../../Figures/Supervised/final_figures/supp_table_regression_model_loss_cv.docx")
    
}

plotMixtureDistributions <- function(input = read_csv("../../Results/Predictions/Supervised/final_pipeline_internal_preds.csv"), target = "final_prediction", pred = ".pred_class", path = "../../Figures/Supervised/full_pipeline_internal_mix_ratios", subtitle_label = "WashU Model Predicting on WashU Data"){
  gt::gt(loss)
  
  
  preds <- makeMixtureRatioPredictions(input, fluid_list = fluid_names) |> mutate(contaminating_fluid_type_pred = ifelse(grepl("D5", contaminating_fluid_pred), "Any D5", "Non-D5"))
  
  gg_mix_input <- 
    preds |> 
      select(matches("pred"), target = !!target, pred = !!pred) |> 
      mutate(tag = factor(paste(pred, target, sep = "_"), levels = c("0_0", "0_1", "1_0", "1_1"), labels = c("TN", "FN", "FP", "TP")))
  
  medians <- gg_mix_input |> summarise(median = median(mix_ratio_pred, na.rm = T), .by = c(tag))
  
  theme_set(theme_ns)
  
  gg_mix_dist <- 
    ggplot() +
      stat_ecdf(data = gg_mix_input |> filter(tag == "TN"), aes(mix_ratio_pred, color = tag), na.rm = T, linewidth = 1, show.legend = F) + 
      stat_ecdf(data = gg_mix_input |> filter(tag %in% c("FP", "TP")), aes(mix_ratio_pred, color = tag), na.rm = T, linewidth = 2, show.legend = F) + 
      geom_text(data = medians |> filter(tag %in% c("TN", "FP", "TP")), aes(x = median + 0.0125, y = 0.5, label = tag, color = tag), hjust = 0, fontface = "bold", show.legend = F) + 
      labs(color = "") +
      scale_x_continuous(name = "Mixture Ratio", limits = c(0, 0.5)) + 
      scale_y_continuous(name = "Cumulative Proportion", limits = c(0, 1), breaks = c(0, 0.5, 1)) +
      scale_color_manual(values = c("darkred", "grey75", "black")) + 
      ggtitle("Distribution of Contamination Severity in Basic Metabolic Panel Results", subtitle = subtitle_label) +
      theme(legend.position = c(0.8, 0.25), legend.background = element_blank(), legend.key = element_blank(), plot.title = element_text(face = "bold", margin = margin(0, 0, 2, 0)))
  ggsave(paste0(path, "_mixture_distributions.png"), gg_mix_dist, width = 8, height = 4, dpi = 600)
  ggsave(paste0(path, "_mixture_distributions.pdf"), gg_mix_dist, width = 8, height = 4, dpi = 600)
  
  list(input = gg_mix_input, plot = gg_mix_dist)
  
}

compareFluidLists <- function(list_1 = list(name = "Final", fluids = c("NS", "SW", "D5W", "D5NS", "D5LR")), list_2 = list(name = "All", fluids = fluid_names)){
  
  pipeline_preds_list_1 <- 
    makePipelinePredictions(
      sim_preds_raw = aggregateTargetObjects(fluid_list = list_1$fluids, step_name = "wf_reviewed", feature_set = "results_with_priors", models_string = models_string), 
      tech_preds_raw = aggregateTargetObjects(fluid_list = c("Tech"), step_name = "wf_calibrated", feature_set = "results_with_priors", models_string = models_string)
    )
  
  pipeline_conf_mat_1 <- getMetricsAndConfusionMatrix(pipeline_preds_list_1)
  
  pipeline_preds_list_2 <- 
    makePipelinePredictions(
      sim_preds_raw = aggregateTargetObjects(fluid_list = list_2$fluids, step_name = "wf_reviewed", feature_set = "results_with_priors", models_string = models_string), 
      tech_preds_raw = aggregateTargetObjects(fluid_list = c("Tech"), step_name = "wf_calibrated", feature_set = "results_with_priors", models_string = models_string)
    )
  
  pipeline_conf_mat_2 <- getMetricsAndConfusionMatrix(pipeline_preds_list_2)
  
  gg_input <- 
    bind_rows(pipeline_conf_mat_1$metrics |> mutate(label = list_1$name), pipeline_conf_mat_2$metrics |> mutate(label = list_2$name))
  
  gg_fluid_comparison <- 
    ggplot(gg_input |> filter(Metric %in% c("sensitivity", "specificity", "mcc")), aes(x = Value, y = Metric, fill = label)) + 
      geom_col(position = "dodge", alpha = 0.75) + 
      geom_text(aes(x = 0.01, color = label, label = label), position = position_dodge(width = 0.9), alpha = 0.75, hjust = 0, fontface = "bold.italic", vjust = -1.28, size = 4) +
      geom_text(aes(x = Value - 0.02, color = label, label = round(Value, digits = 3)), position = position_dodge(width = 0.9), hjust = 1, vjust = 2, fontface = "bold", size = 4) +
      ylab("Prediction Source") + 
      scale_x_continuous(name = "", limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
      scico::scale_fill_scico_d(palette = "davos", begin = 0.2, end = 0.8) + 
      scale_color_manual(values = c(scico::scico(1, palette = "davos", direction = -1), "black", "black")) + 
      facet_wrap(~Metric, scales = "free", labeller = as_labeller(c("mcc" = "MCC", "sensitivity" = "Sensitivity", "specificity" = "Specificity"))) +
      theme(axis.text.y.left = element_blank(), legend.position = "none")
  ggsave("../../Figures/Supervised/final_figures/supp_fig_fluid_list_comparison_metrics_barchart.png", gg_fluid_comparison, width = 8.5, height = 3, dpi = 600)
  ggsave("../../Figures/Supervised/final_figures/supp_fig_fluid_list_comparison_metrics_barchart.pdf", gg_fluid_comparison, width = 8.5, height = 3, dpi = 600)
  
}

compareResultDistributions <- function(internal = pipeline_preds, external = pipeline_preds_yale){
  
  input = 
    bind_rows(internal |> mutate(source = "WashU"), external |> mutate(source = "Yale", anion_gap = sodium - chloride - co2_totl)) |>
      select(any_of(lab_strings_bmp), source)
  
  gg_input = 
     input |> 
      pivot_longer(cols = -source, names_to = "analyte", values_to = "result") 
  
  ggplot(gg_input, aes(result, color = source)) + 
    stat_ecdf(alpha = 0.75, linewidth = 2) + 
    facet_wrap(~analyte, scales = "free", labeller = as_labeller(analyte_labels), nrow = 3) + 
    xlab("Result Value") + ylab("Cumulative Density") + 
    ggtitle("Comparison of Analyte Distributions Across Sites") + 
    theme(axis.text.y.left = element_blank(), legend.position = c(0.95, 0.1), legend.background = element_blank(), legend.title = element_blank(), legend.box = element_blank(), legend.box.background = element_blank(), legend.key = element_blank(), plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(0, 0, 8, 0)))
  ggsave("../../Figures/Supervised/final_figures/supp_fig_4_result_distributions_by_analyte_and_site.png", width = 8, height = 8, dpi = 600)
  ggsave("../../Figures/Supervised/final_figures/supp_fig_4_result_distributions_by_analyte_and_site.pdf", width = 8, height = 8, dpi = 600)
  ggsave("../../Figures/Supervised/final_figures/supp_fig_4_result_distributions_by_analyte_and_site.svg", width = 8, height = 8, dpi = 600)
  
  library(gtsummary)
  tbl_summary(input, by = source, statistic = everything() ~ "{median}") |>
    italicize_levels() |>
    bold_labels() |> 
    add_p() |>  
    bold_p() |> 
    add_ci(method = list(pattern = "{stat} ({ci})")) |>  
    modify_header(all_stat_cols() ~ "**{level}**<br>*N = {n}*")
  
}

compareRCVbyClassLabel <- function(input = pipeline_yale_mixture_ratios |> getReferenceChangeValueFlags(), pred = ".pred_class", target = "final_prediction", out_path = "../../Figures/Supervised/final_figures/prospective_RCV_table"){
  
  library(gtsummary)
  
  gt_out <- 
    input |>
      mutate(tag = factor(paste(get(pred), get(target), sep = "_"), levels = c("0_0", "0_1", "1_0", "1_1"), labels = c("Neither", "Experts Only", "ML Only", "ML + Experts"))) |> 
      tbl_summary(by = tag, include = c("total_exceeding_RCV"), type = list(total_exceeding_RCV ~ "continuous"), statistic = list(all_continuous() ~ "{mean}"), label = list(total_exceeding_RCV ~ "Analytes Exceeding RCV")) |> 
      italicize_levels() |>
      bold_labels() |> 
      add_p() |>  
      bold_p() |> 
      add_ci(pattern = "{stat} {ci}", statistic = list(all_continuous() ~ "[{conf.low} - {conf.high}]")) |> 
      modify_header(all_stat_cols() ~ "**{level}**<br>*N = {n}*", label = " ") |> 
      modify_footnote(everything() ~ NA, abbreviation = T) |> 
      modify_footnote(update = everything() ~ NA) |>
      as_gt()
  gt::gtsave(gt_out, paste0(out_path, ".html"))
  gt::gtsave(gt_out, paste0(out_path, ".docx"))
  
}

compareSimToTechPreds <- function(input = pipeline_preds){
  
  sim_metrics <- metrics_binary(input, truth = "final_prediction", estimate = "sim_pred", event_level = "second") |> transmute(Metric = factor(.metric), value = .estimate, label = "Fluid-Specific")
  tech_metrics <- metrics_binary(input, truth = "final_prediction", estimate = "tech_pred", event_level = "second") |> transmute(Metric = factor(.metric), value = .estimate, label = "General")
  pipeline_metrics <- metrics_binary(input, truth = "final_prediction", estimate = ".pred_class", event_level = "second") |> transmute(Metric = factor(.metric), value = .estimate, label = "Combined")
  
  gg_input <- bind_rows(sim_metrics, tech_metrics, pipeline_metrics)
  
  gg_combo_bar <- 
    ggplot(gg_input |> filter(Metric %in% c("sensitivity", "specificity", "mcc")), aes(x = value, y = Metric, fill = label)) + 
      geom_col(position = "dodge", alpha = 0.75) + 
      geom_text(aes(x = 0.01, color = label, label = label), position = position_dodge(width = 0.9), alpha = 0.75, hjust = 0, fontface = "bold.italic", vjust = -1.28, size = 4) +
      geom_text(aes(x = value - 0.02, color = label, label = round(value, digits = 3)), position = position_dodge(width = 0.9), hjust = 1, vjust = 2, fontface = "bold", size = 4) +
      ylab("Prediction Source") + 
      scale_x_continuous(name = "", limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
      scico::scale_fill_scico_d(palette = "davos", begin = 0.2, end = 0.8) + 
      scale_color_manual(values = c(scico::scico(1, palette = "davos", direction = -1), "black", "black")) + 
      facet_wrap(~Metric, scales = "free", labeller = as_labeller(c("mcc" = "MCC", "sensitivity" = "Sensitivity", "specificity" = "Specificity"))) +
      theme(axis.text.y.left = element_blank(), legend.position = "none")
  ggsave("../../Figures/Supervised/final_figures/sim_vs_tech_vs_combo_metrics_barchart.png", gg_combo_bar, width = 8.5, height = 3, dpi = 600)
  ggsave("../../Figures/Supervised/final_figures/sim_vs_tech_vs_combo_metrics_barchart.pdf", gg_combo_bar, width = 8.5, height = 3, dpi = 600)
  
  
}

compareMethodPerformance <- function(input = pipeline_preds){
  
  ml <- getMetricsAndConfusionMatrix(input)
  current <- comparePredsToCurrentWorkflow(input)
  multivariate_deltas <- comparePredsToYaleRules(input)
  
  gg_input <- bind_rows(ml$metrics |> mutate(label = "ML"), current$metrics |> mutate(label = "Current"), multivariate_deltas$metrics |> mutate(label = "Yale Deltas")) |> mutate(label = factor(label, levels = c("Current", "Yale Deltas", "ML")))
  
  gg_combo_bar <- 
    ggplot(gg_input |> filter(Metric %in% c("sensitivity", "specificity", "mcc", "ppv", "npv")), aes(x = Value, y = Metric, fill = label)) + 
    geom_col(position = "dodge", alpha = 0.75) + 
    geom_text(aes(x = 0.01, color = label, label = label), position = position_dodge(width = 0.9), alpha = 0.75, hjust = 0, fontface = "bold.italic", vjust = -1.28, size = 4) +
    geom_text(aes(x = Value - 0.02, color = label, label = round(Value, digits = 3)), position = position_dodge(width = 0.9), hjust = 1, vjust = 2, fontface = "bold", size = 4) +
    ylab("Prediction Source") + 
    scale_x_continuous(name = "", limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
    scico::scale_fill_scico_d(palette = "davos", begin = 0.2, end = 0.8) + 
    scale_color_manual(values = c(scico::scico(1, palette = "davos", direction = -1), "black", "black")) + 
    facet_wrap(~Metric, scales = "free", labeller = as_labeller(c("mcc" = "MCC", "sensitivity" = "Sensitivity", "specificity" = "Specificity"))) +
    theme(axis.text.y.left = element_blank(), legend.position = "none")
  ggsave("../../Figures/Supervised/final_figures/current_vs_multi_vs_ml_metrics_barchart.png", gg_combo_bar, width = 8.5, height = 3, dpi = 600)
  ggsave("../../Figures/Supervised/final_figures/current_vs_multi_vs_ml_metrics_barchart.pdf", gg_combo_bar, width = 8.5, height = 3, dpi = 600)
  ggsave("../../Figures/Supervised/final_figures/current_vs_multi_vs_ml_metrics_barchart.svg", gg_combo_bar, width = 8.5, height = 3, dpi = 600)
  
}

plotApplicabilityScores <- function(train = tar_read(preprocessed_bmp_inputs)[[1]][["data"]], internal = prospective_preds |> arrange(drawn_dt_tm), external = pipeline_preds_yale, apd_model = tar_read(apd_pca)){
  
  theme_set(theme_ns)
  
  train_apd = getApplicabilityScores(input = train, apd_model) |> bind_cols(train |> select(drawn_dt_tm)) |> drop_na()
  internal_apd = getApplicabilityScores(input = internal, apd_model) |> bind_cols(internal |> select(drawn_dt_tm))
  external_apd = getApplicabilityScores(input = pipeline_preds_yale, apd_model)
  
  train_rolling_sums <- 
    train_apd |> 
      mutate(
        median = zoo::rollmedian(distance, k = 1001, fill = NA)
        )
  
  Kendall::MannKendall(train_rolling_sums$distance)
  
  train_quantiles <- quantile(train_apd$distance, probs = c(0.25, 0.50, 0.75), na.rm = T) %>% bind_cols(., c("25th %ile", "Median", "75th %ile")) |> set_names(c("value", "percentile")) |> mutate(percentile = factor(percentile, levels = c("25th %ile", "Median", "75th %ile")))
  hjusts <- bind_cols(source = c("WashU Train", "WashU Prospective", "Yale Prospective"), hjust = c(0.4, 0.6, 0.8))
  
  gg_input = bind_rows(train_apd |> mutate(source = "WashU Train"), internal_apd |> mutate(source = "WashU Prospective"), external_apd |> mutate(source = "Yale Prospective")) |> left_join(hjusts)

  gg_dens <- 
    gg_input |> 
        ggplot(aes(distance, color = source, label = source, hjust = hjust, linetype = fct_rev(source))) +
        geomtextpath::geom_textdensity(alpha = 0.8, fontface = "bold", size = 3, vjust = -2) + 
        scale_x_continuous(name = "Distance", limits = c(2, 10)) +
        ylab("Proportion") + 
        ggtitle("Distribution of PCA Distances Across Data Sets") + 
        scico::scale_color_scico_d(palette = "batlow", begin = 0.8, end = 0.2) + 
        theme(legend.position = "none", axis.text.y.left = element_blank(), plot.title = element_text(face = "bold.italic"))
  
  gg_time <- 
    gg_input |> 
      filter(source == "WashU Prospective") |> 
      drop_na() |> 
      mutate(
          median = zoo::rollmedian(distance, k = 1001, fill = NA)) |> 
        ggplot(aes(drawn_dt_tm, median)) +
        geom_line(na.rm = T) + 
        geom_hline(data = train_quantiles, aes(yintercept = value, color = percentile), linetype = "dashed") + 
        scico::scale_color_scico_d(palette = "bilbao", begin = 0.8, end = 0.2) + 
        scale_y_continuous(name = "Median Distance", limits = c(2.8, 5)) + 
        xlab("Collection Date") + 
        coord_cartesian(clip = "off") + 
        ggtitle("Rolling Median of PCA Distance Over Time", subtitle = "Internal Prospective Validation Set") + 
        theme(legend.position = "none", plot.title = element_text(face = "bold.italic"))
  
  gg_combo <- ggpubr::ggarrange(gg_dens, gg_time, nrow = 2, ncol = 1, labels = "AUTO")
  ggsave("../../Figures/Supervised/final_figures/fig5_drift_and_applicability_figure.svg", gg_combo, width = 8.5, height = 7)
  ggsave("../../Figures/Supervised/final_figures/fig5_drift_and_applicability_figure.png", gg_combo, width = 8.5, height = 7)
  ggsave("../../Figures/Supervised/final_figures/fig5_drift_and_applicability_figure.pdf", gg_combo, width = 8.5, height = 7)
  
  gg_combo
  
}
