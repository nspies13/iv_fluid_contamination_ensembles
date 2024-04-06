getMetrics <- function(probs, target, threshold = 0.5){
  
  preds = factor(ifelse(probs > threshold, 1, 0), levels = c("0", "1"))
  class = bind_cols(truth = target, estimate = preds)
  num = bind_cols(truth = target, estimate = probs)
  
  conf_mat = conf_mat(class, truth, estimate)
  class_metrics = metrics_binary(class, truth = truth, estimate = estimate, event_level = "second") %>% dplyr::select(.metric, .estimate)
  numeric_metrics = metrics_numeric(num, truth = truth, estimate, event_level = "second") %>% dplyr::select(.metric, .estimate)
  
  list(confusion_matrix = conf_mat, metrics = bind_rows(class_metrics, numeric_metrics))
  
}

plotValidationCurves <- function(val_preds, threshold_value = 0.5, label, out_path = "../../Figures/Supervised/curve_plots/"){

  library(gridExtra)
  metric_list <- getMetrics(val_preds$.pred_1, val_preds$target, threshold_value)
  conf_table <- tableGrob(metric_list[[2]] %>% filter(!.metric %in% c("binary_classification_cost")) %>% transmute(Metric = factor(.metric, levels = names(metric_labels), labels = metric_labels), Value = round(.estimate, digits = 2)), theme = ttheme_minimal(base_size = 8, padding = unit(c(4, 4), units = "pt")), rows = NULL)
  
  pr = val_preds %>% pr_curve(.pred_1, truth = target, event_level = "second") %>% mutate(curve = "pr") %>% bind_rows(tibble_row(.threshold = 0, recall = 1, precision = 0))
  roc = val_preds %>% roc_curve(.pred_1, truth = target, event_level = "second") %>% mutate(curve = "roc")
  
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
  ggarrange(gg_pr, gg_roc, nrow = 1, ncol = 2) %>%
    annotate_figure(top = text_grob(label, face = "bold.italic", size = 18))
  ggsave(paste0(out_path, ".png"), width = 8, height = 4, dpi = 600)
  ggsave(paste0(out_path, ".pdf"), width = 8, height = 4, dpi = 600)

  list(gg_roc, gg_pr)
  
}

optimizeThreshold <- function(probs, target, metric_to_maximize = "mcc", out_dir = "../../Figures/Models/", label = "") {
  
  threshold_dir = paste0(out_dir, "/Thresholds/", label)
  if (!file.exists(threshold_dir)){
    dir.create(threshold_dir)
  }
  
  seqs <- c(0.0001, 0.001, seq(0.01, 0.49, by = 0.005))
  
  data = bind_cols(truth = target, estimate = probs)
  
  thresholds = probably::threshold_perf(data, truth = truth, estimate = estimate, thresholds = c(seqs, 0.5, 1-seqs), metrics = metrics_binary, event_level = "second")
  threshold_value = thresholds[which.max(thresholds %>% filter(.metric == metric_to_maximize) %>% pluck(".estimate")),".threshold"][[1]]
  
  library(geomtextpath)
  metrics_to_plot = c("mcc", "ppv", "npv")
  
  ggplot(thresholds %>% filter(.metric %in% metrics_to_plot), aes(.threshold, .estimate, color = .metric)) + 
    geom_vline(xintercept = threshold_value, linetype = "dashed", alpha = 0.5) +
    geom_textline(aes(label = toupper(.metric)), linewidth = 2, hjust = 0.45, show.legend = F, fontface = "bold") + 
    scale_color_viridis_d() + 
    ylim(0, 1) +
    annotate("text", x = threshold_value + 0.02, y = 0.33, angle = -90, label = "Optimal Threshold", fontface = "italic", color = "grey50") + 
    xlab("Threshold") + ylab("Metric Value")
  ggsave(paste0(threshold_dir, "/threshold_plot.pdf"), width = 8, height = 8)
  
  threshold_value
  
}

explainModel <- function(wf_fit, validation_set = tar_read(wf_prospective_BJH_Sim_BMP_D5halfNSwK_Base_NNet_results)) {
  
  wf_fit = wf_fit %>% bundle::unbundle()
  
  predictors = wf_fit %>% extract_preprocessor() %>% pluck("var_info") %>% filter(role == "predictor") %>% pluck("variable")
  label = wf_fit %>% extract_preprocessor() %>% pluck("label")
  validation_no_outliers <- validation_set %>% filter(if_all(where(is.numeric), ~(. <= quantile(., 0.975, na.rm = TRUE) & . >= quantile(., 0.025, na.rm= T))))
  fluid = wf_fit %>% extract_preprocessor() %>% pluck("fluid_name")
  out_dir = paste0("../../Figures/Models/Explainers/", label)
  
  if(!dir.exists(out_dir)){dir.create(out_dir)}
  
  library(iml)
  
  predict_wrapper <- function(model, newdata){workflows:::predict.workflow(object = model, new_data = newdata, type = "prob")} 
  predictor <- Predictor$new(model = wf_fit, data = validation_no_outliers %>% dplyr::select(-patient_id, -target), y = as.numeric(validation_no_outliers[["target"]]), predict.function = predict_wrapper, type = "prob", class = 2)
  
  plotVarImps(predictor, out_dir = out_dir)
  plotALE(predictor, out_dir = out_dir)

}

plotVarImps <- function(predictor, out_dir = here::here()) {
  
  var_imps = FeatureImp$new(predictor, loss = "logLoss", compare = "difference", n.repetitions = 10)
  
  plot(var_imps) + 
    scale_x_continuous(name = "Change in Log Loss") + 
    scale_y_discrete(name = "Permuted Feature") + 
    theme()
  ggsave(paste0(out_dir, "/variable_importances.pdf"), width = 8, height = 4) 
  
}

plotALE <- function(predictor, out_dir = here::here()){
  
  ale = FeatureEffects$new(predictor, method = "ale")
  
  gg_input = purrr::map(ale$features, ~ale$effects[[.x]]$results %>% as_tibble() %>% 
        dplyr::transmute(ALE = .value, Result = get(.x), analyte = .x)) %>%
  bind_rows()

  ggplot(gg_input %>% filter(analyte != "anion_gap"), aes(Result, ALE)) + 
    geom_line() +
    facet_wrap(~analyte, scales = "free_x", labeller = as_labeller(analyte_labels), nrow = 2) +
    theme()
  ggsave(paste0(out_dir, "/ale_plots.pdf"), width = 8, height = 6)
  
}

plotALE2D <- function(predictor, out_dir = here::here()){
  
  ale = map(lab_strings_bmp, ~FeatureEffect$new(predictor, feature = c(paste0(.x, "_delta_prior"), paste0(.x, "_delta_post")), method = "ale", grid.size = 1000))
  
  plots <- 
    map2(ale, lab_strings_bmp,
         ~.x$plot() + 
            scico::scale_fill_scico(palette = "roma", direction = -1, begin = 0.1, end = 0.9) +
            labs(title = analyte_labels[[.y]], 
                 x = paste0(analyte_labels[[.y]], " Delta from Prior"),
                 y = paste0(analyte_labels[[.y]], " Delta to Post")) +
            theme(axis.title.x.bottom = element_text()))
  
}



assign("custom_costs", tribble(
  ~truth,   ~estimate, ~cost,
  "0", "1",  1,
  "1", "0",  2), envir = globalenv())

classification_cost_2to1 <- metric_tweak("classification_cost_2to1", yardstick::classification_cost, costs = tribble(
  ~truth,   ~estimate, ~cost,
  "0", "1",  1,
  "1", "0",  2))

classification_cost_1to10 <- metric_tweak("classification_cost_1to10", yardstick::classification_cost, costs = tribble(
  ~truth,   ~estimate, ~cost,
  "0", "1",  1,
  "1", "0",  10))

assign("metrics_global", 
       metric_set(mn_log_loss, mcc, pr_auc, pr_auc, roc_auc, sensitivity, specificity, ppv, npv, accuracy), 
       envir = globalenv())

assign("metrics_numeric", 
       metric_set(mn_log_loss, pr_auc, roc_auc), 
       envir = globalenv())

assign("metrics_binary", 
       metric_set(mcc, sensitivity, specificity, ppv, npv, accuracy), 
       envir = globalenv())

assign("metrics_regression", 
       metric_set(huber_loss, huber_loss_pseudo, smape, mape, rmse, rsq, mae, msd), 
       envir = globalenv())
