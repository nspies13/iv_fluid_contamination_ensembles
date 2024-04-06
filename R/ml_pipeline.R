definePredictorSets <- function(input = train){
  
  list(
    results = input %>% select(any_of(!!lab_strings)) %>% names(),
    results_with_priors = input %>% select(any_of(!!lab_strings), matches("_prior"), -matches("post|real|prop|gap")) %>% names(),
    current_with_deltas = input %>% select(any_of(!!lab_strings), matches("_delta_prior"), -matches("post|real|prop|gap")) %>% names(),
    deltas = input %>% select(matches("delta"), -matches("post|prop|gap")) %>% names(),
    current_with_post = input %>% select(any_of(!!lab_strings), matches("_post"), -matches("prior|real|prop|delta|gap")) %>% names(),
    all = input %>% select(any_of(!!lab_strings_bmp), matches("_delta_prior|_delta_post"), -matches("real|prop|gap")) %>% names()
  )
  
}

defineWorkflow <- function(rec, model, train_cohort = "", fluid_name = "", pred_set = "", recipe_name = "", model_name = "", label = "", ...){
  
  wf = workflow(rec, model)

  wf

}

tuneWorkflow <- function(wf = wf, cv = cv_comments, label = ""){
  
  param = wf %>% extract_parameter_set_dials()
  grid = length(param %>% pluck(1)) + 3
  
  if(any(wf %>% extract_parameter_set_dials() %>% pluck("id") == "mtry")){
    param = wf %>% extract_parameter_set_dials() %>% finalize(wf %>% extract_preprocessor() %>% pluck("template"))
  }
  
  library(doParallel)
  library(parallelly)
  cl <- parallelly::makeClusterPSOCK(5)
  doParallel::registerDoParallel(cl)
  
  wf_tuned = 
    wf %>% 
      tune_grid(
        resamples = vfold_cv(wf %>% extract_preprocessor() %>% pluck("template"), v = 5), 
        grid = grid,
        param_info = param,
        metrics = metrics_global, 
        control = control_grid(verbose = F, allow_par = T, parallel_over = "resamples",
                                save_pred = F, save_workflow = F, event_level = "second"))
  
  tune_output = list(wf = wf, wf_tuned = wf_tuned, cv_metrics = wf_tuned %>% collect_metrics() |> mutate(label = paste0(label, "_", Sys.time())), best_params = wf_tuned %>% select_best("mn_log_loss"))
  write_csv(tune_output$cv_metrics, "../../Results/Models/Tuning/tune_outputs.csv", append = T, col_names = T)
  
  tune_output
  
}

tuneRegWorkflow <- function(wf = wf, cv = cv_comments, label = ""){
  
  param = wf %>% extract_parameter_set_dials()
  grid = length(param %>% pluck(1)) + 3
  
  if(any(wf %>% extract_parameter_set_dials() %>% pluck("id") == "mtry")){
    param = wf %>% extract_parameter_set_dials() %>% finalize(wf %>% extract_preprocessor() %>% pluck("template"))
  }
  
  library(doParallel)
  library(parallelly)
  cl <- parallelly::makeClusterPSOCK(5)
  doParallel::registerDoParallel(cl)
  
  wf_tuned = 
    wf %>% 
    tune_grid(
      resamples = vfold_cv(wf %>% extract_preprocessor() %>% pluck("template"), v = 5), 
      grid = grid,
      param_info = param,
      metrics = metrics_regression, 
      control = control_grid(verbose = F, allow_par = T, parallel_over = "resamples",
                             save_pred = F, save_workflow = F, event_level = "second"))
  
  tune_output = list(wf = wf, wf_tuned = wf_tuned, cv_metrics = wf_tuned %>% collect_metrics() |> mutate(label = paste0(label, "_", Sys.time())), best_params = wf_tuned %>% select_best("huber_loss"))
  write_csv(tune_output$cv_metrics, "../../Results/Models/Tuning/tune_outputs.csv", append = T, col_names = T)
  
  tune_output
  
}

tuneWorkflowBayes <- function(wf = wf){
  
  param = wf %>% extract_parameter_set_dials()
  grid = length(param %>% pluck(1)) + 2
  
  if(any(wf %>% extract_parameter_set_dials() %>% pluck("id") == "mtry")){
    param = wf %>% extract_parameter_set_dials() %>% finalize(wf %>% extract_preprocessor() %>% pluck("template"))
  }
  
  library(doParallel)
  library(parallelly)
  cl <- parallelly::makeClusterPSOCK(5)
  doParallel::registerDoParallel(cl)
  
  wf_tuned = 
    wf %>% 
    tune_bayes(
      resamples = group_vfold_cv(wf %>% extract_preprocessor() %>% pluck("template"), group = patient_id, v = 5),
      iter = 20,
      initial = grid, 
      param_info = param,
      metrics = metrics_global, 
      control = control_bayes(verbose_iter = T, verbose = T, no_improve = 7, uncertain = 3, allow_par = T, parallel_over = "resamples",
                             save_pred = F, save_workflow = F, event_level = "second"))
  
  tune_output = list(wf = wf, cv_metrics = wf_tuned %>% collect_metrics(), best_params = wf_tuned %>% select_best("mn_log_loss"))
  
  tune_output
  
}

fitAll <- function(tune_output){
  
  wf = tune_output[["wf"]]
  rec = wf %>% extract_preprocessor()
  
  train = rec %>% pluck("template")

  wf_fit = wf %>% finalize_workflow(tune_output[["best_params"]]) %>% fit(train)
  
  wf_fit = wf_fit %>% butcher::butcher(verbose = T)
  
  wf_fit$metadata = list(
    train_cohort = rec$train_cohort,
    train_type = rec$train_type,
    panel = rec$panel,
    fluid_name = rec$fluid_name,
    pred_set = rec$pred_set,
    recipe_name = rec$recipe_name,
    label = rec$label)
  
  wf_fit %>% bundle::bundle()
  
}

saveVetiverWorkflow <- function(wf_fit = tar_read(wf_fit_BJH_Sim_BMP_NS_Base_XGB_results_with_priors), model_board = board_s3("preanalytical-error-models", versioned = T, region = "us-east-2"), label = "tmp", ...){
  
  library(pins)
  library(vetiver)

  wf_fit = wf_fit %>% bundle::unbundle()
  metadata = wf_fit[["metadata"]]

  v <- 
    vetiver::vetiver_model(
      model = wf_fit, 
      model_name = paste0("wf_vetiver_", label),
      metadata = metadata,
      versioned = T
    )
  
  model_board %>%
    vetiver_pin_write(
      vetiver_model = v,
      check_renv = T
    )
  
  v
  
}

calibratePredictedProbabilities <- function(preds = tar_read(wf_reviewed_BJH_Tech_BMP_Base_XGB_results_with_priors), target = "final_prediction", path =  "../../Figures/Supervised/test_calibration_plot"){
  
  library(probably)
  
  calibrator = probably::cal_estimate_isotonic(preds, truth = !!target)
  
  cal_preds <- cal_apply(preds, calibrator, pred_class = ".pred_class")
  gg_cal_plot <- cal_plot_windowed(cal_preds, truth = !!target)
  ggsave(paste0(path, ".png"), gg_cal_plot, width = 4, height = 4)
  
  cal_preds
  
}

buildApplicabilityModel <- function(train = tar_read(preprocessed_bmp_inputs)[[1]][["data"]], features = tar_read(pred_cols_tar)[["results_with_priors"]]){
  
  library(applicable)
  
  input = train |> select(any_of(features)) |> drop_na()
  apd_model <- apd_pca(input)
  
  apd_model
  
}

getApplicabilityScores <- function(input = tar_read(prospective), apd_model = apd_model){
  
  library(applicable)
  
  score(apd_model, input) |> select(distance, distance_pctl)
  
}