##### Set Configs #####
library(targets)
library(tarchetypes) 
library(tidyverse)
library(tidymodels)
tidymodels_prefer()
conflicted::conflicts_prefer(future::run)

tar_option_set(
  packages = c("tidyverse", "tidymodels"), format = "qs", deployment = "worker", 
  envir = globalenv(), tidy_eval = T, error = "continue", iteration = "list")

future::plan(future.callr::callr, workers = 8)
tar_source()
tar_source("../Helpers/SetGlobals.R")
tar_source("../Simulation/R/simulate_contamination.R")
tar_source("../Retrospective/R/retrospective_contamination.R")
tar_source("../Unsupervised/R/unsupervised_pipeline.R")
theme_set(theme_ns)

##### Prep ML Inputs #####
load_inputs <- list(
  
  # Read preprocessed data from LIS extract board.
  tar_target(preprocessed_bmp_inputs, list(list(data = data_board %>% pin_read("BJH_bmps_with_error_flags") %>% 
                                                  drop_na(any_of(lab_strings_bmp)) %>% 
                                                  rename(patient_id = person_id), label = "BJH"))),
  tar_target(contam_sim_tar, read_csv("../../Data/contam_sim_all.csv") |> mutate(label = ifelse(label == "Water", "SW", label))),
  tar_target(train_cohort, preprocessed_bmp_inputs[["label"]], pattern = map(preprocessed_bmp_inputs)),
  tar_target(panels, c("BMP")),
  tar_target(prospective, data_board %>% pin_read("BJH_bmps_2023_with_error_flags")),
  
  # Filter BMPs with error comments
  tar_target(bmp_uncontaminated, preprocessed_bmp_inputs[["data"]] %>% filter(!code_comment & !unlikely_comment & !contam_comment & !mislabel_comment & !invalid_comment & !credited_comment) %>% select(-matches("prop")) %>% mutate(index = row_number()), pattern = map(preprocessed_bmp_inputs)),

  # Make train-test splits
  tar_target(train_input, list(data = bmp_uncontaminated, train_cohort = train_cohort), pattern = map(bmp_uncontaminated, train_cohort)), 

  # Define fluids
  tar_target(fluids_tar, fluids),
  tar_target(fluid_names_tar, names(fluids_tar)),
  
  # Define mixture parameters
  tar_target(contam_rate, c(0.25)),
  
  # Simulate contamination in training data
  tar_target(train, list(train = makeSimulatedTrainingData(train_input[["data"]], 
                                              makeMixRatios(contam_rate, train_input[["data"]], fluid_names_tar), 
                                              fluids_tar, fluid_names_tar) |> select(-matches("real|anion_gap")), 
                         train_cohort = train_input[["train_cohort"]], 
                         fluid_name = fluid_names_tar), 
             pattern = cross(train_input, contam_rate, map(fluid_names_tar, fluids_tar))),

  # Define predictor sets
  tar_target(pred_cols_tar, definePredictorSets(train[[1]][["train"]])),
  
  # Prepare ML inputs
  tar_target(ml_inputs, tibble(train_cohort = train %>% map("train_cohort") %>% unlist(), 
                               train_type = "Sim",
                               fluid_name = train %>% map("fluid_name") %>% unlist(),
                               train = train %>% map("train"))),
  tar_target(yale_prospective, read_csv("../../Data/yale_ml_input.csv") |> mutate(final_prediction = factor(final_prediction))),
  
  # Build APD model
  tar_target(apd_pca, buildApplicabilityModel(preprocessed_bmp_inputs[[1]][["data"]], features = pred_cols_tar[["current_with_deltas"]]))
  
)

##### Build ML Workflows #####
build_workflows <- tar_map(

  values = expand_grid(cohort_map = c("BJH"),
                       type_map = c("Sim"),
                       panel_map = c("BMP"),
                       recipe_map = getRecipeCalls(),
                       model_map = getModelCalls(),
                       fluid_map = fluid_names,
                       set_map = c("current_with_deltas")) %>%
    mutate(label = str_replace_all(paste(cohort_map, type_map, panel_map, fluid_map, recipe_map, model_map, set_map, sep = "_"), "define|Recipe", "")),
  names = "label",

  tar_target(train, ml_inputs %>%
               filter(train_cohort == cohort_map & train_type == type_map & fluid_name == fluid_map) %>%
               pluck("train", 1) %>%
               drop_na(any_of(pred_cols_tar[[set_map]]), target) %>%
               select(target, any_of(pred_cols_tar[[set_map]]))),
  tar_target(rec, recipe_map(ml_inputs %>%
                               filter(train_cohort == cohort_map & train_type == type_map & fluid_name == fluid_map) %>%
                               pluck("train", 1) %>%
                               drop_na(any_of(pred_cols_tar[[set_map]]), target) %>%
                               select(target, any_of(pred_cols_tar[[set_map]])),
                             train_cohort = cohort_map, panel = panel_map, fluid_name = fluid_map, pred_set = set_map, pred_cols = pred_cols_tar[[set_map]], label = label)),
  tar_target(model, model_map()),
  tar_target(wf, defineWorkflow(rec = rec, model = model[[1]])),
#  tar_target(wf, workflow() |> add_model(model[[1]]) |> add_recipe(recipe(target ~ ., train))),
  tar_target(wf_tuned, tuneWorkflow(wf, label = label)),
  tar_target(wf_cv_metrics, wf_tuned |> pluck("cv_metrics")),
  tar_target(wf_fit, fitAll(wf_tuned)),
#  tar_target(wf_explain, explainModel(wf_fit)),
  tar_target(wf_vetiver, saveVetiverWorkflow(wf_fit, model_board, label = label), deployment = "main"),
  tar_target(wf_sim_sensitivity, augment(bundle::unbundle(wf_fit), contam_sim_tar |> slice_sample(n = 1000, by = c("mix_ratio", "label")))),
  
# Assess Model Performance on Real Expert Labels
    tar_target(wf_prospective, bind_cols(prospective |> drop_na(matches(paste0(lab_strings, collapse = "|"))), predict(bundle::unbundle(wf_fit), prospective |> drop_na(matches(paste0(lab_strings, collapse = "|")))), predict(bundle::unbundle(wf_fit), prospective |> drop_na(matches(paste0(lab_strings, collapse = "|"))), type = "prob"))),
    tar_target(wf_reviewed, getManualReviewLabels(wf_prospective, consensus = readxl::read_xlsx("../../Data/results_to_review_formatted_consensus.xlsx"))),
#    tar_target(wf_calibrated, calibratePredictedProbabilities(wf_reviewed, path = paste0("../../Figures/Supervised/", label, "_cal_plot"))),
    tar_target(wf_test_metrics, getMetrics(wf_reviewed$.pred_1, wf_reviewed$final_prediction)),
    tar_target(wf_test_curves, makeCurves(wf_reviewed, target = "final_prediction", path = paste0("../../Figures/Supervised/curve_plots/", label, "_expert_label_roc_pr_curves"))),
    tar_target(wf_pre_post_scatter, makePrePostScatterplots(wf_reviewed, target = "final_prediction", path = paste0("../../Figures/Supervised/pre_post_plots/", label, "_pre_post_scatter"), analytes_to_include = c("sodium", "chloride", "co2_totl", "calcium", "glucose"))), 
    #tar_target(wf_sim_sensitivity, makeSimSensitivityPlots(wf_fit = wf_fit, path = paste0("../../Figures/Supervised/", label, "_sim_sensitivity"))),
    
    tar_target(wf_yale_reviewed, bind_cols(yale_prospective |> drop_na(matches(paste0(lab_strings, collapse = "|"))), predict(bundle::unbundle(wf_fit), yale_prospective |> drop_na(matches(paste0(lab_strings, collapse = "|")))), predict(bundle::unbundle(wf_fit), yale_prospective |> drop_na(matches(paste0(lab_strings, collapse = "|"))), type = "prob"))),
    #tar_target(wf_yale_calibrated, calibratePredictedProbabilities(wf_yale_reviewed, path = paste0("../../Figures/Supervised/", label, "_cal_plot_yale"))),
    tar_target(wf_yale_validation_metrics, getMetrics(wf_yale_reviewed$.pred_1, wf_yale_reviewed$final_prediction)),
    tar_target(wf_yale_validation_curves, makeCurves(wf_yale_reviewed, target = "final_prediction", path = paste0("../../Figures/Supervised/curve_plots/", label, "_roc_pr_curves_yale"))),
    tar_target(wf_yale_pre_post_scatter, makePrePostScatterplots(wf_yale_reviewed, target = "final_prediction", path = paste0("../../Figures/Supervised/pre_post_plots/", label, "_pre_post_scatter_yale"), analytes_to_include = c("sodium", "chloride", "co2_totl", "calcium", "glucose")))

)

#### Pipeline for single, combined model ####
build_binary_comment_predictors <- list(

  tar_target(input_comments, preprocessed_bmp_inputs[["data"]] %>% filter(!code_comment & !mislabel_comment & !invalid_comment & !credited_comment) %>% select(-matches("prop")) %>% mutate(index = row_number(), target = factor(as.numeric(unlikely_comment | contam_comment))), pattern = map(preprocessed_bmp_inputs)),
  tar_target(smote_rec, recipe(target ~ ., data = input_comments |> select(target, matches(paste(lab_strings, collapse = "|"))) |> drop_na()) |> themis::step_smote(target, over_ratio = contam_rate) |> prep(), pattern = map(input_comments)),
  tar_target(train_comments, bake(smote_rec[[1]], new_data = NULL) |> mutate(patient_id = row_number())),

  tar_map(

    values = expand_grid(cohort_map = c("BJH"), type_map = c("Tech"), panel_map = c("BMP"),
                       recipe_map = getRecipeCalls(), model_map = getModelCalls(),
                       set_map = c("current_with_deltas")) %>%
    mutate(label = str_replace_all(paste(cohort_map, type_map, panel_map, recipe_map, model_map, set_map, sep = "_"), "define|Recipe", "")),
    names = "label",

    tar_target(rec, recipe_map(train_comments %>% drop_na(any_of(pred_cols_tar[[set_map]]), target) %>% select(target, patient_id, any_of(pred_cols_tar[[set_map]])), train_cohort = cohort_map, panel = panel_map, fluid_name = "Tech", pred_set = set_map, pred_cols = pred_cols_tar[[set_map]], label = label)),
    tar_target(model, model_map()),
    tar_target(wf, defineWorkflow(rec = rec, model = model[[1]])),
    tar_target(wf_tuned, tuneWorkflow(wf, label = label)),
    tar_target(wf_cv_metrics, wf_tuned |> pluck("cv_metrics")),
    tar_target(wf_fit, fitAll(wf_tuned)),
    tar_target(wf_sim_sensitivity, augment(bundle::unbundle(wf_fit), contam_sim_tar |> slice_sample(n = 1000, by = c("mix_ratio", "label")))),
#    tar_target(wf_explain, explainModel(wf_fit)),
    tar_target(wf_vetiver, saveVetiverWorkflow(wf_fit, model_board), deployment = "main"),

    tar_target(wf_prospective, bind_cols(prospective |> drop_na(matches(paste0(lab_strings, collapse = "|"))), predict(bundle::unbundle(wf_fit), prospective |> drop_na(matches(paste0(lab_strings, collapse = "|")))), predict(bundle::unbundle(wf_fit), prospective |> drop_na(matches(paste0(lab_strings, collapse = "|"))), type = "prob"))),
    tar_target(wf_reviewed, getManualReviewLabels(wf_prospective, consensus = readxl::read_xlsx("../../Data/results_to_review_formatted_consensus.xlsx"))),
    tar_target(wf_calibrated, calibratePredictedProbabilities(wf_reviewed, path = paste0("../../Figures/Supervised/calibration_plots/", label, "_cal_plot"))),
    tar_target(wf_test_metrics, getMetrics(wf_calibrated$.pred_1, wf_calibrated$final_prediction)),
    tar_target(wf_test_curves, makeCurves(wf_calibrated, target = "final_prediction", path = paste0("../../Figures/Supervised/curve_plots/", label, "_roc_pr_curves"))),
    tar_target(wf_pre_post_scatter, makePrePostScatterplots(wf_calibrated, target = "final_prediction", path = paste0("../../Figures/Supervised/pre_post_plots/", label, "_pre_post_scatter"), analytes_to_include = c("sodium", "chloride", "co2_totl", "calcium", "glucose"))), 
    #tar_target(wf_sim_sensitivity, makeSimSensitivityPlots(wf_fit = wf_fit, path = paste0("../../Figures/Supervised/", label, "_sim_sensitivity"))),

    tar_target(wf_yale_reviewed, bind_cols(yale_prospective |> drop_na(matches(paste0(lab_strings, collapse = "|"))), predict(bundle::unbundle(wf_fit), yale_prospective |> drop_na(matches(paste0(lab_strings, collapse = "|")))), predict(bundle::unbundle(wf_fit), yale_prospective |> drop_na(matches(paste0(lab_strings, collapse = "|"))), type = "prob"))),
    tar_target(wf_yale_calibrated, calibratePredictedProbabilities(wf_yale_reviewed, path = paste0("../../Figures/Supervised/calibration_plots/", label, "_cal_plot_yale"))),
    tar_target(wf_yale_validation_metrics, getMetrics(wf_yale_calibrated$.pred_1, wf_yale_calibrated$final_prediction)),
    tar_target(wf_yale_validation_curves, makeCurves(wf_yale_calibrated, target = "final_prediction", path = paste0("../../Figures/Supervised/curve_plots/", label, "_roc_pr_curves_yale"))),
    tar_target(wf_yale_pre_post_scatter, makePrePostScatterplots(wf_yale_calibrated, target = "final_prediction", path = paste0("../../Figures/Supervised/pre_post_plots/", label, "_pre_post_scatter_yale"), analytes_to_include = c("sodium", "chloride", "co2_totl", "calcium", "glucose"))) 

  )

)

#### Build Mixture Ratio Regression Models #####
build_regression_workflows <- tar_map(
  
  
  values = expand_grid(cohort_map = c("BJH"),
                       type_map = c("Sim"),
                       panel_map = c("BMP"),
                       recipe_map = getRecipeCalls(),
                       model_map = getRegModelCalls(),
                       fluid_map = fluid_names,
                       set_map = c("all")) %>%
    mutate(label = str_replace_all(paste(cohort_map, type_map, panel_map, fluid_map, recipe_map, model_map, set_map, sep = "_"), "define|Recipe", "")),
  names = "label",
  
  tar_target(rec, recipe_map(ml_inputs %>%
                               filter(train_cohort == cohort_map & train_type == type_map & fluid_name == fluid_map) %>%
                               pluck("train", 1) %>%
                               slice_sample(prop = 0.1) |> 
                               bind_rows(contam_sim_tar |> filter(label == fluid_map) |> select(-patient_id)) |> 
                               mutate(target = mix_ratio) |> 
                               drop_na(any_of(pred_cols_tar[[set_map]]), target) %>%
                               select(patient_id, target, any_of(pred_cols_tar[[set_map]])),
                             train_cohort = cohort_map, panel = panel_map, fluid_name = fluid_map, pred_set = set_map, pred_cols = pred_cols_tar[[set_map]], label = label)),
  tar_target(model, model_map()),
  tar_target(wf, defineWorkflow(rec = rec, model = model[[1]])),
  tar_target(wf_regression_tuned, tuneRegWorkflow(wf, label = label)),
  tar_target(wf_regression_fit, fitAll(wf_regression_tuned)),
  #  tar_target(wf_regression_explain, explainModel(wf_regression_fit)),
  tar_target(wf_regression_vetiver, saveVetiverWorkflow(wf_regression_fit, model_board), deployment = "main"),
  
  tar_target(wf_regression_prospective, bind_cols(prospective |> drop_na(matches(paste0(lab_strings, collapse = "|"))), predict(bundle::unbundle(wf_regression_fit), prospective |> drop_na(matches(paste0(lab_strings, collapse = "|")))))),
  
  tar_target(wf_regression_yale_reviewed, bind_cols(yale_prospective |> drop_na(matches(paste0(lab_strings, collapse = "|"))), predict(bundle::unbundle(wf_regression_fit), yale_prospective |> drop_na(matches(paste0(lab_strings, collapse = "|"))))))
  
)

#### Run Pipeline ####

list(load_inputs, 
     build_binary_comment_predictors, 
     build_workflows 
     #build_regression_workflows
     )
