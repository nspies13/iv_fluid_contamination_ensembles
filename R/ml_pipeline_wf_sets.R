
defineWorkflowSets <- function(train){
  
  recipes = getWorkflowSetRecipes(train, train_cohort = "", train_type = "", panel = "", fluid_name = "")
  models = c(defineXGB())
  
  wf_set = workflow_set(recipes, models, cross = T) %>% butcher::butcher(verbose = T)
  
  wf_set

}

finalizeParameters <- function(wf_set = workflow_set){
  
  params = map(wf_set$wflow_id, 
               ~ wf_set %>% extract_workflow(.x) %>% extract_parameter_set_dials() %>% finalize(wf_set %>% extract_workflow(.x) %>% extract_preprocessor() %>% pluck("template")))
  
  library(foreach)
  tmp = foreach(i = 1:length(wf_set$wflow_id)) %do% {
    if(params[i] %>% map(1) %>% pluck(1) %>% length() > 0){
      wf = wf_set %>% extract_workflow(wf_set$wflow_id[i])
      wf_set = wf_set %>% option_add(id = wf_set$wflow_id[i],
                                     grid = ifelse(length(params %>% pluck(1, 1)) > 3, 10, 5),
                                     param_info = wf %>% extract_parameter_set_dials() %>% finalize(wf %>% extract_preprocessor() %>% pluck("template")))
    }
  }
  
  wf_set = tmp[[length(tmp)]]
  
  wf_set
  
}

tuneWorkflowSet <- function(wf_set = workflow_set, cv = cv){
  
  library(doParallel)
  library(parallelly)
  cl <- parallelly::makeClusterPSOCK(32)
  doParallel::registerDoParallel(cl)
  
  wf_set_tuned <- 
    wf_set %>% 
      workflow_map(
        fn = "tune_bayes", 
        iter = 20,
        initial = 5,
        resamples = cv,
        metrics = metrics_global,
        control = control_bayes(verbose = F, allow_par = T, save_pred = F, save_workflow = T, event_level = "second", seed = 12345, uncertain = 3, no_improve = 7, parallel_over = "everything")
      )
  
  wf_set_tuned %>% option_remove(resamples)
  
}

makeFitWorkflowList <- function(wf_set, wf_set_tuned, train = tar_read(train_comments)){
  
  params = 
    map(wf_set$wflow_id, ~wf_set_tuned %>% 
          extract_workflow_set_result(.x) %>% 
          select_best("mn_log_loss"))
  
  wf_fit = map2(wf_set$wflow_id, params, 
                ~extract_workflow(wf_set, .x) %>% 
                  finalize_workflow(.y) %>% 
                  fit(train))
  
  wf_fit %>% 
    set_names(wf_set$wflow_id) %>%
    map(~bundle::bundle(butcher::butcher(.x)))
  
}

predictAll <- function(wf_list = tar_read(wf_fit_list) %>% map(bundle::unbundle), new_data = tar_read(prospective)){
  
  preds <- map2(wf_list, names(wf_list), ~predict(.x, new_data) %>% setNames(paste0(.y, "_pred_class"))) %>% bind_cols()
  
  real_probs <- map2(wf_list, names(wf_list), ~predict(.x, new_data, type = "prob") %>% select(.pred_Physiologic) %>% setNames(paste0(.y, "_real_prob"))) %>% bind_cols()
  real_probs$mean_real_prob <- rowMeans(real_probs)
  real_probs$max_real_prob <- apply(real_probs, 1, function(x) max(x, na.rm = T))
  real_probs$min_real_prob <- apply(real_probs, 1, function(x) min(x, na.rm = T))
  
  bind_cols(new_data, preds, real_probs)
  
}

makeStack <- function(wf_set_tuned){
  
  library(stacks)
  wf_stack <- 
    stacks() %>% 
      add_candidates(wf_set_tuned) 
  
  wf_blend <- 
    wf_stack %>%
      blend_predictions(penalty = 0, mixture = 0, times = 1) 
  
  wf_stack_fit <-
    wf_blend %>% 
      fit_members()
  
  wf_stack_fit %>% 
    butcher::butcher() %>% 
    bundle::bundle()
  
}

saveVetiverStack <- function(wf_stack = wf_stack, proto = proto, metadata = metadata, model_board = model_board, model_name = "Stack", ...){
  
  library(pins)
  library(vetiver)
  
  v <- 
    vetiver::vetiver_model(
      model = wf_stack,
      model_name = model_name,
      save_prototype = proto,
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
