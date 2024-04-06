defineBaseRecipe <- function(train_input, train_cohort = "", train_type = "", panel = "", fluid_name = "", pred_set = "", pred_cols = "", recipe_name = "", label = "", ...){

  out_recipe <- 
    recipe(target ~ ., train_input) |> 
      step_normalize(recipes::all_numeric_predictors())
  
  out_recipe$train_cohort = train_cohort
  out_recipe$train_type = train_type
  out_recipe$panel = panel
  out_recipe$fluid_name = fluid_name
  out_recipe$pred_set = pred_set
  out_recipe$recipe_name = "base"
  out_recipe$label = label
  
  out_recipe
  
}

definePCARecipe <- function(train_input, train_cohort = "", train_type = "", panel = "", fluid_name = "", pred_set = "", pred_cols = "", recipe_name = "", label = "", ...){
  
  train = train_input %>% drop_na()
  
  pca_recipe <- 
    recipe(train) %>% 
    update_role(patient_id, new_role = "metadata") %>%
    update_role_requirements("metadata", bake = F) %>%
    update_role(!!!pred_cols, new_role = "predictor") %>%
    update_role(target, new_role = "outcome") %>%
    step_pca(recipes::all_predictors(), num_comp = 5, options = list(center = T, scale = T), keep_original_cols = T)
  
  pca_recipe$train_cohort = train_cohort
  pca_recipe$train_type = train_type
  pca_recipe$panel = panel
  pca_recipe$fluid_name = fluid_name
  pca_recipe$pred_set = pred_set
  pca_recipe$recipe_name = "pca"
  pca_recipe$label = label
  
  out = pca_recipe %>% butcher::axe_env()
  
  out
  
}

defineVarImpRecipe <- function(train_input, val_set = "", train_cohort = "", train_type = "", panel = "", fluid_name = "", pred_set = "", pred_cols = "", recipe_name = "", label = "", ...){
  
  train_input = train_input %>% drop_na()
  
  train_input = train_input %>% drop_na()
  split = group_initial_split(train_input, group = patient_id, prop = 0.9)
  train = training(split)
  validation_set = testing(split)
  
  rf <- rand_forest(mode = "classification") %>% set_engine("ranger", importance = "permutation")
  rf_fit <- fit(rf, target ~ ., train %>% select(-patient_id) %>% slice_sample(prop = 0.2))
  vi <- vip::vi(rf_fit)
  var_num = ifelse(nrow(vi)*0.35 < 5, 5, nrow(vi)*0.35)
  
  vi_cols <- vi[1:var_num, ][["Variable"]]
  
  vi_recipe <- 
    recipe(train) %>% 
    update_role(patient_id, new_role = "metadata") %>%
    update_role_requirements("metadata", bake = F) %>%
    update_role(!!!vi_cols, new_role = "predictor") %>%
    update_role(target, new_role = "outcome") %>%
    step_normalize(!!!vi_cols)
  
  vi_recipe$train_cohort = train_cohort
  vi_recipe$panel = panel
  vi_recipe$fluid_name = fluid_name
  vi_recipe$pred_set = pred_set
  vi_recipe$recipe_name = "vi"
  vi_recipe$label = label
  vi_recipe$validation_set <- validation_set
  
  out = vi_recipe %>% butcher::axe_env()
  
  out
  
}

defineEmbedRecipe <- function(train_input, val_set = "", train_cohort = "", train_type = "", panel = "", fluid_name = "", pred_set = "", pred_cols = "", recipe_name = "", label = "", ...){
  
  train_input = train_input %>% drop_na()
  
  train_input = train_input %>% drop_na()
  split = group_initial_split(train_input, group = patient_id, prop = 0.9)
  train = training(split)
  validation_set = testing(split)
  
  embed_recipe <- 
    recipe(train) %>% 
      update_role(patient_id, new_role = "metadata") %>%
      update_role_requirements("metadata", bake = F) %>%
      update_role(!!!pred_cols, new_role = "predictor") %>%
      update_role(target, new_role = "outcome") %>%
      step_normalize(!!!pred_cols) %>%
      embed::step_umap(!!!pred_cols, keep_original_cols = T, 
                       metric = "manhattan", num_comp = 3, neighbors = 50, min_dist = 0, 
                       options = list(fast_sgd = T, local_connectivity = 10, bandwidth = 10, 
                                      init = "agspectral", ret_model = T)) %>%
      step_pca(!!!pred_cols, num_comp = 3, keep_original_cols = T)
      
  
  embed_recipe$train_cohort = train_cohort
  embed_recipe$train_type = train_type
  embed_recipe$panel = panel
  embed_recipe$fluid_name = fluid_name
  embed_recipe$pred_set = pred_set
  embed_recipe$recipe_name = "embed"
  embed_recipe$label = label
  embed_recipe$validation_set <- validation_set
  
  embed_recipe
  
}

getRecipeCalls <- function(){
  
  rlang::syms(c("defineBaseRecipe"))
  
}

getWorkflowSetRecipes <- function(train_input, train_cohort = "", train_type = "", panel = "", fluid_name = ""){
  
  base_rec <- 
    recipe(train_input) %>%
      update_role(everything(), new_role = "metadata") %>%
      update_role_requirements("metadata", bake = F) %>%
      update_role(target, new_role = "outcome") %>%
      step_mutate(anion_gap = sodium - chloride - co2_totl,
                  osmolality = 2 * sodium + (bun / 2.8) + (glucose / 18))
  
  results_rec <- 
    base_rec %>% 
      update_role(any_of(!!lab_strings_bmp_no_gap), matches("anion|osmo"), new_role = "predictor") 

  priors_rec <- 
    base_rec %>% 
      update_role(matches(paste0(!!lab_strings_bmp_no_gap, collapse = "|")), -matches("post|target"), new_role = "predictor") %>%
      step_mutate(anion_gap_prior = sodium_prior - chloride_prior - co2_totl_prior,
                  osmolality_prior = 2 * sodium_prior + (bun_prior / 2.8) + (glucose_prior / 18)) %>%
      step_mutate(anion_gap_delta = anion_gap - anion_gap_prior,
                  osmolality_delta = osmolality - osmolality_prior)

  all_rec <- 
    base_rec %>% 
      update_role(matches(paste(!!lab_strings_bmp_no_gap, collapse = "|")), new_role = "predictor") %>%
      step_mutate(anion_gap_prior = sodium_prior - chloride_prior - co2_totl_prior,
                  osmolality_prior = 2 * sodium_prior + (bun_prior / 2.8) + (glucose_prior / 18),
                  anion_gap_post = sodium_post - chloride_post - co2_totl_post,
                  osmolality_post = 2 * sodium_post + (bun_post / 2.8) + (glucose_post / 18)) %>%
      step_mutate(anion_gap_delta = anion_gap - anion_gap_prior,
                  osmolality_delta = osmolality - osmolality_prior,
                  anion_gap_delta = anion_gap - anion_gap_post,
                  osmolality_delta = osmolality - osmolality_post)
  
  list(results = results_rec, priors = priors_rec, all = all_rec) %>% map(., ~butcher::butcher(.x, verbose = T))
  
}

