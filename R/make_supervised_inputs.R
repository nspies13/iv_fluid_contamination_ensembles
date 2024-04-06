simulateContaminationRow <- function(input, mix_ratio, fluid){
  
  cols <- names(fluid)[which(names(fluid) %in% names(input))]
  
  output <- input %>%
    dplyr::mutate(across(all_of(cols), ~(1 - mix_ratio) * . + fluid[cur_column()] * mix_ratio)) %>%
    select(all_of(cols))
  
  output %>% 
    mutate(across(c("sodium", "chloride", "co2_totl", "bun", "glucose"), ~round(.))) %>%
    mutate(across(c("potassium_plas", "calcium"), ~round(., 1))) %>%
    mutate(creatinine = round(creatinine, 2)) %>% 
    mutate(anion_gap = sodium - chloride - co2_totl) %>%
    mutate(mix_ratio = mix_ratio)
  
}

makeMixRatios <- function(contam_rate = 0.25, contam_train_input, fluid_names_tar){
  
  minimum_significant_contamination = ifelse(grepl("D5", fluid_names_tar), 0.05, 0.1)
  
  out = round(rbeta(contam_rate * nrow(contam_train_input), 1, 5), 2) + minimum_significant_contamination
  out[out >= 1] <- 0.99
  
  out
  
}

makeBinarySimTrainSeed <- function(input = tar_read(preprocessed_bmp_inputs)[[1]][["data"]], train = tar_read(train)){
  
  train_all = train %>% map("train") %>% map(., ~filter(., label != "Patient")) %>% bind_rows() %>% distinct() %>% select(-matches("real"))
  
  input = input %>% filter(!code_comment) %>% mutate(target = factor(ifelse(unlikely_comment | contam_comment, 1, 0), levels = c(0, 1)))

  train_common = train_all %>% filter(label %in% c("NS", "SW", "LR")) %>% slice_sample(prop = 0.25)
  train_uncommon = train_all %>% filter(!label %in% c("NS", "SW", "LR")) %>% slice_sample(prop = 0.05)
  
  output = bind_rows(input, train_common, train_uncommon)
  
  arrow::write_feather(output, "../../Data/tabnet_seed.feather")
  
  output
  
}

makeSimulatedTrainingData <- function(input, mix_ratios, fluid, fluid_name, predictors){
  
  input = input %>% mutate(across(any_of(lab_strings), ~.x, .names = "{col}_real"))
  split = initial_split(input, prop = length(mix_ratios)/nrow(input))
  
  unsim_rows <- testing(split)
  sim_rows <- training(split)
  
  sim_rows = sim_rows %>% mutate(mix_ratio = mix_ratios, label = fluid_name)
  unsim_rows = unsim_rows %>% mutate(mix_ratio = 0, label = "Patient")
  
  tmp = simulateContaminationRow(sim_rows, sim_rows$mix_ratio, fluid)
  
  sim_rows[,names(tmp)] <- tmp 
  
  sim_rows= 
    sim_rows %>% 
      mutate(across(any_of(lab_strings), ~. - get(paste0(cur_column(), "_prior")), .names = "{col}_delta_prior")) %>%
      mutate(across(any_of(lab_strings), ~get(paste0(cur_column(), "_post")) - ., .names = "{col}_delta_post"))
  
  output = bind_rows(unsim_rows, sim_rows) %>% mutate(target = factor(ifelse(label == "Patient", 0, 1), levels = c("0", "1"))) %>% select(patient_id, mix_ratio, label, target, matches(paste(lab_strings, collapse = "|")))
  
  output
  
}

compareAbsoluteToProportionalDeltas <- function(){
  
  train = tar_read(train_210c6c17)
  
  pred_cols = train %>% select(matches(paste(lab_strings, collapse = "|")), target) %>% names()
  
  rf <- rand_forest(mode = "classification") %>% set_engine("ranger", importance = "permutation")
  rf_fit <- fit(rf, target ~ ., train %>% select(all_of(pred_cols), target) %>% drop_na())
  vi <- vip::vi(rf_fit)
  pred_cols <- vi[1:5, "Variable"][[1]]
  
  write_delim(vi, "../../Results/Models/VarImps/absolute_vs_proportional_delta_varImps_comparison.txt", delim = '\t')
  
}

makeAnomalyWithResolutionTrainSet <- function(input = bmp_no_comment[[1]], fluid_name = "NS", AwR_negative = 1, AwR_positive = 4){
  
  counts <- 
    input %>% 
      transmute(AwR_count = getAnomalyResolution(., fluid_name = fluid_name, mix_ratio_threshold = -0.05)) %>% 
      rowSums()
  
  output <- 
    input %>% 
      mutate(counts = counts) %>%
      filter(counts <= AwR_negative | counts >= AwR_positive) %>%
      mutate(target = factor(ifelse(counts > AwR_negative, "1", "0"), levels = c("0", "1")),
             label = factor(ifelse(target == 1, fluid_name, "Real"), levels = c("Real", fluid_name)))
  
  output
    
}

plotTrainingBoxplots <- function(train, pred_cols = lab_strings_bmp_no_gap, label_col = "label"){
  
  gg_long <- 
    train %>% 
      select(all_of(pred_cols), all_of(label_col)) %>%
      pivot_longer(pred_cols, names_to = "feature", values_to = "value")
  
  gg_in = 
     gg_long %>%
      group_by(feature, label) %>%
      summarise(p05 = quantile(value, probs = c(0.05)),
                p25 = quantile(value, probs = c(0.25)),
                p50 = quantile(value, probs = c(0.50)),
                p75 = quantile(value, probs = c(0.75)),
                p95 = quantile(value, probs = c(0.95))) %>%
      mutate(axis = as.numeric(label), y = max(p95), label_y = min(p05))
 
  library(ggpubr)
  library(rstatix) 
  means <- 
    map(unique(gg_long$feature), 
      ~tidy(t.test(data = gg_long %>% filter(feature == .x), value ~ label))) %>%
        bind_rows() %>%
        mutate(
          feature = unique(gg_long$feature),
          ci_label = paste0("95% CI: [", round(conf.low, digits = 1), ", ", round(conf.high, digits = 1), "]"))
  
  gg_input = left_join(gg_in, means)
  
  ggplot(data = gg_input, aes(fill = label)) +
    geom_segment(aes(x = axis, y = p05, xend = axis, yend = p95)) +
    geom_rect(aes(xmin = axis - 0.25, ymin = p25, xmax = axis + 0.25, ymax = p75)) +
    geom_segment(aes(x = axis - 0.25, y = p50, xend = axis + 0.25, yend = p50), color = "grey75", linewidth = 1.25) +
    geom_text(aes(label = ci_label, y = label_y), x = 1.5, hjust = 0) +
    scale_x_continuous(limits = c(0.5, 2.5), breaks = c(1,2), labels = c(levels(gg_input$label))) +
    scale_fill_manual(values = c("grey40", "darkred")) +
    facet_wrap(~feature, scales = "free", nrow = 2, labeller = as_labeller(analyte_labels)) +
    coord_flip() + ylab("Result") + xlab("Label") + 
    theme(legend.position = "none", axis.title.y.left = element_blank())
  ggsave(paste0("../../Figures/Models/Pipeline/train_AwR_boxplots_", fluid_name, ".pdf"), width = 8.5, height = 6)
  
}
