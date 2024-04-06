defineXGB <- function(){
  
  list(
    xgb = 
      boost_tree(mode = "classification", trees = tune(), learn_rate = tune(), tree_depth = tune(), loss_reduction = tune()) %>% 
        set_engine("xgboost"))
  
  
}

defineRegXGB <- function(){
  
  list(
    xgb = 
      boost_tree(mode = "regression", trees = tune(), learn_rate = tune(), tree_depth = tune(), loss_reduction = tune()) %>% 
      set_engine("xgboost"))
  
  
}

defineLightGBM <- function(){
  
  library(bonsai)
  library(lightgbm)
  list(
    lgb = boost_tree(mode = "classification", mtry = 5, 
                     trees = tune(), learn_rate = tune(), tree_depth = tune(), 
                     min_n = tune(), loss_reduction = tune()) %>% 
      set_engine("lightgbm")
  )
  
}


defineCART <- function(){
  
  list(
    cart = 
      decision_tree(mode = "classification", 
                    tree_depth = 3, cost_complexity = tune()) %>% 
      set_engine("rpart")
    
  )
  
}

defineLogReg <- function(){
  
  list(
    lr =
      logistic_reg(mode = "classification", 
                   penalty = tune(), mixture = tune()) %>% 
      set_engine("glmnet")
  )
  
}

defineRegGLM <- function(){
  
  list(
    lr =
      linear_reg(mode = "regression", 
                   penalty = tune(), mixture = tune()) %>% 
      set_engine("glmnet")
  )
  
}

defineSVM <- function(){
  
  list(
    svm = 
      svm_linear(mode = "classification", 
                 cost = tune(), margin = tune()) %>% 
      set_engine("kernlab")
  )
  
}

defineRF <- function(){
  
  list(
    rf = 
      rand_forest(mode = "classification", 
                  mtry = tune(), trees = tune()) %>% 
      set_engine("ranger")
  )
  
}

defineNNet <- function(){
  
  list(
    nnet = 
      mlp(mode = "classification", 
          hidden_units = tune(), epochs = tune()) %>% 
      set_engine("nnet")
  )
  
}

defineRegNNet <- function(){
  
  list(
    nnet = 
      mlp(mode = "regression", 
          hidden_units = tune(), epochs = tune()) %>% 
      set_engine("nnet")
  )
  
}

getModelCalls <- function(){
  
  rlang::syms(c("defineXGB", "defineNNet"))
  
}

getRegModelCalls <- function(){
  
  rlang::syms(c("defineRegNNet", "defineRegGLM"))
  
}

defineAllModels <- function(){
  
  xgb <- boost_tree(mode = "classification", mtry = tune(), trees = tune(), learn_rate = tune(), tree_depth = tune(), min_n = tune(), loss_reduction = tune(), sample_size = tune(), stop_iter = tune()) %>% set_engine("xgboost")
  lgb <- boost_tree(mode = "classification", mtry = tune(), trees = tune(), learn_rate = tune(), tree_depth = tune(), min_n = tune(), loss_reduction = tune()) %>% set_engine("lightgbm")
  cart <- decision_tree(mode = "classification", tree_depth = 3, cost_complexity = tune()) %>% set_engine("rpart")
  lr <- logistic_reg(mode = "classification", penalty = tune(), mixture = tune()) %>% set_engine("glmnet")
  svm <- svm_linear(mode = "classification", cost = tune(), margin = tune()) %>% set_engine("kernlab")
  rf <- rand_forest(mode = "classification", mtry = tune(), trees = tune()) %>% set_engine("ranger")
  nnet <- mlp(mode = "classification", hidden_units = tune(), penalty = tune(), epochs = tune()) %>% set_engine("nnet")
  
  list(xgb = xgb, nnet = nnet)
  
}
