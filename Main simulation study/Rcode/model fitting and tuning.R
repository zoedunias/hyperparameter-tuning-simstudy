## Functions performance measures 

# Calculation C-statistic
C_stat <- function(predicts, outcome){
  
  predicts <- as.matrix(predicts)
  
  # Obtain number of categories of the outcome
  cats <- sort(unique(outcome))
  n_cats <- length(cats)
  
  n0 <- sum(outcome == cats[2])
  n1 <- length(outcome) - n0
  
  r <- rank(predicts[,2])
  prob_pos <- as.numeric(r[outcome == cats[2]])

  S0 <- sum(prob_pos)
  c_stat <- (S0 - n0 * (n0 + 1)/2) / (as.numeric(n0) * as.numeric(n1))

  CI <- ci.auc(as.factor(outcome), predicts[,2], method = "delong")
  
  if (CHECK_PM){
    cat("Comparison CStat: ", CI, " and ", c_stat, "\n")
  }
  return(c(c_stat, CI[1], CI[3]))
}

# Calculation square root of mean squared prediction error (rMSPE)
rMSPE <- function(predicts, prob){
  
  rMSPE <- sqrt(mean((prob - predicts[,2])^2))
  SE <- sd((predicts[,2] - prob)^2) / sqrt(length(prob))

  if (CHECK_PM){
    cat("Comparison rMSPE: ", RMSE(pred = predicts[,2], obs = prob), " and ", rMSPE, "\n")
  }
  
  return(c(rMSPE, SE))
}

# Calculation mean absolute prediction error (MAPE)  
MAPE <- function(predicts, prob){
  
  MAPE <- mean(abs(prob - predicts[,2]))
  
  SE <- sd(abs(predicts[,2] - prob)) / sqrt(length(prob))

  if (CHECK_PM){
    cat("Comparison MAPE: ", MAE(prob, predicts[,2]), " and ", MAPE, "\n")
  }
  
  return(c(MAPE, SE))
}

# Calculation calibration intercept and calibration slope
calibration <- function(lp, outcome){
  lp <- as.matrix(lp[,2])
  idx <- (lp == 0)
  lp[idx] <- lp[idx] + 1e-12
  idx <- (lp == 1)
  lp[idx] <- lp[idx] - 1e-12
  
  logitp <- log(lp/(1-lp))
  fit <- glm(outcome ~ logitp, family = "binomial")
  
  slope <- coef(fit)[2]
  if (nrow(summary(fit)$coefficients) < 2) {
    SE <- NA
  } else {
    SE <- summary(fit)$coefficients[2,2]
  }

  if (CHECK_PM){
    cat("Comparison CS: ", slope, "\n")
    valprob <- val.prob(model_predict[,2], test_y_input)
    print(valprob)
  }
  
  return(c(slope, SE))
}

# Calculation calibration in the large (CIL)
CIL <- function(predicts, outcome){
  CIL <- mean(predicts[,2] - outcome)

  SE <- sd(predicts[,2] - outcome) / sqrt(length(outcome))

  return(c(CIL, SE))
}

LL <- function(predicts, outcome){
  idx <- (predicts == 0)
  predicts[idx] <- predicts[idx] + 1e-12
  idx <- (predicts == 1)
  predicts[idx] <- predicts[idx] - 1e-12
  
  sum(outcome*log(predicts)+(1-outcome)*log(1-predicts))
}

# Calculation Cox-Snell, Nagelkerke and Mcfadden Rsquared
pseudo_Rsqrs <- function(predicts, outcome){ 
  LL_fit  <- LL(predicts[,2], outcome) 
  LL_null <- LL(predicts = mean(outcome), outcome)   
  cox <- 1-exp(-(LL_fit-LL_null)*2/nrow(outcome)) 
  cox_max <- 1 - exp(2 * nrow(outcome) ^ (-1) * LL_null)
  c("cox"=cox,"nagelkerke"=cox/cox_max,"mcfadden"=1-LL_fit/LL_null)
}

deviance_optimization <- function(data, lev = NULL, model = NULL){
  outcome <- ifelse(data$obs == lev[1], 1, 0)
  predicts <- data[, lev[1]]

  idx <- (predicts == 0)
  predicts[idx] <- predicts[idx] + 1e-12
  idx <- (predicts == 1)
  predicts[idx] <- predicts[idx] - 1e-12
  #cat("Ones: " , sum((predicts==1)))
  #cat("Zeros: " , sum((predicts==0)))
  c(deviance = -2*sum(outcome*log(predicts) + 
                        ((1-outcome)*log(1-predicts))),
    twoClassSummary(data, lev = lev))
}


## Functions model fitting and tuning 

# Prepare tuneGrid
grid <- function(model_label, train_x, lambda_input) {
  
  # grid parameter space
  # lambda_input = seq(0, 2, length=201)
  alpha_input = seq(0, 1, length= 11)
  
  if (model_label == "standard") {
    tuneGrid = expand.grid(alpha = 1, 
                           lambda = lambda_input)
    model = "glm"
    } 
  if (model_label == "ridge") { 
    tuneGrid = expand.grid(alpha = 0, 
                           lambda = lambda_input)
    model = "glmnet"
  }
  
  if (model_label == "lasso") { 
    tuneGrid = expand.grid(alpha = 1, 
                           lambda = lambda_input)
    model = "glmnet"
  }
  
  if (model_label == "elastic net") {
    tuneGrid = expand.grid(alpha = alpha_input, 
                           lambda = lambda_input)
    model = "glmnet"
  }
  
  if (model_label == "random forest") {
    tuneGrid = expand.grid(mtry = array(1:ncol(train_x)),
                           splitrule = "extratrees",
                           #splitrule = c("extratrees", "gini"),
                           min.node.size = array(1:10))
    model = "ranger"
  }
  
  if (model_label == "svmRadial") {
    tuneGrid = expand.grid(sigma = 2^seq(-15, 3, length = 10), 
                           C = 2^seq(-5, 15, length = 11))
    model = "svmRadialSigma"
  }
  
  output <- list("tuneGrid" = tuneGrid, "model" = model)
  return(output)
}

# Tuning, fit model and predict
fit_model <- function(train_x_input, train_y_input, test_x_input, test_y_input, prob, model_input, method_input, number_input, repeats_inputs, selection_input, grid, model_label){

  train_y_factor_input <- factor(train_y_input, labels = c('benign','malignant'))
  # Tuning
  train_control <- trainControl(method = method_input,
                                number = number_input,
                                search = "grid",
                                #summaryFunction = twoClassSummary,
                                summaryFunction = deviance_optimization,
                                classProbs = TRUE,
                                verboseIter = FALSE,
                                selectionFunction = selection_input)
  if (method_input == "repeatedcv"){
    repeats_input <- repeats_inputs[number_input/5]
    
    train_control <- trainControl(method = method_input,
                                  number = number_input,
                                  repeats = repeats_input,
                                  search = "grid",
                                  #summaryFunction = twoClassSummary,
                                  summaryFunction = deviance_optimization,
                                  classProbs = TRUE,
                                  verboseIter = FALSE,
                                  selectionFunction = selection_input)
  }
  
  
  if (model_input == "ranger"){
    # Train the model
    model_fit <- train(x = train_x_input, y = train_y_factor_input,   
                       method = model_input,
                       tuneGrid = grid,
                       #metric = "ROC",
                       metric = "deviance",
                       maximize = 0,
                       #num.trees = 128,
                       num.trees = 500,
                       trControl = train_control)
    # Predict using the fitted model
    model_predict <- predict.train(model_fit, test_x_input, type = 'prob')
    
  }else if (model_label == "standard"){
    dat_train <- as.data.frame(cbind(train_y_input, train_x_input))
    dat_test <- as.data.frame(cbind(test_y_input, test_x_input))
    model_fit <- glm(formula = train_y_input ~ ., family = binomial, data = dat_train)
    coefs <- summary(model_fit)$coefficients
    SEs <- coefs[,2]
    if (sum(as.numeric(SEs > 70))) {
      cat("MYWARNING: SEPARATION (SE of COEFS LARGER THAN 70), SE coefs: ", SEs, "\n")      
    }
    
    parestimates <- coefs[-1,1]
    if (sum(parestimates) < 1e-6) {
      cat("MYWARNING: (possible) degenerate linear predictors (COEFS SMALLER THAN 1e-6: CHECK WHETHER 0), estimates coefs: ", parestimates, "\n")      
    }
    
    model_predict <- predict(model_fit, dat_test[-1], type = 'response')
    model_predict <- cbind(model_predict, new_col = model_predict)
    model_fit$bestTune <- list(NA, NA)
    
  }else {
    
      # Train the model
      model_fit <- train(x = train_x_input, y = train_y_factor_input,   
                         method = model_input,
                         tuneGrid = grid,
                         #metric = "ROC",
                         metric = "deviance",
                         maximize = 0,
                         trControl = train_control)
      
      parestimates <- coef(model_fit$finalModel, (model_fit$finalModel)$lambdaOpt)
      parestimates <- parestimates[-1]
      if (sum(parestimates) < 1e-6) {
        cat("MYWARNING: (possible) degenerate linear predictors (COEFS SMALLER THAN 1e-6: CHECK WHETHER 0), estimates coefs: ", parestimates, "\n")      
      }
      
      # Predict using the fitted model
      model_predict <- predict.train(model_fit, test_x_input, type = 'prob')
  }
  
  # Calculate performance measures
  C_stat <- C_stat(model_predict, test_y_input) # C-statistic
  #valprob <- val.prob(model_predict[,2], test_y_input)
  #print(valprob)
  calib <- calibration(model_predict, test_y_input) # Calibration intercept and slope
  #valprob <- 0
  CIL <- CIL(model_predict, test_y_input) # Calibration in the large
  rMSPE <- rMSPE(model_predict, prob) # Root mean squared prediction error
  #caret <- RMSE(pred = model_predict[,2], obs = prob)
  MAPE <- MAPE(model_predict, prob) # Mean absolute prediction error
  #inbuiltmae <- MAE(prob, model_predict[,2])
  #cox <- pseudo_Rsqrs(model_predict, test_y_input)
  
  outputtuning <- model_fit$bestTune
  if (length(outputtuning)==2) {
    outputtuning[[3]] <- NA
  }
  output <- c(C_stat, calib, CIL, rMSPE, MAPE, outputtuning[[1]], outputtuning[[2]], outputtuning[[3]], model_label, model_input, method_input, number_input, selection_input, paste(method_input, number_input, selection_input, sep = "_"))
  
  #return(output)
  return(list(output = output, outputtuning = outputtuning))
}

modeling_tuning <- function(train_x, train_y, test_x, test_y, prob, model_labels, method_inputs, number_inputs, selection_inputs){
  
  combination_index <- 1

  ## Bootstrap
  number_input_boot = 500
  selection_input_boot = "best"
  
  ## LOOCV
  method_input_LOOCV = "LOOCV"
  
  besttunes <- list()
  
  cat(c("model label = standard", "\n", "model = glm", "\n", "tuning method = n/a", "\n", "number = n/a", "\n", "selection_input = n/a", "\n", "\n"))
  if (sum(train_y) < 8){
    results <- c(rep(NA, times = nr_of_perf_measures), "standard", "glm", "n/a", "n/a", "n/a", "n/a")
    cat("MYWARNING: degenerate outcome distributions (ONES OCCURS LESS THAN 8 TIMES IN OUTCOME), number of ones: ", sum(train_y), "\n")
    cat("Results are: ", results, "\n")
  } else {
    r <- fit_model(train_x_input = train_x, train_y_input = train_y, test_x_input = test_x, test_y_input = test_y, prob = prob, model_input = "glm", method_input = "", number_input = 0, repeats_inputs = repeats_inputs, selection_input = "", grid = 0, model_label = "standard")
    results <- r$output
  }
  combination_index <- combination_index + 1
  
  lambda_input <- seq(0, 2, length=201)
  
  for (model_label in model_labels){
    model_vars <- grid(model_label, train_x, lambda_input = lambda_input)
    grid <- model_vars$tuneGrid
    model_input <- model_vars$model
    
    
    for (method_input in method_inputs){
      if ((method_input == "boot") | (method_input == "LOOCV")){
        
        cat(c("model label =", model_label, "\n", "model =", model_input, "\n", "tuning method =", method_input, "\n", "number =", number_input_boot, "\n", "selection_input = best", "\n", "\n"))
        if (sum(train_y) < 8){
          r <- c(rep(NA, times = nr_of_perf_measures), model_label, model_input, method_input, number_input, selection_input, paste(method_input, number_input, selection_input, sep = "_"))
          results <- rbind(results, r)
          besttunes[[combination_index]] <- NULL 
          cat("MYWARNING: degenerate outcome distributions (ONES OCCURS LESS THAN 8 TIMES IN OUTCOME), number of ones: ", sum(train_y), "\n")
          cat("Results are: ", r, "\n")
          next  
        }
        r <- fit_model(train_x, train_y, test_x, test_y, prob, model_input, method_input, number_input_boot, repeats_inputs, "best", grid, model_label)
        results <- rbind(results, r$output)
        besttunes[[combination_index]] <- r$outputtuning 
        
        combination_index <- combination_index + 1
      }
      
      if (method_input == "cv" | method_input == "repeatedcv"){ 
        for (number_input in number_inputs){
          for (selection_input in selection_inputs){  
            cat(c("model label =", model_label, "\n", "model =", model_input, "\n", "tuning method =", method_input, "\n", "number =", number_input, "\n", "selection_input", selection_input, "\n", "\n"))
            if (sum(train_y) < 8){
              r <- c(rep(NA, times = nr_of_perf_measures), model_label, model_input, method_input, number_input, selection_input, paste(method_input, number_input, selection_input, sep = "_"))
              results <- rbind(results, r)
              besttunes[[combination_index]] <- NULL 
              cat("MYWARNING: degenerate outcome distributions (ONES OCCURS LESS THAN 8 TIMES IN OUTCOME), number of ones: ", sum(train_y), "\n")
              cat("Results are: ", r, "\n")
              next  
            }
            
            r <- fit_model(train_x, train_y, test_x, test_y, prob, model_input, method_input, number_input, repeats_inputs, selection_input, grid, model_label)
            results <- rbind(results, r$output)
            besttunes[[combination_index]] <- r$outputtuning 
            
            combination_index <- combination_index + 1
          }
        }
      }
    }
  }
  colnames(results) <- c("C-statistic", "C-statistic CI_lower", "C-statistic CI_upper", "CS", "CS SE","CIL", "CIL SE", "rMSPE", "rMSPE SE", "MAPE", "MAPE SE", "hyperparameter 1","hyperparameter 2", "hyperparameter 3", "model label", "package", "tuning method", "number folds", "selection criterium", "combination label")
  return(as.data.frame(results))
  #return(list('performance measures' = results, 'best tuning parameter' = besttunes))
}

plot_function <- function(data){
  
  data_df <- data.frame(data)
  data_df$CS <- as.numeric(data_df$CS)
  #data_df$CS <- round(data_df$CS,3)
  print(data_df$CS)
  
  ggplot(data_df, mapping = aes(x = combination.label, y = CS, group = model.label))+
    geom_line()+
    #geom_point(aes(shape = model.label), size = 2)+
    geom_point(shape = 2, size = 2)+
    facet_wrap(vars(model.label))+
    #ylim(0.5,1.5)+
    ylim(0.5,2.0)+
    theme_bw(base_size = 13)+
    #theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust = 0.5))+
    labs(x = "Tuning procedures", y = "Calibration slope", title = 'Predictive performance of tuning procedures for 4 models')+
    
    scale_x_discrete(labels = c('Bootstrap','10-fold CV','1SE 10-fold CV','5-fold CV','1SE 5-fold CV','LOOCV'))
  
    #theme_light()
    #scale_x_discrete(, expand=c(0.05, 0.05))
    #scale_y_continuous(breaks=seq(0,2,0.2))
    
}
#pm <- test$'performance measures'
#plot_function(data = pm)










