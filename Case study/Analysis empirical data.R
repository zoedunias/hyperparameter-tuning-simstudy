## Load required packages ##
library(glmnet)
library(rms)
library(caret)
library(ranger)
library(dplyr)
library(kernlab)
library(pROC)
#library(pmsampsize)

## Define constants ##
HPC <- TRUE

# Tune and modeling
model_labels <- c("ridge", "lasso", "elastic net", "random forest")
method_inputs <- c("cv", "boot") 
number_inputs <- c(5, 10)
selection_inputs <- c("best", "oneSE")

if (!HPC){
  # TEST SETTINGS
  model_labels <- c("random forest")
  method_inputs <- c("boot") 
  number_inputs <- c(10)
  selection_inputs <- c("best", "oneSE")
}

## Functions performance measures (and sample size calculation) ##

# Approximate R2
approximate_R2 <- function(auc, prev, n = 1000000){
  # define mu as a function of the C statistic
  mu <- sqrt(2) * qnorm(auc)
  # simulate large sample linear prediction based on two normals
  # for non-eventsN(0, 1), events and N(mu, 1)
  LP <- c(rnorm(prev*n, mean=0, sd=1), rnorm((1-prev)*n, mean=mu, sd=1))
  y <- c(rep(0, prev*n), rep(1, (1-prev)*n))
  # Fit a logistic regression with LP as covariate;
  # this is essentially a calibration model, and the intercept and
  # slope estimate will ensure the outcome proportion is accounted
  # for, without changing C statistic
  fit <- lrm(y~LP)
  max_R2 <- function(prev){
    1-(prev^prev*(1-prev)^(1-prev))^2
  }
  return(list(R2.nagelkerke = as.numeric(fit$stats['R2']),
              R2.coxsnell = as.numeric(fit$stats['R2']) * max_R2(prev)))
}

# Calculation C-statistic
C_stat <- function(predicts, outcome){
  
  predicts <- as.matrix(predicts)
  
  # Obtain number of categories of the outcome
  cats <- sort(unique(outcome))
  n_cats <- length(cats)
  
  # Number of category ...
  n0 <- sum(outcome == cats[2])
  # Number of category ...
  n1 <- length(outcome) - n0
  
  # Rank probabilities of ... 
  r <- rank(predicts[,2])
  # Area
  prob_pos <- as.numeric(r[outcome == cats[2]])
  #plot(prob_pos)
  
  S0 <- sum(prob_pos)
  #print(S0)
  #print(n0)
  #print(n1)
  c_stat <- (S0 - n0 * (n0 + 1)/2) / (as.numeric(n0) * as.numeric(n1))
  
  SE <- sd(prob_pos/ (as.numeric(n0) * as.numeric(n1)) / sqrt(length(prob_pos)) )  
  
  CI_lower <- c_stat - 1.96*SE
  CI_upper <- c_stat + 1.96*SE
  
  CI <- ci.auc(as.factor(outcome), predicts[,2], boot.n = 500)
  #return(c(c_stat, CI_lower, CI_upper))
  return(c(c_stat, CI[1], CI[3]))
}

# Calculation square root of mean squared prediction error (rMSPE)
rMSPE <- function(predicts, outcome){
  rMSPE <- sqrt(mean((outcome - predicts[,2])^2))
  return(rMSPE)
}

# Calculation mean absolute prediction error (MAPE)  
MAPE <- function(predicts, outcome){
  MAPE <- mean(abs(outcome - predicts[,2]))
  return(MAPE)
}

# Calculation calibration intercept and calibration slope
calibration <- function(lp, outcome){
  lp <- as.matrix(lp[,2])
  
  idx <- (lp == 0)
  lp[idx] <- lp[idx] + 1e-12
  idx <- (lp == 1)
  lp[idx] <- lp[idx] - 1e-12
  
  #p <- predict.lrm(lp)
  logitp <- log(lp/(1-lp))
  fit <- glm(outcome ~ logitp, family = "binomial")
  
  slope <- coef(fit)[2]
  
  CI_lower <- slope - 1.96*summary(fit)$coefficients[2,2]
  CI_upper <- slope + 1.96*summary(fit)$coefficients[2,2]
  
  return(c(slope, CI_lower, CI_upper))
}

# Calculation calibration in the large (CIL)
CIL <- function(predicts, outcome){
  CIL <- mean(predicts[,2] - outcome)
  
  SE <- sd(predicts[,2] - outcome) / sqrt(nrow(outcome))
  CI_lower <- CIL - 1.96*SE
  CI_upper <- CIL + 1.96*SE
  
  return(c(CIL, CI_lower, CI_upper))
}

LL <- function(predicts, outcome){
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


## Functions model fitting and tuning 

# Prepare tuneGrid
grid <- function(model_label, train_x) {
  
  # grid parameter space
  lambda_input <- seq(0, 2, length=201)
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
                           min.node.size = array(1:10),
                           splitrule = c("extratrees"))
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
fit_model <- function(train_x_input, train_y_input, test_x_input, test_y_input, model_input, method_input, number_input, selection_input, grid, model_label){
  
  train_y_factor_input <- factor(train_y_input, labels = c('benign','malignant'))

  # Tuning
  train_control <- trainControl(method = method_input,
                                number = number_input,
                                search = "grid",
                                classProbs = TRUE,
                                verboseIter = FALSE,
                                selectionFunction = selection_input)
  
  
  if (model_input == "ranger"){
    # Train the model
    model_fit <- train(x = train_x_input, y = train_y_factor_input,   
                       method = model_input,
                       tuneGrid = grid,
                       #num.trees = 128,
                       num.trees = 500,
                       trControl = train_control)
    # Predict using the fitted model
    model_predict <- predict.train(model_fit, test_x_input, type = 'prob')
  }else if (model_label == "standard"){
    dat <- as.data.frame(cbind(train_y_input, train_x_input))
    colnames(dat) <- c('outcome_binary', 'Age', 'lesdmax', 'papnr')
    model_fit <- glm(formula = outcome_binary ~ Age + lesdmax + as.factor(papnr), family = binomial, data = dat)
    model_predict <- predict(model_fit, data.frame(test_x_input), type = 'response')
    model_predict <- cbind(model_predict, new_col = model_predict)

  }else {
    
    # Create dummy variables for categorical papnr variable
    train_x_factor_input <- train_x_input[,1:2]
    test_x_factor_input <- test_x_input[,1:2]
    for (i in 1:4){
      train_x_factor_input <- cbind(train_x_factor_input, as.integer(train_x_input[,3]==i))
      test_x_factor_input <- cbind(test_x_factor_input, as.integer(test_x_input[,3]==i))
    }
    
    # Train the model
    model_fit <- train(x = train_x_factor_input, y = train_y_factor_input,   
                       method = model_input,
                       tuneGrid = grid,
                       trControl = train_control)
    # Predict using the fitted model
    model_predict <- predict.train(model_fit, test_x_factor_input, type = 'prob')
    
    #print(coef(model_fit$finalModel,unlist(model_fit$bestTune))) # Check the parameter coefficients
  }
  
  
  #plot(test_y_input, model_predict[,2])
  
  # Calculate performance measures
  C_stat <- C_stat(model_predict, test_y_input) # C-statistic
  calib <- calibration(model_predict, test_y_input) # Calibration intercept and slope
  #valprob <- val.prob(model_predict[,2], test_y_input)
  #valprob
  #valprob <- 0
  CIL <- CIL(model_predict, test_y_input) # Calibration in the large
  #rMSPE <- rMSPE(model_predict, test_y_input) # Root mean squared prediction error
  #rmse <- sqrt(mean((test_y_input - model_predict[,2])^2))
  #caret <- RMSE(pred = model_predict[,2], obs = test_y_input)
  #MAPE <- MAPE(model_predict, test_y_input) # Mean absolute prediction error
  #inbuiltmae <- MAE(test_y_input, model_predict[,2])
  cox <- pseudo_Rsqrs(model_predict, test_y_input)
  
  #output <- c(C_stat, calib[1], calib[2], valprob, CIL, rMSPE, rmse, caret, MAPE, inbuiltmae, cox[2], model_label, model_input, method_input, number_input, selection_input)
  #output <- c(C_stat, calib[1], calib[2], CIL, rMSPE, MAPE, cox[2], model_label, model_input, method_input, number_input, selection_input, paste(method_input, number_input, selection_input, sep = "_"))
  output <- c(C_stat[1], C_stat[2], C_stat[3], calib[1], calib[2], calib[3], CIL[1], CIL[2], CIL[3], cox[1], model_label, model_input, method_input, number_input, selection_input, paste(method_input, number_input, selection_input, sep = "_"))
  
  outputtuning <- model_fit$bestTune
  
  #return(output)
  return(list(output = output, outputtuning = outputtuning))
}

modeling_tuning <- function(model_labels, method_inputs, number_inputs, selection_inputs){
  
  combination_index <- 1
  model_label_index <- 1
  
  ## Bootstrap
  number_input_boot = 500
  selection_input_boot = "best"
  
  ## LOOCV
  method_input_LOOCV = "LOOCV"
  
  
  #results <- matrix(0, (length(model_labels)*length(number_inputs)*length(selection_inputs) + length(model_labels) + length(model_labels)), 16)
  #colnames(results) <- c("C-statistic", "calibration intercept", "CS", "valprob", "CIL", "rMSPE", "rmse", "caret", "MAPE", "inbuiltmae", "Cox-Snell Rsqr", "model label", "package", "tuning method", "number folds", "selection criterium")
  #results <- matrix(0, (length(model_labels)*length(number_inputs)*length(selection_inputs) + length(model_labels) + length(model_labels)), 13)
  #colnames(results) <- c("C-statistic", "calibration intercept", "CS", "CIL", "rMSPE", "MAPE", "Cox-Snell Rsqr", "model label", "package", "tuning method", "number folds", "selection criterium")
  
  results <- data.frame(matrix(nrow = (length(model_labels)*length(number_inputs)*length(selection_inputs) + length(model_labels) + length(model_labels)+1), ncol = 16))
  #colnames(results) <- c("C-statistic", "CS", "CIL", "Cox-Snell Rsqr", "model label", "package", "tuning method", "number folds", "selection criterium", "combination label")
  
  besttunes <- list()
  
  #fit_slr <- glm.fit(x = train_x, y = train_y, family = binomial(link=logit))
  #predict_slr <- predict(fit_slr, test_x, type = "response")
  #plot(test_y, predict_slr[,2])
  #plot(train_y, fit_slr$fitted.values)
  #return()
  
  
  cat(c("model label = standard", "\n", "model = glm", "\n", "tuning method = n/a", "\n", "number = n/a", "\n", "selection_input = n/a", "\n", "\n"))
  #results[combination_index,] <- fit_model(train_x, train_y_factor, test_x, test_y, model_input, method_input, number_input_boot, "best", grid, model_label)
  #results[combination_index,]<- fit_model(train_x, train_y_factor, test_x, test_y, model_input, method_input, number_input_boot, "best", grid, model_label)
  result <- fit_model(train_x, train_y, test_x, test_y, "glm", "", 0, "", 0, "standard")
  results[combination_index,] <- result$output
  
  combination_index <- combination_index + 1
  
  
  for (model_label in model_labels){
    model_vars <- grid(model_label, train_x)
    grid <- model_vars$tuneGrid
    model_input <- model_vars$model
    
    
    for (method_input in method_inputs){
      if ((method_input == "boot") | (method_input == "LOOCV")){
        
        cat(c("model label =", model_label, "\n", "model =", model_input, "\n", "tuning method =", method_input, "\n", "number =", number_input_boot, "\n", "selection_input = best", "\n", "\n"))
        #results[combination_index,] <- fit_model(train_x, train_y_factor, test_x, test_y, model_input, method_input, number_input_boot, "best", grid, model_label)
        #results[combination_index,]<- fit_model(train_x, train_y_factor, test_x, test_y, model_input, method_input, number_input_boot, "best", grid, model_label)
        result <- fit_model(train_x, train_y, test_x, test_y, model_input, method_input, number_input_boot, "best", grid, model_label)
        results[combination_index,] <- result$output
        besttunes[[combination_index]] <- result$outputtuning 
        
        combination_index <- combination_index + 1
      }
      
      if (method_input == "cv"){ 
        for (number_input in number_inputs){
          for (selection_input in selection_inputs){  
            cat(c("model label =", model_label, "\n", "model =", model_input, "\n", "tuning method =", method_input, "\n", "number =", number_input, "\n", "selection_input", selection_input, "\n", "\n"))
            
            #results[combination_index,] <- fit_model(train_x, train_y_factor, test_x, test_y, model_input, method_input, number_input, selection_input, grid, model_label)
            #results[combination_index,] <- fit_model(train_x, train_y_factor, test_x, test_y, model_input, method_input, number_input, selection_input, grid, model_label)
            result <- fit_model(train_x, train_y, test_x, test_y, model_input, method_input, number_input, selection_input, grid, model_label)
            results[combination_index,] <- result$output
            besttunes[[combination_index]] <- result$outputtuning 
            
            combination_index <- combination_index + 1
          }
        }
        
      }
      
      #save(results, file = "results_RP_C75_sz355.RData")
    }
  }
  #colnames(results) <- c("C-statistic", "calibration intercept", "CS", "CIL", "rMSPE", "MAPE", "Cox-Snell Rsqr", "model label", "package", "tuning method", "number folds", "selection criterium", "combination tuning method")
  colnames(results) <- c("C-statistic", "C-statistic CI_lower", "C-statistic CI_upper", "CS", "CS CI_lower", "CS CI_upper","CIL", "CIL CI_lower", "CIL CI_upper", "Cox-Snell Rsqr", "model label", "package", "tuning method", "number folds", "selection criterium", "combination label")
  
  #save(results, file = "randomforest_pm_c75_sz338.RData")
  #save(besttunes, file = "randomforest_tp_C75_sz355.RData")
  #return(results)
  return(list('performance measures' = results, 'best tuning parameter' = besttunes))
}



## Data set up ##
# Load data
if (!HPC){
  setwd("~/Master jaar 2/Master thesis/Research report") # Set working directory
}
full_data <- read.table("iota_utrecht.txt", header=T, sep="\t") # Load data
if (!HPC){
  setwd("~/Phd project Julius Centre/Project hyperparameter tuning/Case study") # Set working directory
}
# Outcome
data_y <- as.data.frame(full_data$outcome1)
colnames(data_y) <- c('outcome_binary') 
# Predictors
data_x <- as.data.frame(cbind(full_data$Age, full_data$lesdmax, full_data$papnr))
colnames(data_x) <- c('Age', 'lesdmax', 'papnr')

n <- nrow(full_data) # Number of observations

## Sample size calculations ##
outcomep <- mean(full_data$outcome1) # Calculate proportion of outcome
set.seed(100)
R2 <- approximate_R2(auc = 0.75, prev = outcomep, n=nrow(full_data)) # Calculate R^2 for dataset
#pmsampsize(type = "b", rsquared = as.numeric(R2[2]), parameters = 6, shrinkage = 0.9, prevalence = outcomep)
# Minimum sample size = 338

## Sample data set for reduced sample size
set.seed(100) # Set seed
train_index <- sample(1:n, 338)

## Split data into train and test set
# Create training dataset
train_x <- as.matrix(data_x[train_index, ])
test_x <- as.matrix(data_x[-train_index, ])

# Create test dataset
train_y <- as.matrix(data_y[train_index, ])
test_y <- as.matrix(data_y[-train_index, ])

## Execute analysis ##
set.seed(100)
results <- modeling_tuning(model_labels, method_inputs, number_inputs, selection_inputs)

save(results, file = "case_study_results.RData")
