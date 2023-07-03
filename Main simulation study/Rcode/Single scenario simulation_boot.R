

generate_x_y <- function(n_samples, current_predictor, betas){
  # Specify means of the predictors
  mus <- c(rep(0, times = current_predictor)) 
  
  # Create covariance matrix of the predictors
  covariance_row <- rep(0.2, times = current_predictor*current_predictor) # Covariance of length n_predictors * n_predictors
  covariance_matrix <- matrix(covariance_row, ncol = current_predictor) # Matrix of covariances of n_predictors by n_predictors  
  diag(covariance_matrix) <- rep(1, times = current_predictor)  # Replace diagonal of matrix with 1's (variances)
  
  # Generate predictor variables vector
  X <- mvrnorm(n = n_samples, mu = mus, Sigma = covariance_matrix)
  
  # Calculate probabilities
  linear <- betas[1] + betas[2]*rowSums(X[,1:(0.5*current_predictor)]) + 3*betas[2]*rowSums(X[,(0.5*current_predictor+1):(0.75*current_predictor)])
  prob <- 1 / (1 + exp(-linear))
  # Generate outcome variable
  y <- rbinom(n = n_samples, size = 1, prob = prob)
  
  return(cbind(y, X))
}

single_scenario_simulation_tuning <- function(start_n_predictors, start_prevs, start_shrinkages, N_train, n_valid, betas_matrix, betas_matrix_valid, shrinkages, prevs, n_predictors, start_simulation, end_simulation){
  
  if (begin_simulation_flag > 0) {
    repeat {
      file_name_recovery_bk <- paste(paste(path_name, "data/sim_results/simulation_data", sep = "/"), "i_simulation", i_simulation, "shrinkage", start_shrinkages, "prevalence", start_prevs, "n_predictors", start_n_predictors, "seed_valid", seed_valid, "seed_train", seed_train, "batch", batch_number, "boot_s.RData", sep="_")
      file_exist <- dir.exists(file_name_recovery_bk)
      if (!file_exist) {
        cat("MYWARNING: System recovered from CRASH at: ", start_simulation, "\n")
        break
      } else {
        start_simulation <- i_simulation + 1
      }
    }
  }
  
  length_simulation <- length(seq(start_simulation, end_simulation))
  current_predictor <- n_predictors[start_n_predictors]
  cat("Predictor = ", current_predictor, "\n")
  
  prev_target <- prevs[start_prevs]
  cat("Prevalence = ", prev_target, "\n")
  
  seed_valid <- scenario_n + 100
 
  
  cat("Seed number validation = ", seed_valid, "\n")   
  set.seed(seed_valid)
  
  betas_valid <- betas_matrix_valid[((start_n_predictors-1)*length(n_predictors) + start_prevs), 3:4]
  n_repeats <- 0 
  repeat{ 
    data_valid <- generate_x_y(n_valid, current_predictor, betas_valid)
    if (sum(data_valid[,1]) > 0) {
      cat("Generate outcome validation data succesfull", "\n")
      break 
    }
    n_repeats <- n_repeats + 1
    if (n_repeats > 10) {
      cat("MYWARNING: After 10 attemps, outcome validation data has zero events", "\n")
      return(NULL) 
    }
  }
  
  cat("Shrinkage = ", shrinkages[start_shrinkages], "\n")
  n_samples <- N_train[[start_shrinkages]]
  
  n <- n_samples[start_prevs,start_n_predictors]
  cat("Sample size = ", n, "\n")
  
  
  #beta_ext <- c(betas[1], rep(betas[2], times = 0.5*current_predictor), rep(3*betas[2], times = 0.25*current_predictor), rep(0, times = 0.25*current_predictor))
  y_valid <- data_valid[,1]
  X_valid <- data_valid[,-1]

  for (i_simulation in start_simulation:end_simulation){
    cat("Simulation run = ", i_simulation, "\n")
    
    seed_train <- i_simulation + 200 + max_simulations*(scenario_n - 1 + batch_number*length(n_predictors)*length(prevs)*length(shrinkages) )
    set.seed(seed_train)
    cat("Seed number train = ", seed_train, "\n")   
    
    betas <- betas_matrix[((start_n_predictors-1)*length(n_predictors) + start_prevs), 3:4]
    
    # Calculate probabilities
    linear <- betas[1] + betas[2]*rowSums(X_valid[,1:(0.5*current_predictor)]) + 3*betas[2]*rowSums(X_valid[,(0.5*current_predictor+1):(0.75*current_predictor)])
    prob <- 1 / (1 + exp(-linear))
    
    data <- generate_x_y(n, current_predictor, betas)
    #beta_ext <- c(betas[1], rep(betas[2], times = 0.5*current_predictor), rep(3*betas[2], times = 0.25*current_predictor), rep(0, times = 0.25*current_predictor))
    y <- data[,1]
    X <- data[,-1]
    actual_prev <- mean(y)
    cat("True prevalence = ", actual_prev, "\n")
    
    tuning <- modeling_tuning(train_y = y, train_x = X, test_x = X_valid, test_y = y_valid, prob = prob, model_labels, method_inputs, number_inputs, selection_inputs)
    nrow_tuning <- nrow(tuning)
    tuning$shrinkage <- rep(shrinkages[start_shrinkages], times = nrow_tuning)
    tuning$n_predictors <- rep(current_predictor, times = nrow_tuning)
    tuning$prevalence <- rep(prev_target, times = nrow_tuning)
    tuning$actual_prev <- rep(actual_prev, times = nrow_tuning)
    tuning$sample_size <- rep(n, times = nrow_tuning)
    tuning$n_simulation <- rep(i_simulation, times = nrow_tuning)
    tuning$seed_train <- rep(seed_train, times = nrow_tuning)
    tuning$seed_valid <- rep(seed_valid, times = nrow_tuning)
    #if (i_simulation == start_simulation){
    #  tuning_sum <- matrix(0, nrow = nrow_tuning, ncol = nr_of_perf_measures)
    #  tuning_sum2 <- matrix(0, nrow = nrow_tuning, ncol = nr_of_perf_measures)
    #  slope_vector <- matrix(0, nrow = nrow_tuning, ncol = length_simulation)
    #  hyperpar_vector_1 <- matrix(0, nrow = nrow_tuning, ncol = length_simulation)
    #  hyperpar_vector_2 <- matrix(0, nrow = nrow_tuning, ncol = length_simulation)
    #}
    #tuning_sum <- tuning_sum + as.numeric(as.character(unlist(tuning[,1:nr_of_perf_measures])))
    #tuning_sum2 <- tuning_sum2 + (as.numeric(as.character(unlist(tuning[,1:nr_of_perf_measures]))))^2
    #slope_vector[,i_simulation-start_simulation+1] <- as.numeric(as.character(unlist(tuning$CS)))
    #hyperpar_vector_1[,i_simulation-start_simulation+1] <- as.numeric(as.character(unlist(tuning$"hyperparameter 1")))
    #hyperpar_vector_2[,i_simulation-start_simulation+1] <- as.numeric(as.character(unlist(tuning$"hyperparameter 2")))
    
    file_name_recovery_bk <- paste(paste(path_name, "data/sim_results/simulation_data", sep = "/"), "i_simulation", i_simulation, "shrinkage", start_shrinkages, "prevalence", start_prevs, "n_predictors", start_n_predictors, "seed_valid", seed_valid, "seed_train", seed_train, "batch", batch_number, "boot_s.RData", sep="_")
    save(list = c("tuning", "i_simulation", "start_shrinkages", "start_prevs", "start_n_predictors"), file = file_name_recovery_bk)
    cat("Backup file successfully saved: ", file_name_recovery_bk, "\n\n\n")

    # save: tuning_sum, tuning_sum2, slope_vector, i_simulation, start_shrinkages, start_prevs, start_n_predictors, tuning_total
    #save(list = c("tuning", "tuning_sum", "tuning_sum2", "slope_vector", "hyperpar_vector_1", "hyperpar_vector_2","i_simulation", "start_shrinkages", "start_prevs", "start_n_predictors"), file = file_name_recovery)
    #cat("Data file successfully saved: ", file_name_recovery, "\n\n\n")
  }

  return(tuning)
}


RMSD_added <- function(slope_vector, tuning){
  
  diff_slope <- log(1)-log(slope_vector)
  
  RMSD <- sqrt(rowMeans(diff_slope^2, na.rm=T))
  
  SE <- rowSds(diff_slope, na.rm=T) / sqrt(ncol(slope_vector))
  CI_lower <- RMSD - 1.96*SE
  CI_upper <- RMSD + 1.96*SE
  
  tuning$RMSD <- RMSD
  tuning$"RMSD CI_lower" <- CI_lower
  tuning$"RMSD CI_upper" <- CI_upper
  
  
  return(tuning)
}

pivot_columns <- function(tuning, tuning_mean, tuning_MC_SE){
  tuning_total <- NULL
  names_tuning <- names(tuning)
  new_names <- NULL
  for (i in 1:nr_of_perf_measures){
    i_name <- names_tuning[i]
    new_names <- c(new_names, i_name)
    new_names <- c(new_names, paste("MC SE ", i_name))
    tuning_total <- cbind(tuning_total, tuning_mean[,i])
    tuning_total <- cbind(tuning_total, tuning_MC_SE[,i])
  }
  tuning_total <- cbind(tuning_total, tuning[,-c(1:nr_of_perf_measures)])
  new_names <- c(new_names, names_tuning[-c(1:nr_of_perf_measures)])
  names(tuning_total) <- new_names
  return(tuning_total)
}

single_scenario_post_proccessing <- function(start_n_predictors, start_prevs, start_shrinkages, N_train, n_valid, betas_matrix, betas_matrix_valid, shrinkages, prevs, n_predictors, start_simulation, end_simulation){
  
  tuning_mean <- tuning_sum / length_simulation
  tuning_MC_SE <- sqrt((tuning_sum2 / length_simulation) - tuning_mean^2) / sqrt(length_simulation)
  tuning_extended <- pivot_columns(tuning, tuning_mean, tuning_MC_SE)
  
  tuning_extended$"C-statistic" <- tuning_extended$"C-statistic" - AUC_target 
  tuning_extended$"C-statistic CI_lower" <- tuning_extended$"C-statistic CI_lower" - AUC_target 
  tuning_extended$"C-statistic CI_upper" <- tuning_extended$"C-statistic CI_upper" - AUC_target
  names(tuning_extended)[names(tuning_extended)=="C-statistic"] <- "delta C-statistic"
  tuning_extended$CS <- rowMedians(slope_vector)
  
  tuning_extended <- RMSD_added(slope_vector = slope_vector, tuning = tuning_extended)
  tuning_extended <- tuning_extended[, c(1:40, 55:57, 41:54)]
  
  return(tuning_extended)
}
