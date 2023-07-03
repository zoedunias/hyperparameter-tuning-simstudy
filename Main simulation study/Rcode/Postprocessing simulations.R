RMSD_added <- function(slope_vector, tuning){
  
  slope_vector[slope_vector <= 0] <- 1e-8
  
  diff_slope <- log(1)-log(slope_vector)
  
  RMSD <- sqrt(rowMeans(diff_slope^2, na.rm=T))
  
  SE <- rowSds(diff_slope, na.rm=T) / sqrt(ncol(slope_vector))
  #CI_lower <- RMSD - 1.96*SE
  #CI_upper <- RMSD + 1.96*SE
  
  tuning$RMSD <- RMSD
  #  tuning$"RMSD CI_lower" <- CI_lower
  #  tuning$"RMSD CI_upper" <- CI_upper
  tuning$RMSD_SE <- SE
  
  
  return(tuning)
}

pivot_columns <- function(tuning, tuning_mean, tuning_MC_SE){
  tuning_total <- NULL
  names_tuning <- names(tuning)
  new_names <- NULL
  for (i in 1:nr_of_perf_measures){
    i_name <- names_tuning[i]
    new_names <- c(new_names, i_name)
    new_names <- c(new_names, paste("Empirical SE", i_name))
    tuning_total <- cbind(tuning_total, tuning_mean[,i])
    tuning_total <- cbind(tuning_total, tuning_MC_SE[,i])
  }
  tuning_total <- cbind(tuning_total, tuning[,-c(1:nr_of_perf_measures)])
  new_names <- c(new_names, names_tuning[-c(1:nr_of_perf_measures)])
  names(tuning_total) <- new_names
  return(tuning_total)
}

single_scenario_post_proccessing <- function(start_n_predictors, start_prevs, start_shrinkages){
  number_of_files_exist <- 0
  
  seed_valid <- scenario_n + 100
  file_name_all_tuning <- paste(paste(path_name, "data/post_processed_results/all_simulation_data", sep = "/"), "shrinkage", start_shrinkages, "prevalence", start_prevs, "n_predictors", start_n_predictors, "boot", onlyBoot, "s.RData", sep="_")
  file_name_tuning_ext <- paste(paste(path_name, "data/post_processed_results/post_processed_results_scenario", sep = "/"), "shrinkage", start_shrinkages, "prevalence", start_prevs, "n_predictors", start_n_predictors, "boot", onlyBoot, "s.RData", sep="_")
  #cat("File name for all scenarios: ", file_name_all_tuning, "\n")  
  
  for (batch_number in 0:9) { 
    cat("Batch number: ", batch_number, "\n")
    start_simulation <- batch_number*50 + 1
    end_simulation <- start_simulation + 49
    for (i_simulation in start_simulation:end_simulation){
      seed_train <- i_simulation + 200 + max_simulations*(scenario_n - 1 + batch_number*length(n_predictors)*length(prevs)*length(shrinkages) )
      # Load simulation data files
      if (onlyBoot == 0) {
        file_name_recovery_bk <- paste(paste(path_name, "data/all_sim_results_no_boot/simulation_data", sep = "/"), "i_simulation", i_simulation, "shrinkage", start_shrinkages, "prevalence", start_prevs, "n_predictors", start_n_predictors, "seed_valid", seed_valid, "seed_train", seed_train, "batch", batch_number, "s.RData", sep="_")
      } else {
        file_name_recovery_bk <- paste(paste(path_name, "data/all_sim_results_boot/simulation_data", sep = "/"), "i_simulation", i_simulation, "shrinkage", start_shrinkages, "prevalence", start_prevs, "n_predictors", start_n_predictors, "seed_valid", seed_valid, "seed_train", seed_train, "batch", batch_number, "boot_s.RData", sep="_")
      }
      #cat("Simulation run = ", i_simulation, "\n")
      
      if (file.exists(file_name_recovery_bk)){
        load(file = file_name_recovery_bk)
        #cat("File does exist: ", file_name_recovery_bk, "\n")
        number_of_files_exist <- number_of_files_exist + 1
      }
      else{
        cat("File does not exist: ", file_name_recovery_bk, "\n")
        next
      }
      if (number_of_files_exist == 1){
        # Memory allocation
        sim_seq <- i_simulation

	idx_NA_CS <- is.na(tuning$CS)
	idx_neg_CS <- (tuning$CS[!idx_NA_CS] < 0)
	if (sum(idx_neg_CS) > 0) {
		tuning$"C-statistic"[idx_neg_CS] <- 1 - as.numeric(tuning$"C-statistic"[idx_neg_CS])
		tuning$CS[idx_neg_CS] <- -as.numeric(tuning$CS[idx_neg_CS])  
	}
        tuning_all <- tuning
        tuning_sum <- matrix(0, nrow = nrow(tuning), ncol = nr_of_perf_measures)
        tuning_sum <- tuning_sum + as.numeric(as.character(unlist(tuning[,1:nr_of_perf_measures])))
        tuning_sum2 <- matrix(0, nrow = nrow(tuning), ncol = nr_of_perf_measures)
        tuning_sum2 <- tuning_sum2 + as.numeric(as.character(unlist(tuning[,1:nr_of_perf_measures])))^2
        
        slope_vector <- as.numeric(as.character(unlist(tuning$CS)))
        hyperpar_vector_1 <- as.numeric(as.character(unlist(tuning$"hyperparameter 1")))
        hyperpar_vector_2 <- as.numeric(as.character(unlist(tuning$"hyperparameter 2")))
        hyperpar_vector_3 <- as.numeric(as.character(unlist(tuning$"hyperparameter 3")))
      } else {
        sim_seq <- c(sim_seq, i_simulation)
	idx_NA_CS <- is.na(tuning$CS)
	idx_neg_CS <- (tuning$CS[!idx_NA_CS] < 0)
	if (sum(idx_neg_CS) > 0) {
		tuning$"C-statistic"[idx_neg_CS] <- 1 - as.numeric(tuning$"C-statistic"[idx_neg_CS])
		tuning$CS[idx_neg_CS] <- -as.numeric(tuning$CS[idx_neg_CS])  
	}
        slope_vector <- cbind(slope_vector, as.numeric(as.character(unlist(tuning$CS))))
        hyperpar_vector_1 <- cbind(hyperpar_vector_1, as.numeric(as.character(unlist(tuning$"hyperparameter 1"))))
        hyperpar_vector_2 <- cbind(hyperpar_vector_2, as.numeric(as.character(unlist(tuning$"hyperparameter 2"))))
        hyperpar_vector_3 <- cbind(hyperpar_vector_3, as.numeric(as.character(unlist(tuning$"hyperparameter 3"))))
        # Aggregation of data over all simulations
        tuning_sum <- tuning_sum + as.numeric(as.character(unlist(tuning[,1:nr_of_perf_measures])))
        tuning_sum2 <- tuning_sum2 + as.numeric(as.character(unlist(tuning[,1:nr_of_perf_measures])))^2
        tuning_all <- rbind(tuning_all, tuning)
      }
    }
    #cat("Number of files exist in batch ", batch_number, " is: ", number_of_files_exist, "\n")
  }
  if (exists("tuning") & (number_of_files_exist > 1))
  {
    cat("Total number of files exist: ", number_of_files_exist, "\n")

    tuning_mean <- tuning_sum / number_of_files_exist
    tuning_MC_SE <- sqrt((tuning_sum2 / number_of_files_exist) - tuning_mean^2) / sqrt(number_of_files_exist)
    #tuning_MC_SE <- tuning_mean
    tuning_extended <- pivot_columns(tuning, tuning_mean, tuning_MC_SE)
    
    max_slopes <- rowMaxs(slope_vector, na.rm = TRUE)
    slope_nas <- is.na(slope_vector)
    for (k in 1:nrow(slope_vector)){
      #cat("Before correction: ", slope_vector[k, slope_nas[k,]], "\n")
      slope_vector[k, slope_nas[k,]] <- max_slopes[k]
      #cat("After correction: ", slope_vector[k, slope_nas[k,]], "\n")

      tuning_extended$"Empirical SE CS"[k] <- IQR(slope_vector[k,])
    }
    names(tuning_extended)[names(tuning_extended)=="Empirical SE CS"] <- "Empirical IQR CS"
    cat("anyNA in slope_vector: ", anyNA(slope_vector), "\n")
    tuning_extended$CS <- rowMedians(slope_vector)
    
    tuning_extended <- RMSD_added(slope_vector = slope_vector, tuning = tuning_extended)
    #Next line is not correct partially due to less columns in RMSD_added
    #tuning_extended <- tuning_extended[, c(1:2*nr_of_perf_measures, (2*nr_of_perf_measures+15):(2*nr_of_perf_measures+17), (2*nr_of_perf_measures+1):(2*nr_of_perf_measures+14))]
    save(list = c("tuning_all", "start_shrinkages", "start_prevs", "start_n_predictors"), file = file_name_all_tuning)
    save(list = c("tuning_extended", "hyperpar_vector_1", "hyperpar_vector_2", "slope_vector"), file = file_name_tuning_ext)
  } else {
    tuning_extended <- NA
    hyperpar_vector_1 <- NA
    hyperpar_vector_2 <- NA
  }
  
  return(list("tuning_result" = tuning_extended, "hyperpar_vector_1" = hyperpar_vector_1, "hyperpar_vector_2" = hyperpar_vector_2))
}
