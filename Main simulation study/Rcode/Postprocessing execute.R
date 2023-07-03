start_time <- proc.time()

CHECK_PM <- FALSE  # Flag representing check for performance measures (TRUE) or not (FALSE)
nr_of_perf_measures <- 14

HPC <- TRUE
if (HPC) {
  path_name <- ".." 
} else {
  setwd("~/Master jaar 2/Master thesis")
  path_name <- "." 
}

source(paste(path_name, "Rcode/libraries_file.R", sep = "/"))
source(paste(path_name, "Rcode/Postprocessing simulations.R", sep = "/"))


max_simulations <- 50

AUC_target <<- 0.75
n_predictors <- c(8, 16, 32) # Specify vector of number of predictors
prevs <- c(0.1, 0.3, 0.5)
shrinkages <- c(0.6, 0.9, 0.97) # Define calibration slopes (shrinkage factors)
scenarios <- cbind(c(rep(1, times = 9), 
                     rep(2, times = 9)), 
                   rep(c(rep(1, times = 3), 
                         rep(2, times = 3),
                         rep(3, times = 3)), times = 2),
                   rep(c(1,2,3), times = 6))

for (scenario_n in 1:nrow(scenarios)){
    start_n_predictors <- scenarios[scenario_n, 1] 
    start_prevs <- scenarios[scenario_n, 2]
    start_shrinkages <- scenarios[scenario_n, 3]
    cat("Predictor: ", start_n_predictors, "\n")
    cat("Prevalence: ", start_prevs, "\n")
    cat("Shrinkage: ", start_shrinkages, "\n")
    file_name_all_tuning_boot_added <- paste(paste(path_name, "data/post_processed_results/all_simulation_data", sep = "/"), "shrinkage", start_shrinkages, "prevalence", start_prevs, "n_predictors", start_n_predictors, "s.RData", sep="_")
    file_name_tuning_ext_boot_added <- paste(paste(path_name, "data/post_processed_results/post_processed_results_scenario", sep = "/"), "shrinkage", start_shrinkages, "prevalence", start_prevs, "n_predictors", start_n_predictors, "s.RData", sep="_")

  for (onlyBoot in 0:1) {
    post_processed_results <- single_scenario_post_proccessing(start_n_predictors, start_prevs, start_shrinkages)
    file_name_all_tuning <- paste(paste(path_name, "data/post_processed_results/all_simulation_data", sep = "/"), "shrinkage", start_shrinkages, "prevalence", start_prevs, "n_predictors", start_n_predictors, "boot", onlyBoot, "s.RData", sep="_")
    file_name_tuning_ext <- paste(paste(path_name, "data/post_processed_results/post_processed_results_scenario", sep = "/"), "shrinkage", start_shrinkages, "prevalence", start_prevs, "n_predictors", start_n_predictors, "boot", onlyBoot, "s.RData", sep="_")
    file_exists <- 0
    if (file.exists(file_name_all_tuning)){
      load(file = file_name_all_tuning)
      file_exists <- 1
      cat("Tuning without boot exist " , scenario_n,  "\n")
    }
    if (file.exists(file_name_tuning_ext)){
      load(file = file_name_tuning_ext)
    }
    
    if ( (file_exists==1) & onlyBoot == 0) {
      tuning_all_boot_added <- tuning_all
      tuning_extended_boot_added <- tuning_extended
      hyperpar_vector_1_boot_added <- hyperpar_vector_1
      hyperpar_vector_2_boot_added <- hyperpar_vector_2
      slope_vector_boot_added <- slope_vector
      file_exists <- 0
    } else if ( (file_exists==1) & onlyBoot ==1) {
      tuning_all_boot_added <- rbind(tuning_all_boot_added, tuning_all)
      tuning_extended_boot_added <- rbind(tuning_extended_boot_added, tuning_extended)
      hyperpar_vector_1_boot_added <- rbind(hyperpar_vector_1_boot_added, hyperpar_vector_1)
      hyperpar_vector_2_boot_added <- rbind(hyperpar_vector_2_boot_added, hyperpar_vector_2)
      slope_vector_boot_added <- rbind(slope_vector_boot_added, slope_vector)
      cat("Tuning only boot exist \n")
    } else {
      #cat("Tuning doesnt exist \n")
    }
  }
  file_name_all_tuning_boot_added <- paste(paste(path_name, "data/post_processed_results/all_simulation_data", sep = "/"), "shrinkage", start_shrinkages, "prevalence", start_prevs, "n_predictors", start_n_predictors, "s.RData", sep="_")
  file_name_tuning_ext_boot_added <- paste(paste(path_name, "data/post_processed_results/post_processed_results_scenario", sep = "/"), "shrinkage", start_shrinkages, "prevalence", start_prevs, "n_predictors", start_n_predictors, "s.RData", sep="_")

  save(list = c("tuning_all_boot_added"), file = file_name_all_tuning_boot_added)
  save(list = c("tuning_extended_boot_added", "hyperpar_vector_1_boot_added", "hyperpar_vector_2_boot_added", "slope_vector_boot_added"), file = file_name_tuning_ext_boot_added)
}

