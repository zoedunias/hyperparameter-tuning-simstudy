start_time <- proc.time()

HPC <- TRUE # Flag representing whether script runs on HPC (TRUE) or on local computer (FALSE)
CHECK_PM <- FALSE  # Flag representing check for performance measures (TRUE) or not (FALSE)
nr_of_perf_measures <- 14

if (HPC){ # When ran on HPC
  job_args <- commandArgs(trailingOnly = T)
  if (length(job_args) < 5){
    print("No arguments given")
    return()
  }
  path_name <- "/home/julius_bs/zdunias/Thesis" 
}else{ # When ran locally
  setwd("~/Master jaar 2/Master thesis")
  path_name <- "." 
  job_args <- c(1, 0, 50, 0, 0) # 55 = modulo 50 is the scenario number (remainder is the piece of the simulation), 0 = start simulation run from beginning
}
source(paste(path_name, "Rcode/libraries_file.R", sep = "/"))
source(paste(path_name, "Rcode/Single scenario simulation_boot.R", sep = "/"))
source(paste(path_name, "Rcode/model fitting and tuning.R", sep = "/"))
load(file = paste(path_name, "data/generated_betas.RData", sep = "/"))
load(file = paste(path_name, "data/train_sample_sizes_doubled.RData", sep = "/"))
load(file = paste(path_name, "data/generated_betas_valid.RData", sep = "/"))


TASK_idx <- as.numeric(job_args[1])  # 1-10 1st scenario, 11-20 2nd scenario, ... , 261-270 27th scenario
begin_simulation_flag <- as.numeric(job_args[2])
n_data_window <- as.numeric(job_args[3])
batch_number <- as.numeric(job_args[4]) #0 -> 1 to 50 simulations, 1 -> 51 to 100 simulations, ... 
onlyBoot <- as.numeric(job_args[5])

max_simulations <- 50
length_data_window <- max_simulations/n_data_window
scenario_n <- floor((TASK_idx-1)/n_data_window) + 1 # 10 times 1, 10 times 2, ...
start_simulation <- mod(TASK_idx, n_data_window) * length_data_window + 1 + batch_number * max_simulations # 1-5, 6-10, ..., 46-50, 51-55, ..., 96-100, ..., 266-270 
end_simulation <- start_simulation + (length_data_window-1) # start_simulation + 4
cat("Start simulation from: ", start_simulation, " to: ", end_simulation, "\n")

AUC_target <<- 0.75
n_predictors <- c(8, 16, 32) # Specify vector of number of predictors
prevs <- c(0.1, 0.3, 0.5)
shrinkages <- c(0.6, 0.9, 0.97) # Define calibration slopes (shrinkage factors)
scenarios <- cbind(c(rep(1, times = 9), 
                     rep(2, times = 9),
                     rep(3, times = 9)), 
                   rep(c(rep(1, times = 3), 
                         rep(2, times = 3),
                         rep(3, times = 3)), times = 3),
                   rep(c(1,2,3), times = 9))
start_n_predictors <- scenarios[scenario_n, 1] 
start_prevs <- scenarios[scenario_n, 2]
start_shrinkages <- scenarios[scenario_n, 3]

## Tune and modeling
model_labels <- c("ridge", "lasso", "elastic net", "random forest")
if (onlyBoot) {
  method_inputs <- c("boot") 
} else {
  method_inputs <- c("cv", "repeatedcv") 
}
number_inputs <- c(5, 10)
selection_inputs <- c("best", "oneSE")
repeats_inputs <- c(20, 10)

# TEST SETTINGS
#model_labels <- c("ridge")
#method_inputs <- c("cv") 
#number_inputs <- c(5)
#selection_inputs <- c("oneSE")


performance_measures <- single_scenario_simulation_tuning(start_n_predictors = start_n_predictors, start_prevs = start_prevs, start_shrinkages = start_shrinkages, N_train = N_train, n_valid = 1e5, betas_matrix = betas_matrix, betas_matrix_valid = betas_matrix_valid, shrinkages = shrinkages, prevs = prevs, n_predictors = n_predictors, start_simulation = start_simulation, end_simulation = end_simulation)

cat("process time = ", (proc.time()-start_time))
