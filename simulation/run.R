library(dplyr);library(purrr);library(MASS);library(stringr);library(Hmisc);library(survey);library(tidyr);library(WeightIt)
source("sim_functions.R")
source("utils.R")
source("config.R")

### Calulate the start time
start_time = Sys.time()

### Set the seed (from command line argument)
temp = commandArgs(TRUE)
f = as.numeric(temp[1])
set.seed(f)

pop_data_params = readRDS('pop_data/pop_data_params.rds')

results = seq_along(pop_data_params) %>% #[scenarios_to_use] %>% 
  map_df(~run_pops(., sample_data_params, fit_params, n_samples))

saveRDS(results, paste0("sim_results/results_",f,".rds"))

### calculate the end time
end_time = Sys.time()
print(difftime(end_time,start_time,units="mins"))