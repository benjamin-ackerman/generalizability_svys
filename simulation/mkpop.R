library(dplyr);library(purrr);library(MASS);library(stringr)
source("sim_functions.R")
source("utils.R")
source("config.R")

start_time = Sys.time()

### Set the seed
set.seed(13783)

### Make the population parameters file if this is the first time you're running this
if(!file.exists("pop_data/pop_data_params.rds")){
  saveRDS(pop_data_params, "pop_data/pop_data_params.rds")
  
  seq_along(pop_data_params) %>% 
    map(function(i){
      pop = data_pop(N, pop_data_params[[i]]$x_cor, pop_data_params[[i]]$scale_survey, pop_data_params[[i]]$scale_rct, 
                     pop_data_params[[i]]$coefs_survey, pop_data_params[[i]]$coefs_rct, 
                     pop_data_params[[i]]$scale_tx_het)
      saveRDS(pop, paste0("pop_data/pop_",i,".rds"))
    })
}

### Make the population data files
if(!identical(pop_data_params, readRDS("pop_data/pop_data_params.rds"))){
  
  unlink("pop_data/*")
  
  seq_along(pop_data_params) %>% 
    map(function(i){
      pop = data_pop(N, pop_data_params[[i]]$x_cor, pop_data_params[[i]]$scale_survey, pop_data_params[[i]]$scale_rct, 
                     pop_data_params[[i]]$coefs_survey, pop_data_params[[i]]$coefs_rct, 
                     pop_data_params[[i]]$scale_tx_het)
      saveRDS(pop, paste0("pop_data/pop_",i,".rds"))
    })
  saveRDS(pop_data_params,"pop_data/pop_data_params.rds")
  
} else{
  cat("Population data in config are already made! \n")
}

end_time = Sys.time()
print(difftime(end_time,start_time,units="mins"))