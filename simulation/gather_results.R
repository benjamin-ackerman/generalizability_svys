library(purrr);library(dplyr)
source("utils.R")
source("config.R")

if(fit_params[[1]]$bootstrap){
  list.files("sim_results/") %>% 
    map_df(~readRDS(paste0("sim_results/",.))) %T>%
    saveRDS(., paste0("combined_results/",format(Sys.Date(),"%m%d%Y"),"_survey_sim_results_full.rds")) %>%
    group_by(survey_weighted, transportability_weighted, weights_multiplied_method, weighting_method, variables_excluded, x_cor, scale_survey, scale_rct, scale_tx_het, coefs_rct, coefs_survey) %>% 
    mutate(covered = coverage_bootstrap(CI_l, CI_u, mean(PATE_true))) %>% 
    summarise(PATE_true = mean(PATE_true),
              ATE = mean(ATE),
              bias = mean(bias),
              coverage = mean(covered),
              #g_index = mean(g_index),
              n_overlap = mean(n_overlap),
              n_rct = mean(n_rct),
              n_survey = mean(n_survey),
              asmd_survey = mean(asmd_survey),
              asmd_rct = mean(asmd_rct)) %>% 
    ungroup() %>% 
    saveRDS(paste0("combined_results/",format(Sys.Date(),"%m%d%Y"),"_survey_sim_results.rds"))
} else{
  list.files("sim_results/") %>% 
    map_df(~readRDS(paste0("sim_results/",.))) %T>%
    saveRDS(., paste0("combined_results/",format(Sys.Date(),"%m%d%Y"),"_survey_sim_results_full.rds")) %>%
    group_by(survey_weighted, transportability_weighted, weights_multiplied_method, weighting_method, variables_excluded, x_cor, scale_survey, scale_rct, scale_tx_het, coefs_rct, coefs_survey) %>% 
    mutate(covered = coverage(ATE, ATE_se, mean(PATE_true))) %>% 
    summarise(PATE_true = mean(PATE_true),
              ATE = mean(ATE),
              bias = mean(bias),
              coverage = mean(covered),
              #g_index = mean(g_index),
              n_overlap = mean(n_overlap),
              n_rct = mean(n_rct),
              n_survey = mean(n_survey),
              asmd_survey = mean(asmd_survey),
              asmd_rct = mean(asmd_rct)) %>% 
    ungroup() %>% 
    saveRDS(paste0("combined_results/",format(Sys.Date(),"%m%d%Y"),"_survey_sim_results.rds"))
}

