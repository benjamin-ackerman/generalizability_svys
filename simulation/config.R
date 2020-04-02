library(dplyr);library(purrr);library(stringr)

################################
###### Parameters to Vary ######
################################

# Population size:
N = 1000000

# Number of samples to draw from each population
n_samples = 1

# Parameters for population data: scale and true S models
pop_data_params = list(x_cor = c(0, 0.3, 0.6),
                       scale_survey = seq(0,1,by = 0.3),
                       scale_rct = seq(0,1,by = 0.3),
                       coefs_survey = list(list(b1 = 2,b2 = 0, b3 = 1, b4 = 1, b5 = 1, b6 = 0)),
                       coefs_rct = list(list(b1 = 1,b2 = 1, b3 = 2, b4 = 0, b5 = 1, b6 = 0)),
                       scale_tx_het = c(0,0.3)) %>% 
  cross()

# Sample sizes for the survey and trial:
sample_data_params = list(prop_rct = 0.0006,
                          prop_survey = 0.004)

# Parameters for fitting models and estimating the ATE:
fit_params = list(survey_weights = c(FALSE, TRUE),
                  transport_weights = c(FALSE,TRUE),
                  multiply_weights = FALSE,
                  weight_method = c("ps","super_moodie","super_luedtke"),
                  transport_formula = list(list(),
                                           list("x1"),
                                           list("x1","x3","x5")),
                  bootstrap = FALSE,
                  n_bootstrap = 1000) %>% 
  cross()

### Only keep scenarios for the three different couplings of coefficients
scenarios = seq_along(pop_data_params) %>% 
  map(~(((identical(pop_data_params[[.]]$coefs_survey,list(b1 = 2,b2 = 0, b3 = 1, b4 = 1, b5 = 1, b6 = 0)) &
            identical(pop_data_params[[.]]$coefs_rct, list(b1 = 1,b2 = 1, b3 = 2, b4 = 0, b5 = 1, b6 = 0))) |
           
           (identical(pop_data_params[[.]]$coefs_survey,list(b1 = 20,b2 = 0, b3 = 1, b4 = 1, b5 = 1, b6 = 0)) &
              identical(pop_data_params[[.]]$coefs_rct, list(b1 = 1,b2 = 0, b3 = 20, b4 = 1, b5 = 1, b6 = 0))) |
           
           (identical(pop_data_params[[.]]$coefs_survey,list(b1 = -2,b2 = 0, b3 = -1, b4 = -1, b5 = -1, b6 = 0)) &
              identical(pop_data_params[[.]]$coefs_rct, list(b1 = 2,b2 = 0, b3 = 1, b4 = 1, b5 = 1, b6 = 0)))
  ))) %>% unlist() %>% which()
pop_data_params = pop_data_params[scenarios]

### Get rid of scenarios where transport weights aren't there but survey weights are
nonscenarios = seq_along(fit_params) %>%
  map(~((fit_params[[.]]$transport_weights == FALSE & fit_params[[.]]$survey_weights == TRUE) | 
          (fit_params[[.]]$multiply_weights == TRUE & fit_params[[.]]$survey_weights == TRUE) | 
          (fit_params[[.]]$multiply_weights == TRUE & fit_params[[.]]$transport_weights == FALSE))) %>%
  unlist() %>% which()
if(length(nonscenarios) != 0){
  fit_params = fit_params[-nonscenarios]
}

rm(scenarios)
rm(nonscenarios)
