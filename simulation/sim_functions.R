source("buchanan_variance.R")

### Assign the true models of survey and trial selection to the x's
data_pop = function(N, x_cor, scale_survey, scale_rct, coefs_survey, coefs_rct, scale_tx_het){
  x = mvrnorm(N, mu = c(0,0,0,0,0,0),Sigma = matrix(c(1,x_cor,0,0,0,0,
                                                      x_cor,1,0,0,0,0,
                                                      0,0,1,x_cor,0,0,
                                                      0,0,x_cor,1,0,0,
                                                      0,0,0,0,1,x_cor,
                                                      0,0,0,0,x_cor,1),
                                                    nrow=6))

  if(ncol(x) != 6){stop("x must have six columns!")}
  
  colnames(x) = paste0("x",1:6)
  
  # # Specify the true propensity score model for survey selection
  # p_survey = expit(scale_survey*(coefs_survey$b1*x[,1] + coefs_survey$b2*x[,2] + coefs_survey$b3*x[,3] + 
  #                                  coefs_survey$b4*x[,4] + coefs_survey$b5*x[,5] + coefs_survey$b6*x[,6]))
  # 
  # # Specify the true propensity score model for the trial selection
  # p_rct = expit(scale_rct*(coefs_rct$b1*x[,1] + coefs_rct$b2*x[,2] + coefs_rct$b3*x[,3] + 
  #                                coefs_rct$b4*x[,4] + coefs_rct$b5*x[,5] + coefs_rct$b6*x[,6]))
  
  # Specify the true propensity score model for survey selection
  survey = expit(scale_survey*(coefs_survey$b1*x[,1] + coefs_survey$b2*x[,2] + coefs_survey$b4*x[,4] + 
                                 coefs_survey$b3*x[,3] + coefs_survey$b5*x[,5] + coefs_survey$b6*x[,6]))
  
  # Specify the true propensity score model for the trial selection
  rct = expit(scale_rct*(coefs_rct$b1*x[,1] + coefs_rct$b2*x[,2] + coefs_rct$b3*x[,3] + 
                           coefs_rct$b5*x[,5] + coefs_rct$b4*x[,4] + coefs_rct$b6*x[,6]))

  # Gather the data together
  data = data.frame(p_survey = survey, p_rct = rct,
                    x1 = x[,1], x2=x[,2], x3=x[,3],x4=x[,4], x5=x[,5], x6 = x[,6],
                    stringsAsFactors = FALSE)

  # Generate the potential outcomes for everyone in the population
  data$Y0 = rnorm(nrow(data),mean = 0, sd = 1)
  data$Y1 = rnorm(nrow(data),mean = 0 + 2 + scale_tx_het*(data$x1 + data$x2 +
                    data$x3 + data$x4 + data$x5 + data$x6), 1)
  
  return(list(data = data,
              PATE_true = mean(data$Y1 - data$Y0),
              scale_survey = scale_survey,
              scale_rct = scale_rct,
              coefs_survey = coefs_survey,
              coefs_rct = coefs_rct,
              scale_tx_het = scale_tx_het,
              x_cor = x_cor))
}

### Randomly assign S = survey or S = rct for all individuals
assign_samplings = function(df_pop_obj, prop_rct, prop_survey){
  df_pop = df_pop_obj$data
  
  N = nrow(df_pop)
  gc()
  # Generate S from the sample selection probabilities and then randomly sample the trial and validation data
  # df_pop$S_survey = rbinom(N, 1, prob = df_pop$p_survey)
  # df_pop$S_rct = rbinom(N, 1, prob = df_pop$p_rct)
  df_pop$rownum = rownames(df_pop)
  gc()
  df_pop$rct = df_pop$p_rct*prop_rct/mean(df_pop$p_rct)
  gc()
  df_pop$survey = df_pop$p_survey*prop_survey/mean(df_pop$p_survey)
  gc()
  df_pop$p0 = 1 - df_pop$rct - df_pop$survey
  gc()
  df_pop$S = as.vector(rMultinom(df_pop[,c("survey","rct","p0")],1))
  gc()
  df_pop_obj$data = df_pop
  gc()
  
  # Save the sample sizes for that simulation run
  raw_n_rct = sum(df_pop$S == "rct")
  raw_n_survey = sum(df_pop$S == "survey")
  
  asmd_survey = asmd(df_pop$survey[which(df_pop$S == "survey")], df_pop$survey)#[which(df_pop$S != "survey")])
  asmd_rct = asmd(df_pop$rct[which(df_pop$S == "rct")], df_pop$rct)#[which(df_pop$S != "rct")])
  
  # Save 
  df_pop_obj$revision_output = data.frame(scale_survey = df_pop_obj$scale_survey, 
                                          scale_rct = df_pop_obj$scale_rct,
                                          prop_survey, prop_rct, 
                                          raw_n_survey, raw_n_rct, 
                                          asmd_survey, asmd_rct)
  
  return(df_pop_obj)
}

### Modify this to sample the trial and population separately based on the two models
data_samples = function(df_pop_obj, seed = NULL){
  if(!is.null(seed)){set.seed(seed)}
  
  df_pop = df_pop_obj$data
  
  df_pop_obj$asmd_survey = asmd(df_pop$survey[which(df_pop$S == "survey")], df_pop$survey[which(df_pop$S != "survey")])
  df_pop_obj$asmd_rct = asmd(df_pop$rct[which(df_pop$S == "rct")], df_pop$rct[which(df_pop$S != "rct")])
  
  # Create the survey sample
  survey = df_pop %>% 
    filter(S == "survey") %>% 
    mutate(S = 0,
           A = 0,
           Y = Y1*A + Y0*(1-A))
  
  # Create the rct sample
  rct = df_pop %>% 
    filter(S == "rct") %>% 
    mutate(S = 1)
  
  n_rct = nrow(rct)
  
  rct$A = rbinom(n_rct,1,0.5)
  rct$Y = rct$Y1*rct$A + rct$Y0*(1-rct$A)

  # Gather the data together
  data = rct %>%
    bind_rows(survey)
  
  df_pop_obj$data = data

  return(df_pop_obj)
}

data_gen = function(pop, #population data created from data_pop (saved as RDS)
                    prop_rct, # proportion of pop for rct
                    prop_survey, # proportion of pop for survey
                    n_samples = 1 # Number of times to randomly sample from each population
                    ){ 
  
  pop_S = assign_samplings(pop, prop_rct, prop_survey)
  
  if(n_samples == 1){
    samples = data_samples(pop_S,seed = NULL)
  } else{
    samples = 1:n_samples %>% 
      map(~data_samples(pop_S, .))
  }
  
  return(samples)
}

### Turn list of parameters to include in model to a formula
s_formula = function(variables_to_exclude){
  vars = paste0("x",1:6)
  
  vars_exclude = unlist(variables_to_exclude)
  # Get rid of any fake variables that shouldn't belong
  vars_exclude = vars_exclude[stringr::str_detect(vars_exclude,"^x[1-6]")]
  vars_exclude = stringr::str_replace(vars_exclude,"_","^")
  vars_exclude = stringr::str_replace(vars_exclude,"-","*")
  
  if(length(vars) == 0){
    model = paste0("S~", paste(vars,collapse="+"))
  }
  else{
    model = paste0("S~", paste(setdiff(vars,vars_exclude),collapse="+"))
  }
  return(model)
}

# Look at the mean difference between the validation and trial
x_mean_diffs = function(data){
  means = data %>%
    group_by(S) %>% 
    summarise_at(.vars=c("x1","x2","x3","x4","x5","x6"), .funs=mean) %>% 
    dplyr::select(-S) %>% 
    t() %>% 
    as.data.frame() %>% 
    mutate(diff = V1-V2,
           covariate = row.names(.)) %>% 
    dplyr::select(covariate,diff)
  
  return(means)
}

# Calculate the weighted means of the covariates in the trial
x_means = function(data){
  trial_covariate_means = data %>% 
    filter(S == 1) %>% 
    summarise_at(paste0("x",1:6),funs(weighted.mean(.,transport_weights)))
  
  return(trial_covariate_means)
}

### Estimate the survey weights
estimate_survey_weights = function(data_obj, survey_weights, weight_method = "ps"){
  dat = data_obj$data
  
  if(survey_weights == FALSE){
    dat$s_weights = rep(1,nrow(dat))
  }
  
  if(survey_weights == TRUE & weight_method != "rf"){
    # Weight the sample membership model using the inverse probability of survey weights in the population data
    dat$s_weights = ifelse(dat$S == 0, 1/dat$survey, 1)
    
    # Try normalizing the survey weights by the sum of the weights!!! Dec. 4, 2019
    # dat = dat %>%
    #   group_by(S) %>%
    #   mutate(s_weights = s_weights/mean(s_weights)) %>%
    #   as.data.frame()
  }
  
  if(survey_weights == TRUE & weight_method == "rf"){
    dat = dat %>% 
      filter(S == 0) %>% 
      mutate(n_persons = ceiling(1/survey)) %>% 
      tidyr::uncount(n_persons) %>% 
      bind_rows(data %>% filter(S == 1))
  }
  gc()
  data_obj$data = dat
  return(data_obj)
}

### Estimate the transportability weights
estimate_transport_weights = function(data_obj, transport_weights,transport_formula,multiply_weights, weight_method = "ps"){
  if(weight_method == "super_luedtke"){
    weight_method = "super"
    sl_library = c("SL.gam","SL.glm","SL.glm.interaction","SL.nnet","SL.rpart")
  } else if(weight_method == "super_moodie"){
    weight_method = "super"
    sl_library = c("SL.glmnet","SL.knn","SL.mean","SL.randomForest")
  } else{
    sl_library = NULL
  }
  
  dat = data_obj$data
  
  if(transport_weights == FALSE){
    ps = rep(0.5, nrow(dat))
    dat$ps = ps
    dat$transport_weights = ifelse(dat$S == 1, (1-ps)/ps, 0)
  }
  if(transport_weights == TRUE & weight_method != "rf"){
    #transport_design = svydesign(id = ~1, data = dat, weights = dat$s_weights)
    #ps = predict(svyglm(as.formula(s_formula(transport_formula)), design = transport_design, family = 'quasibinomial'), type='response')
    
    ps = weightit(as.formula(s_formula(transport_formula)), data = dat, method = weight_method, 
                       estimand="ATT",focal="0",s.weights = dat$s_weights,
                       stop.method="ks.mean",#gbm parameters
                       SL.library = sl_library)$ps #superlearner parameters
    dat$ps = ps
    dat$transport_weights = ifelse(dat$S == 1, (1-ps)/ps, 0)
  }
  
  if(transport_weights == TRUE & weight_method == "rf"){
    ps = predict(randomForest::randomForest(as.formula(s_formula(transport_formula)), 
                                            data = dat, na.action = na.omit, sampsize = 454, ntree = 1500), 
                 type = "prob")[, 2]
    dat$ps = ps
    dat$transport_weights = ifelse(dat$S == 1, (1-ps)/ps, 0)
  }
  gc()
  if(multiply_weights == TRUE){
    dat$transport_weights = dat$transport_weights*(1/dat$survey)
  }
  
  data_obj$data = dat
  gc()
  return(data_obj)
}

### Function to multiply x coefficients by the x means (weighted by the transportability weights)
cov_effects = function(model_coef, covariate_mean){
  return(model_coef * covariate_mean)
}

### Estimate the PATE
estimate_ATE = function(data_obj, survey_weights, transport_weights,transport_formula,multiply_weights, weight_method){
  # Extract data frame from data generating object
  estimate_weights = estimate_survey_weights(data_obj, survey_weights, weight_method) %>% 
    estimate_transport_weights(transport_weights,transport_formula,multiply_weights, weight_method) 
  gc()
  dat = estimate_weights$data
  
  ATE_design = svydesign(id = ~1, data = dat %>% filter(S == 1), weights = dat$transport_weights[which(dat$S == 1)])
  ATE_model = svyglm(Y ~ A, design = ATE_design, family='gaussian')
  #ATE_model = lm(Y ~ A, data = dat %>% filter(S == 1), weights = transport_weights)
  gc()
  ATE = summary(ATE_model)$coefficients["A",1]
  ATE_se = summary(ATE_model)$coefficients["A",2]
  
  ### SIM REVISION JULY 2: ESTIMATE EFFECTS USING BUCHANAN CODE FOR VARIANCE AS WELL
  vars_included = setdiff(paste0("x",1:6), unlist(transport_formula))
  
  buchanan_results <- ipsw_var(dat, vars_included, "ps", "s_weights", "A", "Y")
  
  return(list(data = dat,
              ATE = ATE,
              ATE_se = ATE_se,
              ATE_2 = buchanan_results$PATE2,
              ATE_se_2 = buchanan_results$se2))
}

### Bootstrap the survey sample
bootstrap_survey = function(data_obj, survey_weights, transport_weights, transport_formula,multiply_weights, weight_method){
  dat = data_obj$data
  
  # Sample the survey with replacement
  dat = dat %>%
    filter(S == 1) %>% 
    # Bootstrap the trial too
    sample_n(nrow(dat %>% filter(S == 1)), replace = TRUE) %>% 
    # Bootstrap the survey using stratified sampling
    bind_rows(dat %>% 
                filter(S == 0) %>% 
                {if(data_obj$scale_survey == 0)
                  sample_n(., nrow(dat %>% filter(S == 1)), replace = TRUE) 
                  else
                    mutate(., strata = cut(survey, breaks = c(0,quantile(survey, probs = seq(.1, 1, 0.1))),
                                              labels=paste0("decile_",1:10))) %>% 
                    group_by(strata) %>% 
                    mutate(n_h=n(),
                           ID = paste0(strata,"_",row_number())) %>% 
                    sample_n(n_h - 1, replace = TRUE) %>% 
                    ungroup() %>% 
                    group_by(ID) %>% 
                    mutate(m_hi_star = n()) %>% 
                    ungroup() %>% 
                    mutate(survey = survey*(n_h-1)/n_h *(1/m_hi_star)) %>%
                    dplyr::select(-strata,-n_h,-ID,-m_hi_star)
                })
  
  data_obj$data = dat
  
  estimate = estimate_ATE(data_obj, survey_weights, transport_weights, transport_formula, multiply_weights, weight_method)$ATE
  
  return(estimate)
}

# Gather and return all results
get_results = function(data_obj, survey_weights, transport_weights, transport_formula, bootstrap, n_bootstrap,multiply_weights, weight_method){
  gc()
  if(bootstrap){
    estimates = 1:n_bootstrap %>% 
      map_dbl(~bootstrap_survey(data_obj, survey_weights, transport_weights, transport_formula, multiply_weights))
    ATE = mean(estimates)
    ATE_se = sd(estimates)
    CI_l = quantile(estimates,0.025)
    CI_u = quantile(estimates,0.975)
  } else{
    ATE_obj = estimate_ATE(data_obj, survey_weights, transport_weights, transport_formula, multiply_weights, weight_method)
    ATE = ATE_obj$ATE
    ATE_se = ATE_obj$ATE_se
    
    ### Buchanan results
    ATE_2 = ATE_obj$ATE_2
    ATE_se_2 = ATE_obj$ATE_se_2
  }
  gc()
  bias = ATE - data_obj$PATE_true
  ### Buchanan results
  bias_2 = ATE_2 - data_obj$PATE_true
  
  # Calculate generalizability index as metric of similarity between trial and survey
  #dat = estimate_ATE(data_obj, survey_weights, TRUE, list(),FALSE)$data
  #g_index = gen_index(dat$ps[which(dat$S == 1)],dat$ps[which(dat$S == 0)])
  
  # Find the number of people overlapping in both samples
  n_overlap = sum(table(data_obj$data$rownum) > 1)

  # Save which variables were excluded
  if(length(transport_formula) == 0){
    vars_excluded = "none"
  } else{
    vars_excluded = paste(transport_formula,collapse = "_")
  }
  
  if(bootstrap){
  results = data.frame(survey_weighted = survey_weights,
                       transportability_weighted = transport_weights,
                       weights_multiplied_method = multiply_weights,
                       weighting_method = weight_method,
                       variables_excluded = vars_excluded,
                       x_cor = data_obj$x_cor,
                       scale_survey = data_obj$scale_survey,
                       scale_rct = data_obj$scale_rct,
                       scale_tx_het = data_obj$scale_tx_het,
                       coefs_rct = paste(unlist(data_obj$coefs_rct),collapse="_"),
                       coefs_survey = paste(unlist(data_obj$coefs_survey),collapse="_"),
                       PATE_true = data_obj$PATE_true,
                       ATE = ATE,
                       CI_l = CI_l,
                       CI_u = CI_u,
                       bias = bias,
                       #g_index = g_index,
                       n_overlap = n_overlap,
                       n_rct = sum(data_obj$data$S),
                       n_survey = nrow(data_obj$data) - sum(data_obj$data$S),
                       asmd_survey = data_obj$asmd_survey,
                       asmd_rct = data_obj$asmd_rct,
                       stringsAsFactors = FALSE)
  row.names(results) = NULL
  } else{
    results = data.frame(survey_weighted = survey_weights,
                         transportability_weighted = transport_weights,
                         weights_multiplied_method = multiply_weights,
                         weighting_method = weight_method,
                         variables_excluded = vars_excluded,
                         x_cor = data_obj$x_cor,
                         scale_survey = data_obj$scale_survey,
                         scale_rct = data_obj$scale_rct,
                         scale_tx_het = data_obj$scale_tx_het,
                         coefs_rct = paste(unlist(data_obj$coefs_rct),collapse="_"),
                         coefs_survey = paste(unlist(data_obj$coefs_survey),collapse="_"),
                         PATE_true = data_obj$PATE_true,
                         ATE = ATE,
                         ATE_se = ATE_se,
                         bias = bias, 
                         ### ADD RESULTS FROM BUCHANAN VARIANCE ESTIMATOR
                         ATE_2 = ATE_2,
                         ATE_se_2 = ATE_se_2,
                         bias_2 = bias_2,
                         #g_index = g_index, 
                         n_overlap = n_overlap,
                         n_rct = sum(data_obj$data$S),
                         n_survey = nrow(data_obj$data) - sum(data_obj$data$S),
                         asmd_survey = data_obj$asmd_survey,
                         asmd_rct = data_obj$asmd_rct,
                         stringsAsFactors = FALSE)
    row.names(results) = NULL
  }
  
  return(results)
}


### Run the simulations using the population data by specifying the seed (to sample 1000 times...)
run_samples = function(df_pop_obj, dat_params, fit_params, seed){
  sample_data_obj = data_samples(df_pop_obj, seed)
  
  gc()

  results = seq_along(fit_params) %>% 
    map_df(~get_results(sample_data_obj, fit_params[[.]]$survey_weights, fit_params[[.]]$transport_weights, 
                        fit_params[[.]]$transport_formula, fit_params[[.]]$bootstrap, fit_params[[.]]$n_bootstrap,
                        fit_params[[.]]$multiply_weights, fit_params[[.]]$weight_method))
  return(results)
}

run_pops = function(population, dat_params, fit_params, n_samples){
  gc()
  pop = assign_samplings(readRDS(paste0("pop_data/pop_",population,".rds")), dat_params$prop_rct, dat_params$prop_survey)
  gc()
  if(n_samples == 1){
    results = run_samples(pop, dat_params, fit_params, NULL)
  } else{
    results = 1:n_samples %>% 
      map_df(~run_samples(pop, dat_params, fit_params, .)) %>% 
      group_by(survey_weighted, transportability_weighted, variables_excluded, scale_survey, scale_rct, x_cor, scale_tx_het) %>% 
      summarise(PATE_true = mean(PATE_true),
                ATE = mean(ATE),
                bias = mean(bias),
                g_index = mean(g_index)) %>% 
      ungroup()
  }
  
  return(results)
}

getSampleSizeASMD <- function(single_population, dat_params){
  gc()
  test_obj <- seq_along(dat_params) %>% 
    map_df(~assign_samplings(readRDS(paste0("pop_data/pop_",single_population,".rds")), dat_params[[.]]$prop_rct, dat_params[[.]]$prop_survey)$revision_output)
  gc()
  return(test_obj)
}

##### PLOTTING FUNCTIONS #####
plot_bias = function(data, x, tx_het, coefs_s, coefs_r, weight_grp, path=NULL){
  d = data %>% 
    filter(coefs_rct == coefs_r, coefs_survey == coefs_s) %>% 
    ungroup()
  
  ### Set the colors to be the same as the data example
  myColors <- brewer.pal(4,"Set1")[c(2,3,4,1)]
  names(myColors) <- c("hat(Delta)[naive]","hat(Delta)[transport]","hat(Delta)[svy.wtd]","hat(Delta)[svy.wtd2]")

  lines = c("none" = 1,"lr" = 2, "gbm" = 3, "super_luedtke" = 5, "super_moodie" = 4)
  
  p = data %>% 
    filter(x_cor == x, scale_tx_het == tx_het, coefs_rct == coefs_r, coefs_survey == coefs_s,
           #variables_excluded == "paste(\"None omitted\")",
           weighting_group %in% weight_grp) %>% #,
           #weighting_method %in% c("gbm","none")) %>%
    ungroup() %>% 
    mutate(weighting_method = factor(weighting_method, levels = c("none","lr","gbm","super_luedtke","super_moodie")),
           weighting_group = factor(weighting_group, levels = c(0,1,2,3), labels = c("hat(Delta)[naive]","hat(Delta)[transport]","hat(Delta)[svy.wtd]","hat(Delta)[svy.wtd2]"))) %>% 
    ggplot() +
    geom_line(aes(x = asmd_survey, y = bias, linetype=weighting_method, col = weighting_group)) +
    facet_grid(scale_rct ~ variables_excluded, label = 'label_parsed') +
    labs(color = "Weighting Method",
         #x = expression(paste(gamma[2],": Scale of Difference between Survey and Population")),
         x = "ASMD of survey selection probabilities between survey and population",
         y = "Bias") +
    ggtitle("PATE Bias") +
    theme_bw() +
    ylim(min(d$bias), max(d$bias)) +
    scale_colour_manual(name = "Estimator",
                        labels = label_parse(),
                        values = myColors) +   
    scale_linetype_manual(name = "Weight Estimation Method",
                          labels = c("Unweighted","Logistic Regression","GBM","SuperLearner (Luedtke)","SuperLearner (Moodie)"),
                          values = lines) +
    guides(color = guide_legend(order=1, title.position = "top"),
           linetype = guide_legend(order=2, title.position = "top")) +
    theme(legend.text = element_text(hjust = 0),
          legend.position = "bottom", legend.box="vertical", legend.title.align = 0.5)
  
  # Automatically save if path is present
  if(!is.null(path)){
    ggsave(paste0(path,"bias_",paste(weight_grp,collapse="_"),".png"),p,width=8,height=8,units="in",scale=1)
  }
  
  return(p)
}

plot_coverage = function(data, x, tx_het, coefs_s, coefs_r, weight_grp,path=NULL){
  d = data %>% 
    filter(coefs_rct == coefs_r, coefs_survey == coefs_s) %>% 
    ungroup()
  
  ### Set the colors to be the same as the data example
  myColors <- brewer.pal(4,"Set1")[c(2,3,4,1)]
  names(myColors) <- c("hat(Delta)[naive]","hat(Delta)[transport]","hat(Delta)[svy.wtd]","hat(Delta)[svy.wtd2]")
  
  lines = c("none" = 1,"lr" = 2, "gbm" = 3, "super_luedtke" = 5, "super_moodie" = 4)
  
  p = data %>% 
    filter(x_cor == x, scale_tx_het == tx_het, coefs_rct == coefs_r, coefs_survey == coefs_s,
           #variables_excluded == "paste(\"None omitted\")",
           weighting_group %in% weight_grp) %>% #,
           #weighting_method %in% c("gbm","none")) %>%
    ungroup() %>% 
    mutate(weighting_method = factor(weighting_method, levels = c("none","lr","gbm","super_luedtke","super_moodie")),
           weighting_group = factor(weighting_group, levels = c(0,1,2,3), labels = c("hat(Delta)[naive]","hat(Delta)[transport]","hat(Delta)[svy.wtd]","hat(Delta)[svy.wtd2]"))) %>% 
    ggplot() +
    geom_line(aes(x = asmd_survey, y = coverage*100, linetype=weighting_method, col = weighting_group)) +
    facet_grid(scale_rct ~ variables_excluded, label = 'label_parsed') +
    labs(color = "Weighting Method",
         #x = expression(paste(gamma[2],": Scale of Difference between Survey and Population")),
         x = "ASMD of survey selection probabilities between survey and population",
         y = "95% CI Coverage (%)") +
    ggtitle("Empirical 95% CI Coverage") +
    theme_bw() +
    ylim(0,100)  +
    scale_colour_manual(name = "Estimator",
                        labels = label_parse(),
                        values = myColors) +   
    scale_linetype_manual(name = "Weight Estimation Method",
                          labels = c("Unweighted","Logistic Regression","GBM","SuperLearner (Luedtke)","SuperLearner (Moodie)"),
                          values = lines) +
    guides(color = guide_legend(order=1, title.position = "top"),
           linetype = guide_legend(order=2, title.position = "top")) +
    theme(legend.text = element_text(hjust = 0),
          legend.position = "bottom", legend.box="vertical", legend.title.align = 0.5)
  
  # Automatically save if path is present
  if(!is.null(path)){
    ggsave(paste0(path,"coverage_",paste(weight_grp,collapse="_"),".png"),p,width=8,height=8,units="in",scale=1)
  }
  
  return(p)
}

plot_coverage_lines = function(data, w, x, tx_het, coefs_s, coefs_r,
                               s_scale, r_scale, vars){
  weight_lab = ifelse(w == 0,"Naive",ifelse(
    w == 1,"Transported w/o Survey Weights","Transported w/ Survey Weights"
  ))
  
  x_limits = data %>% 
    filter(scale_tx_het == tx_het, scale_rct == r_scale, scale_survey == s_scale, 
           coefs_rct == coefs_r, coefs_survey == coefs_s, x_cor == x, variables_excluded == vars) %>% 
    summarise(min_ATE = min(ATE - 1.96*ATE_se),
              max_ATE = max(ATE + 1.96*ATE_se))
  
  d = data %>% 
    filter(weighting_group == w, scale_tx_het == tx_het, scale_rct == r_scale, scale_survey == s_scale, 
           coefs_rct == coefs_r, coefs_survey == coefs_s, x_cor == x, variables_excluded == vars) %>% 
    arrange(desc(ATE)) %>% 
    mutate(id = row_number())
  
  d$test = coverage(d$ATE, d$ATE_se, 2)
  c = mean(d$test)
  
  p = d %>% 
    ggplot(aes(x = as.factor(id), y = ATE)) +
    geom_errorbar(aes(ymin=as.numeric(ATE-1.96*ATE_se), ymax=as.numeric(ATE+1.96*ATE_se), colour = as.factor(test)), width=.1) +
    #geom_point() + 
    geom_hline(yintercept = 2, color="green",linetype="dashed")+
    geom_hline(yintercept = mean(d$ATE), color="light blue",linetype="dashed")+
    theme_classic()+theme(text = element_text(family="Times New Roman"),
                          axis.title.y=element_blank(),
                          axis.text.y=element_blank(),
                          axis.ticks.y=element_blank())+
    labs(y = "mu_v 95% CI")+coord_flip() +
    #ylim(x_limits$min_ATE, x_limits$max_ATE) +
    ylim(0, 3.5) +
    ggtitle(paste0("Empirical Coverage: ",signif(c*100,3),"% (",weight_lab,")")) +
    theme(legend.position = "none",
          plot.title = element_text(size = 10)) +
    scale_color_manual(breaks = c("0","1"),
                       values=c("red", "black"))
  
  return(p)
}

### Make all of the plots at once
make_all_plots = function(data_name){
  results = readRDS(paste0("simulation/results/",data_name))
  
  results = results %>% 
    mutate(weighting_group = survey_weighted + transportability_weighted + weights_multiplied_method,
           weighting_group = ifelse(weights_multiplied_method, weighting_group + 1, weighting_group),
           asmd_survey = ifelse(is.nan(asmd_survey),0, asmd_survey),
           asmd_rct = ifelse(is.nan(asmd_rct),0, asmd_rct),
           numeric_scale_rct = scale_rct) %>% 
    dplyr::select(-survey_weighted, -transportability_weighted, -weights_multiplied_method) %>% 
    filter(variables_excluded %in% c("none","x1","x1_x3_x5"))
  results$scale_rct = paste0("gamma[1] == ",results$scale_rct)
  results$variables_excluded = as.factor(results$variables_excluded)
  levels(results$variables_excluded) = c(expression(paste("None omitted")), 
                                         expression(paste(X[1], " omitted")), 
                                         expression(paste(X[1],", ",X[3]," & ",X[5]," omitted")))
  
  path = paste0("simulation/figures/",str_split(data_name,"_",simplify = TRUE)[1],"_")
  
  plot_bias(results,x = 0.3, tx_het = 0.3, coefs_s = "2_0_1_1_1_0", coefs_r = "1_1_2_0_1_0", path)
  
  ### ASMD to scale
  asmd_to_scale = results %>% 
    ggplot() +
    geom_line(aes(x = scale_survey, y = asmd_survey, colour = as.factor(x_cor))) +
    labs(color = "X Correlation",
         x = expression(paste(gamma[2],": Scale of Difference between Sample and Population")),
         y = "ASMD of selection probabilities between sample and population") 
  ggsave(paste0(path,"asmd_to_scale.png"),asmd_to_scale, width = 11, height = 7, units = "in")
  
  plot_coverage(results,x = 0.3, tx_het = 0.3, coefs_s = "2_0_1_1_1_0", coefs_r = "1_1_2_0_1_0", path)
}


plot_coverage_bstrap = function(data, x, tx_het, coefs_s, coefs_r, weight_grp,path=NULL){
  d = data %>% 
    filter(coefs_rct == coefs_r, coefs_survey == coefs_s) %>% 
    ungroup()
  
  ### Set the colors to be the same as the data example
  myColors <- brewer.pal(4,"Set1")[c(2,3,4,1)]
  names(myColors) <- c("hat(Delta)[naive]","hat(Delta)[transport]","hat(Delta)[svy.wtd]","hat(Delta)[svy.wtd2]")
  
  lines = c("FALSE" = 1,"TRUE" = 2)
  
  p = data %>% 
    filter(x_cor == x, scale_tx_het == tx_het, coefs_rct == coefs_r, coefs_survey == coefs_s,
           #variables_excluded == "paste(\"None omitted\")",
           weighting_group %in% weight_grp) %>% #,
    #weighting_method %in% c("gbm","none")) %>%
    ungroup() %>% 
    mutate(weighting_group = factor(weighting_group, levels = c(0,1,2,3), labels = c("hat(Delta)[naive]","hat(Delta)[transport]","hat(Delta)[svy.wtd]","hat(Delta)[svy.wtd2]"))) %>% 
    ggplot() +
    geom_line(aes(x = scale_survey, y = coverage*100, linetype=bootstrap, col = weighting_group)) +
    facet_grid(scale_rct ~ variables_excluded, label = 'label_parsed') +
    labs(color = "Weighting Method",
         x = expression(paste(gamma[2],": Scale of Difference between Survey and Population")),
         #x = "ASMD of survey selection probabilities between survey and population",
         y = "95% CI Coverage (%)") +
    ggtitle("Empirical 95% CI Coverage") +
    theme_bw() +
    ylim(0,100)  +
    scale_colour_manual(name = "Estimator",
                        labels = label_parse(),
                        values = myColors) +   
    scale_linetype_manual(name = "Variance Estimation",
                          labels = c("Sandwich Estimator","Double Bootstrap"),
                          values = lines) +
    guides(color = guide_legend(order=1, title.position = "top"),
           linetype = guide_legend(order=2, title.position = "top")) +
    theme(legend.text = element_text(hjust = 0),
          legend.position = "bottom", legend.box="vertical", legend.title.align = 0.5)
  
  # Automatically save if path is present
  if(!is.null(path)){
    ggsave(paste0(path,"coverage_",paste(weight_grp,collapse="_"),".png"),p,width=8,height=7,units="in",scale=1)
  }
  
  return(p)
}
