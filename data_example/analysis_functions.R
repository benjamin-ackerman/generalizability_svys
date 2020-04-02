library(tidyr)
library(purrr)

compareNA <- function(v1,v2) {
  same <- (v1 == v2) | (is.na(v1) & is.na(v2))
  same[is.na(same)] <- FALSE
  return(same)
}

get_pscores = function(data, trial_var, method){
  weight_obj = weighting(outcome = NULL, treatment = NULL, trial = trial_var, 
                         selection_covariates = covariates, data = data,
                         selection_method = method)
  probs = data.frame(prob = unlist(weight_obj$participation_probs),
                     Sample = str_match(names(unlist(weight_obj$participation_probs)),or("population","trial")),
                     Method = method)
  row.names(probs)=NULL
  return(probs)
}

covtab = function(data, trial, covariates, method, svy_wt_var,normalize_method = "none",transport_method="ps",trim_weights=FALSE, trim_pctile = .97){
  if(!(method %in% c("naive","svy","transport","both"))){
    stop("Must choose method from naive, svy, transport, or both", call. = FALSE)
  }
  
  if(!(transport_method %in% c("ps","gbm","super_luedtke","super_moodie"))){
    stop("Must choose transport model fitting method between ps, gbm, super_luedtke and super_moodie", call. = FALSE)
  }
  
  if(!(normalize_method %in% c("none","sum_wt","mean_wt","trial_mean_wt"))){
    stop("Must choose survey weight normalizing method between none, sum_wt, mean_wt and trial_mean_wt", call. = FALSE)
  }
  
  if(normalize_method == "none"){
    data = data %>%
      mutate(svy_wt_normalized = get(svy_wt_var)) %>% 
      ungroup() %>% 
      as.data.frame()
  }else if(normalize_method == "sum_wt"){
    data = data %>%
      group_by(get(trial)) %>% 
      mutate(svy_wt_normalized = get(svy_wt_var)/sum(get(svy_wt_var))) %>% 
      ungroup() %>% 
      as.data.frame()
  }else if(normalize_method == "mean_wt"){
    data = data %>%
      group_by(get(trial)) %>% 
      mutate(svy_wt_normalized = get(svy_wt_var)/mean(get(svy_wt_var))) %>% 
      ungroup() %>% 
      as.data.frame()
  }else if(normalize_method == "trial_mean_wt"){
    n_trial = data %>% 
      filter(get(trial) == 1) %>% 
      summarise(n_row = n()) %>% 
      pull()
    
    data = data %>%
      group_by(get(trial)) %>% 
      mutate(svy_wt_normalized = get(svy_wt_var)*n_trial/sum(get(svy_wt_var))) %>% 
      ungroup() %>% 
      as.data.frame()
  }
  
  if(method %in% c("svy","both")){
    data$svy = ifelse(data[,trial] == 1, 1, data$svy_wt_normalized)
  } else{
    data$svy = 1
  }

  if(transport_method == "super_luedtke"){
    transport_method = "super"
    sl_library = c("SL.gam","SL.glm","SL.glm.interaction","SL.nnet","SL.rpart")
  } else if(transport_method == "super_moodie"){
    transport_method = "super"
    sl_library = c("SL.glmnet","SL.knn","SL.mean","SL.randomForest")
  } else{
    sl_library = NULL
  }
  
  if(method == "naive"){
    data$w = 1
  } else if(method == "svy"){
    data$w = ifelse(data[,trial] == 1, 1, data$svy)
  } else if(method %in% c("transport","both")){
    data$p = weightit(as.formula(paste0(trial," ~",paste(covariates,collapse="+"))), 
                      data = data, method = transport_method, 
                      estimand="ATT",focal="0",s.weights = data$svy,
                      stop.method="ks.mean",n.trees=20000,
                      SL.library = sl_library)$ps
    data$w = ifelse(data[,trial] == 1, (1-data$p)/data$p, data$svy)
    if(trim_weights){
      cutoff = as.numeric(quantile(data$w[which(data[,trial] == 1)], trim_pctile))
      data$w[which(data[,trial] == 1)] = ifelse(data$w[which(data[,trial] == 1)] > cutoff, cutoff, data$w[which(data[,trial] == 1)])
    }
  } 
  
  dmy <- caret::dummyVars(" ~ .", data = data[, covariates], 
                          sep = "_")
  expanded.data = data.frame(trial = data[, trial], 
                             predict(dmy, newdata = data[, covariates]), 
                             weights = data$w)
  names(expanded.data) = stringr::str_replace_all(names(expanded.data), 
                                                  "\\.", "-")
  means.tab = expanded.data %>% 
    dplyr::group_by(trial) %>% 
    dplyr::summarise_at(names(expanded.data)[-1], funs(weighted.mean(., weights))) %>%
    dplyr::select(-weights) %>% 
    t() %>% 
    as.data.frame()
  
  names(means.tab)[which(means.tab["trial", ] == "1")] = "trial"
  names(means.tab)[which(means.tab["trial", ] == "0")] = "population"
  means.tab = means.tab[!rownames(means.tab) %in% c("trial", "X.Intercept."), ]
  
  n_trial = as.numeric(table(expanded.data$trial)["1"])
  n_pop = as.numeric(table(expanded.data$trial)["0"])
  
  sd.tab = expanded.data %>% 
    dplyr::group_by(trial) %>% 
    # dplyr::summarise_at(names(expanded.data)[which(!names(expanded.data) %in% c("trial", "X.Intercept."))], 
    #                     funs(sum(weights *(. - weighted.mean(., weights))^2)/sum(weights))) %>% 
    dplyr::summarise_at(names(expanded.data)[which(!names(expanded.data) %in% c("trial", "X.Intercept."))], 
                       funs(sum(weights)/(sum(weights)^2 - sum(weights^2))*sum(weights*(. - weighted.mean(., weights))^2))) %>% 
    dplyr::select(-weights) %>% 
    t() %>% 
    as.data.frame()
  
  sd.tab$pooled_sd = sqrt(((n_trial - 1) * sd.tab$V1 + (n_pop - 1) * sd.tab$V2)/(n_trial + n_pop - 2))
  
  names(sd.tab)[which(sd.tab["trial", ] == "1")] = "trial_var"
  names(sd.tab)[which(sd.tab["trial", ] == "0")] = "population_var"
  names(sd.tab)[3] = "pooled_sd"
  sd.tab = sd.tab[!rownames(sd.tab) %in% c("trial", "X.Intercept."),]
  
  covariate_table = means.tab %>% 
    dplyr::bind_cols(sd.tab) %>% 
    dplyr::mutate(ASMD = round(abs((trial - population)/pooled_sd), 3)) %>% 
    dplyr::select(trial, population, ASMD)
  
  if(method != "naive"){
    names(covariate_table)[which(names(covariate_table) == 
                                   "population")] = "population (weighted)"
  }
  if(method %in% c("transport","both")){
    names(covariate_table)[which(names(covariate_table) == 
                                   "trial")] = "trial (weighted)"
  }
  
  row.names(covariate_table) = setdiff(names(expanded.data), 
                                       c("trial", "X.Intercept.", "weights"))
  
  ASMD = data.frame(variables = row.names(covariate_table), covariate_table$ASMD)
  names(ASMD)[2] = method
  
  data$trim_weights = trim_weights
  data$outcome_method = method
  data$normalize_method = normalize_method
  data$transport_method = transport_method
  
  return(list(ASMD = ASMD,
              covtab = covariate_table,
              data = as.data.frame(data)))
}

cov_asmds = function(data, trial, covariates, methods, svy_wt_var, normalize_method = "none",transport_method = "ps", trim_weights=FALSE,trim_pctile = .97){
  ASMDs = methods %>% 
    map(~covtab(data, trial, covariates, ., svy_wt_var,normalize_method,transport_method,trim_weights, trim_pctile)$ASMD) %>% 
    reduce(left_join, by = "variables") %>% 
    gather(method, asmd, !!methods) %>% 
    mutate(method = factor(method, levels=c("svy","transport","both")),
           normalize_svywt = normalize_method,
           trim_weights = trim_weights,
           transport_method = transport_method)

  return(ASMDs)
}

loveplot_svytransport = function(asmds,methods){
  d = asmds %>% filter(method %in% methods)
  
  loveplot = d %>% 
    ggplot() +
    geom_point(aes(x = as.numeric(asmd), y = variables, colour = method2), size = 3)+
    labs(x = "ASMD", y = "") +
    theme(axis.text.y = element_text(size=12))+
    geom_vline(xintercept = 0) +
    xlim(0,max(as.numeric(asmds$asmd)))#+
    #facet_grid(normalize_svywt~trim_weights)
  
  return(loveplot)
}

outcome = function(data, trial, covariates, outcome_method, svy_wt_var, outcome_var, treatment_var, normalize_method = "none",transport_method = "ps", trim_weights=FALSE,trim_pctile = .97){
  d = covtab(data, trial, covariates, outcome_method, svy_wt_var, normalize_method, transport_method,trim_weights,trim_pctile)$data
  
  if(length(table(d[,outcome_var]))!=2){
    ATE_design = svydesign(id = ~1, data = d %>% filter(get(trial) == 1), weights = d$w[which(d[,trial] == 1)])
    model = svyglm(as.formula(paste0(outcome_var,"~",treatment_var)), design = ATE_design, family='gaussian')
    PATE = summary(model)$coefficients[treatment_var,"Estimate"]
    CI = as.numeric(confint(model)[treatment_var,])
  } else{
    ATE_design = svydesign(id = ~1, data = d %>% filter(get(trial) == 1), weights = d$w[which(d[,trial] == 1)])
    model = svyglm(as.formula(paste0(outcome_var,"~",treatment_var)), design = ATE_design, family='quasibinomial')
    PATE = exp(summary(model)$coefficients[treatment_var,"Estimate"])
    CI = as.numeric(exp(confint(model)[treatment_var,]))
  }
  outcome_df = data.frame(PATE = PATE, CI_l = CI[1], CI_u = CI[2], 
                          method = outcome_method, 
                          normalized_svywt = normalize_method,
                          trim_weights=trim_weights,
                          transport_method = transport_method)
  return(outcome_df)
}



