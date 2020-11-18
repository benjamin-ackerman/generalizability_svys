library(cowplot)
library(RColorBrewer)
library(caret)
library(generalize)
library(dplyr)
library(purrr)
library(stringr)
library(rebus)
library(WeightIt)
library(forcats)
source("data_example/code/analysis_functions.R")

data = readRDS("data_example/data/2019-12-04_merged_data_sud.rds") %>% 
  filter(!is.na(female),!is.na(race_eth)) %>%
  as.data.frame()

covariates = c("female","age_cat","race_eth","education","marital_status","employment","iv_drug","prior_tx")

covtab_unweighted = covtab(data, "S", covariates, "naive", "svy_wt",trim_weights=FALSE)$covtab
covtab_svy_wtd = covtab(data, "S", covariates, "svy", "svy_wt",trim_weights=FALSE)$covtab
covtab_transport = covtab(data, "S", covariates, "transport", "svy_wt",trim_weights=FALSE)$covtab
covtab_transport_svy_wtd = covtab(data, "S", covariates, "both", "svy_wt",trim_weights=FALSE)$covtab

########## CHECK BALANCE ##############
params = list(normalize_method = c("none","sum_wt","trial_mean_wt","mean_wt"),
              trim_weights = c(TRUE,FALSE)) %>% 
  cross()

# ASMD boxplots
asmds = seq_along(params) %>% 
  map_df(~cov_asmds(data,"S",covariates,c("svy","transport","both"),"svy_wt",params[[.]]$normalize_method,trim_weights = params[[.]]$trim_weights)) %>%
  mutate(normalize_svywt = factor(normalize_svywt, levels = c("none","sum_wt","trial_mean_wt","mean_wt")),
         trim_weights = ifelse(trim_weights,"Trimmed Weights","Untrimmed Weights")) 

ggsave("data_example/figures/loveplots_sud.png",loveplot_svytransport(asmds), width=11,height=7,units = 'in')

ggsave("data_example/figures/asmds_sud_allmethods.png",
       ggplot(asmds) +
         geom_boxplot(aes(x = method, y = as.numeric(asmd),fill = normalize_svywt))+
         labs(x = "Method", y = "ASMD") +
         facet_grid(.~trim_weights), 
       width=11,height=7,units = 'in')

# Read in ASMDs from cluster
asmds = readRDS("data_example/code/asmds_sud.rds")
asmds2 = cov_asmds(data,"S",covariates,"naive","svy_wt") %>% 
  mutate(transport_method = "gbm", method = "naive", trim_weights = "Untrimmed Weights") %>% 
  bind_rows(asmds) %>% 
  filter(normalize_svywt == "none",trim_weights == "Untrimmed Weights",
         transport_method %in% c("gbm")) %>% 
  mutate(variables = fct_recode(variables, `Age 18-25` = "age_cat_18-25", `Age 26-34` = "age_cat_26-34",
                                `Age 35-49` = "age_cat_35-49", `Age 50+` = "age_cat_50-", `Less than HS` = "educationLess-than-high-school",
                                `College or more` = "educationCollege-or-more", `High School` = "educationHigh-school",
                                `Employed full time` = "employmentFull-time",`Other Job` = "employmentOther", `Employed part-time` = "employmentPart-time",
                                `Unemployed` = 'employmentUnemployed', `Female` = "female", `Past IV drug use` = "iv_drug", `Married` = "marital_statusMarried",
                                `Separated/Divorced` = "marital_statusSeparated-divorced-widowed",`Single` = "marital_statusSingle--never-married",
                                `Prior SUD tx` = "prior_tx",`Black` = "race_ethBlack",`Hispanic` = "race_ethHispanic",`Other Race/Ethnicity` = "race_ethOther",
                                `White` = "race_ethWhite"),
         method = factor(method, levels = c("naive","svy","transport","both")),
         method2 = fct_recode(method, `Naive` = "naive", `Naive \n(vs. pop)`="svy",`Transported` = "transport",
                             `Transported + \nSvy Weights` = "both"))

myColors <- brewer.pal(4,"Set1")
names(myColors) <- levels(asmds2$method2)
colScale <- scale_colour_manual(name = "Method",values = myColors)

ggsave("data_example/figures/asmds_sud_manuscript.png",
  loveplot_svytransport(asmds2,c("svy","transport","both")) +
  scale_y_discrete(limits = rev(levels(asmds2$variables))) +
  colScale,
  width=11,height=7,units='in')


asmds = asmds %>% 
  mutate(method = fct_recode(method,`Unweighted`="svy",`Transported`="transport",`Transported + \nSvy Weights` = "both"))
ggsave("data_example/figures/loveplot_sud_jobtalk.png", 
       loveplot_svytransport(asmds %>% 
                        filter(normalize_svywt == "none",trim_weights == "Untrimmed Weights",
                               transport_method %in% c("gbm"))),
       width=11,height=7,units='in')


ggsave("data_example/figures/asmds_sud_jobtalk.png",
       asmds %>%
         filter(normalize_svywt == "none",trim_weights == "Untrimmed Weights",
                transport_method %in% c("gbm")) %>% 
       ggplot() +
         geom_boxplot(aes(x = method, y = as.numeric(asmd),fill=method))+
         labs(x = "Method", y = "ASMD"),
       width=7,height=7,units = 'in')

# Weights Distribution by Normalization
params = list(normalize_method = c("none","sum_wt","trial_mean_wt","mean_wt"),
              trim_weights = c(TRUE,FALSE),
              outcome_method = "both") %>% 
  cross()

d = seq_along(params) %>% 
  map_df(~covtab(data, "S", covariates, params[[.]]$outcome_method, "svy_wt", 
                 params[[.]]$normalize_method, "ps",params[[.]]$trim_weights)$data)

# Weights Distribution by Transport Method
params = list(trim_weights = c(TRUE,FALSE),
              outcome_method = "both",
              transport_method = c("ps","gbm","super_luedtke","super_moodie")) %>% 
  cross()

d = seq_along(params) %>% 
  map_df(~covtab(data, "S", covariates, params[[.]]$outcome_method, "svy_wt", 
                 normalize_method = "none", transport_method = params[[.]]$transport_method,
                 params[[.]]$trim_weights)$data)

ggsave("data_example/figures/weights_sud.png",
       d %>% 
         filter(S == 1) %>% 
         ggplot() +
         geom_histogram(aes(w)) +
         facet_wrap(trim_weights ~ normalize_method,scales="free",nrow=2),
       width=11,height=7,units = 'in') 

## Weights by method from cluster
d = readRDS("data_example/code/weights_sud.rds")
ggsave("data_example/figures/weights_sud_methods.png",
       d %>% 
         filter(S == 1) %>% 
         ggplot() +
         geom_histogram(aes(w)) +
         facet_wrap(trim_weights ~ transport_method,scales="free",nrow=2),
       width=11,height=7,units = 'in') 

########## OUTCOME ANALYSIS ###########
params = list(method = c("naive","transport","both"),
              normalize_method = "none",
              trim_weights = FALSE,
              transport_method = "gbm") %>% 
              #normalize_method = c("none","sum_wt","trial_mean_wt","mean_wt"),
              #trim_weights = c(TRUE,FALSE),
              #transport_method = c("ps","gbm","super_luedtke","super_moodie")) %>%
  cross()

results = seq_along(params) %>%
  map_df(~outcome(data,"S",covariates,params[[.]]$method,"svy_wt","abstinence","treatment",
                  params[[.]]$normalize_method,trim_weights = params[[.]]$trim_weights,
                  transport_method = params[[.]]$transport_method))

saveRDS(results, "results_sud.rds")
print("SUD results done!")


### Read results from the cluster:
results = readRDS("data_example/code/results_sud.rds")
ggsave("data_example/figures/PATE_sud_jobtalk.png", 
       results %>% 
         mutate(method = factor(method,levels = c("naive","transport","both")),#, labels = c("Unweighted","Transported","Transported + \nSvy Weights")),
                method2 = fct_recode(method, `Naive \n(vs. pop)`="naive",`Transported` = "transport",
                                     `Transported + \nSvy Weights` = "both"),
                normalized_svywt = factor(normalized_svywt, levels = c("none","sum_wt","trial_mean_wt","mean_wt")),
                trim_weights = ifelse(trim_weights,"Trimmed Weights","Untrimmed Weights")) %>%
         filter(normalized_svywt == "none",
                trim_weights == "Untrimmed Weights",
                transport_method == "gbm") %>% 
         ggplot() +
         geom_point(aes(x = method2, y = PATE, colour = method2), position = position_dodge(0.3),size=3)+
         geom_errorbar(aes(x = method2, ymin = CI_l,ymax=CI_u, colour = method2), position = position_dodge(0.3),width = .2) +
         geom_hline(yintercept = 1, col = 'black',lty=2)+
         labs(x = "Estimator",y = "Odds Ratio of Substance Abstinence")+
         theme(axis.text.y = element_text(size=12),
               axis.text.x = element_text(size=12))+
         colScale,
       width=11,height=7,units = 'in')

mutate(method = fct_recode(method,`Unweighted`="svy",`Transported`="transport",`Transported + \nSvy Weights` = "both"))



###################### FIGURES FOR PAPER: ######################
asmds = readRDS("data_example/code/asmds_sud.rds")
asmds2 = cov_asmds(data,"S",covariates,"naive","svy_wt") %>% 
  mutate(transport_method = "gbm", method = "naive", trim_weights = "Untrimmed Weights") %>% 
  bind_rows(asmds) %>% 
  filter(normalize_svywt == "none",trim_weights == "Untrimmed Weights",
         transport_method %in% c("gbm")) %>% 
  mutate(variables = factor(variables, levels = levels(asmds$variables)[c(12,1:4,18:21,5:7,14:16,8:11,13,17)]),
         variables = fct_recode(variables, `Age 18-25` = "age_cat_18-25", `Age 26-34` = "age_cat_26-34",
                                `Age 35-49` = "age_cat_35-49", `Age 50+` = "age_cat_50-", `Less than HS` = "educationLess-than-high-school",
                                `College or more` = "educationCollege-or-more", `High School` = "educationHigh-school",
                                `Employed full time` = "employmentFull-time",`Other Job` = "employmentOther", `Employed part-time` = "employmentPart-time",
                                `Unemployed` = 'employmentUnemployed', `Female` = "female", `Past IV drug use` = "iv_drug", `Married` = "marital_statusMarried",
                                `Separated/Divorced` = "marital_statusSeparated-divorced-widowed",`Single` = "marital_statusSingle--never-married",
                                `Prior SUD tx` = "prior_tx",`Black` = "race_ethBlack",`Hispanic` = "race_ethHispanic",`Other Race/Ethnicity` = "race_ethOther",
                                `White` = "race_ethWhite"),
         method = factor(method, levels = c("naive","svy","transport","both")),
         method2 = fct_recode(method, `Naive` = "naive", `Naive \n(vs. pop)`="svy",`Transported` = "transport",
                              `Transported + \nSvy Weights` = "both"))

myColors <- brewer.pal(4,"Set1")
names(myColors) <- levels(asmds2$method2)
colScale <- scale_colour_manual(name = "Method",values = myColors)

loveplot_manuscript = loveplot_svytransport(asmds2,c("svy","transport","both")) +
  scale_y_discrete(limits = rev(levels(asmds2$variables))) +
  colScale + 
  theme_light() #+
  #theme(legend.position = "none")

##### Make table and asmd plot
tab = covtab_unweighted %>% 
  dplyr::select(-ASMD) %>% 
  left_join(covtab_svy_wtd %>% dplyr::select(-ASMD), by = "trial")

tab = tab*100

rownames(tab) = rownames(covtab_svy_wtd) =  paste0(c("Female","Age 18-25","Age 26-34","Age 35-49","Age 50+","Black","Hispanic","Other Race","White","College or More","High School","Less than HS","Married","Separated, Divorced, Widowed","Single","Employed Full-time","Employed Other","Employed Part-time","Unemployed","Past IV drug use","Prior SUD treatment")," (%)")
colnames(tab)[1] = "CTN-0044 \n (trial)"
colnames(tab)[2] = "NSDUH \n (survey)"
colnames(tab)[3] = "Population \n (weighted survey)"

tt <- ttheme_default(base_family ="Arial")
tbl <- tableGrob(round(tab,1), theme=tt)

obj = plot_grid(ggdraw() + draw_grob(tbl),loveplot_manuscript,labels = c('A', 'B'), label_size = 12,
                rel_heights = c(5,1))

ggsave("data_example/figures/manuscript_sud_covariates.png",
       obj,
       width=13,height=7,units = 'in')

### Make result plot
results = readRDS("data_example/code/results_sud.rds")

results_plot = results %>% 
  mutate(method = factor(method,levels = c("naive","transport","both")),#, labels = c("Unweighted","Transported","Transported + \nSvy Weights")),
         method2 = fct_recode(method, `Naive \n(vs. pop)`="naive",`Transported` = "transport",
                              `Transported + \nSvy Weights` = "both"),
         normalized_svywt = factor(normalized_svywt, levels = c("none","sum_wt","trial_mean_wt","mean_wt")),
         trim_weights = ifelse(trim_weights,"Trimmed Weights","Untrimmed Weights")) %>%
  filter(normalized_svywt == "none",
         trim_weights == "Untrimmed Weights",
         transport_method == "gbm") %>% 
  ggplot() +
  geom_point(aes(x = method2, y = PATE, colour = method2), position = position_dodge(0.3),size=3)+
  geom_errorbar(aes(x = method2, ymin = CI_l,ymax=CI_u, colour = method2), position = position_dodge(0.3),width = .2) +
  geom_hline(yintercept = 0, col = 'black',lty=2)+
  labs(x = "Estimator",y = "Risk Difference of Substance Abstinence")+
  # theme(axis.text.y = element_text(size=12),
  #       axis.text.x = element_text(size=12))+
  colScale +
  theme_light()

ggsave("data_example/figures/manuscript_SUD_results.png",
       results_plot)
