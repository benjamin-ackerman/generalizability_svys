library(cowplot)
library(WeightIt)
library(caret)
library(generalize)
library(dplyr)
library(purrr)
library(stringr)
library(rebus)
source("data_example/code/analysis_functions.R")

data = readRDS("data_example/data/2019-12-03_nutrition_data.rds") %>% 
  filter(is.na(VISIT) | VISIT == "6mo") %>% 
  mutate(svy_weights = ifelse(is.na(svy_weights),1,svy_weights)) %>%
  filter(!is.na(sex),!is.na(age_cat),!is.na(BMI),!is.na(black),!is.na(education),!is.na(SBP),!is.na(DBP)) %>% 
  as.data.frame()

covariates = c("sex","age_cat","BMI","black","education")

########## Calculate ASMDs ##########
params = list(normalize_method = c("none","sum_wt","trial_mean_wt","mean_wt"),
              trim_weights = c(TRUE,FALSE),
              transport_method = c("ps","gbm","super_luedtke","super_moodie")) %>%
  cross()

asmds = seq_along(params) %>%
  map_df(~cov_asmds(data,"S",covariates,c("svy","transport","both"),"svy_wt",
                    params[[.]]$normalize_method,
                    trim_weights = params[[.]]$trim_weights,
                    transport_method = params[[.]]$transport_method)) %>%
  mutate(normalize_svywt = factor(normalize_svywt, levels = c("none","sum_wt","trial_mean_wt","mean_wt")),
         transport_method = factor(transport_method, levels = c("ps","gbm","super_luedtke","super_moodie")),
         trim_weights = ifelse(trim_weights,"Trimmed Weights","Untrimmed Weights"))

saveRDS(asmds, "asmds_nutrition.rds")
print("Nutrition ASMDs done!")

########## OUTCOME ANALYSIS ###########
params = list(method = c("naive","transport","both"),
              normalize_method = c("none","sum_wt","trial_mean_wt","mean_wt"),
              trim_weights = c(TRUE,FALSE),
              transport_method = c("ps","gbm","super_luedtke","super_moodie")) %>%
  cross()

results = seq_along(params) %>%
  map_df(~outcome(data,"S",covariates,params[[.]]$method,"svy_wt","SBP_CHANGE","T",
                  params[[.]]$normalize_method,trim_weights = params[[.]]$trim_weights,
                  transport_method = params[[.]]$transport_method))

saveRDS(results, "results_nutrition.rds")
print("Nutrition results done!")


###################### FIGURES FOR PAPER: ######################
asmds = readRDS("data_example/code/asmds_nutrition.rds")
asmds2 = cov_asmds(data,"S",covariates,"naive","svy_weights") %>% 
  mutate(transport_method = "gbm", method = "naive", trim_weights = "Untrimmed Weights") %>% 
  bind_rows(asmds) %>% 
  filter(normalize_svywt == "none",trim_weights == "Untrimmed Weights",
         transport_method %in% c("gbm")) %>% 
  mutate(variables = factor(variables, levels = c("sex","age_cat_1","age_cat_2","age_cat_3",
                                                  "age_cat_4","age_cat_5","age_cat_6",
                                                  "BMI","black","education")),
         variables = fct_recode(variables, `Male` = "sex", `Age 25-40` = "age_cat_1", `Age 41-45` = "age_cat_2",
                                `Age 46-50` = "age_cat_3", `Age 51-55` = "age_cat_4", `Age 56-60` = "age_cat_5",
                                `Age 60+` = "age_cat_6", `BMI` = "BMI", 
                                 `Black` = "black",`College or more` = "education"),
         method = factor(method, levels = c("naive","svy","transport","both")),
         method2 = fct_recode(method, `Naive` = "naive", `Naive \n(vs. pop)`="svy",`Transported` = "transport",
                              `Transported + \nSvy Weights` = "both"))

myColors <- brewer.pal(4,"Set1")
names(myColors) <- levels(asmds2$method2)
colScale <- scale_colour_manual(name = "Method",values = myColors)

loveplot_manuscript = loveplot_svytransport(asmds2,c("svy","transport","both")) +
  scale_y_discrete(limits = rev(levels(asmds2$variables))) +
  colScale + 
  theme_light()

### Panel the covariate table and loveplot
covtab_unweighted = covtab(data, "S", covariates, "naive", "svy_weights",trim_weights=FALSE)$covtab
covtab_svy_wtd = covtab(data, "S", covariates, "svy", "svy_weights",trim_weights=FALSE)$covtab

tab = covtab_unweighted %>% 
  dplyr::select(-ASMD) %>% 
  left_join(covtab_svy_wtd %>% dplyr::select(-ASMD), by = "trial")

tab[-8,] = tab[-8,]*100

rownames(tab) = rownames(covtab_svy_wtd) =  c("Male (%)","Age 25-40 (%)","Age 41-45 (%)","Age 46-50 (%)",
                                              "Age 51-55 (%)","Age 56-60 (%)","Age 60+ (%)" ,"BMI (mean)",
                                              "Black (%)","College or more (%)")
colnames(tab)[1] = "PREMIER \n (trial)"
colnames(tab)[2] = "NHANES \n (survey)"
colnames(tab)[3] = "Population \n (weighted survey)"

tt <- ttheme_default(base_family ="Arial")
tbl <- tableGrob(round(tab,1), theme=tt)

# Plot chart and table into one object
obj = plot_grid(ggdraw() + draw_grob(tbl),loveplot_manuscript,labels = c('A', 'B'), label_size = 12,
                rel_heights = c(5,1))

ggsave("data_example/figures/manuscript_nutrition_covariates.png",
       obj,
       width=11,height=4,units = 'in')

### Make results plot
results = readRDS("data_example/code/results_nutrition.rds")

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
  labs(x = "Estimator",y = "PATE: Change in Systolic Blood Pressure")+
  colScale +
  theme_light()

ggsave("data_example/figures/manuscript_nutrition_results.png",
       results_plot)
