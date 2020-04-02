### Setup ###
# Load Packages
library(haven);library(dplyr);library(stringr);library(RNHANES)

# Path to PREMIER data
premier_path = "data_example/data/PREMIER/SASv8_FILES/"

########## NHANES ##########
# Here is the link to the NHANES data description from 2003-2004 (year that PREMIER was conducted)
# https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx?BeginYear=2003

# Download dietary data and demographics from 2003-2004
nhanes_diet = nhanes_load_data("DR1TOT_C", "2003-2004", demographics = TRUE)
# weight variables: WTINT2YR - full sample 2 year interview weight, WTDRD1 - Dietary day one sample weight
nhanes_bmi = nhanes_load_data("BMX_C", "2003-2004", demographics = TRUE)

# blood pressure
nhanes_bp = nhanes_load_data("BPX_C", "2003-2004")

### Just select what I need from NHANES and recode the variables:
  ### ID, age, sex, BMI, race (black/no), education, survey weights, urinary sodium?
### PROBLEM: how can we simultaneously check the positivity assumption (i.e. restrict the target population) and retain survey weights??
nhanes = nhanes_bmi %>% 
  left_join(nhanes_bp, by = "SEQN") %>% 
  dplyr::select(ID = SEQN, age = RIDAGEYR, sex = RIAGENDR, race = RIDRETH2, education = DMDEDUC,svy_weights = WTINT2YR,BMI = BMXBMI, SBP = BPXSY1, DBP = BPXDI1) %>% 
  filter(age >= 25, BMI >= 18.5 & BMI <= 45, SBP >= 120|DBP >= 80) %>% 
  mutate(sex = abs(sex - 2),
         black = ifelse(race == 2, 1, 0),
         age_cat = cut(age, breaks = c(0, 40, 45,50,55,60,85), labels = 1:6),
         education = ifelse(education %in% c(1,2), 0, ifelse(education == 3, 1, NA)),
         DBP = ifelse(DBP < 40, NA, DBP),
         SBP = ifelse(SBP < 70, NA, SBP)) %>% 
  dplyr::select(-race,-age)

########## PREMIER ##########
# Note that original TX variable was: 1 - control, 2 - Tx, 3 - DASH + Tx
# Demographics
premier_dem = read_sas(paste0(premier_path,"masterp.sas7bdat"))

# Education variable
premier_ed = read_sas(paste0(premier_path,"pthistp.sas7bdat")) %>% 
  mutate(more_than_hs = ifelse(EDU_REL == 1, 0, 1))

# Urinary sodium (True outcome - FOR TESTING METHOD ONLY)
premier_urine = read_sas(paste0(premier_path,"labp.sas7bdat"))

# Self-reported sodium (Mismeasured outcome)
premier_selfreport = read_sas(paste0(premier_path,"foodanalp.sas7bdat"))

## Blood pressure outcomes
premier_outcomes = read_sas(paste0(premier_path,"outcomep.sas7bdat")) %>% 
  filter(VISIT %in% c(8,10)) %>% 
  mutate(VISIT = ifelse(VISIT == 8, "6mo","18mo"))

### Merge PREMIER self-report and sodium data
premier_sodium = premier_selfreport %>% 
  left_join(premier_urine, by = c("STUDY_ID","VISIT","TX","COHORT")) %>% 
  dplyr::select(STUDY_ID, VISIT, TX, COHORT, sodium_self = RINA, sodium_urine = MGURNA) %>% 
  mutate(log_sodium_self = log(sodium_self), log_sodium_urine = log(sodium_urine)) %>% 
  filter(VISIT == 8) %>%
  mutate(VISIT = "6mo") %>% 
  bind_rows(premier_selfreport %>%
              left_join(premier_urine, by = c("STUDY_ID","VISIT","TX","COHORT")) %>% 
              dplyr::select(STUDY_ID, VISIT, TX, COHORT, sodium_self = RINA, sodium_urine = MGURNA) %>% 
              mutate(log_sodium_self = log(sodium_self), log_sodium_urine = log(sodium_urine)) %>% 
              filter(VISIT == 10) %>% 
              mutate(VISIT = "18mo"))

### Merge all of the PREMIER data together
premier = premier_dem %>% 
  left_join(premier_ed, by = c("STUDY_ID", "COHORT","TX")) %>% 
  left_join(premier_sodium, by = c("STUDY_ID","COHORT","TX")) %>%
  left_join(premier_outcomes, by = c("STUDY_ID","COHORT","TX", "VISIT")) %>%
  mutate(ID = str_replace(STUDY_ID, "PREMIER",""),
         S = 1, T = ifelse(TX == 1, 0, 1),
         BMI = WEIGHT_BASE/((HEIGHT/100)^2),
         sex = abs(SEX - 2),
         age_cat = forcats::fct_collapse(as.factor(AGE_REL),
                                         `1` = c("1","2","3"),
                                         `2` = "4",
                                         `3` = "5",
                                         `4` = "6",
                                         `5` = "7",
                                         `6` = c("8","9","10"))) %>% 
  dplyr::select(VISIT, ID, S, T, sex, age_cat, BMI, black = RACE_REL, education = more_than_hs, SBP=SBP_BASE, DBP=DBP_BASE,
                SBP_CHANGE)#, 
         #log_sodium_self, log_sodium_urine, SBP_BASE:WT_CHANGE)

### Combined PREMIER and OPEN data, and refactor age ###
data = premier %>%
  bind_rows(nhanes %>% 
              mutate(ID = as.character(ID),
                     VISIT = NA,
                     SBP_CHANGE = NA,
                     T = NA,
                     S = 0))

### Save merged data ###
saveRDS(data,paste0("data_example/data/",Sys.Date(),"_merged_data_nutrition.rds"))
