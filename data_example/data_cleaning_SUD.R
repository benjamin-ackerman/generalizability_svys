library(stringr)
library(dplyr)
library(readr)

### Load NSDUH data
# Can be found here: https://www.datafiles.samhsa.gov/study-dataset/national-survey-drug-use-and-health-2011-nsduh-2011-ds0001-nid13707
# Survey weights/design variables: weight - ANALWT_C, strata - vestr, replicate - verep
load("data_example/data/NSDUH/NSDUH_2011.RData")

# Variable names for variables indicating past 30 day illicit substance use: (just use the variable "SUMMON")
#substance_use = names(PUF2011_090718)[str_detect(names(PUF2011_090718),"rec")]

## Limit the target population to 18 and older, and reported any illicit drug use in past 30 days
nsduh = PUF2011_090718 %>% 
  filter(AGE2 > 6, summon == 1) %>%  #Filter to 18+ who used any illicit drugs in the past 30 days
  dplyr::select(age_cat = AGE2, female = irsex, marital_status = irmarit, race_eth = NEWRACE2, 
                education = EDUCCAT2, employment = empstaty, prior_tx = txilalev,
                #sud_tx_last_year = illtrmt,
                iv_drug = NEDFLG01,
                svy_wt = ANALWT_C, strata = vestr, replication = verep) %>% 
  mutate(age_cat = cut(age_cat, breaks = c(0,12,14,15,17), labels = c("18-25","26-34","35-49","50+")),
         female = female - 1,
         race_eth = cut(race_eth, breaks = c(0,1,2,6,7), labels = c("White","Black","Other","Hispanic")),
         education = cut(education, breaks = c(0, 1, 2, 4), labels = c("Less than high school","High school", "College or more")),
         employment = cut(employment, breaks = c(0,1,2,3,4), labels = c("Full-time","Part-time","Unemployed","Other")),
         marital_status = cut(marital_status, breaks = c(0,1,3,4), labels = c("Married","Separated/divorced/widowed","Single, never married")))

### Load CTN-44 data:
ctn_path = "data_example/data/CTN_44/data/"

# Get age, sex and race/ethnicity from demographics
ctn_demographics = read_csv(paste0(ctn_path,"DEMOGRAPHICS/CSV/CTN44_DEM.csv")) %>% 
  mutate(age_cat = cut(age, breaks = c(0, 25, 34, 49, 10000), labels = c("18-25","26-34","35-49","50+")),
         female = ifelse(sex == "Female",1, ifelse(sex == "Male", 0, NA)),
         race_eth = ifelse(ethnic == "Hispanic or Latino","Hispanic",
                           ifelse(amindian == "Yes" | asian == "Yes" | hawaiian == "Yes" | raceoth == "Yes", "Other",
                                  ifelse(white == "Yes" & black == "Yes", "Other",
                                         ifelse(white == "Yes", "White",
                                                ifelse(black == "Yes", "Black", NA)))))) %>% 
  dplyr::select(usubjid, arm, female, age_cat, race_eth)

# Get education, marital status and employment from screening
ctn_screening = read_csv(paste0(ctn_path,"ENROLLMENT-RANDOMIZATION/CSV/CTN44_SCRN.csv")) %>% 
  mutate(education = ifelse(sclvledu == "High school not complete", "Less than high school",
                            ifelse(sclvledu == "High school graduate or GED", "High school","College or more")),
         marital_status = ifelse(scmarita == "Married/remarried","Married",scmarita),
         employment = ifelse(str_detect(scemploy, "Full time"), "Full-time",
                             ifelse(str_detect(scemploy, "Part time"), "Part-time",
                                    ifelse(scemploy == "Unemployed","Unemployed","Other")))) %>% 
  dplyr::select(usubjid, arm, education, marital_status, employment)

# Past IV drug use
ctn_past_drugs = read.csv(paste0(ctn_path,"DRUG USE-SELF REPORT/CSV/CTN44_TFB.csv"))
use_routes = names(ctn_past_drugs)[str_detect(names(ctn_past_drugs),"rou")]

ctn_iv_ever = ctn_past_drugs %>% 
  filter(assessdays <= 0) %>% 
  filter_at(vars(use_routes), any_vars(. == "IV Injection")) %>% 
  group_by(usubjid) %>% 
  summarise(iv_drug = 1)

# Past treatment
ctn_prior_tx = read_csv(paste0(ctn_path,"ADHERENCE-COMPLIANCE/CSV/CTN44_NMS.csv")) %>% 
  filter(visno == "00") %>% 
  dplyr::select(usubjid, arm, prior_tx = nmdrugal) %>% 
  mutate(prior_tx = ifelse(prior_tx == "Yes", 1, 0))
  

## Outcomes: Abstinence at 12 weeks (shown to be significant in study)
####### My analysis #####
ctn_urine = read_csv(paste0(ctn_path,"DRUG USE-BIOLOGIC/CSV/CTN44_UDS.csv"))

maxday = ctn_urine %>% 
  filter(visno %in% c("23","24"),!is.na(assessdays)) %>%
  group_by(usubjid,arm) %>% 
  mutate(nrows = n()) %>% 
  ungroup() %>% 
  mutate(assessdays = ifelse(visno == "23" & nrows == 1, assessdays + 4, assessdays)) %>% 
  group_by(usubjid,arm) %>% 
  summarise(maxday = max(assessdays, na.rm=TRUE)) %>%
  ungroup() %>%
  dplyr::select(usubjid, arm, maxday)

ctn_urine=ctn_urine %>% 
  left_join(maxday, by = c("usubjid","arm")) %>%
  group_by(usubjid) %>% 
  filter(visno %in% c("23","24"))
  
positive_urine = ctn_urine %>% 
  group_by(usubjid,arm) %>% 
  filter_at(vars(BZD:BARB), any_vars(. == "Positive")) %>% 
  summarise(pos_test = 1)

negative_urine = ctn_urine %>% 
  group_by(usubjid,arm) %>% 
  filter_at(vars(BZD:BARB), all_vars(. == "Negative")) %>% 
  ungroup() %>% 
  filter(!usubjid %in% positive_urine$usubjid) %>% 
  group_by(usubjid,arm) %>% 
  summarise(pos_test = 0)

# Allow urine test to be missing
ctn_urine = positive_urine %>% 
  bind_rows(negative_urine)

# Self-reported abstinence:
ctn_selfreport = read.csv(paste0(ctn_path,"DRUG USE-SELF REPORT/CSV/CTN44_TFB.csv")) %>%
  left_join(maxday, by = c("usubjid","arm")) %>% 
  # Select last four visits from maxvisit
  filter(assessdays > maxday-4, assessdays <= maxday,!is.na(tfsubalc)) %>% 
  ungroup() %>% 
  # Merge gender variable
  left_join(ctn_demographics %>% dplyr::select(usubjid, female), by = "usubjid") %>% 
  # Recode alcohol use and recode any substance use/alcohol use based on it
  mutate(tfalcohl2 = ifelse(tfalcohl == "No", as.character(tfalcohl), 
                            ifelse(female == 0 & tfnmdrnk < 5, "No",
                            ifelse(female == 1 & tfnmdrnk < 4, "No",as.character(tfalcohl)))),
         tfsubalc2 = ifelse(tfsubalc == 0, tfsubalc, 
                            ifelse(tfcannab == "Yes"| tfcocain== "Yes"| tfamphet== "Yes"| 
                                     tfmetamp== "Yes"|tfoxycod== "Yes"|tfmethad== "Yes"|tfopiate== "Yes"|
                                     tfecstas== "Yes"|tfbarbit== "Yes"|tfbenzod== "Yes"|tfothdrg== "Yes"|
                                     tfalcohl2== "Yes", 1, 0))) %>% 
  group_by(usubjid, arm) %>% 
  # Self report is 1 if any reported use in those four days, 0 otherwise
  summarise(self_report = ifelse(sum(tfsubalc2, na.rm=TRUE)>=1, 1, 0))

# Merge the CTN data together
ctn = ctn_demographics %>% 
  full_join(ctn_screening, by = c("usubjid","arm")) %>% 
  left_join(ctn_iv_ever, by = c("usubjid")) %>% 
  left_join(ctn_prior_tx, by = c("usubjid","arm")) %>% 
  left_join(ctn_selfreport, by = c("usubjid","arm")) %>% 
  left_join(ctn_urine, by = c("usubjid","arm")) %>% 
  filter(!is.na(arm)) %>% 
  ungroup() %>% 
  mutate(iv_drug = ifelse(is.na(iv_drug), 0, iv_drug),
         treatment = ifelse(str_detect(arm,"TES"), 1, 0),
         # Following guidelines, if any are missing, abstinence is 0, otherwise if no positive test AND no self reported drug use, abstinence is 1
         abstinence = ifelse(is.na(pos_test) & is.na(self_report), NA,
                             ifelse(compareNA(pos_test,0) & compareNA(self_report,0), 1, 0))) %>% 
  dplyr::select(-arm,-usubjid)#, -pos_test, -self_report)

### MERGE TRIAL AND SURVEY TOGETHER
data = ctn %>% 
  mutate(S = 1) %>% 
  bind_rows(nsduh %>% 
              mutate(S = 0))

### Save the data
saveRDS(data,paste0("data_example/data/",Sys.Date(),"_merged_data_sud.rds"))
