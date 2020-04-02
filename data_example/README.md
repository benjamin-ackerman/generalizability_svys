# Data Examples Code for "Generalizing Randomized Trial Findings to a Target Population using Complex Survey Population Data"

This repo contains all of the necessary code to replicate the data examples in the manuscript titled "Generalizing Randomized Trial Findings to a Target Population using Complex Survey Population Data" 

### R scripts

- `data_cleaning_nutrition.R` and `data_cleaning_SUD.R` contain all code for pre-processing the data from the four studies across the two data examples: PREMIER and NHANES (Nutrition) and CTN-0044 and NSDUH (SUD). This script also merges the respective two studies for each example, harmonizing and dichotomizing covariates for consistency.
- `analysis_functions.R` contains the functions used to run the analyses.
- `analysis_nutrition.R` and `analysis_sud.R` execute the functions in `analysis_functions.R` and create the tables and figures of results that are displayed in the manuscript.

### Data Availability Statement
Data from the PREMIER trial that support the findings of this study are available upon request through NHLBI. Restrictions apply to the availability of these data, which were used under license for this study. Data are available at https://biolincc.nhlbi.nih.gov/studies/premier/ with the permission of NHLBI.

Data from the CTN-0044 trial that support the findings of this study are available upon request through the NIDA CTN. Restrictions apply to the availability of these data, which were used under license for this study. Data are available at https://datashare.nida.nih.gov/ with the permission of the NIDA CTN.

Data from NHANES that support the findings of this study were derived from the following resources available in the public domain: https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx?BeginYear=2003

Data from NSDUH that support the findings of this study were derived from the following resources available in the public domain: https://www.datafiles.samhsa.gov/study/national-survey-drug-use-and-health-nsduh-2011-nid13563