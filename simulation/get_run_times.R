library(readr)
library(dplyr)

read.csv("all.txt",header=FALSE, stringsAsFactors = FALSE) %>% 
  mutate(run_times = parse_number(V1)) %>% 
  dplyr::select(run_times) %>% 
  saveRDS(paste0("combined_results/",format(Sys.Date(),"%m%d%Y"),"_sim_run_times.rds"))

