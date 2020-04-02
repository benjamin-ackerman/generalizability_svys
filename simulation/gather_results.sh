module load conda_R
cat survey_sim.o* >> all.txt
Rscript get_run_times.R
rm all.txt
rm survey_sim.o*
rm survey_sim.e*

Rscript gather_results.R
