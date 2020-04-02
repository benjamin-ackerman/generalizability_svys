for f in {10001..11000}
do
qsub -N survey_sim -cwd -l mem_free=1.5G,h_vmem=1.5G run.sh  $f
done