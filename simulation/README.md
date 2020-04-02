# Simulation Code for "Generalizing Randomized Trial Findings to a Target Population using Complex Survey Population Data"

This repo contains all of the necessary code to run the simulations for the manuscript titled "Generalizing Randomized Trial Findings to a Target Population using Complex Survey Population Data." The code is set up to run on a remote cluster (originally run on the JHPCE cluster at Johns Hopkins).

### Setup

Before running the simulations on the cluster, make sure you have all of the files from this repo's `code/` directory in this `simulation/` folder, along with the following sub-directories set up:

- `pop_data/` (for the created population data files)
- `sim_results/` (for the created simulation results)
- `combined_results/` (for the combined results)

You can do this by running the following code in the command line:
```
mkdir pop_data
mkdir sim_results
mkdir combined_results
```

### config.R
The `config.R` file defines *all* of the parameters to run the simulations, including parameters you wish to vary. There are four major groups of objects to define in the config file:

- `N`: The size of the target population
- `pop_data_params`: This list specifies all parameters to create the population data, including
    - `x_cor`: the pairwise correlation between the 6 covariates
    - `scale_survey`: a scaling parameter dictating how much the covariates influence survey selection
    - `scale_rct`: a scaling parameter dictating how much the covariates influence trial selection
    - `coefs_survey`: coefficients for each covariate's relative influence on survey selection
    - `coefs_rct`: coefficients for each covariate's relative influence on trial selection	
    - `scale_tx_het`: a scaling parameter dictating the amount of treatment effect heterogeneity in the outcome due to the covariates
- `sample_data_params`: This list specifies the parameters to generate the trial and complex survey sample sizes
    - `prop_rct`: overall proportion of target population to sample for the trial
    - `prop_survey`: overall proportion of the target population to sample for the complex survey
- `fit_params`: parameters to specify what method is applied to estimate the population average treatment effect
	- `survey_weights`: logical indicator of whether or not to use the survey weights
	- `transport_weights`: logical indicator of whether or not to estimate transportability weights
	- `multiply_weights`: logical indicator to try multiplying survey weights and transportability weights
	- `weight_method`: specify estimation method for the sample membership model
	- `transport_formula`: specify which variables to omit when fitting the sample membership model
	- `bootstrap`: logical indicator of whether or not to obtain boostrapped variance estimate
	- `n_bootstrap`: number of bootstrap samples
   
### Running the simulations on the JHPCE cluster

1. Update and upload your `config.R` file to the cluster directory.
2. Create the population data files (this also checks if the population data already exist, in which case they're not re-generated)
```
qsub -N mkpop -cwd -l mem_free=1G,h_vmem=1G mkpop.sh 
```
3. Run the simulations (`sim.sh` submits 1000 jobs with 1000 different seeds, and each creates a results file in `sim_results/`):
```
sh sim.sh
```
4. Gather the simulation run times and combine the results across the 1000 runs (this is done in an active session. Make sure to navigate to the directory after typing `qrsh`):
```
qrsh
sh gather_results.sh
exit
```
