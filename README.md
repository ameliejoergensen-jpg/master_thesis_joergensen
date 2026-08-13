# master_thesis_joergensen
All code required to reproduce the results of the Master's thesis submitted by Amelie Jörgensen in partial fulfilment of the requirements for the degree M.Sc. Statistics.


The folder simulation_code contains all code required to run the simulation.
- git_master_functions.R is the source script with all the functions
- git_master_store_scenarios_complex.R stores the emission and dwell time specification for the simulation scenarios
- git_master_make_scenarios_complex.R sets the remaining simulation settings and creates master_scenarios_complex.csv
- master_scenarios_complex.csv gives an overview of the simulation scenarios
- git_run_array.R (scenario 1-4, 6-9) and git_run_array_spike.R (scenario 5, 10) simulate data according to the scenario settings, fit the 5 models, and store the results
- collect_runs.R collects the files produced by the run_array scripts into a single RDA per simulation scenario
- submit_collect.sh submits collect_runs.R on an HPC
- submit_array.sh submits git_run_array.R (scenario 1-4, 6-9) on an HPC
- submit_array_spike.sh submits git_run_array_spike.R (scenario 5, 10) on an HPC
- submit.sh, submit_spike.sh submit both the submit_array files and the collect job on an HPC
- git_master_final_evaluation.Rmd is the script for evaluating the simulation results and producing the final figures

The folder simulation_results_rda contains the results of the 10 simulation scenarios (100 runs each).

The folder actigraph_code contains the scripts for analysis of the dataset DEPRESJON (see https://dl.acm.org/doi/10.1145/3204949.3208125).
- git_master_act_functions.R is the source script, a lot of overlap with git_master_functions.R but adjusted to dataset
- git_master_actigraph.Rmd contains the traditional metrics, sensitivity analysis, fitting of Markov models on entire dataset, plots
- git_leave_one_out_classification.Rmd contains the leave-one-out cross-validation on the DEPRESJON dataset
