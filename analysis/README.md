## Data analysis steps

This directory contains a set of scripts to analyze the simulation results and the radiocarbon dataset. These scripts assume that the simulations have been run according to the script [../neolithic_cpp/simulation_runs.fsh](../neolithic_cpp/simulation_runs.fsh); they can be adapted to other parameter combinations by changing the input files and parameter ranges.

### 1. SPD of radiocarbon data

The script [c14_data_spd.r](c14_data_spd.r) performs an aggregation of the radiocarbon dataset into the overlapping tilings and calculated summed probability distributions (SPDs) in each tile for the study period.

Inputs:
```
population/2021-09-08-dataset-nodup-no-sd-filter.csv
new_grids_all.RData
```

Output:
```
population/c14_spd_new_10000_5000.RData
```

Runtime: up to a few hours; multiple CPU cores are used if available.


### 2. ACF of radiocarbon data

The script [c14_data_acf.r](c14_data_acf.r) calculates autocorrelation function (ACFs) of the SPDs calculated in the previous step after detrending (using either an exponential or logistic function). It identifies the first minimum in the ACFs as the typical scale of periodic patterns. Additionally, it calculates the coefficient of variation (CV) in each case as well. ACF minimum locations and CV values are saved into CSV files to be used as "baseline" distributions later.

Input: results from the previous step

Outputs:
```
population/acf_result/acf_peaks_min_{exp,log}_r{5,7,10}.csv
population/acf_result/all_res_cv_{exp,log}_r{5,7,10}.csv
```
Additionally, figures in the same directory summarizing the results.

Runtime: a few minutes


### 3. ACFs of simulation results

The script [acf_res2.r](acf_res2.r) analysis the regionally aggregated population numbers of the simulation results directly, calculating ACFs and creating a set of plots of the distribution of ACF minima.

Inputs: <br> results from the previous step
```
simulation/run*spaggr.out.gz
```

Outputs: 
```
simulation/run*acf.png
simulation/acf_cmb_{c,n}.png
```
Plots of ACF minimum distributions.

Runtime: a few minutes


### 4. Sampled simulation results

The script [sim_dates.r](sim_dates.r) reads the simulation results after a sampling process (as done by the code [sample_res.cpp](../neolithic_cpp/sample_res.cpp)), adds noise and simulates the aggregation and analysis steps as done for the radiocarbon dataset. The results are aggregate distributions of ACF minima and CV values (as in the previous step) along with plots of individual example time series and ACFs.

Inputs: <br> results from processing the C14 data
```
simulation/sample_pop*.dat
```

Output: plots of ACF and CV distributions after the sampling and aggregation process.

Runtime: up to 6 hours (when using 4 cores for the main analysis)





