## Data analysis steps

This directory contains a set of scripts to analyze the simulation results and the radiocarbon dataset. These scripts assume that the simulations have been run according to the script [../neolithic_cpp/simulation_runs.fsh](../neolithic_cpp/simulation_runs.fsh); they can be adapted to other parameter combinations by changing the input files and parameter ranges.


### 1. Preprocessing of radiocarbon data

The script [c14_data_preprocess_new.r](c14_data_preprocess_new.r) performs a binning on the radiocarbon dataset using the [rcarbon](https://github.com/ahb108/rcarbon) package as a preprocessing step.

Inputs:
```
population/2021-09-08-dataset-nodup-no-sd-filter.csv
```

Output:
```
population/c14_data_bins1.csv
```

Runtime: less than a minute.


### 2. SPD and ACF of radiocarbon data

The script [c14_model_new.py](c14_model_new.py) performs the main processing necessary on the radiocarbon dataset. This includes aggregating in space (into overlapping grid cells), calculation of SPDs, fitting of a logistic trend, detrending by generating 1,000 synthetic datasets in each region and identification of ACF minima. The output also includes coefficients of variation and p-values that characterize the deviation from the logistic trend in each region.

This step uses the output of the previous step, along with the grids created during the preprocessing.

Inputs:
```
population/c14_data_bins1.csv
new_grids_all.csv
```

Outputs:
```
population/res_min.csv
population/res_pval.csv
population/res_cv.csv
```

Runtime: up to 9-10 hours.


### 3. Plot of radiocarbon data results

The script [c14_model_plots1.r](c14_model_plots1.r) reads the results of the previous steps and produces plots of regional time series, ACFs and distribution of ACF minima, CV values and p-values among the regions.

Inputs: output of the previous step

Outputs: panels for Fig. 2 in the main text and SI Figs. S6-S9.


### 4. ACFs of simulation results

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



