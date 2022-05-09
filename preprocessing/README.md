## Preprocessing steps

This directory contains a set of scripts that should be run to preprocess data necessary for the simulation. Before running, each script should be edited and the `base_dir` variable set to base directory for the data files (see [here](../data)). Scripts can be run in full or step by step in an interactive environment. Scripts should be run in the order given here. Note: runtimes are given as informative and correspond to a mid-range laptop from 2021.

### 1. Preprocess climate datasets

This is done by [climate_aggregate.r](climate_aggregate.r), which includes extracting the data used in the next steps and doing some preliminary aggregation.

Inputs:
```
climate/AgMERRA_*_prate.nc4
climate/AgMERRA_*_tavg.nc4
climate/bias_regrid_{tas,pr}_*.nc
```

Outputs:
```
climate/climate_processed_2.5_5kyr.RData
climate/climate_processed_5_7.5kyr.RData
climate/climate_processed_7.5_10kyr.RData
climate/AgMERRA_months_avg_eu_df.RData
```

Runtime: ~ 1 day


### 2. Estimate relative variation in crop yields over the study period

This is done by [agri_model_estimate.r](agri_model_estimate.r).

Inputs: <br> results from the previous step
```
climate/LPJmL_spring_wheat_ggcmi_phase2_emulator_A0.nc4
climate/LPJmL_winter_wheat_ggcmi_phase2_emulator_A0.nc4
```

Outputs:
```
climate/crop_data_cmb.{dat, RData}
```

Runtime: several hours


### 3. Preprocess GAEZ datasets

This is done by [gaez_convert.py](gaez_convert.py), converting the GAEZ datasets into a more convenient format, giving estimates of potential yield and usable land share in the rectangular grid used by GAEZ.

Inputs:
```
gaez_crop/gaez_yield_spring_wheat_low.zip
gaez_crop/gaez_yield_winter_wheat_low.zip
gaez_crop/gaez_water.zip
gaez_crop/gaez_slope_30p.zip
gaez_crop/ne_10m_land.shp
```

Outputs:
```
climate/eu_combined_estimate.csv
climate/eu_land_share.csv
```

Runtime: ~1-2 hours


### 4. Convert results into the hexagon cells from DGGRID

This is done by [dggrid_process.py](dggrid_process.py), which calculates overlaps between the rectangular grid used by GAEZ and the climate data and the hexagon cells, and calculates weighted averages of the relevant variables in the DGGRID cells. The study area is then further cut to the shape established manually and included among the data files.

Inputs: <br> results from the previous steps
```
gaez_crop/ne_50m_land.shp
neolithic_4000bc_cut_ne.shp
```

Outputs:
```
gaez_crop/eu_dggrid_coords_K_land4a.csv
gaez_crop/eu_dggrid_land_share.csv
climate/cell_ids_dggrid.csv
```

Runtime: ~1 hour


### 5. Create a series of regions used for spatial aggregation

This step creates a set of quasi-rectangular tilings of the study area that is used later for aggregating the simulation results and the radiocarbon dates. This is done by [sp_aggr_new.r](sp_aggr_new.r).

Inputs: none

Outputs:
```
new_grids_all.{RData,csv}
grid_map{1,2}.png
```

Note: the PNG images are examples of the tilings, created only if the variable `do_plot` is set to `TRUE`.


### 6. Match spatial aggregation tiles with DGGRID hexagons

This is necessary for processing simulation outputs later, done by [region_dggrid_match.py](region_dggrid_match.py).

Inputs:
```
gaez_crop/eu_dggrid_coords_K_land4a.csv
climate/crop_data_cmb.dat
new_grids_all.csv
```

Outputs:
```
new_grid_aggr.csv
new_grid_ids.csv
new_grid_aggr_flat.csv
new_dgid_aggr_filter.dat
dggrid/dggrid_poly.csv
```

Runtime: a few minutes









