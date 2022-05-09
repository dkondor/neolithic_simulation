## Data sources used

This directory contains some auxilliary datasets used in the simulation and analysis, along with the instructions to download and generate the additional datasets needed for the simulation.

All analysis steps (i.e. all scripts in this repository) assume a structure where all data is stored under a base directory that is set in the `base_dir` variable at the beginning of each script. This should be adjusted to the actual location that is suitable for data storage. The contents of this directory should be copied there as a first step. Typically, at least 100 GiB free space is needed for downloading all necessary data sets, performing the preprocessing and preparing and running the simulation. Having the data files on an SSD is recommended.

### 1. Simulation space

[DGGRID](https://github.com/sahrk/DGGRID) is used for creating the simulation space. The script [create_dggrid.fsh](create_dggrid.fsh) contains the example commands to create this (assuming that DGGRID has been compiled from source with the binary under `$dggrid_bin_dir`). This creates a subdirectory `dggrid` under the base data directory for the output.

Alternatively, the resulting grid can be downloaded from [here](https://www.dropbox.com/s/tnsbiar9ha0vqkm/isea3h12eu.zip?dl=0); extract the contents to a `dggrid` subdirectory under the base data directory, e.g. using the following commands (run from the data base directory):
```
mkdir dggrid
cd dggrid
wget https://www.dropbox.com/s/tnsbiar9ha0vqkm/isea3h12eu.zip
unzip isea3h12eu.zip
```

### 2. Climate and land area datasets

These include a past climate dataset from (Armstrong et al., 2019), the crop yield emulators from (Franke et al., 2020) and the baseline climate dataset (NASA). These can be downloaded in an automated way by the commands in the script [agri_data_download.fsh](agri_data_download.fsh) (or manually from the URLs in it). All of these should be saved in a subdirectory `climate`.

The download script also include shapefiles with the land areas of Earth from the [Natural Earth](https://www.naturalearthdata.com/) dataset. These should be saved in a subdirectory `gaez_crop` and extracted there.

References:
 - Edward Armstrong, Peter O. Hopcroft, and Paul J. Valdes (2019). <br> A simulated Northern Hemisphere terrestrial climate dataset for the past 60,000 years. Scientific Data, 6. doi: [10.1038/s41597-019-0277-1](http://dx.doi.org/10.1038/s41597-019-0277-1).
 - James A. Franke, Christoph Müller, Joshua Elliott, Alex C. Ruane, Jonas Jägermeyr, Juraj Balkovic, Philippe Ciais, Marie Dury, Pete D. Falloon, Christian Folberth, Louis François, Tobias Hank, Munir Hoffmann, R. Cesar Izaurralde, Ingrid Jacquemin, Curtis Jones, Nikolay Khabarov, Marian Koch, Michelle Li, Wenfeng Liu, Stefan Olin, Meridel Phillips, Thomas A.M. Pugh, Ashwan Reddy, Xuhui Wang, Karina Williams, Florian Zabel, and Elisabeth J. Moyer (2020). <br> The GGCMI Phase 2 emulators: global gridded crop model responses to changes in CO 2 , temperature, water, and nitrogen (version 1.0). Geoscientific Model Development, 13(5):3995–4018. doi: [10.5194/gmd-13-3995-2020](https://doi.org/10.5194/gmd-13-3995-2020).
 - NASA AgMERRA dataset: https://data.giss.nasa.gov/impacts/agmipcf/


### 3. Agricultural productivity data (GAEZ)

Data can be accessed (after registration which is free and available to anyone) at https://www.gaez.iiasa.ac.at/. Data files need to be downloaded manually, using the selections in the menu system.

All data should be saved in a subdirectory called `gaez_crop` (in the data base directory).

Data files used in the current simulation:

#### 3.1. Yield estimates

Select the following main category:
```Suitability and Potential Yield -> Agro-climatic yield -> Agro-climatically attainable yield```

Additional selections (using the menu on the left):
```
Crop: Winter wheat
Water supply: rain-fed
Input level: low
Time period: Baseline (1961-990)
```

Select the `Visualization and download` option from the menu and the option `Standard bundle archive as ZIP file`. Save the resulting file as `gaez_yield_winter_wheat_low.zip` (note: no need to extract it).

Do the same for spring wheat, by changing the `Crop` selection, and save it as `gaez_yield_spring_wheat_low.zip`


#### 3.2. Slope dataset

Share of land with slope above 30%. This is available under:
```
Land Resources -> Terrain Resources -> Terrain slope > 30%
```
Again, this should be saved using the `Standard bundle archive as ZIP file` option in the `Visualization and download` menu, under the file name `gaez_slope_30p.zip`.


#### 3.3. Area covered by water

Share of grid cells covered by water, available as:

```
Land Resources -> Land cover -> Water bodies
```

Save this file as `gaez_water.zip`.


### 4. Radiocarbon data

Most of the data is available using the [c14bazAAR](https://github.com/ropensci/c14bazAAR) R package. The data file used in our analysis is available by request from us. It should be saved in a subdirectory `population`.





