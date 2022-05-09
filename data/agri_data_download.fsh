#!/usr/bin/fish
# Shell script for downloading the climate and agriculture model datasets
# This uses the Fish shell (https://fishshell.com/)

# base data directory -- change this to whatever is appropriate!
set base_dir ~/CSH/HoloSim/data

# directory to save the datasets to (climate subdir is used)
set outdir $base_dir/climate
# ensure that the download directory exists
mkdir -p $outdir


# 1. Baseline climate dataset
for i in (seq 1980 2010)
wget -P $outdir https://data.giss.nasa.gov/impacts/agmipcf/agmerra/AgMERRA_"$i"_prate.nc4
wget -P $outdir https://data.giss.nasa.gov/impacts/agmipcf/agmerra/AgMERRA_"$i"_tavg.nc4
end


# 2. Past climate dataset

# which years to download -- adjust this to what's needed
# data used in the current study: 2.5k to 10k BP (550 BCE to 8050 BCE)
set years 2.5_5 5_7.5 7.5_10

for y in $years
wget -P $outdir https://dap.ceda.ac.uk/badc/deposited2018/HadCM3B_60Kyr_Climate/data/temp/bias_regrid_tas_"$y"kyr.nc
wget -P $outdir https://dap.ceda.ac.uk/badc/deposited2018/HadCM3B_60Kyr_Climate/data/precip/bias_regrid_pr_"$y"kyr.nc
end


# 3. Agricultural model
# (LPJmL, spring and winter wheat without adaptation)
wget -P $outdir https://zenodo.org/record/3592453/files/LPJmL_spring_wheat_ggcmi_phase2_emulator_A0.nc4
wget -P $outdir https://zenodo.org/record/3592453/files/LPJmL_winter_wheat_ggcmi_phase2_emulator_A0.nc4


# 4. Coastlines
set outdir $base_dir/gaez_crop
mkdir -p $outdir
wget -P $outdir https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_land.zip
wget -P $outdir https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/50m/physical/ne_50m_land.zip
unzip $outdir/ne_10m_land.zip -d $outdir
unzip $outdir/ne_50m_land.zip -d $outdir






