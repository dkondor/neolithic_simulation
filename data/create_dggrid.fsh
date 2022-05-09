#!/usr/bin/fish
# Shell commands to run DGGRID to create the hexagon grid used in the simulations
# This script uses the syntax of the Fish shell (https://fishshell.com/)
# (the main difference is in setting the shell variables)

# base data directory -- change this to whatever is appropriate!
set base_dir ~/CSH/HoloSim/data
# note: the base directory is assumed to contain the input files:
#  dggrid_eu.meta -- main file describing the grid creation parameters
#  eu_boundary2.{cpg,dbf,prj,shp,shx} -- bounding box or the study area (grid is generated inside this box)

# directory where DGGRID (binaries) are installed -- change this to whatever is appropriate!
set dggrid_bin_dir ~/CSH/HoloSim/DGGRID/build/src/apps/dggrid
# note: download DGGRID from here: https://github.com/sahrk/DGGRID
# and follow installation instructions here: https://github.com/sahrk/DGGRID/blob/master/INSTALL.md


# 1. create the output directory
mkdir -p $base_dir/dggrid

# 2. switch to the base data directory -- this is important since the metafile contains relative paths
cd $base_dir

# 3. run DGGRID
$dggrid_bin_dir/dggrid dggrid_eu.meta



 
