# Neolithic population model simulation and scripts

This repository contains the simulation code and related scripts for our paper:

Dániel Kondor, James S. Bennett, Detlef Gronenborn, Nicolas Antunes, Daniel Hoyer, Peter Turchin (2022). <br>
Explaining population booms and busts in Mid-Holocene Europe. SocArXiv preprint: https://osf.io/preprints/socarxiv/c32up


## Content of this repository

This repository is organized into the following subdirectories:
 - [analysis](analysis): scripts for evaluating temporal patterns in the simulation results and in radiocarbon data
 - [data](data): scripts and instructions for getting the necessary data along with some auxilliary data files used
 - [neolithic_cpp](neolithic_cpp): simulation code, C++ version
 - [neolithic_py](neolithic_py): simulation code, Python version
 - [preprocessing](preprocessing): scripts doing preprocessing steps necessary before running the simulation


## Requirements

Scripts in this repository use R and Python 3, while the simulation code has two versions, in C++ and Python 3 (either version can be used independently of the other). Currently, code for spatial aggregation of simulation results, and a sampling of population is only available in the C++ version. Running the scripts and the simulation typically assumes a UNIX (e.g. Linux, OSX, WSL) environment, but all code should be easily adaptable to run on Windows as well. Having at least 8 GiB of RAM is recommended for both the preprocessing and the simulation steps.

All scripts assume that all data is stored within one base directory, according to the structure defined in the successive steps. All scripts contain a variable (`base_dir`) that is set in the beginning to the root of this directory. This should be adjusted to the actual data location, and the data files under the [data](data) directory of this repository should be copied there (of course, it is possible to use the data subdirectory as well).

Main dependencies:

 - [R](https://www.r-project.org/), tested on version 4.0.4, should work on any recent version
 - [Python 3](https://www.python.org/), tested on version 3.9.7, should work on any recent version
 - a C++17 capable compiler, tested on [GCC](https://www.gnu.org/software/gcc/), version 11.2.0 (for compiling the C++ version of the simulation)

R packages used:

 - ggplot2
 - ggmap
 - rworldmap
 - foreach
 - parallel
 - doParallel
 - [rcarbon](https://github.com/ahb108/rcarbon)
 - reshape2
 - ncdf4

All packages can be installed with the built-in package manager (the `install.packages()` function).

Python libraries used (beyond the standard library):

 - [pandas](https://pandas.pydata.org/)
 - [geopandas](https://geopandas.org/)
 - [numpy](https://numpy.org/)
 - [shapely](https://github.com/shapely/shapely)
 - [networkx](https://networkx.org/)
 - [tqdm](https://tqdm.github.io/)
 - [scipy](https://www.scipy.org)

All of these should be possible to install via `pip / pip3`.

Additional software dependencies:

 - [DGGRID](https://github.com/sahrk/DGGRID)
 - [Fish shell](https://github.com/fish-shell/fish-shell) (syntax used in some scripts; should be easy to adapt to other shells as well)
 - [Cairo](https://cairographics.org/) (optional, used for creating visualizations and videos of simulation outputs)
 - [ffmpeg](https://ffmpeg.org/) (optional, used in examples for creating videos of simulation outputs)


## Steps to run the simulation

1. Look at the [datasets used](data) and steps necessary to download and generate them.

2. Perform the necessary [preprocessing](preprocessing) steps.

3. Run the simulation; the [C++](neolithic_cpp) is suitable for reproducing the results in our paper.

4. Perform the [data analysis](analysis) steps on the simulation outputs and the radiocarbon data.



