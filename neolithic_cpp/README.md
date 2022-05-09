# Neolithic simulation

Agent-based simulation of interactions among groups of farmers and raiders in a neolithic setting, C++ version.

## Examples for the simulation

Examples of a few cases are given below. A more complete set of examples that are used for the main analysis of our paper are given in the script [simulation_runs.fsh] -- it is recommended to run that script in full before progressing to the data processing steps.

The below commands are provided to be run in a standard UNIX (Linux, OSX) shell (such as `bash`). They require a `C++17` compiler such as a recent version of GCC or CLang (to use CLang, replace the `g++` commands with `clang++`). The code for creating images and videos (`spout2png.cpp`) uses the [Cairo](https://www.cairographics.org/) graphics libraries. [FFMPEG](https://ffmpeg.org/) is used for creating videos. Running the below code will require at least 8GB free memory.

All code is assumed to be run from the current directory, and assumes a shell variable `$base_dir` set to the base directory where all data files are stored after the preprocessing steps. Files needed for the simulation and results are also written here, under the subdirectory `simulation`; if this does not exist, it should be created before running any of the below commands:
```
mkdir -p $base_dir/simulation
```


### 1. Compile the simulation code

Manual compilation from the command line with GCC:
```
g++ -o cm create_helper_matrix.cpp read_table_cpp.cpp -std=gnu++17 -O3 -march=native -lm -fopenmp
g++ -o n2 neolithic2.cpp read_table_cpp.cpp -std=gnu++17 -O3 -march=native -lm -fopenmp -pthread
g++ -o n2w neolithic2w.cpp read_table_cpp.cpp -std=gnu++17 -O3 -march=native -lm -fopenmp -pthread
g++ -o spas spaggr_simple.cpp read_table_cpp.cpp -O3 -march=native
g++ -o sr sample_res.cpp read_table_cpp.cpp -O3 -march=native -lm
g++ -o sp2png spout2png.cpp read_table_cpp.cpp cellnet_img.cpp -O3 -march=native -std=gnu++17 -I/usr/include/cairo -I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include -I/usr/include/pixman-1 -I/usr/include/uuid -I/usr/include/freetype2 -I/usr/include/libpng16 -lcairo
```

Notes:
 - The `-fopenmp` switch is optional and will enable parallelization of some parts of the code; it requires support for OpenMP.
 - Compiling the `sp2png` program (the last command) above requires options locating the Cairo graphics libraries. The development headers for this needs to be installed, e.g. on Ubuntu, the `libcairo2-dev` package. The command line above includes default locations on typical Linux distributions; this likely requires adjustments on other systems. On many systems, the `pkg-config --libs --cflags cairo` command can be used in determining the options needed.
 - If there are compilation issues, the [Python version](../neolithic_py) of the simulation code might be easier to use, especially for the video creation code.


Alternatively, the [Meson](https://mesonbuild.com/) build system can be used to compile the codes:
```
meson -Dbuildtype=release build
ninja -C build
```
This will place the executables under the newly created `build` subdirectory. Note that the following examples assumes the output of manual compilation, i.e. that all executables are under the current directory. Copy them or create symbolic links as follows:
```
ln -s build/cm
ln -s build/n2
ln -s build/n2w
ln -s build/spas
ln -s build/sr
ln -s build/sp2png
```


### 2. Create the base probability distributions for the migrations

These are used by the simulation later. The main parameter is given with the `-G` switch that determines the characteristic distance in kilometers (the parameter of the exponential function that gives the decrease of migration probability with distance). In the current exapmles, G = 40 km is used for the simulation variant without conflict, and G = 80 km is used for the case with conflict.


```
./cm -o $base_dir/simulation/matrix_dg_G40zC10.bin -G 40 -k -K 150 -Kt 50 -cc $base_dir/gaez_crop/eu_dggrid_coords_K_land4a.csv -ct -n $base_dir/dggrid/isea3h12eun.nbr -nr -nC 0.1 -nE $base_dir/eu_dggrid_extra_edges.dat -Ks $base_dir/gaez_crop/eu_dggrid_land_share.csv -z -t 4
./cm -o $base_dir/simulation/matrix_dg_G80zC10.bin -G 80 -k -K 150 -Kt 50 -cc $base_dir/gaez_crop/eu_dggrid_coords_K_land4a.csv -ct -n $base_dir/dggrid/isea3h12eun.nbr -nr -nC 0.1 -nE $base_dir/eu_dggrid_extra_edges.dat -Ks $base_dir/gaez_crop/eu_dggrid_land_share.csv -z -t 4
```

These result in output files with an approximate size of 6.2 GB (assuming a study area with 36400 cells)



### 3. Run the simulation

Two examples are given here, for the cases of parameter combinations in the main text of our paper.


3.1. Model without conflict (i.e. p_E = 0), G = 40 km, s = 4, time range: 7000 BCE to 2000 BCE:
```
./n2 -G 40 -H2 -Hf0 $base_dir/simulation/matrix_dg_G40zC10.bin -p -d 1 -s 1 -S 5000 -o $base_dir/simulation/runG40C10ws4sp.out.gz -oz -op 1 -i 1607919 -I 100000 -If 0.5 -Ip 1 -Ce 0 -Cb 0 -De 1 -Dm 0.5 -Db 0.25 -Dl 200 -T -t 4 -w $base_dir/climate/cell_ids_dggrid.csv $base_dir/climate/crop_data_cmb.dat -wm -wC -wf -wF -7000 -ws 4 > $base_dir/simulation/runG40C10ws4.out
```

3.2. Model with conflict, G = 80 km, p_E = 0.05, p_A = 0.2, s = 4:
```
./n2w -Hf $base_dir/simulation/matrix_dg_G80zC10.bin -Rr -Rp 1 -RP -RA 0.2 -R 8 -m -Rs 0 -d 1 -s 1 -S 5000 -E 20 -a 1.0 -o $base_dir/simulation/runmr_G80_C10_E20_A5_s4sp.out.gz -oz -op 1 -i 1596250 -I 5000 -If 0.5 -Ip 0.5 -Ce 0 -Cb 0 -De 1 -Dm 0.5 -Db 0.25 -Dl 200 -w $base_dir/climate/cell_ids_dggrid.csv $base_dir/climate/crop_data_cmb.dat -wm -wC -wf -wF -7000 -ws 4 > $base_dir/simulation/runmr_G80_C10_E20_A5_s4.out
```


### 4. Spatial aggregation

This step aggregates population numbers in the tiling that was created for spatial aggregation in the preprocessing step (examples are given for the above two cases of the simulation:
```
zcat $base_dir/simulation/runG40C10ws4sp.out.gz | ./spas -m $base_dir/new_grid_aggr_flat.csv -f $base_dir/new_dgid_aggr_filter.dat | gzip -c > $base_dir/simulation/runG40C10ws4spaggr.out.gz
zcat $base_dir/simulation/runmr_G80_C10_E20_A5_s4sp.out.gz | ./spas -m $base_dir/new_grid_aggr_flat.csv -f $base_dir/new_dgid_aggr_filter.dat | gzip -c > $base_dir/simulation/runmr_G80_C10_E20_A5_s4spaggr.out.gz
```



### 5. Compile videos

This uses the `sp2png` program that has two modes of operation. The first one outputs a series of PNG images that can be later used by any video editing software to create video files:
```
mkdir $base_dir/simulation/video_example
zcat $base_dir/simulation/runG40C10ws4sp.out.gz | ./sp2png -oy -7000 -op 10 -om 230 -ow 1608 -oh 900 -o $base_dir/simulation/video_example/runG40C10ws4 $base_dir/gaez_crop/eu_dggrid_coords_K_land4a.csv -n $base_dir/dggrid/isea3h12eun.nbr -cp $base_dir/dggrid/dggrid_poly.csv
```
Note that the input is always read from the standard input; this should be redirected from the previously created output file. This will create a series of PNG images under `$base_dir/simulation/video_example` that can be used e.g. as an input to `ffmpeg` or similar software.

An alternative way is to run ffmpeg directly while generating the images. This can be achieved with the `-of` option that will try to run ffmpeg with the correct options and the given file as the final output, supplying a set of images to encode as a video directly:

```
zcat $base_dir/simulation/runG40C10ws4sp.out.gz | ./sp2png -oy -7000 -op 10 -om 230 -ow 1608 -oh 900 -of $base_dir/simulation/runG40C10ws4.mkv -oc h264 $base_dir/gaez_crop/eu_dggrid_coords_K_land4a.csv -n $base_dir/dggrid/isea3h12eun.nbr -cp $base_dir/dggrid/dggrid_poly.csv
```
Note that the output container format will depend on the filename given (multiple formats are supported by ffmpeg; the above example uses the MKV container format). The actual video codec can be selected with the `-oc` parameter; in the above example, H264 is used which works reasonably fast.

Second example for the simulation version with conflict:
```
zcat $base_dir/simulation/runmr_G80_C10_E20_A5_s4sp.out.gz | ./sp2png -oy -7000 -op 10 -om 200 -ow 1608 -oh 900 -of $base_dir/simulation/runmr_G80_C10_E20_A5_s4.mkv -oc h264 $base_dir/gaez_crop/eu_dggrid_coords_K_land4a.csv -n $base_dir/dggrid/isea3h12eun.nbr -cp $base_dir/dggrid/dggrid_poly.csv
```



## Program parameters

A more detailed description of the command line arguments of the individual programs.

## cm / create_helper_matrix.cpp

This program creates a binary file that contains the cumulative distribution function of migration probabilities to be used during the simulations. Essentially, the distributions form an NxN matrix (for N cells); given that any one probability can be calculated on the fly, it is sufficient to store NxN/2 values. Optionally, a matrix of distances can be stored as well; this is useful if non-trivial adjustments such as scaling factors or long-range connections are used.

Command line parameters used:

 - `-o fn` output file name
 - `-K num` scale carrying capacities so that the average capacity of a cell is the value given here (recommended to use, since the input files contain arbitrary values)
 - `-Kg` additionally scale the capacity in each cell by the cosine of the latitude (to account for cells whose size depends on the latitude -- this is not used with DGGRID, where size of cells is constant)
 - `-Ks fn` read additional scaling factors for each cell from the given file (this is used for taking into account the share of actual usable soil in each cell); this option can be used multiple times, resulting in the product of all scaling factors used
 - `-c fn` or `-cc fn` name of the file with a list of cell IDs, coordinates and (optionally) base carrying capacities; by default, this file is expected in a TSV format without header; use the `-cc` form if the input is a CSV with header
 - `-k` has to be given if the above file contains carrying capacities (e.g. GAEZ base estimates); otherwise a constant value is assumed (note: this is used in all cases)
 - `-n fn` or `-ne fn` name of the file containing the neighbor relations among cells; use `-n` if this file is in the DGGRID output format (all neighbors of a cell are listed on a line), and the `-ne` variant if this file is an edgelist (each line contains two node IDs); the input file should be symmetric, i.e. it should contain each edge both ways
 - `-nC num` scaling factor for travel distance between coastal cells; calculated distance is multiplied by this, so giving a number < 1 allows faster / further travel along coastal cells
 - `-nE fn` name of a file containing additional edges (as an edgelist) that should be added to the cell neighbor relations; this is currently used to add links between non-adjacent cells
 - `-nr` if given, distances are calculated along neighbor edges taking into account any scaling factor 
 - `-z` calculate distances along shortest paths in the cell network, taking into account any scaling, and store a matrix of such distances in the output
 - `-C` instead of creating a binary file, read one and check that the probability distributions are correct in it
 - `-G dist` characteristic distance used in the migration probabilities (exponential decay)
 - `-P exp` if given, a power-law is used for distance dependence with the exponent given here instead of an exponential
 - `-t num` number of threads to use for distance and probability calculations

## Main simulation: n2, n2w (neolithic2.cpp, neolithic2w.cpp)

### Common options:

These program runs the simulation in one realization; `n2` does not include conflicts, `n2w` allows the creation of aggressors and conflicts (the main difference is that migration probabilities are kept track of in different ways).

Initial conditions and basic settings:
 - `-i ID` ID of the cell where to place a starting population (a location in Anatolia is used in the examples: cells 1596250 and 1607919 for the cases with and without conflict respectively)
 - `-I pop` size of the starting population to use (default: 100); if this is larger than the carrying capacity of the starting cell, the additional population is distributed in neighboring cells
 - `-If factor` when distributing the initial population, cells are only filled up to this factor (e.g. 0.5); this option is useful to avoid an initial "explosion" of migrations that would happen if the initial conditions include cells close to their carrying capacity
 - `-Ip factor` when distributing the initial population, any neighbor cell is chosen with this probability (e.g. 0.5); this can be used to start from a less tightly packed configuration of occupied cells
 - `-S n` number of steps (years) to run the simulation for
 - `-K num` scale carrying capacities so that the average over all cells is this value (can be adjusted regardless of the value given when creating the migration distributions, since this means scaling the whole distribution with the same number)
 - `-T` if given, the average distance of newly settled cells (from the starting cell) is tracked and output
 - `-o fn` detailed output filename; the population in each cell is written here every 10 years
 - `-oz` if present, the above output file is compressed with `/bin/gzip`
 - `-s num` number to use to seed the random number generator (current time is used by default)
 - `-t num` number of parallel OpenMP threads to use for evaluating population growth and split-offs (default: 1, i.e. no parallelism as it matters very little in this case)
 - `-r r` base population growth rate (per year)
 - `-d r` population collapse factor (if population is above carrying capacity, it collapses to this ratio of it (values below one can model a larger collapse)

Farmers' dynamics:
 - `-Dl num` Dunbar number to use (default: 150)
 - `-Db p` probability of a group splitting off if the population of a cell is equal to the Dunbar number
 - `-Dm r` minimum ratio of population to the Dunbar number to consider a possibility of split-off
 - `-De e` exponent used to calculate the probability of split-off due to approaching the Dunbar limit
 - `-Cb p` probability of a group splitting off if the population of a cell is equal to the currnt local carrying capacity (note: this is set to zero in all simulations)
 - `-Cm r` minimum ratio of population to the carrying capacity to consider a possibility of split-off
 - `-Ce e` exponent used to calculate the probability of split-off due to approaching the carrying capacity limit

Note: each year, the possiblity of a group splitting off (and migrating elsewhere) is considered for each cell. This is done by first estimating a probability for this to happen based on the Dunbar limit and carrying capacity limit. In each case, this is done in the following way:
```
P = p * ((x - r) / (1 - r))**e
```
where `x` is the current population relative to the limit considered (either the Dunbar number, given by the `-Dl` parameter, or the local carrying capacity), and `p`, `r` and `e` are the parameters given above. By default, both mechanisms are considered, and a migration occurs if either of them yields it by random choice (given the probability computed as above). Setting the `-Db` or `-Cb` value to zero disables that mechanism.

Climate data (for variation in agriculture):

 - `-w fn1 fn2` specifies the names of two files with the data to use; `fn1` should be the mapping between weather cells and simulation cells (i.e. `cell_ids_dggrid.csv` in this case), `fn2` should be the file with the actual variations (`crop_data_cmb.dat` in this case)
 - `-wm` if given, missing cells are allowed; without this, every cell has to appear in all years in the input file
 - `-wc` if given, cut the simulation area to cells that have climatic variation data
 - `-wf` if given, the mapping between weather cells and simulation cells contains weighting factors; this means that multiple weather cells overlap with one simulation cell and a weighted average is used to calculate any actual effect of climate variation; this option is necessary when using the hexagon grid (maybe this should be made the default)
 - `-wC` if given the mapping file is in CSV format (this is necessary for the above files)
 - `-wF year` if given, assume that the start year of the simulation is the one given here (as an absolute number, in calendar years, i.e. BCE / CE); data before this date in the input file is skipped; if not given, the first year is the input file is assumed to be the beginning of the simulation
 - `-ws factor` scale any variation in yield by this factor; this can be used to exaggerate the effect of climatic variability
 


### n2 (neolithic2.cpp)

Command line parameters specific to the simulation version without conflict:
 - `-Hf fn` name of the binary file containing the pre-computed probability distributions (created with `cm`)
 - `-Hf0 fn` name of the binary file as before, but do not actually use the distribution in it, only the network and distance matrix (this is recommended since re-creating the distributions as needed can be more efficient in this case)
 - `-H2` use a more space-efficient version of storing the migration probabilities where probabilities can be recalculated as needed


### n2w (neolithic2w.cpp)

Command line parameters specific to the simulation version with conflict:

 - `-Hf fn` name of the binary file containing the pre-computed probability distributions (created with `cm`)
 - `-R dist` characteristic distance in the probability of choosing a target for aggressors (this can be varied independently from the migration probability distribution parameter)
 - `-Rc` if a cell is succesfully defended from an attack, it turns into aggressors with a probability of 50% (not used in the main simulation)
 - `-Rr` if given, aggressors cannot revert back to farmers (they die out instead if they cannot find a target and win; this leaves their cell empty; this is used in the simulation for the main results)
 - `-RP` if given, aggressors take into account the current farmer population when choosing a target (they prefer cells with larger populations)
 - `-Rm dist` maximum distance that aggressors are able to attack (default is 10x the distance given by the `-R` option)
 - `-m` aggressors are considered mobile: after each successfull attack, they take over the attacked cell
 - `-a prob` probability of an attack (by aggressors) being successful (e.g. setting this to 1 means that they always win; this is used in the default version of the simulation)
 - `-Ra dist` if given, the probability of successful attack (by aggressors) depends on the distance exponentially, using the characteristic distance given here (i.e. as `a * exp(- d / dist)`, where `a` is the base success rate given by the `-a` parameter, `d` is the distance of raiders to the target cell, and `dist` is the paremeter value given here; not used in the main analysis)
 - `-E num` preference to choose an empty cell for group migrations (typical value: 10, default value is 1); i.e. a larger value will result in less aggressors to be created until the landscape is more saturated; note: this is the inverse of the p_E parameter in the paper



## sp2pgm / spout2pgm.cpp

This is a simple program that converts the detailed output of the simulation to [PGM](https://en.wikipedia.org/wiki/Netpbm#File_formats) images that can be used to construct a video. Of course, this requires a rectangular grid, with IDs assigned as done for the GAEZ-based grid. Each cell is represented by a single pixel in the output (scaling up can be done later with image / video processing tools). Values are read from the standard input, and the filenames are generated by appending the year to a base name given as a command line parameter.

Parameters:

 - `-x num` size of the grid on the x axis (note: cell IDs are expected to be `ID = y * xsize + x`, so it is important to use the correct size)
 - `-y num` size of the grid on the y axis (need to be given as the output is allocated in advance)
 - `-X num` maximum x dimensions of output; if given, the output is cut to this size from the right (otherwise, the full grid is output)
 - `-Y num` maximum y dimansions of output; if given, the output is cut to this size from the top (otherwise, the full grid is output)
 - `-o fn` base filename of the output; the year is appended to this
 - `-m num` maximum allowed value in the input; values larger than this will cause errors, unless the `-s` option is used (see below)
 - `-M num` maximum value used in the output, typically 255 (a value larger than 255 might cause issues with some image viewers)
 - `-s` saturate input values above the given maximum


## sp2png / spout2png.cpp

This program converts the detailed output of the simulation run among arbitrary regions to a series of PNG images that can be used to create a video. It also supports directly running ffmpeg to create a video output file.

Parameters:

 - `-ow num` width of the output images (note that images are scaled in a way so that the whole area fits in the requested size)
 - `-oh num` height of the output images
 - `-op num` output images in this intervals (e.g. 10 means every 10 years)
 - `-c fn` or `-cc fn` name of the file with a list of cell IDs, coordinates and (optionally) base carrying capacities; by default, this file is expected in a TSV format without header; use the `-cc` form if the input is a CSV with header
 - `-n fn` or `-ne fn` name of the file containing the neighbor relations among cells; use `-n` if this file is in the DGGRID output format (all neighbors of a cell are listed on a line), and the `-ne` variant if this file is an edgelist (each line contains two node IDs); the input file should be symmetric, i.e. it should contain an edge both ways
 - `-cp fn` name of the file containing the actual cell polygons in a simple CSV format (note: needs to be ordered in a way that corresponds to how they should be displayed)
 - `-om num` maximum population value expected in the input (needs to be given since the data is not read in advance)
 - `-og num` gamma correction to apply to the output color scale (essentially scales output values when selecting the color of cells)
 - `-o fn` base output filename when creating a set of PNG images; image numbers are appended based on time step
 - `-of fn` output video filename; if this option is used, this program will try to run `ffmpeg` directly to create a video
 - `-oc codec` change the video codec used (only matters when using `ffmpeg`, i.e. the `-of` option); this should be something that `ffmpeg` supports, the default is `vp9` (setting this to `h264` can often result in faster encoding); additional options can also be given as part of this
 - `oq num` change the video encoding quality (only matters when using `ffmpeg`, i.e. the `-of` option), `num` should be a small integer (default: 15); lower numbers result in better quality and larger file sizes; setting a too low quality (a high number) will result in noticable issues, e.g. blurry or choppy videos
 

Currently, the color scale used is hardcoded (goes from dark to light blue, bandits are displayed in red), but this should be made configurable in the future. If the `-of` option is used, no PNG images are used. The `ffmpeg` executable should be in PATH to run it successfully. It is given the following options:
```
ffmpeg -y -f rawvideo -pix_fmt bgra -s %ux%u -r 24 -i - -codec:v vp9 -crf 15 %s
```
(where the size and the output filename are filled in based on the options)


