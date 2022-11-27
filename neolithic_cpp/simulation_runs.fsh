#!/usr/bin/fish

# simulation_runs.fsh -- run the simulations for a set of parameter values (examples for videos and regional time series)
# Notes:
#  - this script uses the syntax of the fish shell (https://github.com/fish-shell/fish-shell or https://fishshell.com/)
#  - instead of running the whole script, it can make sense to run each section separately as needed by copying it to a shell window



# base data directory -- change this to whatever is appropriate!
set base_dir ~/CSH/HoloSim/data

# directory of the simulation code (by default it is the current directory) -- set this to whatever is appropriate
set program_dir .


# create the data subdirectory for the simulations
mkdir -p $base_dir/simulation


###########################################################################################################
# 1. Create the matrices representing the simulation space for the cases used in this script
# note: each matrix is ~6.2 GiB in size -- the total size of four matrices is ~25 GiB
# note: the last parameter ('-t 4') determines how many CPU cores should be used -- adjust this as appropriate
for G1 in 20 30 40 80
$program_dir/cm -o $base_dir/simulation/matrix_dg_G$G1"z"C10.bin -k -K 150 -Kt 50 -cc $base_dir/gaez_crop/eu_dggrid_coords_K_land4a.csv -ct -n $base_dir/dggrid/isea3h12eun.nbr -nr -nC 0.1 -nE $base_dir/eu_dggrid_extra_edges.dat -G $G1 -Ks $base_dir/gaez_crop/eu_dggrid_land_share.csv -z -t 4
echo $G1
end


###########################################################################################################
# 2. Simulation variant including conflict
set G1 80 # only G = 80 km case is used
set s 2 # s = 2 (scale of climate variability)
set y1 7000 # simulation start year (BCE)
set E1 1 5 10 20 100 200 # 1 / p_E parameter
set A1 1 2 5 10 20 50 # 1 / p_a parameter

# 2.1. run the simulation
for E in $E1
for A in $A1
set pa2 (math 1/$A)
$program_dir/n2w -Hf $base_dir/simulation/matrix_dg_G"$G"zC10.bin -Rr -RP -RA $pa2 -R 8 -R2 -RR 0.01 -Rc -Rp 1 -Rs 0 -RS -d 1 -s 1 -S 5500 -E $E -a 0.5 -o $base_dir/simulation/runs_G"$G1"_C10_E"$E"_A"$A"_s"$s"sp.out.gz -oz -op 1 -i 1596250 -I 5000 -If 0.5 -Ip 0.5 -Ce 0 -Cb 0 -De 1 -Dm 0.5 -Db 0.25 -Dl 200 -w $base_dir/climate/cell_ids_dggrid.csv $base_dir/climate/crop_data_cmb.dat -wm -wC -wf -wF -$y1 -ws $s > $base_dir/simulation/runs_G"$G1"_C10_E"$E"_A"$A"_s$s.out
echo $E $A
end
end

# 2.2. aggregate results in space for further analysis
for E in $E1
for A in $A1
zcat $base_dir/simulation/runs_G"$G1"_C10_E"$E"_A"$A"_s"$s"sp.out.gz | $program_dir/spas -m $base_dir/new_grid_aggr_flat.csv -f $base_dir/new_dgid_aggr_filter.dat | gzip -c > $base_dir/simulation/runs_G"$G1"_C10_E"$E"_A"$A"_s"$s"spaggr.out.gz
echo $E $A
end
end

# 2.3. create videos of the results with FFMPEG
for E in $E1
for A in $A1
zcat $base_dir/simulation/runs_G"$G1"_C10_E"$E"_A"$A"_s"$s"sp.out.gz | $program_dir/sp2png -oy -$y1 -op 10 -om 200 -ow 1608 -oh 900 -of $base_dir/simulation/runs_G"$G1"_C10_E"$E"_A"$A"_s"$s".mp4 -oc "h264 -pix_fmt yuv420p" -n $base_dir/dggrid/isea3h12eun.nbr -cp $base_dir/dggrid/dggrid_poly.csv
end
end



##############################################################################################################
# 3. Simulation variant without conflict, vary the speed of spread and the strength of climate variations
set G1 20 30 40
set s1 1 2 3 4
set y1 7000

# 3.1. run the simulation
for G in $G1
for s in $s1
$program_dir/n2 -G $G -H2 -Hf0 $base_dir/simulation/matrix_dg_G"$G"zC10.bin -p -d 1 -s 1 -S 5500 -o $base_dir/simulation/runG"$G"C10ws"$s"sp.out.gz -oz -op 1 -i 1607919 -I 100000 -If 0.5 -Ip 1 -Ce 0 -Cb 0 -De 1 -Dm 0.5 -Db 0.25 -Dl 200 -T -t 4 -w $base_dir/climate/cell_ids_dggrid.csv $base_dir/climate/crop_data_cmb.dat -wm -wC -wf -wF -$y1 -ws $s > $base_dir/simulation/runG"$G"C10ws"$s".out
echo $G $s
end
end

# 3.2. aggregate results in space for further analysis
for G in $G1
for s in $s1
zcat $base_dir/simulation/runG"$G"C10ws"$s"sp.out.gz | $program_dir/spas -m $base_dir/new_grid_aggr_flat.csv -f $base_dir/new_dgid_aggr_filter.dat | gzip -c > $base_dir/simulation/runG"$G"C10ws"$s"spaggr.out.gz
echo $G $s
end
end

# 3.3. create videos of the results with FFMPEG
for G in $G1
for s in $s1
zcat $base_dir/simulation/runG"$G"C10ws"$s"sp.out.gz | $program_dir/sp2png -oy -$y1 -op 10 -om 200 -ow 1608 -oh 900 -of $base_dir/simulation/runG"$G"C10ws"$s".mp4 -oc "h264 -pix_fmt yuv420p" -n $base_dir/dggrid/isea3h12eun.nbr -cp $base_dir/dggrid/dggrid_poly.csv
end
end



