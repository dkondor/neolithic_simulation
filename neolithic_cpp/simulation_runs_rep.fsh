#!/usr/bin/fish

# simulation_runs_rep.fsh -- run all simulations with 100 realizations
# Notes:
#  - this script uses the syntax of the fish shell (https://github.com/fish-shell/fish-shell or https://fishshell.com/)
#  - instead of running the whole script, it can make sense to run each section separately as needed by copying it to a shell window
#  - the total runtime can be 1-3 weeks
#  - this script was originally run on a machine with an Intel Xeon Phi CPU that has a very large number of low performance cores; in accordance with this, it utilizes 50-100 cores when running the simulation and processing steps
#  - it is recommended to change the "ncores" variable, both below and in the R script acf_new_aggr_rep.r as well to the actual available number of CPU cores

# base data directory -- change this to whatever is appropriate!
set base_dir ~/CSH/HoloSim/data
set data_dir $base_dir

# directory of the simulation code (by default it is the current directory) -- set this to whatever is appropriate
set code_dir .

set N 100 # number of repeated realizations (in each simulation case)
set y1 7000 # simulation start year
set ncores 50 # number of cores to use when running the simulation
# note: this also needs to be adjusted in the R analysis script (acf_new_aggr_rep.r)!!

# output directory
set outdir $base_dir/simulation/runs_rep
mkdir -p $outdir




begin
set d1 (date)
echo $d1


# 1. main model variant (stationary aggressors)
echo G,s,A,E,i,seed > $outdir/seeds_run1.csv
for G in 50 60 70 80
for s in 0 1 2 3 4
for A in 1 2 5 10 20 50
set pa2 (math 1/$A)
for E in 1 5 10 20 100 200
set fnbase $outdir/res_G"$G"_C10_E"$E"_A"$A"_s"$s"
for i in (seq $N)
set seed (head -c 4 /dev/urandom | hexdump -e "\"%u\"")
echo $G,$s,$A,$E,$i,$seed >> $outdir/seeds_run1.csv
echo $code_dir/n2w -Hf $data_dir/matrix_dg_G"$G"zC10.bin -Rr -Rp 1 -RP -RA $pa2 -R 8 -R2 -RR 0.01 -Rc -Rp 1 -Rs 0 -RS -d 1 -s $seed -S 6000 -E $E -a 0.5 -o - -op 1 -i 1596250 -I 5000 -If 0.5 -Ip 0.5 -Ce 0 -Cb 0 -De 1 -Dm 0.5 -Db 0.25 -Dl 200 -w $data_dir/climate/cell_ids_dggrid.csv $data_dir/climate/crop_data_cmb.dat -wm -wC -wf -wF -$y1 -ws $s \| $code_dir/spas -m new_grid_aggr_flat_r500.csv -f new_dgid_aggr_filter.dat \| gzip -c \> $fnbase"_"$i.out.gz
end
end | $code_dir/rm -t $ncores

echo \n
for E in 1 5 10 20 100 200
set fnbase $outdir/res_G"$G"_C10_E"$E"_A"$A"_s"$s"
echo Rscript --verbose --no-save $code_dir/acf_new_aggr_rep.r $fnbase $N
end | $code_dir/rm -t 2

rm $outdir/*out.gz

echo $G $A $s
date

end
end
end


# 2. model with conflict, but no aggressors (first order model)
echo G,s,E,i,seed > $outdir/seeds_run1nr.csv
for G in 50 60 70 80
for s in 0 1 2 3 4
for E in 1 5 10 20 100 200
set fnbase $outdir/resnr_G"$G"_C10_E"$E"_s"$s"
for i in (seq $N)
set seed (head -c 4 /dev/urandom | hexdump -e "\"%u\"")
echo $G,$s,$E,$i,$seed >> $outdir/seeds_run1nr.csv
echo $code_dir/n2w -Hf $data_dir/matrix_dg_G"$G"zC10.bin -Rr -Rp 1 -RP -RA 0 -R 8 -R2 -RR 0.01 -Rc -Rp 0 -Rs 0 -RS -d 1 -s $seed -S 6000 -E $E -a 0.5 -o - -op 1 -i 1596250 -I 5000 -If 0.5 -Ip 0.5 -Ce 0 -Cb 0 -De 1 -Dm 0.5 -Db 0.25 -Dl 200 -w $data_dir/climate/cell_ids_dggrid.csv $data_dir/climate/crop_data_cmb.dat -wm -wC -wf -wF -$y1 -ws $s \| $code_dir/spas -m new_grid_aggr_flat_r500.csv -f new_dgid_aggr_filter.dat \| gzip -c \> $fnbase"_"$i.out.gz
end
end | $code_dir/rm -t $ncores

for E in 1 5 10 20 100 200
set fnbase $outdir/resnr_G"$G"_C10_E"$E"_s"$s"
echo Rscript --verbose --no-save $code_dir/acf_new_aggr_rep.r $fnbase $N
end | $code_dir/rm -t 2
rm $outdir/*out.gz

echo $G,$s
date

end
end


date


# 3. model without conflict (split-off groups can only target empty cells)
echo G,s,i,seed > $outdir/seeds_run1n.csv
for s in 1 2 3 4
for G in 20 30 40
set fnbase $outdir/resn_G"$G"_s"$s"
for i in (seq $N)
set seed (head -c 4 /dev/urandom | hexdump -e "\"%u\"")
echo $G,$s,$A,$E,$i,$seed >> $outdir/seeds_run1n.csv
echo $code_dir/n2w -Hf $data_dir/matrix_dg_G"$G"zC10.bin -Rp 0 -a 0 -d 1 -s $seed -S 6000 -E 100 -o - -op 1 -i 1607919 -I 100000 -If 0.5 -Ip 1 -Ce 0 -Cb 0 -De 1 -Dm 0.5 -Db 0.25 -Dl 200 -w $data_dir/climate/cell_ids_dggrid.csv $data_dir/climate/crop_data_cmb.dat -wm -wC -wf -wF -7000 -ws $s \| $code_dir/spas -m new_grid_aggr_flat_r500.csv -f new_dgid_aggr_filter.dat \| gzip -c \> $fnbase"_"$i.out.gz
end
end | $code_dir/rm -t $ncores

for G in 20 30 40
set fnbase $outdir/resn_G"$G"_s"$s"
echo Rscript --verbose --no-save $code_dir/acf_new_aggr_rep.r $fnbase $N
end | $code_dir/rm -t 2

rm $outdir/*out.gz

echo $s
date
end


set d2 (date)
echo $d1\n$d2
end



# 4. additional simulation variants:
# 1st case: stationary aggressors, but they move when created (original fission still happens): no -RS option
# 2nd case: aggressors do not consider cells' population when choosing a target: no -RP option
# 3rd case: logistic model for fission probabilities (Alberti, 2014): new -DL option
# -> 8 total combinations, run all of these for two parameter combinations


set outdir $base_dir/simulation/runs_rep_new
mkdir -p $outdir

# simulation variants
set ppar "" "-RP"
set pnames 0 P
set spar "" "-RS"
set snames 0 S
set lpar "" "-DL"
set lnames 0 L

# all combinations -- only G = 80, s = 2, A = 10 (i.e. p_A = 0.1); + E = 1 and 5 (i.e. p_E = 1 and 0.2)
echo G,s,A,E,P,S,L,i,seed > $outdir/seeds_run1.csv
for G in 80
for s in 2
for A in 10
set pa2 (math 1/$A)
for E in 1 5
for i0 in 1 2
for j in 1 2
for k in 1 2
set fnbase $outdir/res$pnames[$i0]$snames[$j]$lnames[$k]_G"$G"_C10_E"$E"_A"$A"_s"$s"
for i in (seq $N)
set seed (head -c 4 /dev/urandom | hexdump -e "\"%u\"")
echo $G,$s,$A,$E,$pnames[$i0],$snames[$j],$lnames[$k],$i,$seed >> $outdir/seeds_run1.csv
echo $code_dir/n2w -Hf $data_dir/matrix_dg_G"$G"zC10.bin -Rp 1 $ppar[$i0] -RA $pa2 -R 8 -R2 -RR 0.01 -Rc -Rp 1 -Rs 0 $spar[$j] -d 1 -s $seed -S 6000 -E $E -a 0.5 -o - -op 1 -i 1596250 -I 5000 -If 0.5 -Ip 0.5 -Ce 0 -Cb 0 -De 1 -Dm 0.5 -Db 0.25 -Dl 200 $lpar[$k] -w $data_dir/climate/cell_ids_dggrid.csv $data_dir/climate/crop_data_cmb.dat -wm -wC -wf -wF -$y1 -ws $s \| $code_dir/spas -m new_grid_aggr_flat_r500.csv -f new_dgid_aggr_filter.dat \| gzip -c \> $fnbase"_"$i.out.gz
end
end
end
end
end | $code_dir/rm -t $ncores

for E in 1 5
for i0 in 1 2
for j in 1 2
for k in 1 2
set fnbase $outdir/res$pnames[$i0]$snames[$j]$lnames[$k]_G"$G"_C10_E"$E"_A"$A"_s"$s"
echo Rscript --verbose --no-save acf_new_aggr_rep.r $fnbase $N
end
end
end
end | $code_dir/rm -t 2

end







