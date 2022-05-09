#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  nwrun1.py -- run the neolithic agent-based simulation
#  
#  Copyright 2021 Daniel Kondor <kondor@csh.ac.at>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

import sys
import neolithic

def main(args):
	nn = neolithic.neolithic()
	
	# main parameters for the simulation
	helper_matrix_fn = None # previously created helper matrix to load
	no_helpers = False # if true, the actual probabilities are not used, only the cells and distances
	steps = 2000 # number of steps to take
	K1 = None # average carrying capacity to set after loading
	Kthresh = None # threshold for cells' carrying capacity
	
	out_base = None # detailed output file
	out_period = 10 # period for writing detailed output
	
	start_cell = None # cell where to seed the simulation
	starting_pop = None # if given, try to distribute this amount of population in the beginning
	start_pop_fill_factor = 1.0 # starting population should fill ratio to carrying capacity
	start_pop_pchoice = 1.0 # probability of each cell included in the start population

	
	weather_mapping = None # use a mapping from cell IDs to weather cell IDs
	weather_mapping_csv = False # true if weather mapping file is in CSV format
	weather_file = None # climate-based variability data (multiplicative factor for carrying capacities)
	wupd = neolithic.weather_update()
	
	net_fn = None # filename with the connections among cells
	net_edges = False # if true, the above file is an edgelist (otherwise it is the neighbors written by DGGRID)
	net_extra_edges = list() # read extra connections (edges) from the files given in this list
	net_recalculate_distances = False # set to true to (re-)calculate distances along all edges after reading the input
	
	coords_fn = None # filename with the cell coordinates and carrying capacities
	coords_csv = False # set to true if the above file is in CSV format
	coords_have_types = False # set to true if the above file includes cell "types"
		# type is expected as an integer ID; currently, the only supported type is 1, which indicates coastal cells
	coords_have_factors = False # set to true if the above file includes additional scaling factors for travel distances
	create_matrix = False # set to true to create a probability helper matrix instead of running the simulation
	create_distances = False # set to true if distances among all cells should be calculated and saved as a matrix as well
	
	scale_K_geo = False
	scale_K_fn = list()
	
	i = 1
	while i < len(args):
		r1 = nn.parse_arg(args, i)
		if r1 < -1:
			raise BaseException('Invalid argument: {}!\n'.format(args[i]))
		elif r1 > 0:
			i += r1
		# here, r1 == 0
		elif args[i][1] == 'w':
			if len(args[i]) == 2:
				if i + 2 >= len(args):
					raise BaseException('-w needs two arguments!\n')
				weather_mapping = args[i+1];
				weather_file = args[i+2];
				i += 2
			elif args[i][2] == 'm':
				wupd.weather_allow_missing = True
			elif args[i][2] == 'c':
				wupd.weather_cut = True
			elif args[i][2] == 'f':
				wupd.weather_factors = True
			elif args[i][2] == 'C':
				weather_mapping_csv = True
			elif args[i][2] == 'F':
				wupd.first_year = int(args[i+1])
				i += 1
			elif args[i][2] == 's':
				wupd.sd_factor = float(args[i+1])
				i += 1
			else:
				raise BaseException('Invalid argument: {}!\n'.format(args[i]))
			i += 1
		elif args[i][1] == 'i':
			start_cell = int(args[i+1])
			i += 2
		elif args[i][1] == 'I':
			if len(args[i]) == 2:
				starting_pop = int(args[i+1])
				if starting_pop <= 0:
					raise BaseException('Starting population must be > 0!\n')
			elif args[i][2] == 'f':
				start_pop_fill_factor = float(args[i+1])
				if start_pop_fill_factor < 0.0 or start_pop_fill_factor > 1.0:
					raise BaseException('Invalid argument: {} {} (must be between 0 and 1)!\n'.format(args[i], args[i+1]))
			elif args[i][2] == 'p':
				start_pop_pchoice = float(args[i+1])
				if start_pop_pchoice < 0.0 or start_pop_pchoice > 1.0:
					raise BaseException('Invalid argument: {} {} (must be between 0 and 1)!\n'.format(args[i], args[i+1]))
			else:
				raise BaseException('Invalid argument: {}!\n'.format(args[i]))
			i += 2
		elif args[i][1] == 'o':
			if len(args[i]) > 2 and args[i][2] == 'p':
				out_period = int(args[i+1])
				if out_period <= 0:
					raise BaseException('Invalid argument: {} {} (must be positive)!\n'.format(args[i], args[i+1]))
			else:
				out_base = args[i+1]
			i += 2
		elif args[i][1] == 'S':
			steps = int(args[i+1])
			i += 2
		elif args[i][1] == 'K':
			if len(args[i]) > 2:
				if args[i][2] == 'g':
					scale_K_geo = True
					i += 1
				elif args[i][2] == 's':
					scale_K_fn.append(args[i+1])
					i += 2
				elif args[i][2] == 't':
					Kthresh = float(args[i+1])
					i += 2
				else:
					raise BaseException('Invalid argument: {}!\n'.format(args[i]))
			else:
				K1 = float(args[i+1])
				i += 2
		elif args[i][1] == 'H':
			helper_matrix_fn = args[i+1]
			if len(args[i]) > 2 and args[i][2] == '1':
				no_helpers = True
			i += 2
		elif args[i][1] == 'c':
			if len(args[i]) == 2:
				coords_fn = args[i+1]
				i += 2
			elif args[i][2] == 'c':
				coords_csv = True
				coords_fn = args[i+1]
				i += 2
			elif args[i][2] == 't':
				coords_have_types = True
				i += 1
			elif args[i][2] == 'l':
				coords_have_factors = True
				i += 1
			else:
				raise BaseException('Invalid argument: {}!\n'.format(args[i]))
		elif args[i][1] == 'n':
			if len(args[i]) == 2:
				net_fn = args[i+1]
				i += 2
			elif args[i][2] == 'e':
				net_edges = True
				net_fn = args[i+1]
				i += 2
			elif args[i][2] == 'r':
				net_recalculate_distances = True
				i += 1
			elif args[i][2] == 'C':
				nn.g.coastal_travel_factor = float(args[i+1])
				if nn.g.coastal_travel_factor < 0.0 or nn.g.coastal_travel_factor > 1.0:
					raise BaseException('Invalid argument: {} {} (must be between 0 and 1)!\n'.format(args[i], args[i+1]))
				i += 2
			elif args[i][2] == 'E':
				net_extra_edges.append(args[i+1])
				i += 2
			else:
				raise BaseException('Invalid argument: {}!\n'.format(args[i]))
		elif args[i][1] == 'M':
			create_matrix = True
			i += 1
		elif args[i][1] == 'z':
			create_distances = True
			i += 1
		else:
			raise BaseException('Invalid argument: {}!\n'.format(args[i]))
		
	if (helper_matrix_fn is None or create_matrix) and (net_fn is None or coords_fn is None):
		raise BaseException('No input data given (use the -Hf parameter or the -n and -c parameters)!\n')
	
	if helper_matrix_fn is not None and not create_matrix:
		nn.read(helper_matrix_fn, no_helpers)
		for pt in nn.g:
			nn.g[pt].K = nn.orig_K[pt]
		Kthresh = 0.0 # not supported in this case
	else:
		nn.g.load_net(net_fn, coords_fn, coords_csv, net_edges, read_edge_dists = False,
			read_types = coords_have_types, read_travel_factors = coords_have_factors)
		for fn in net_extra_edges:
			nn.g.read_edges(fn)
		if net_recalculate_distances:
			nn.g.calculate_edge_dists()
			nn.g.use_net_dist = True
		if create_distances:
			nn.g.create_distance_matrix()
		nn.create_global_ids()
		if scale_K_geo:
			nn.g.scale_K_lat()
		for fn in scale_K_fn:
			nn.g.scale_K(fn, True)
		if Kthresh is None:
			Kthresh = 0.0
	
	if K1 is not None or Kthresh > 0.0:
		nn.g.adjust_K(K1, Kthresh)
		
	if create_matrix:
		nn.write(helper_matrix_fn)
		return 0
	
	if weather_mapping:
		wupd.read_mapping(weather_mapping, weather_mapping_csv, nn)
		wupd.wf = open(weather_file, 'r')
	
	if start_cell is None:
		for pt in nn.g:
			start_cell = pt
			break
	if starting_pop is not None:
		for pt in nn.g.neighbors(start_cell, 1e9, True):
			if start_pop_pchoice < 1.0:
				if nn.rng.random() > start_pop_pchoice:
					continue
			c1 = nn.g[pt[0]]
			tmp1 = min(starting_pop, int(c1.K * start_pop_fill_factor))
			if tmp1 <= 0:
				continue
			c1.N = tmp1
			c1.st = neolithic.cell.state.FARMING
			starting_pop -= tmp1
			if starting_pop <= 0:
				break
	else:
		c1 = nn.g[start_cell]
		c1.N = 100
		c1.st = neolithic.cell.state.FARMING
	
	outf1 = None
	if out_base is not None:
		outf1 = open(out_base, 'w')
	
	for i in range(steps):
		# process climate data, update carrying capacities
		# (note: this is a no-op if no climate files were opened)
		wupd.step_one_year(nn, i)
		
		# step the simulation
		nn.step(1.0)
		nn.create_group_orders()
		nn.apply_group_mig(i)
		
		if outf1 is not None and i % out_period == 0:
			# write detailed output: each cell that is not empty
			for pt in nn.g:
				c1 = nn.g[pt]
				if c1.st == neolithic.cell.state.FARMING:
					outf1.write('{}\t{}\t{}\t0\n'.format(i, pt, c1.N))
				elif c1.st == neolithic.cell.state.RAIDERS:
					outf1.write('{}\t{}\t{}\t1\n'.format(i, pt, c1.N))
		
		
		# output: count cells of different type
		raider_cells = 0
		raider_pop = 0
		farmer_cells = 0
		farmer_pop = 0
		for pt in nn.g:
			c1 = nn.g[pt]
			if c1.st == neolithic.cell.state.FARMING:
				farmer_cells += 1
				farmer_pop += c1.N
			elif c1.st == neolithic.cell.state.RAIDERS:
				raider_cells += 1
				raider_pop += c1.N	
		sys.stdout.write('{}\t{}\t{}\t{}\t{}\n'.format(i, farmer_cells, farmer_pop, raider_cells, raider_pop))
		sys.stderr.write('\r{} / {} steps complete'.format(i, steps))
	
	sys.stderr.write('\n')
	
	if outf1 is not None:
		outf1.close()
	
	return 0

if __name__ == '__main__':
#    import os
#    os.environ["OPENBLAS_NUM_THREADS"] = "1" # needed before importing numpy to avoid starting excessive extra threads
	sys.exit(main(sys.argv))

