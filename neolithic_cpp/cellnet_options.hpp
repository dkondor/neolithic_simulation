/*
 * cellnet_options.hpp -- options only relevant when running the
 * 	simulation with a network of cells
 * 
 * Copyright 2021 Daniel Kondor <kondor@csh.ac.at>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */

#ifndef CELLNET_OPTIONS_HPP
#define CELLNET_OPTIONS_HPP
#include "neolithic_base.hpp"
#include "cellnet.hpp"

struct cellnet_options : public common_options<cellnet> {	
	char* net_file = nullptr; /* read a network from this file */
	bool net_file_edges = false; /* set to true if the file with the network is an edgelist */
	bool net_file_bin = false; /* set to true if the cellnet object should be read from a binary file (without the probability helpers) */
	bool net_file_write = false; /* set to true if the cellnet object should be written to a binary file */
	bool net_recalculate_distances = false; /* set to true to recalculate distances along edges after reading the network */
	
	std::vector<const char*> extra_edges_files; /* files containing additional edges */
	
	char* coords_file = nullptr; /* read cell (center) coordinates from this file */
	bool coords_file_csv = false;
	bool coords_file_has_types = false; /* set to true if the above file includes information about coastal cells (as cell type) */
	bool coords_file_has_factors = false; /* set to true if the above file includes land travel scaling factors */
	double coastal_travel_factor = 1.0; /* factor used for coastal travel */
	
	bool have_K = false; /* set to true if we should load carrying capacities from the above input file */
	std::vector<const char*> scale_K_fn; /* scale productivity values by these given as inputs */
	bool scale_K_geo = false; /* scale productivity values by a factor corresponding to the latitude */
	
	/* use a distance matrix */
	bool use_distance_matrix = false;
	
	/* parse specialized options; note: this is NOT a virtual function,
	 * so this needs to be called on an instance of this object and not
	 * on the base class */
	int parse_option(int argc, char** argv, int i, uint32_t flags = 0) {
		int res = common_options<cellnet>::parse_option(argc, argv, i, flags);
		if(res != 0) return res;
		switch(argv[i][1]) {
			case 'K':
				switch(argv[i][2]) {
					case 'g':
						scale_K_geo = true;
						res = 1;
						break;
					case 's':
						scale_K_fn.push_back(argv[i+1]);
						res = 2;
						break;
				}
				break;
			case 'k':
				have_K = true;
				res = 1;
				break;
			case 'c':
				switch(argv[i][2]) {
					case 't':
						coords_file_has_types = true;
						res = 1;
						break;
					case 'l':
						coords_file_has_factors = true;
						res = 1;
						break;
					case 'c':
						coords_file_csv = true;
						/* fallthrough*/
					case 0:
						coords_file = argv[i+1];
						res = 2;
						break;
					default:
						res = -1;
						break;
				}
				break;
			case 'n':
				switch(argv[i][2]) {
					case 'r':
						net_recalculate_distances = true;
						res = 1;
						break;
					case 'C':
						coastal_travel_factor = atof(argv[i+1]);
						res = 2;
						break;
					case 'E':
						extra_edges_files.push_back(argv[i+1]);
						res = 2;
						break;
					case 'w':
						net_file_write = true;
						res = 1;
						break;
					case 'b':
						net_file_bin = true;
						/* fallthrough */
					case 'e':
						net_file_edges = true;
						/* fallthrough */
					case 0:
						net_file = argv[i+1];
						res = 2;
						break;
					default:
						res = -1;
						break;
				}
				break;
			case 'z':
				/* use a distance matrix */
				use_distance_matrix = true;
				res = 1;
				break;
		}
		return res;
	}
	
	/* load a network data (only if cell_collection supports it) */
	void read_net_data(nbase<cellnet>& nn) {
		if(helper_matrix_fn) throw std::runtime_error("common_options::read_net_data() should only be called if no external helper was given!\n");
		if(!(net_file && coords_file)) throw std::runtime_error("common_options::read_net_data(): missing input filenames!\n");
		
		uint32_t flags = 0;
		if(have_K) flags |= cellnet::load_net_flags::read_K;
		if(net_file_edges) flags |= cellnet::load_net_flags::read_edges;
		if(coords_file_has_types) flags |= cellnet::load_net_flags::read_types;
		if(coords_file_has_factors) flags |= cellnet::load_net_flags::read_travel_factors;
		
		nn.g.load_net(net_file, coords_file, flags, coords_file_csv);
		for(auto fn : extra_edges_files) nn.g.read_edges(fn);
		if(net_recalculate_distances) nn.g.calculate_edge_distances(coastal_travel_factor);
		
		if(have_K) {
			if(scale_K_geo) nn.g.scale_K_lat();
			for(const char* fn1 : scale_K_fn) nn.g.scale_K(fn1, true);
		}
		else for(const auto& pt : nn.g) nn.g.at(pt).K = K1;
	}
};

#endif

