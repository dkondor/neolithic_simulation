/*
 * create_helper_matrix.cpp -- create a binary matrix with pre-computed
 * 	probability distributions used for migrations
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
 */


#include <stdio.h>
#include <string.h>
#include <random>
#include <optional>
#include <vector>
#include <unordered_map>
#include <iostream>
#include "cellnet.hpp"
#include "cellnet_options.hpp"
#include "neolithic_w.hpp"
#include "neolithic.hpp"

int main(int argc, char **argv)
{
	neolithic_w<cellnet> nnw;
	neolithic<cellnet> nn0;
	cellnet_options opts;
	
	bool check_matrix = false; /* instead of creating a matrix, open one and check it for correctness */
	bool write_combined_yields = false; /* instead of creating the matrix, just write out the calculated (scaled) yields for each cell */
	
	bool have_raiders = true; /* flag to determine which model to use */
	/* 1st pass among the options to determine the above flag */
	for(int i = 1; i < argc; i++) if(argv[i][0] == '-' && argv[i][1] == 'W' && argv[i][2] == '0')
		have_raiders = false;
	
	nbase<cellnet>& nn = have_raiders ? dynamic_cast<nbase<cellnet>&>(nnw) :
		dynamic_cast<nbase<cellnet>&>(nn0);
	
	bool use_global_neighbors = have_raiders; /* if we have raiders, we use global neighbors by default */
	bool write_all_neighbors = false; /* instead of creating / checking the matrix, only write out all of the neighbors in all helpers */
	bool write_helpers = true; /* whether we should write the helper distributions */
	
	for(int i = 1; i < argc; i++) {
		if(argv[i][0] == '-') {
			bool found = true;
			switch(argv[i][1]) {
				case 'N':
					write_all_neighbors = true;
					/* fallthrough */
				case 'C':
					check_matrix = true;
					break;
				case 'W':
					switch(argv[i][2]) {
						case 0:
							write_combined_yields = true;
							break;
						case '0':
							break;
						case 'H':
							write_helpers = false;
							break;
						default:
							found = false;
							break;
					}
					break;
				case 'g':
					if(argv[i][2] == '0') use_global_neighbors = false;
					else found = false;
					break;
				default:
					found = false;
					break;
			}
			if(!found) {
				int x = nn.parse_option(argc, argv, i);
				if(x == 0) x = opts.parse_option(argc, argv, i, 0);
				if(x > 0) i += (x - 1);
				else {
					fprintf(stderr, "Unknown parameter: %s!\n", argv[i]);
					return 1;
				}
			}
		}
		else {
			fprintf(stderr, "Unknown parameter: %s!\n", argv[i]);
			return 1;
		}
	}
	
	if(!(opts.out_base || write_combined_yields)) {
		fprintf(stderr, "No output file name given!\n");
		return 1;
	}
	
	if(check_matrix) {
		nn.read_combined_matrix(opts.out_base);
		if(write_all_neighbors) nn.write_all_neighbors(stdout);
		else if(write_combined_yields) {
			FILE* out1 = stdout;
			for(const auto& pt : nn.g) {
				unsigned int id = pt;
				fprintf(out1, "%u\t%g\n", id, nn.g.at(pt).K);
			}
		}
		else {
			double s1 = nn.check_helper();
			fprintf(stderr, "%f\n", s1);
		}
	}
	else {
		opts.read_net_data(nn);
		if(opts.use_distance_matrix) nn.g.create_distance_matrix();
		opts.do_adjust_K(nn, false);
		
		if(write_combined_yields) {
			FILE* out1 = stdout;
			for(const auto& pt : nn.g) {
				/* note: we know that cell IDs are integers */
				unsigned int id = pt;
				fprintf(out1, "%u\t%g\n", id, nn.g.at(pt).K);
			}
		}
		else nn.write_combined_matrix(opts.out_base, false, use_global_neighbors, !write_helpers);
	}
	
	return 0;
}

