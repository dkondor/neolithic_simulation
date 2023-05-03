/*
 * neolithic_common.hpp -- common helper to run the simulation with or
 * 	without raiders
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


#ifndef NEOLITHIC_COMMON_HPP
#define NEOLITHIC_COMMON_HPP

#include "neolithic_base.hpp"
#include "neolithic.hpp"
#include "neolithic_w.hpp"
#include "cellnet.hpp"
#include "cellnet_options.hpp"
#include <type_traits>
#include <stdio.h>


/* output to files -- need to be specialized to be used! */
template<bool have_raiders>
struct output_files {
	FILE* summary_out = stdout;
	FILE* detail_out = nullptr;
	FILE* split_pop_dist_out = nullptr;
	unsigned int out_period = 10;
	bool track_avg_distance = false;
	bool write_progress = true;
	unsigned int steps = 0; /* total number of steps (only used for writing the progress) */
	
	using nclass = typename std::conditional<have_raiders, neolithic_w<cellnet>, neolithic<cellnet> >::type;
	
	void operator()(const nclass& nn, unsigned int i, const step_res<nclass>& sr) {
		FILE* out1 = summary_out;
		FILE* outf = detail_out;
		
		if(out1) {
			if constexpr (have_raiders) {
				fprintf(out1, "%u\t%u\t%" PRIu64 "\t%u\t%" PRIu64 "\t%u\t%u", i, sr.t.farming_cells, sr.t.farmers,
					sr.t.raider_cells, sr.t.raiders, nn.raiders_created, nn.raiders_reverted);
				if(track_avg_distance) {
					fprintf(out1, "\t%f\t%u", sr.sum_dist / (double)sr.cnt_dist, sr.total_mig);
					fprintf(out1, "\t%f", sr.sum_dist_new / (double)sr.cnt_dist_new);
				}
				
				fprintf(out1, "\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n", sr.farmer_targets[0], sr.farmer_targets[1], sr.farmer_targets[2],
					sr.raider_targets[0], sr.raider_targets[1], sr.raider_targets[2], sr.group_mig_rejected);	
			}
			else {
				fprintf(out1, "%u\t%u\t%" PRIu64, i, sr.t.farming_cells, sr.t.farmers);
				if(track_avg_distance) {
					fprintf(out1, "\t%f\t%u\t%u\t0\t0", sr.sum_dist / (double)sr.cnt_dist, sr.cnt_dist, sr.total_mig);
					fprintf(out1, "\t%f", sr.sum_dist_new / (double)sr.cnt_dist_new);
				}
				fprintf(out1, "\t%u\n", sr.group_mig_rejected); /* number of attempted split-offs that failed due to no available target */
			}
		}
		
		if(split_pop_dist_out) {
			for(unsigned int N : sr.farmer_split_pop) fprintf(split_pop_dist_out, "%u\t%u\n", i, N);
		}
		
		if(outf && (i % out_period == 0)) nn.g.write_pop(outf, i);
		
		if(write_progress) {
			fprintf(stderr, "%u / %u steps complete\r", i, steps);
			fflush(stderr);
		}
	}
};


template<bool have_raiders>
class neolithic_sim {
	protected:
		using nclass = typename std::conditional<have_raiders, neolithic_w<cellnet>, neolithic<cellnet> >::type;
		
		nclass nn;
		cellnet_options opts;
		
		unsigned int nrepeat = 0; /* if set, do repeated runs and output only the average result */
		unsigned int nthreads = 1; /* number of threads to use for the above */
		
		/* match between cell IDs and region IDs for output -- only used if we are running repeated realizations */
		char* region_match_fn = nullptr;
	
		/* calculate the average and standard deviation based on the sum and sum of squares */
		static void do_avg_sd(double sum, double sum2, double N, double& avg, double& sd) {
			avg = sum / N;
			sd = sqrt(sum2 / N - avg * avg);
		}
	
	public:
		int run(int argc, char** argv) {
			for(int i = 1; i < argc;) {
				int x = nn.parse_option(argc, argv, i);
				if(x == 0) x = opts.parse_option(argc, argv, i, common_options<cellnet>::option_flags::SIMULATION |
					common_options<cellnet>::option_flags::WEATHER | common_options<cellnet>::option_flags::HELPER);
				if(x < 0) {
					/* argv[i][0] != '-' or invalid combination */
					fprintf(stderr, "Unknown parameter: %s!\n", argv[i]);
					return 1;
				}
				if(x > 0) {
					/* found a parameter that belongs to the common options or the simulation class */
					i += x;
					continue;
				}
				
				switch(argv[i][1]) {
					case 'o':
						switch(argv[i][2]) {
							case 'R':
								region_match_fn = argv[i+1];
								break;
							default:
								fprintf(stderr, "Invalid parameter: %s!\n", argv[i]);
								return 1;
						}
						i++; /* all cases have one argument */
						break;
					case 'N':
						if(argv[i][2] == 't') nthreads = atoi(argv[i+1]);
						else nrepeat = atoi(argv[i+1]);
						i++;
						break;
					default:
						fprintf(stderr, "Unknown parameter: %s!\n", argv[i]);
						return 1;
				}
				i++; /* need to step i */
			}
			
			if(opts.weather_allow_missing && opts.weather_cut) {
				fprintf(stderr, "Incompatible options!\n");
				return 1;
			}
			
			/* read the simulation space */
			if(opts.helper_matrix_fn) nn.read_combined_matrix(opts.helper_matrix_fn, opts.helper_only_net);
			else {
				opts.read_net_data(nn);
				if(opts.use_distance_matrix) nn.g.create_distance_matrix();
			}
			opts.do_adjust_K(nn, (opts.helper_matrix_fn != nullptr));
			opts.ensure_start_cell(nn);
			
			/* add the starting population */
			nn.set_start_pop(opts.start_cell, opts.starting_pop, false, opts.start_pop_fill_factor, opts.start_pop_pchoice);

			/* open the detailed output file (if needed) */
			FILE* outf = nullptr;
			if(opts.out_base) {
				/* optionally allow writing detailed results to the standard output */
				if(opts.out_base[0] == '-' && opts.out_base[1] == 0) outf = stdout;
				else {
					if(opts.out_base_zip) {
						char* tmp1 = new char[strlen(opts.out_base) + 30];
						sprintf(tmp1, "/bin/gzip -c -1 > %s", opts.out_base);
						outf = popen(tmp1, "w");
						delete[]tmp1;
					}
					else outf = fopen(opts.out_base, "w");
					if(!outf) {
						fprintf(stderr, "Error opening output file %s!\n", opts.out_base);
						return 1;
					}
				}
			}
			
			/* run the simulation */
			if(nrepeat > 1) {
				/* do repeated realizations, save the average */
				const unsigned int steps = opts.steps;
				const unsigned int out_period = opts.out_period;
				
				/* 1. read all the weather data */
				weather_update_all wd;
				opts.read_weather_mapping(wd, nn);
				if(opts.weather_file) wd.read_data(opts.weather_file, nn.g.size());
				
				/* 2. define the output aggregation */
				output_avg<nclass> oa;
				oa.out_period = out_period;
				
				/* 3. read cell ID -- region ID matches if needed */
				std::vector<unsigned int> region_ids_r;
				unsigned int nregions = 0;
				if(region_match_fn && outf) {
					read_table2 rtr(region_match_fn);
					rtr.set_delim(',');
					rtr.read_line(); /* skip header */
					region_ids_r = oa.read_region_ids(rtr);
					nregions = region_ids_r.size();
				}
				
				/* 4. allocate space for the output */
				oa.allocate(steps, outf ? (nregions ?: (unsigned int)nn.g.size()) : 0, have_raiders);
				
				/* run nrepeat simulations */
				if(nthreads > 1) oa.runmultiple(nn, opts, wd, nrepeat, nthreads);
				else for(unsigned int j = 0; j < nrepeat; j++) {
					if(j) nn.set_start_pop(opts.start_cell, opts.starting_pop, true, opts.start_pop_fill_factor, opts.start_pop_pchoice);
					do_one_run(nn, steps, opts.start_cell, wd, oa);
				}
				
				/* write the results */
				if(outf != stdout) for(unsigned int j = 0; j < steps; j++) {
					typename output_avg<nclass>::record avg;
					typename output_avg<nclass>::record sd;
#define AVGSD(name) do_avg_sd(oa.sums[j].name, oa.sum2[j].name, nrepeat, avg.name, sd.name)
					AVGSD(nfarmers);
					AVGSD(pfarmers);
					if(have_raiders) {
						AVGSD(nraiders);
						AVGSD(praiders);
					}
					AVGSD(avg_dist);
#undef AVGSD
					fprintf(stdout, "%u\t%f\t%f\t", j, avg.nfarmers, avg.pfarmers);
					if(have_raiders) fprintf(stdout, "%f\t%f\t", avg.nraiders, avg.praiders);
					fprintf(stdout, "%f\t", avg.avg_dist);
					fprintf(stdout, "%f\t%f\t", sd.nfarmers, sd.pfarmers);
					if(have_raiders) fprintf(stdout, "%f\t%f\t", sd.nraiders, sd.praiders);
					fprintf(stdout, "%f\n", sd.avg_dist);
				}
				
				if(outf) for(unsigned int j1 = 0; j1 < steps / out_period + 1; j1++) {
					unsigned int j = j1 * out_period;
					if(j >= steps) break;
					const auto& tmpf = oa.cells_farmers[j1];
					const auto& tmpr = oa.cells_raiders[j1];
					double nr2 = nrepeat;
					size_t size1 = nregions ? (size_t)nregions : nn.g.size();
					for(size_t k = 0; k < size1; k++) {
						if(tmpf[k].nonzero() || (have_raiders && tmpr[k].nonzero())) {
							unsigned int gid = nregions ? region_ids_r[k] : nn.g.ptid(k);
							fprintf(outf, "%u\t%u\t%f\t%f", j, gid, ((double)tmpf[k].n) / nr2, ((double)tmpf[k].pop) / nr2);
							if(have_raiders) fprintf(outf, "\t%f\t%f\n", ((double)tmpr[k].n) / nr2, ((double)tmpr[k].pop) / nr2);
							else fputc('\n', outf);
						}
					}
				}
			}
			else {
				output_files<have_raiders> of;
				of.detail_out = outf;
				/* note: if we are writing the detailed results to stdout,
				 * we need to supress the summary output that would normally
				 * go there */
				if(outf == stdout) of.summary_out = nullptr;
				of.steps = opts.steps;
				of.write_progress = true;
				of.out_period = opts.out_period;
				of.track_avg_distance = opts.track_avg_distance;
				if(opts.out_split_pop) of.split_pop_dist_out = fopen(opts.out_split_pop, "w");
				
				if(opts.weather_file) {
					weather_update_file wf;
					opts.read_weather_mapping(wf, nn);
					wf.open_file(opts.weather_file);
					do_one_run(nn, opts.steps, opts.start_cell, wf, of, true);
				}
				else do_one_run(nn, opts.steps, opts.start_cell, weather_update_noop(), of, true);
				if(of.split_pop_dist_out) fclose(of.split_pop_dist_out);
			}
			
			if(outf && outf != stdout) {
				if(opts.out_base_zip) pclose(outf);
				else fclose(outf);
			}
			return 0;
		}
};

#endif

