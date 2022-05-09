/*
 * sample_res.cpp -- generate a sample of "dates" as simulated
 * 	C14 dates from a simulation result
 * 
 * Copyright 2022 Daniel Kondor <kondor@csh.ac.at>
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


#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include <cmath>
#include <utility>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <random>
#include <time.h>
#include "read_table_cpp.h"

int main(int argc, char **argv)
{
	double p1 = 1e-6; /* expected rate of generating one sample date (1 / year / km2) */
	double cellsize = 95.978; /* size of cells (in input), in km2 -- needed for normalization */
	unsigned int tmin = 0; /* minimum (simulation) date to use in the output */
	unsigned int tmax = 0; /* maximum simulation date to use (0: use all data) */
	bool norm_pop = true; /* whether to take into account the population numbers or only whether a cell is occupied or empty */
	unsigned int total_dates = 0; /* can give directly the number of dates to generate */
	uint64_t seed = time(0);
	char* fcnts = nullptr; /* file with per-cell counts of the number of events to generate */
	double cnts_factor = 1.0; /* scale the counts from the above file with this factor (if != 1.0) */
	char* aggr_file = nullptr; /* instead of one-to-one match, use an aggregation given by this file */
	
	for(int i = 1; i < argc; i++) if(argv[i][0] == '-') switch(argv[i][1]) {
		case 'p':
			p1 = atof(argv[i+1]);
			i++;
			break;
		case 'c':
			cellsize = atof(argv[i+1]);
			i++;
			break;
		case 'm':
			tmax = atoi(argv[i+1]);
			i++;
			break;
		case 'M':
			tmin = atoi(argv[i+1]);
			i++;
			break;
		case 'n':
			norm_pop = false;
			break;
		case 's':
			seed = strtoul(argv[i+1], 0, 10);
			i++;
			break;
		case 't':
			total_dates = atoi(argv[i+1]);
			i++;
			break;
		case 'C':
			fcnts = argv[i+1];
			i++;
			break;
		case 'f':
			cnts_factor = atof(argv[i+1]);
			i++;
			break;
		case 'a':
			aggr_file = argv[i+1];
			i++;
			break;
		default:
			fprintf(stderr, "Unknown parameter: %s!\n", argv[i]);
			return 1;
	}
	else {
		fprintf(stderr, "Unknown parameter: %s!\n", argv[i]);
		return 1;
	}
	
	
	uint64_t sum1 = 0; // sum of weights (occupied cells or population)
	
	std::unordered_map<uint64_t, unsigned int> cnts; // counts of events to generate
	if(fcnts) {
		read_table2 rt(fcnts);
		while(rt.read_line()) {
			uint64_t id;
			unsigned int cnt;
			if(!rt.read(id, cnt)) break;
			if(!cnt) continue;
			if(cnts_factor != 1.0) cnt = (unsigned int)std::round(cnts_factor * cnt);
			if(!cnt) cnt = 1;
			cnts[id] = cnt;
		}
		
		if(rt.get_last_error() != T_EOF || cnts.size() == 0) {
			fprintf(stderr, "Error reading input:\n");
			rt.write_error(std::cerr);
			return 1;
		}
	}
	
	std::vector<std::pair<uint64_t, std::vector<uint64_t> > > aggr;
	if(aggr_file) {
		if(cnts.size() == 0) {
			fprintf(stderr, "Cannot use aggregation without reading counts first!\n");
			return 1;
		}
		read_table2 rt(aggr_file);
		while(rt.read_line()) {
			uint64_t id;
			if(!rt.read(id)) break;
			std::vector<uint64_t> ids;
			if(!cnts.count(id)) continue; // no need to read aggregation for cells that do not have any events
			while(true) {
				uint64_t tmp;
				if(rt.read(tmp)) ids.push_back(tmp);
				else {
					if(rt.get_last_error() == T_EOL) break;
					fprintf(stderr, "Error reading input:\n");
					rt.write_error(std::cerr);
					return 1;
				}
			}
			aggr.push_back(std::make_pair(id, std::move(ids)));
		}
	}
	
	
	std::unordered_map<uint64_t, std::vector<unsigned int> > pop; // all population
	unsigned int last_year0 = 0;
	
	{
		read_table2 rt(std::cin);
		bool end_year = false;
		while(rt.read_line()) {
			unsigned int year, pop1;
			uint64_t id;
			if(!rt.read(year, id, pop1)) break;
			
			if(aggr.size() == 0 && cnts.size() > 0 && cnts.count(id) == 0) continue; // we only care about cells where we need to generate events
			
			if(year < last_year0) throw std::runtime_error("Input is not ordered by year!\n");
			last_year0 = year;
			
			if(year < tmin) continue;
			if(tmax > 0 && year >= tmax) { end_year = true; break; }
			year -= tmin; // we index years starting from zero
			
			auto& pop0 = pop[id];
			if(pop0.size() > year) throw std::runtime_error("Inconsistent population data!\n");
			if(year && pop0.size() < year) pop0.resize(year, 0); // add zeros to "missing" years
			if(!norm_pop) pop1 = pop1 ? 1 : 0; // in this case, only record if a cell is occupied / empty
			pop0.push_back(pop1);
			sum1 += pop1;
		}
		if( !(rt.get_last_error() == T_EOF || (rt.get_last_error() == T_OK && end_year)) ) {
			fprintf(stderr, "Error reading input:\n");
			rt.write_error(std::cerr);
			return 1;
		}
		
		last_year0++;
		if(tmax && last_year0 < tmax) fprintf(stderr, "Maximum year (%u) is less than expected!\n", last_year0);
	}
	
	double total_size = cellsize * (double)pop.size(); /* total simulation area (i.e. all cells to consider) */
	double total_events = 0.0;
	if(cnts.size()) {
		if(aggr.size()) for(const auto& x : aggr) total_events += cnts.at(x.first);
		else for(const auto& x : cnts) if(pop.count(x.first)) total_events += x.second;
		
		
		if(total_dates) {
			/* scale the counts to get the number of desired dates */
			double factor1 = ((double)total_dates) / ((double)total_events);
			double sum2 = 0.0;
			total_events = 0.0;
			
			if(aggr.size()) for(const auto& x : aggr) {
				unsigned int& y = cnts.at(x.first);
				sum2 += factor1 * (double)y;
				double tmp1 = std::round(sum2);
				y = (unsigned int)(tmp1 - total_events);
				total_events = tmp1;
			}
			else for(auto& x : cnts) if(pop.count(x.first)) {
				sum2 += factor1 * (double) x.second;
				double tmp1 = std::round(sum2);
				x.second = (unsigned int)(tmp1 - total_events);
				total_events = tmp1;
			}
		}
	}
	else if(total_dates) total_events = total_dates;
	else total_events = p1 * total_size * (double)(last_year0 - tmin); /* total number of events we expect to generate */
	double p2 = total_events / (double)sum1; /* probability of generating one event for each occurrance in the input data */
	
	fprintf(stderr, "Number of cells: %lu, total area: %f km2, total number of events expected: %f\n",
		pop.size(), total_size, total_events);
	if(!cnts.size()) fprintf(stderr, "Total weight events in the input: %lu, chance of output per input: %g\n", sum1, p2);
	
	// generate events
	std::uniform_real_distribution<double> d1;
	std::mt19937_64 rng(seed);
	
	FILE* fout = stdout;
	uint64_t n = 0;
	if(aggr.size()) {
		for(const auto& x : aggr) {
			double cnt1 = cnts.at(x.first); // total number of events in this area
			uint64_t sum_pop = 0;
			for(uint64_t id : x.second) {
				auto it = pop.find(id);
				if(it != pop.end()) for(unsigned int y : it->second) sum_pop += (uint64_t)y;
			}
			double p2 = cnt1 / (double)sum_pop; // base probability of creating an event
			for(uint64_t id : x.second) {
				auto it = pop.find(id);
				if(it != pop.end()) for(size_t i = 0; i < it->second.size(); i++) if(it->second[i]) {
					double p3 = it->second[i];
					p3 *= p2; // probability of generating one event (note: it should be << 1 for the next step to make sense)
					double y = d1(rng);
					if(y < p3) {
						fprintf(fout, "%lu\t%lu\n", id, tmin + i); // output is cell ID -- year (as per the simulation dates)
						n++;
					}
				}
			}
		}
	}
	else for(const auto& x : pop) {
		uint64_t id = x.first; // cell ID -- used in the output
		if(cnts.size()) {
			// we sample based on a prescribed number of events
			unsigned int cnt1 = cnts.at(id);
			std::discrete_distribution<unsigned int> dd(x.second.begin(), x.second.end());
			for(unsigned int i = 0; i < cnt1; i++) {
				unsigned int year = dd(rng);
				fprintf(fout, "%lu\t%u\n", id, tmin + year);
				n++;
			}
		}
		else for(size_t i = 0; i < x.second.size(); i++) if(x.second[i]) {
			double p3 = x.second[i];
			p3 *= p2; // probability of generating one event (note: it should be << 1 for the next step to make sense)
			double y = d1(rng);
			if(y < p3) {
				fprintf(fout, "%lu\t%lu\n", id, tmin + i); // output is cell ID -- year (as per the simulation dates)
				n++;
			}
		}
	}
	
	fprintf(stderr, "Total events generated: %lu\n", n);
	
	return 0;
}

