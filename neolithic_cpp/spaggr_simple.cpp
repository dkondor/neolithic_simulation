/*
 * spaggr_simple.cpp -- aggregate simulation results based on a
 * many-to-many matching from DGGRID IDs to region IDs; one grid cell
 * can belong to multiple regions, but region IDs need to be sequential
 * integers
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
#include "read_table_cpp.h"
#include <vector>
#include <utility>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>



int main(int argc, char **argv)
{
	const char* matches_fn = nullptr; /* list of dgid -- region ID pairs */
	const char* filter_fn = nullptr; /* additional filter of hexagons to include */
	std::unordered_map<unsigned int, unsigned int> ptix; /* store a matching from dgids to numeric indices */
	std::vector<unsigned int> region_ids; /* region IDs to aggregate to */
	std::vector<unsigned int> region_ix; /* indexes into the above, based on dgid
		i.e. for cell x, the corresponding regions are region_ids[region_ix[ptix[x]]] -- region_ids[region_ix[ptix[x]+1]] */
	std::vector<double> res; /* result of the aggregation in each time step (sum of the input) */
	std::vector<unsigned int> cnts; /* count of nonempty cell in each time step) */
	bool separate_raiders = false;
	std::vector<double> res2; /* raider population separately (only if the above is true) */
	std::vector<unsigned int> cnts2; /* raider cells separately */
	
	for(int i = 1; i < argc; i++) if(argv[i][0] == '-') switch(argv[i][1]) {
		case 'm':
			matches_fn = argv[i+1];
			i++;
			break;
		case 'f':
			filter_fn = argv[i+1];
			i++;
			break;
		case 'r':
			separate_raiders = true;
			break;
		default:
			fprintf(stderr, "Unknown argument: %s!\n", argv[i]);
			break;
	}
	else fprintf(stderr, "Unknown argument: %s!\n", argv[i]);
	
	if(!matches_fn) {
		fprintf(stderr, "No matches given!\n");
		return 1;
	}
	
	unsigned int maxr = 0; /* maximum region ID */
	unsigned int ng = 0; /* number of grid cells */
	bool have_filter = false;
	{
		std::unordered_set<unsigned int> dgfilter;
		if(filter_fn) {
			have_filter = true;
			read_table2 rtf(filter_fn);
			while(rtf.read_line()) {
				unsigned int id;
				if(!rtf.read(id)) break;
				dgfilter.insert(id);
			}
			if(rtf.get_last_error() != T_EOF) {
				fprintf(stderr, "Error reading filter IDs:\n");
				rtf.write_error(std::cerr);
				return 1;
			}
		}
		
		read_table2 rtm(matches_fn);
		rtm.set_delim(','); /* CSV input */
		rtm.read_line(); /* header */
		std::vector<std::pair<unsigned int, unsigned int> > m;
		
		while(rtm.read_line()) {
			unsigned int x, y;
			if(!rtm.read(x, y)) break;
			if(have_filter && !dgfilter.count(x)) continue;
			if(y > maxr) maxr = y;
			m.push_back(std::make_pair(x, y));
		}
		if(rtm.get_last_error() != T_EOF) {
			fprintf(stderr, "Error reading matches:\n");
			rtm.write_error(std::cerr);
			return 1;
		}
		
		std::sort(m.begin(), m.end()); /* sort by dgid */
		
		/* copy the results */
		region_ids.resize(m.size(), 0);
		res.resize(maxr + 1, 0.0);
		cnts.resize(maxr + 1, 0);
		if(separate_raiders) {
			res2.resize(maxr + 1, 0.0);
			cnts2.resize(maxr + 1, 0);
		}
		for(size_t i = 0; i < m.size(); i++) {
			if(i == 0 || m[i].first != m[i-1].first) {
				/* new dgid */
				ptix[m[i].first] = ng;
				region_ix.push_back(i);
				ng++;
			}
			region_ids[i] = m[i].second;
		}
		region_ix.push_back(m.size()); /* sentinel element */
	}
	
	
	/* read the input, perform the aggregation */
	unsigned int t = 0; /* current year -- we assume only nonnegative numbers */
	read_table2 rt(std::cin);
	FILE* fout = stdout;
	
	while(true) {
		unsigned int y, dgid; /* year, dgid */
		unsigned int N; /* population to aggregate */
		int raiders = 0; /* whether this cell is occupied by raiders */
		
		bool have_line = rt.read_line();
		if(have_line) {
			if(!rt.read(y, dgid, N)) break;
			if(separate_raiders && !rt.read(raiders)) break;
		}
		
		if(!have_line || y != t) {
			/* end of input or new year -- write out current results */
			for(unsigned int j = 0; j <= maxr; j++) if(cnts[j] || (separate_raiders && cnts2[j])) {
				fprintf(fout, "%u\t%u\t%u\t%f", t, j, cnts[j], res[j]);
				if(separate_raiders) {
					fprintf(fout, "\t%u\t%f\n", cnts2[j], res2[j]);
					cnts2[j] = 0;
					res2[j] = 0.0;
				}
				else fputc('\n', fout);
				cnts[j] = 0;
				res[j] = 0.0;
			}
			if(have_line) t = y;
			else break;
		}
		
		auto it = ptix.find(dgid);
		if(it == ptix.end()) continue;
		unsigned int pt = it->second; /* index of the current point */
		for(unsigned int j = region_ix[pt]; j < region_ix[pt+1]; j++) {
			unsigned int r = region_ids[j];
			if(separate_raiders && raiders) {
				cnts2[r]++;
				res2[r] += N;
			}
			else {
				cnts[r]++;
				res[r] += N;
			}
		}
	}
	
	if(rt.get_last_error() != T_EOF) {
		fprintf(stderr, "Error reading input:\n");
		rt.write_error(std::cerr);
		return 1;
	}
	
	return 0;
}



