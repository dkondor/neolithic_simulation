/*
 * runmultiple.cpp -- Execute a set of shell commands in parallel
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


#include <iostream>
#include <fstream>
#include <istream>
#include <thread>
#include <atomic>
#include <string>
#include <vector>
#include <stdlib.h>

int main(int argc, char **argv)
{
	char* fn = 0;
	int nthreads = -1;
	for(int i = 1; i < argc; i++) {
		if(argv[i][0] == '-') switch(argv[i][1]) {
			case 'f':
				fn = argv[i+1];
				i++;
				break;
			case 't':
				nthreads = atoi(argv[i+1]);
				if(!nthreads) {
					fprintf(stderr, "Invalid value for the number of threads!\n");
					return 1;
				}
				i++;
				break;
			default:
				fprintf(stderr, "Unknown parameter: %s!\n", argv[i]);
				break;
		}
		else fprintf(stderr, "Unknown parameter: %s!\n", argv[i]);
	}
	
	if(nthreads < 0) {
		nthreads = std::thread::hardware_concurrency();
		if(nthreads <= 0) nthreads = 1;
	}
	
	std::ifstream fs;
	if(fn) fs.open(fn);
	
	std::istream& is = fn ? fs : std::cin;
	
	std::vector<std::string> cmds;
	for(std::string l; std::getline(is, l);) cmds.push_back(std::move(l));
	
	if(nthreads == 1) {
		for(const auto& cmd : cmds)
			if(system(cmd.c_str()) < 0) throw std::runtime_error("Error running the child process!\n");
	}
	else {
		std::atomic<size_t> j{0};
		std::vector<std::thread> t;
		t.reserve(nthreads);
		auto fn = [&cmds, &j]() {
			while(true) {
				size_t x = j++;
				if(x >= cmds.size()) break;
				if(system(cmds[x].c_str()) < 0) throw std::runtime_error("Error running the child process!\n");
			}
		};
		for(int i = 0; i < nthreads; i++) t.emplace_back(fn);
		for(int i = 0; i < nthreads; i++) t[i].join();
	}
	
	return 0;
}

