/*
 * cellnet.hpp -- representation of a simulation space as arbitrary regions with neighbor relations
 * 	(main target is a hexagonal grid)
 * TODO: it is implicitely assumed that all regions have the same size
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


#ifndef CELLNET_HPP
#define CELLNET_HPP

#include "cell.hpp"
#include "cell_it_helper.hpp"
#include "read_table_cpp.h"
#include "spdists.hpp"
#include <iostream>
#include <vector>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <utility>
#include <random>
#include <math.h>
#include <stdexcept>
#include <stdint.h>
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/**
 * cellnet: organized a set of cells to a network where edges are present among neighbors
 * 
 * Edges can be arbitrary (it is not checked that they represent a planar graph; so they
 * could be used to represent more abstract connection, e.g. along a sea coast).
 * 
 * This class contains functionality related to storing a set of cells, reading network
 * connections and calculating distances among cells (either as great circle distance, or
 * with a Dijkstra-search along the network connections). Currently, arbitrary distance
 * metrics are not supported, but these could be added later.
 * 
 * Basic properties of a cell network (coordinates, connections, etc.) can be saved in a
 * binary file for easier later use.
 */
class cellnet {
public:
	typedef uint32_t point; // points are just IDs that represent regions (we use this typedef to be more generic; see cellgrid.hpp as well)

	// hash function to be used when storing a set of cells (e.g. in an unordered_map)
	struct pointhash {
		size_t operator() (const point& p) const {
			return p; // points are 32-bit integers, identity should work
		}
	};
	
	// parse a given point given as command line arguments;
	// returns the number of arguments consumed or -1 in case of an error
	static int try_parse_point(int argc, char** argv, int i, point& p) {
		if(i >= argc) return -1;
		p = atoi(argv[i]);
		return 1;
	}
	
	// convert a point to an std::tuple that the caller can handle in a more generic way 
	static std::tuple<uint32_t> get_point_as_tuple(const point& p) {
		return std::make_tuple(p);
	}
	
protected:
	
	// store cells in a vector -- this is the main container
	std::vector<cell> cells;
	// separate vector for neighbors -- optionally store distances along the edges as well
	std::vector<std::unordered_map<point, double> > n;
	// flag indicating whether the network has valid edge distances
	bool have_edge_dist = false;
	// separate vector for coordinates
	std::vector<std::pair<double, double> > coords;
	// separate vector for polygons (used when creating images); this is not required to run the simulation
	std::vector<std::vector<std::pair<double, double> > > poly;
	// store a sorted list of valid point IDs -- this is used to easily access points by index
	std::vector<point> ptids;
	// optional reverse hashmap for getting simple IDs
	std::unordered_map<point, unsigned int> ptidsr;
	// optionally a collection of coastal cells -- travel between two coastal cells is faster by the factor given below
	std::unordered_set<point> coastal_cells;
	double coastal_travel_factor = 1.0;
	// optionally a collection of land travel weights (by cell); this can be used to model some terrain features
	std::unordered_map<point, double> land_travel_factors;
	double get_land_travel_factor(const point& pt) const {
		auto it = land_travel_factors.find(pt);
		if(it != land_travel_factors.end()) return it->second;
		else return 1.0;
	}
	
	// distance matrix (optional, if great circle distances are not used)
	spdists<uint16_t> spd;
	bool use_spd = false; // set to true if the above should be used
	
	static void rt_error(read_table2& rt) { rt.write_error(std::cerr); }
	
	/* Helper function to sort the ptids, i.e. sort the ptids array, update the ptidsr mapping,
	 * sort the coordinates and optionally the cell carrying capacities (if sort_K == true).
	 * Note: this does not sort the cell contents (e.g. population, state); cells are assumed to
	 * be empty. This function should be run before starting the simulation, typically as part
	 * of reading the network. */
	void sort_by_ptid(const bool sort_K = false) {
		// sort the IDs so that the order is predictable
		std::sort(ptids.begin(), ptids.end());
		
		// rearrange the coords and cells arrays according to the IDs location as well
		size_t N = ptids.size();
		for(size_t i = 0; i < N; i++) {
			size_t j = ptidsr.at(ptids[i]); // j is the old index of this point
			if(i == j) continue; // if already in the right location, skip
			size_t k = i;
			double tmpK = sort_K ? cells[k].K : 0.0;
			auto tmpc = coords[k];
			do {
				if(j == k) throw std::runtime_error("cellnet::load_net(): inconsistent state while sorting points!\n");
				// update in a loop: we have to move data from position j to k
				if(sort_K) cells[k].K = cells[j].K;
				coords[k] = coords[j];
				ptidsr.at(ptids[k]) = k; // store the new index for this point
				k = j; // advance the iteration
				j = ptidsr.at(ptids[k]); // j is the old index of point at k
			} while(j != i); // end of a loop
			if(sort_K) cells[k].K = tmpK; // here i is the old index of the point at k
			coords[k] = tmpc;
			ptidsr.at(ptids[k]) = k; // store the new index for this point
		}
	}
	
public:
	
	// access a cell with a given ID -- will throw an exception if no such cell exists
	cell& at(const point& p) { return cells.at(ptix(p)); }
	const cell& at(const point& p) const { return cells.at(ptix(p)); }
	// access a cell based on its index
	cell& at_ix(size_t i) { return cells.at(i); }
	const cell& at_ix(size_t i) const { return cells.at(i); }
	// check if a cell exists
	bool cell_exists(const point& p) const {
		if(ptidsr.size()) return ptidsr.count(p);
		auto it = std::lower_bound(ptids.begin(), ptids.end(), p);
		return it != ptids.end() && *it == p;
	}
	// return the total number of cells
	size_t size() const { return cells.size(); }
	// return the centre coordinates of a cell
	const std::pair<double, double>& get_coords(const point& p) const { return coords.at(ptix(p)); }
	/* return the ID of the ith point -- it is assumed that there
	 * exists a linear ordering of all points */
	point ptid(size_t i) const {
		if(i >= size()) throw std::runtime_error("cellnet::ptid(): requested index is out of range!\n");
		return ptids[i];
	}
	/* returns the index of the given point */
	size_t ptix(const point& pt) const {
		if(ptidsr.size()) return ptidsr.at(pt);
		auto it = std::lower_bound(ptids.begin(), ptids.end(), pt);
		if(it == ptids.end() || *it != pt) throw std::runtime_error("cellned::ptid() invalid cell ID!\n");
		return (it - ptids.begin());
	}
	
	// empty sentinel struct for end of iteration
	struct sentinel { };
	
	// iterator among the neighbors of one cell, ordered by distance;
	// this performs a Dijkstra-search on the network, until the given maximum distance is reached
	struct neighbor_iterator {
		protected:
			const cellnet& parent;
			const point& target;
			const double max_dist;
			std::pair<point, double> current;
			bool is_end1 = false;
			
			struct node {
				point pt;
				double dist;
				
				bool operator < (const node& n) const {
					return (dist < n.dist) || (dist == n.dist && pt < n.pt);
				}
			};
			
			std::unordered_set<point> visited; /* already visited points */
			std::set<node> queue; /* queue of nodes to visit next */
			std::unordered_map<point, double> dist_map; /* map of currently assigned distances */
			
			neighbor_iterator(const cellnet& parent_, const point& target_, double max_dist_, bool include_self = false) :
					parent(parent_), target(target_), max_dist(max_dist_) {
				current.first = target;
				current.second = 0.0;
				if(!include_self) advance();
			}
			
			void advance() {
				if(is_end1) return;
				
				visited.insert(current.first);
				for(const auto& y : parent.n.at(parent.ptix(current.first))) {
					point x = y.first;
					if(visited.count(x)) continue;
					double dist1 = parent.have_edge_dist ? current.second + y.second : parent.point_dist(target, x);
					auto it = dist_map.find(x);
					bool queue_insert = false;
					if(it != dist_map.end()) {
						if(dist1 < it->second) {
							queue.erase({x, it->second});
							queue_insert = true;
						}
					}
					else queue_insert = true;
					
					if(queue_insert) {
						queue.insert({x, dist1});
						dist_map[x] = dist1;
					}
				}
				
				if(queue.empty()) {
					is_end1 = true;
					return;
				}
				
				auto it2 = queue.begin();
				current.first = it2->pt;
				current.second = it2->dist;
				queue.erase(it2);
				
				if(current.second > max_dist) is_end1 = true;
			}
			
			friend class cellnet;
		
		public:
			const std::pair<point, double>& operator *() { return current; }
			const std::pair<point, double>* operator ->() { return &current; }
			neighbor_iterator& operator ++() { advance(); return *this; }
			
			/* iterators should not be copied (it is an expensive operation) */
			neighbor_iterator(const neighbor_iterator& it) = delete;
			/* move should be possible, but non-trivial and is not needed */
			neighbor_iterator(neighbor_iterator&& it) = delete;
			
			bool is_end() {
				return is_end1;
			}
/*			bool operator == (const neighbor_iterator& it) { return target == it.target; }
			bool operator != (const neighbor_iterator& it) { return target != it.target; } */
			bool operator == (const sentinel&) { return is_end(); }
			bool operator != (const sentinel&) { return !is_end(); }
			
	};
	
	friend struct neighbor_iterator;
	
	/* start an iteration from the given cell (p) that will run until the given maximum distance; */
	neighbor_iterator neighbors_start(const point& p, double max_dist, bool include_self = false) const {
		return neighbor_iterator(*this, p, max_dist, include_self);
	}
	sentinel neighbors_end() const { return sentinel(); }
	
	/* helper struct to allow C++11 style range for loops to iterate over cells from a given start node */
	struct neighbor_range {
		protected:
			const cellnet& g;
			const point p;
			const double max_dist;
			const bool include_self;
		public:
			neighbor_range(const cellnet& g_, const point& p_, double max_dist_, bool include_self_) :
				g(g_), p(p_), max_dist(max_dist_), include_self(include_self_) { }
			cellnet::neighbor_iterator begin() const { return g.neighbors_start(p, max_dist, include_self); }
			cellnet::sentinel end() const { return g.neighbors_end(); }
	};
	
	/* get all cells in the neighborhood of p, sorted by distance, up to the given maximum distance;
	 * this returns a "range" that can be iterated over with a C++11 range for loop */
	neighbor_range neighbors(const point& p, double max_dist, bool include_self = false) const {
		return neighbor_range(*this, p, max_dist, include_self);
	}
	
	// step the demographic simulation in each cell, according to the given parameters
	// the given callback object is called for each case when a cell reverts to empty (due to its population falling to zero)
	// sequential version, requiring one random number generator
	template<class RNG, class FUN>
	void step_cb(double t, RNG& rng, FUN&& cb, const cell::cell_pars& pars = cell::default_pars())
	{
		for(size_t i = 0; i < cells.size(); i++) {
			cell& c = cells[i];
			if(c.st != cell::state::EMPTY) {
				c.step(t, rng, pars);
				if(c.N == 0) {
					c.st = cell::state::EMPTY;
					cb(ptids[i]);
				}
			}
		}
	}
	
	struct noop_cb {
		void operator()(const point&) const { }
	};
	
	template<class RNG>
	void step(double t, RNG& rng, const cell::cell_pars& pars = cell::default_pars()) {
		step_cb(t, rng, noop_cb(), pars);
	}
	
	// adjust the carrying capacity of each cell randomly based on a given fraction
	// (not used anymore in the simulation)
	template<class rng_type>
	void K_rw(rng_type& rng, double r = 0.05, double Kmin = 1000.0, double Kmax = 1000000.0) {
		const double rp = 1.0 + r;
		const double rm = 1.0 - r;
		std::uniform_int_distribution<int> dst = std::uniform_int_distribution<int>(0, 1);
		for(auto& c : cells) {
			if(dst(rng)) c.K = std::min(Kmax, c.K * rp);
			else c.K = std::max(Kmin, c.K * rm);
		}
	}
	
	/* iterate over all cells in the collection; it is more efficient to use this interface if
	 * all points should be iterated and order is not relevant (i.e. there is no need to order
	 * them based on their distance from any other point)
	 * Note that iterators are proxies and dereference to an rvalue point (not possible to modify their value) */
	using iterator = it_base<cellnet>;
	iterator begin() const { return iterator(*this, 0); }
	iterator cbegin() const { return iterator(*this, 0); }
	iterator end() const { return iterator(*this, size()); }
	iterator cend() const { return iterator(*this, size()); }
	
	// helper function to calculate the distance (in kilometers) between two points given by their geographic coordinates
	static double coords_dist(const std::pair<double, double>& c1, const std::pair<double, double>& c2) {
		const double RADIUS = 6371.0;
		double lon1 = M_PI * c1.first / 180.0;
		double lat1 = M_PI * c1.second / 180.0;
		double lon2 = M_PI * c2.first / 180.0;
		double lat2 = M_PI * c2.second / 180.0;

		double s1 = std::sin ((lat1 - lat2) / 2.0);
		double s2 = std::sin ((lon1 - lon2) / 2.0);
		double r1 = s1 * s1 + s2 * s2 * std::cos (lat1) * std::cos (lat2);
		return RADIUS * 2.0 * std::asin (std::sqrt (r1));
	}
	
	/* Get the geographic distance between two cells with the given IDs.
	 * This always uses the great circle distance */
	double coords_dist(point p1, point p2) const {
		size_t ix1 = ptix(p1);
		size_t ix2 = ptix(p2);
		const auto& c1 = coords[ix1];
		const auto& c2 = coords[ix2];
		return coords_dist(c1, c2);
	}
	
	/* Get the distance among two cells (with the given IDs: p1 and p2).
	 * This uses the previously calculated network-based distances only if create_distance_matrix() was called before
	 * (or the distance matrix was read from a binary file as well). Otherwise it uses the great circle distance. */
	double point_dist(point p1, point p2) const {
		size_t ix1 = ptix(p1);
		size_t ix2 = ptix(p2);
		if(use_spd) return (double)spd.get_dist(ix1, ix2);
		
		const auto& c1 = coords[ix1];
		const auto& c2 = coords[ix2];
		return coords_dist(c1, c2);
	}
	
	struct load_net_flags {
		constexpr static uint32_t read_K = 1; // read carrying capacities along with cell coordinates
		constexpr static uint32_t read_edges = 2; // network is provided as an edgelist (two nodes per line; if not given, use the DGGRID neighbors format)
		constexpr static uint32_t read_edge_dists = 4; // read distances along edges (needs to set read_edges as well)
		constexpr static uint32_t read_types = 8; // read cell "types"; currently only suppored type is 1, which means coastal cells
		constexpr static uint32_t read_travel_factors = 16; // read land travel weighting factors for cells as well
	};
	
	/* Read a set of cells (neighbor connections and center coordinates) from the given files opened with
	 * the read_table2 helper class.
	 * flags can be the combination of the values defined above as part of load_net_flags */
	void load_net(read_table2& rt_net, read_table2& rt_coords, uint32_t flags) {
		// clear existing data
		cells.clear();
		n.clear();
		ptids.clear();
		have_edge_dist = flags & load_net_flags::read_edge_dists;
		
		// load the coordinates first
		size_t N = 0;
		bool need_sort = false; // set to true if the input is not sorted by cell ID
		while(rt_coords.read_line()) {
			uint32_t id;
			std::pair<double, double> c1;
			double K = 0.0;
			int type = 0;
			double land_factor = 1.0;
			if(!rt_coords.read(id, read_bounds_coords(c1))) break;
			if(flags & load_net_flags::read_K && !rt_coords.read(K)) break;
			if(flags & load_net_flags::read_types && !rt_coords.read(type)) break;
			if(flags & load_net_flags::read_travel_factors && !rt_coords.read(land_factor)) break;
			
			if(ptidsr.count(id)) {
				rt_error(rt_coords);
				throw std::runtime_error("cellnet::load_net(): duplicate coordinates!\n");
			}
			
			coords.push_back(c1);
			cells.emplace_back(K);
			ptids.push_back(id);
			ptidsr[id] = N;
			if(N > 0 && ptids[N] < ptids[N-1]) need_sort = true;
			
			if(flags & load_net_flags::read_types && type == 1) coastal_cells.insert(id);
			if(flags & load_net_flags::read_travel_factors) land_travel_factors[id] = land_factor;
			
			N++;
		}
		if(rt_coords.get_last_error() != T_EOF) {
			rt_error(rt_coords);
			throw std::runtime_error("cellnet::load_net(): error reading coordinates!\n");
		}
		
		coords.shrink_to_fit();
		ptids.shrink_to_fit();
		cells.shrink_to_fit();
		
		if(need_sort) sort_by_ptid(flags & load_net_flags::read_K);
		
		// load the network, only include cells that have coordinates
		n.resize(N);
		while(rt_net.read_line()) {
			if(flags & load_net_flags::read_edges) {
				uint32_t id1, id2;
				double dist = -1.0;
				if(!rt_net.read(id1, id2)) break;
				auto it1 = ptidsr.find(id1);
				if(it1 == ptidsr.end() || !ptidsr.count(id2)) continue;
				if(flags & load_net_flags::read_edge_dists && !rt_net.read(read_bounds(dist, 0.0, 40000.0))) break;
				n[it1->second][id2] = dist;
			}
			else {
				uint32_t id;
				if(!rt_net.read(id)) break;
				
				// ignore cells with no coordinates
				auto it1 = ptidsr.find(id);
				if(it1 == ptidsr.end()) continue;
				
				auto& tmp = n.at(it1->second);
				if(!N) throw std::runtime_error("cellnet::load_net(): too many lines!\n");
				N--;
				
				while(true) {
					uint32_t id2;
					if(!rt_net.read(id2)) {
						if(rt_net.get_last_error() != T_EOL) {
							rt_error(rt_net);
							throw std::runtime_error("cellnet::load_net(): error reading the network!\n");
						}
						break;
					}
					if(!ptidsr.count(id2)) continue;
					
					// TODO: should we try to make the network symmetric or check for duplicates here?
					tmp[id2] = -1.0;
				}
				if(tmp.size() == 0) {
					rt_error(rt_net);
					throw std::runtime_error("cellnet::load_net(): cell with no valid neighbors!\n");
				}
			}
		}
		if(rt_net.get_last_error() != T_EOF) {
			rt_error(rt_net);
			throw std::runtime_error("cellnet::load_net(): error reading neighbor connections!\n");
		}
		// TODO: check result when reading an edgelist!
		if(N && !(flags & load_net_flags::read_edges)) throw std::runtime_error("cellnet::load_net(): not all cells have neighbors!\n");
	}
	
	void load_net(read_table2&& rt_net, read_table2& rt_coords, uint32_t flags) { load_net(rt_net, rt_coords, flags); }
	void load_net(read_table2& rt_net, read_table2&& rt_coords, uint32_t flags) { load_net(rt_net, rt_coords, flags); }
	void load_net(read_table2&& rt_net, read_table2&& rt_coords, uint32_t flags) { load_net(rt_net, rt_coords, flags); }
	
	void load_net(const char* net_fn, const char* coords_fn, uint32_t flags, bool coords_file_csv = false) {
		read_table2 rt(coords_fn);
		if(coords_file_csv) {
			rt.set_delim(',');
			rt.read_line();
		}
		load_net(read_table2(net_fn), rt, flags);
	}
	
	
	/* calculate distances along the edges and store them along the network
	 * distance is based on the great circle distance for all neighbor cells,
	 * with scaling applied for coastal cells and based on land travel factors (if provided) */
protected:
	double edge_dist_helper(size_t i, size_t j) const {
		point pt1 = ptids[i];
		point pt2 = ptids[j];
		double dist = coords_dist(coords[i], coords[j]);
		if(coastal_cells.count(pt1) && coastal_cells.count(pt2)) dist *= coastal_travel_factor;
		else dist *= (get_land_travel_factor(pt1) + get_land_travel_factor(pt2)) / 2.0;
		return dist;
	}

public:
	void calculate_edge_distances(double coastal_factor = 1.0) {
		if(use_spd) throw std::runtime_error("Cannot recalculare edge distances if a distance matrix is used!\n");
		coastal_travel_factor = coastal_factor;
		for(size_t i = 0; i < ptids.size(); i++) {
			for(auto& x : n[i]) {
				point pt2 = x.first;
				size_t j = ptix(pt2);
				double dist = edge_dist_helper(i, j);
				x.second = dist;
			}
		}
		have_edge_dist = true;
	}
	
	/* read a set of additional edges from a given file, supplied as an edgelist */
	void read_edges(read_table2& rt, bool read_dist = false, bool make_symmetric = true) {
		while(rt.read_line()) {
			uint32_t id1, id2;
			double dist = -1.0;
			if(!rt.read(id1, id2)) break;
			if(read_dist && !rt.read(read_bounds(dist, 0.0, 40000.0))) break;
			
			size_t i = ptix(id1);
			size_t j = ptix(id2);
			
			if(!read_dist && have_edge_dist) dist = edge_dist_helper(i, j);
			n[i][id2] = dist;
			if(make_symmetric) n[j][id1] = dist;
		}
		if(rt.get_last_error() != T_EOF) {
			rt_error(rt);
			throw std::runtime_error("cellnet::read_edges(): error reading input!\n");
		}
	}
	void read_edges(read_table2&& rt, bool read_dist = false, bool make_symmetric = true) {
		return read_edges(rt, read_dist, make_symmetric);
	}
	void read_edges(const char* fn, bool read_dist = false, bool make_symmetric = true) {
		return read_edges(read_table2(fn), read_dist, make_symmetric);
	}
	
	/* apply a scaling factor to all data, read from the given file
	 * (format is ID, value)
	 * if zero_missing == true, values for missing cells (in the input file) are replaced with zeros */
	void scale_K(read_table2& rt, bool zero_missing = false) {
		std::unordered_set<point> all_pts;
		if(zero_missing) for(const auto& x : ptids) all_pts.insert(x);
		while(rt.read_line()) {
			uint32_t id;
			double v;
			if(!rt.read(id, v)) break;
			auto it = ptidsr.find(id);
			if(it != ptidsr.end()) {
				cells[it->second].K *= v;
				if(zero_missing) all_pts.erase(id);
			}
		}
		
		if(rt.get_last_error() != T_EOF) {
			rt_error(rt);
			throw std::runtime_error("cellned::scale_K(): error reading input!\n");
		}
		
		if(zero_missing) for(point pt : all_pts) cells[ptidsr.at(pt)].K = 0.0;
	}
	
	void scale_K(read_table2&& rt, bool zero_missing = false) { scale_K(rt, zero_missing); }
	
	/* same, but open the input file ourselves */
	void scale_K(const char* fn, bool zero_missing = false, bool is_csv = true, bool header = true) {
		read_table2 rt(fn);
		if(is_csv) rt.set_delim(',');
		if(header) rt.read_line();
		scale_K(rt, zero_missing);
	}
	
	/* scale values based on map projection distortion */
	void scale_K_lat() {
		for(size_t i = 0; i < cells.size(); i++) {
			const auto& co = coords[i];
			double factor = cos(co.second * M_PI / 180.0);
			cells[i].K *= factor;
		}
	}
	
	
	
	/* save / load functionality */
	constexpr static uint64_t bin_file_header = 0x0b5c22f6d3deb068UL;
	constexpr static uint64_t bin_file_have_spd = 1UL;
	constexpr static uint64_t bin_file_edge_dist = 2UL;
	constexpr static uint64_t bin_file_all_flags = 3UL;
	
	/* save the IDs and coordinates of grid cells and the neighbor connections
	 * to the given binary file; returns the number of bytes written, 0 on error
	 * (note that this does not save the carrying capacities) */
	size_t write(FILE* f) const {
		size_t ret = 0;
		uint64_t header1 = bin_file_header;
		if(spd.has_data()) header1 |= bin_file_have_spd;
		if(have_edge_dist) header1 |= bin_file_edge_dist;
		size_t r1 = fwrite(&header1, sizeof(uint64_t), 1, f);
		if(r1 != 1) return 0;
		ret += sizeof(uint64_t);
		
		/* number of cells and number of edges */
		size_t nnodes = cells.size();
		size_t nedges = 0;
		for(const auto& x : n) nedges += x.size();
		r1 = fwrite(&nnodes, sizeof(size_t), 1, f) + fwrite(&nedges, sizeof(size_t), 1, f);
		if(r1 != 2) return 0;
		ret += 2*sizeof(size_t);
		
		/* write all node IDs */
		for(uint32_t id : ptids) {
			r1 = fwrite(&id, sizeof(uint32_t), 1, f);
			if(r1 != 1) return 0;
			ret += sizeof(uint32_t);
		}
		
		/* pad to 8 byte boundary */
		if(nnodes % 2) {
			uint32_t tmp = 0;
			r1 = fwrite(&tmp, sizeof(uint32_t), 1, f);
			if(r1 != 1) return 0;
			ret += sizeof(uint32_t);
		}
		
		/* write the coordinates -- note that coords is in the same order as the node IDs*/
		for(const auto& y : coords) {
			double c1[2] = {y.first, y.second};
			r1 = fwrite(c1, sizeof(double), 2, f);
			if(r1 != 2) return 0;
			ret += 2*sizeof(double);
		}
		
		/* write the network */
		for(size_t i = 0; i < ptids.size(); i++) {
			uint32_t n1[2] = {ptids[i], 0};
			for(const auto& y : n[i]) {
				n1[1] = y.first;
				r1 = fwrite(n1, sizeof(uint32_t), 2, f);
				if(r1 != 2) return 0;
				ret += 2*sizeof(uint32_t);
				if(have_edge_dist) {
					/* write the distance along this edge */
					double d = y.second;
					if(fwrite(&d, sizeof(double), 1, f) != 1) return 0;
					ret += sizeof(double);
				}
			}
		}
		
		/* optionally write the distance matrix */
		if(spd.has_data()) {
			size_t ret2 = spd.write(f);
			if(!ret2) return 0;
			ret += ret2;
		}
		
		return ret;
	}
	
	/* load the parameters from the given memory area (read from the file written by write())
	 * and initialize this class; note that this does not read carrying capacities (they are initialized to zero)
	 * returns the number of bytes read, or 0 on error */
	size_t read(const void* base, size_t base_size) {
		cells.clear();
		n.clear();
		coords.clear();
		ptids.clear();
		ptidsr.clear();
		
		/* validate the header and set some options from it */
		if(base_size < sizeof(uint64_t) + 2*sizeof(size_t)) return 0;
		uint64_t tmp = *((const uint64_t*)base);
		if((tmp & ~bin_file_all_flags) != bin_file_header) return 0;
		use_spd = tmp & bin_file_have_spd;
		have_edge_dist = tmp & bin_file_edge_dist;
		
		/* read the number of nodes and edges */
		const uint8_t* base2 = (const uint8_t*)base;
		base2 += sizeof(uint64_t);
		size_t nnodes, nedges;
		{
			const size_t* tmp2 = (const size_t*)base2;
			nnodes = tmp2[0];
			nedges = tmp2[1];
			base2 += 2*sizeof(size_t);
		}
		
		/* calculate the total size */
		size_t size1 = sizeof(uint64_t) + 2*sizeof(size_t);
		size1 += nnodes * sizeof(uint32_t);
		if(nnodes % 2) size1 += sizeof(uint32_t); /* padding after the node IDs */
		size1 += nnodes * sizeof(double) * 2UL; /* coordinates */
		size1 += nedges * sizeof(uint32_t) * 2UL; /* edges */
		if(have_edge_dist) size1 += nedges * sizeof(double); /* distances along the edges */
		if(base_size < size1) return 0;
		
		/* read the node IDs */
		ptids.reserve(nnodes);
		coords.reserve(nnodes);
		cells.resize(nnodes);
		n.resize(nnodes);
		
		bool need_sort = false;
		{
			/* read the node IDs -- note: node IDs are not necessarily sorted */
			const uint32_t* tmp2 = (const uint32_t*)base2;
			for(size_t i = 0; i < nnodes; i++) {
				uint32_t id = tmp2[i];
				ptids.push_back(id);
				ptidsr[id] = i;
				if(i > 0 && ptids[i] < ptids[i-1]) need_sort = true;
			}
			base2 += nnodes * sizeof(uint32_t);
			if(nnodes % 2) base2 += sizeof(uint32_t);
		}
		
		{
			/* read the coordinates and initialize the cells */
			const double* tmp2 = (const double*)base2;
			for(size_t i = 0; i < nnodes; i++) {
				double lon = tmp2[2*i];
				double lat = tmp2[2*i + 1];
				coords.push_back(std::make_pair(lon, lat));
			}
			base2 += nnodes * sizeof(double) * 2UL;
		}
		
		/* sort the node IDs and coordinates */
		if(need_sort) sort_by_ptid(false);
		
		{
			/* read the network connections */
			for(size_t i = 0; i < nedges; i++) {
				const uint32_t* tmp2 = (const uint32_t*)base2;
				uint32_t n1 = tmp2[0];
				uint32_t n2 = tmp2[1];
				base2 += 2*sizeof(uint32_t);
				
				double d = -1.0;
				if(have_edge_dist) {
					d = *(const double*)base2;
					base2 += sizeof(double);
				}
				
				n.at(ptidsr.at(n1))[n2] = d;
			}
		}
		
		// note: spd.read() will always return nonzero or throw and exception on error
		// (TOOD: handling of errors should be more uniform!)
		if(use_spd) size1 += spd.read(base2, base_size - size1);
		else spd.reset();
		
		return size1;
	}
	
	void create_distance_matrix() {
		use_spd = false;
		spd.allocate(ptids.size(), 0.0);
		int nt1 = 1;
#ifdef _OPENMP
#pragma omp parallel
		{
			nt1 = omp_get_num_threads();
#pragma omp for
#endif
			for(size_t i = 0; i < ptids.size(); i++) {
				for(const auto& x : neighbors(ptids[i], 1e9, false)) {
					double d1 = std::round(x.second);
					if(d1 > 65535.0 || d1 < 0.0) throw std::runtime_error("cellnet::create_distance_matrix(): invalid distance!\n");
					spd.set_dist(i, ptix(x.first), (uint16_t)d1);
				}
				if(nt1 == 1) fprintf(stderr, "%lu / %lu nodes processed\r", i, ptids.size());
			}
#ifdef _OPENMP
		}
#endif
		use_spd = true;
	}
	
	/* convenience function to write the population in all cells where it is nonzero */
	
	/* write the population in all cells where it is nonzero */
	void write_pop(FILE* f, unsigned int step) const {
		for(size_t i = 0; i < cells.size(); i++) {
			const auto& c = cells[i];
			unsigned int pt = ptids.at(i);
			/* to limit the file size, only write out population if
			 * it is nonzero */
			if(c.N > 0) fprintf(f, "%u\t%u\t%u\t%u\n", step, pt, c.N, c.st == cell::state::RAIDERS ? 1 : 0);
		}
	}
};

#endif

