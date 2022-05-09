/*
 * cellgrid.hpp -- representation of a simulation space as a rectangular grid of cells
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


#ifndef CELLGRID_HPP
#define CELLGRID_HPP

#include "cell.hpp"
#include "cell_it_helper.hpp"
#include <vector>
#include <utility>
#include <random>
#include <algorithm>
#include <queue>
#include <unordered_set>
#include <stdexcept>
#include <stdio.h>


class cellgrid {
public:
	struct neighbor_iterator;
	
	struct point {
		unsigned int x;
		unsigned int y;
		point(unsigned int x_, unsigned int y_) : x(x_), y(y_) { }
		point() : x(0), y(0) { }
		point(const point& p) = default;
		bool operator == (const point& p) const { return (x == p.x) && (y == p.y); }
		bool operator != (const point& p) const { return (x != p.x) || (y != p.y); }
		
		// comparison operator for the possibility to include this in an ordered set
		// (used for efficient selection of random points from a set)
		// note: order has no physical meaning and should not be used in the simulation
		bool operator < (const point& p) const {
			return (x < p.x) || (x == p.x && y < p.y);
		}
		
	protected:
		// calculate the distance from a given point
		double dist(double x1, double y1) const {
			double dx = x1 - x;
			double dy = y1 - y;
			// if(dx == 0 || dy == 0) return std::abs(dx + dy);
			return std::sqrt(dx*dx + dy*dy);
		}
		
		double dist(const point& p) const {
			return this->dist(p.x, p.y);
		}
		
		static double dist(const point& p1, const point& p2) { return p1.dist(p2); }
		
		friend class cellgrid;
		friend class cellgrid::neighbor_iterator;
	};
	
	struct pointhash {
		size_t operator() (const point& p) const {
			size_t r = p.x;
			// see e.g. https://stackoverflow.com/questions/2590677/how-do-i-combine-hash-values-in-c0x
			r ^= p.y + 0x9e3779b9UL + (r<<6) + (r>>2);
			return r;
		}
	};
	
	// parse a given point given as command line arguments;
	// returns the number of arguments consumed or -1 in case of an error
	static int try_parse_point(int argc, char** argv, int i, point& p) {
		if(i + 2 > argc) return -1;
		p.x = atoi(argv[i]);
		p.y = atoi(argv[i+1]);
		return 2;
	}
	
	// convert a point to an std::tuple that the caller can handle in a more generic way 
	static std::tuple<unsigned int, unsigned int> get_point_as_tuple(const point& p) {
		return std::make_tuple(p.x, p.y);
	}
	
	double point_dist(const point& p1, const point& p2) const {
		return p1.dist(p2);
	}
	
	point dim; // size of the grid (should be given when created)
	
	std::vector<cell> cells; // all the cells
	cell& at(unsigned int x, unsigned int y) {
		if(x >= dim.x || y >= dim.y) throw std::runtime_error("cellgrid::at(): index out of range!\n");
		return cells[x * dim.y + y];
	}
	const cell& at(unsigned int x, unsigned int y) const {
		if(x >= dim.x || y >= dim.y) throw std::runtime_error("cellgrid::at(): index out of range!\n");
		return cells[x * dim.y + y];
	}
	cell& at(const point& p) { return at(p.x, p.y); }
	const cell& at(const point& p) const { return at(p.x, p.y); }
	cell& at_ix(size_t i) { return cells.at(i); }
	const cell& at_ix(size_t i) const { return cells.at(i); }
	bool cell_exists(const point& p) const { return p.x < dim.x && p.y < dim.y; }
	size_t size() const {
		size_t ret = dim.x;
		ret *= dim.y;
		return ret;
	}
	/* return the ID of the ith point -- it is assumed that there
	 * exists a linear ordering of all points */
	point ptid(size_t i) const {
		if(i >= size()) throw std::runtime_error("cellgrid::ptid(): requested index is out of range!\n");
		point p;
		p.x = i / dim.y;
		p.y = i % dim.y;
		return p;
	}
	/* returns the index of the given point */
	size_t ptix(const point& pt) const {
		if(!(pt < dim)) throw std::runtime_error("cellgrid::ptix(): invalid point!\n");
		size_t i = pt.x;
		i *= dim.y;
		i += pt.y;
		return i;
	}
	
	explicit cellgrid(const point& dim_) : dim(dim_) { init(); }
	explicit cellgrid(unsigned int X = 1, unsigned int Y = 1) : dim(X, Y) { init(); }
	
	// empty sentinel struct for end of iteration
	struct sentinel { };
	
	// iterator over a set of cells within a given distance from a center
	struct neighbor_iterator {
		protected:
			point target;
			/* note: distance is stored along with the point */
			std::pair<point, double> current;
			point dim;
			bool is_end1 = false;
			const double max_dist;
			const bool include_self;
			
			struct node {
				point pt;
				double dist;
				
				bool operator < (const node& n) const {
					return (dist > n.dist) || (dist == n.dist && n.pt < pt);
				}
				
				operator std::pair<point, double> () const { return {pt, dist}; }
			};
			
			std::unordered_set<point, pointhash> visited; /* already visited points */
			std::priority_queue<node> queue; /* queue of nodes to visit next */
			std::unordered_set<point, pointhash> queued; /* points already in the queue */
			
			friend class cellgrid;
			
			neighbor_iterator(const point& dim_, const point& target_, double max_dist_, bool include_self_)
					: target(target_), dim(dim_), max_dist(max_dist_), include_self(include_self_) {
				current.first = target;
				current.second = 0.0;
				if(!include_self) advance();
			}
			
			
			void advance() {
				if(is_end1) return;
				visited.insert(current.first);
				for(unsigned int i = 0; i < 4; i++) {
					point next = current.first;
					switch(i) {
						case 0:
							// left neighbor
							if(!next.x) continue;
							next.x--;
							break;
						case 1:
							// top neighbor
							if(!next.y) continue;
							next.y--;
							break;
						case 2:
							// right neighbor
							next.x++;
							if(next.x >= dim.x) continue;
							break;
						case 3:
							// bottom neighbor
							next.y++;
							if(next.y >= dim.y) continue;
							break;
					}
					if(visited.count(next) || queued.count(next)) continue;
					double dist = target.dist(next);
					if(dist > max_dist) continue;
					queue.push({next, dist});
					queued.insert(next);
				}
				
				if(queue.empty()) {
					is_end1 = true;
					return;
				}
				
				current = queue.top();
				queue.pop();
				queued.erase(current.first);
			}
			
		public:
			const std::pair<point, double>& operator *() { return current; }
			const std::pair<point, double>* operator ->() { return &current; }
			neighbor_iterator& operator ++() { advance(); return *this; }
			// neighbor_iterator operator ++(int) { neighbor_iterator tmp = *this; advance(); return tmp; }
			
			/* iterators should not be copied -- it is an expensive operation */
			neighbor_iterator(const neighbor_iterator& it) = delete;
			/* move should be possible, but non-trivial and is not needed */
			neighbor_iterator(neighbor_iterator&& it) = delete;
			
			bool is_end() {
				return is_end1;
			}
/*			bool operator == (const neighbor_iterator& it) { return target == it.target; }
			bool operator != (const neighbor_iterator& it) { return target != it.target; } */
			bool operator == (const sentinel& s) { return is_end(); }
			bool operator != (const sentinel& s) { return !is_end(); }
	};
	
	neighbor_iterator neighbors_start(const point& p, double max_dist, bool include_self = false) const {
		return neighbor_iterator(dim, p, max_dist, include_self);
	}
	sentinel neighbors_end() const {
		return sentinel();
	}
	
	struct neighbor_range {
		protected:
			const cellgrid& g;
			const point p;
			const double max_dist;
			const bool include_self;
		public:
			neighbor_range(const cellgrid& g_, const point& p_, double max_dist_, bool include_self_)
				: g(g_), p(p_), max_dist(max_dist_), include_self(include_self_) { }
			cellgrid::neighbor_iterator begin() { return g.neighbors_start(p, max_dist); }
			cellgrid::sentinel end() { return g.neighbors_end(); }
	};
	
	neighbor_range neighbors(const point& p, double max_dist, bool include_self = false) const {
		return neighbor_range(*this, p, max_dist, include_self);
	}
	
	// step the demographic simulation in each cell, according to the given parameters
	template<class RNG>
	void step(double t, RNG& rng, const cell::cell_pars& pars = cell::default_pars())
	{
		for(cell& c : cells) if(c.st != cell::state::EMPTY) {
			c.step(t, rng, pars);
			if(c.N == 0) c.st = cell::state::EMPTY;
		}
	}
	
	// adjust the carrying capacity of each cell randomly based on a given fraction
	template<class rng_type>
	void K_rw(rng_type& rng, double r = 0.05, double Kmin = 1000.0, double Kmax = 1000000.0) {
		const double rp = 1.0 + r;
		const double rm = 1.0 - r;
		std::uniform_int_distribution<int> dst = std::uniform_int_distribution<int>(0, 1);
		for(cell& c : cells) {
			if(dst(rng)) c.K = std::min(Kmax, c.K * rp);
			else c.K = std::max(Kmin, c.K * rm);
		}
	}
	
	// iterate over all cells -- iterators dereference to a const point& -- not possible to modify
	using iterator = it_base<cellgrid>;
	iterator begin() const { return iterator(*this, 0); }
	iterator cbegin() const { return iterator(*this, 0); }
	iterator end() const { return iterator(*this, size()); }
	iterator cend() const { return iterator(*this, size()); }
	
	/* Output in PGM format in the given file */
	void write_pgm(const char* fn, unsigned int max_pix_val) {
		FILE* f = fopen(fn, "w");
		if(max_pix_val > 65535) throw std::runtime_error("cellgrid::write_pgm(): Invalid maximum pixel value!\n");
		fprintf(f, "P5\n%u %u\n%u\n", dim.x, dim.y, max_pix_val);
		
		if(max_pix_val < 256) write_rows<uint8_t>(f);
		else write_rows<uint16_t>(f);
		fclose(f);
	}
	
	/* output in PPM (RGB) format, with the color corresponding to cell type
	 * (farmer: red, raider: blue)
	 * max_pop is the possible maximum population value in the grid
	 * if raw == true, do not write a header, only the pixel data */
	void write_ppm(FILE* f, double max_pop, bool raw = false, double gamma = 1.0, bool raider_pop = true) const {
		const unsigned int max_pix_val = 255;
		if(!raw) fprintf(f, "P6\n%u %u\n%u\n", dim.x, dim.y, max_pix_val);
		
		std::vector<uint8_t> row;
		row.resize(dim.x * 3, 0);
		for(unsigned int y = 0; y < dim.y; y++) {
			for(unsigned int x = 0; x < dim.x; x++) {
				double N = at(x, y).N;
				N /= max_pop;
				if(gamma != 1.0) N = pow(N, gamma);
				unsigned int val = (unsigned int)std::round(max_pix_val * N);
				if(val > max_pix_val) val = max_pix_val;
				row[3*x] = 0;
				row[3*x + 1] = 0;
				row[3*x + 2] = 0;
				switch(at(x, y).st) {
					case cell::state::EMPTY:
						if(N > 0) throw std::runtime_error("cellgrid::write_ppm(): Empty cell with nonzero population!\n");
						break;
					case cell::state::FARMING:
						row[3*x + 2] = val; // blue
						break;
					case cell::state::RAIDERS:
						if(!raider_pop) val = 255;
						row[3*x] = val; // red
						break;
				}
			}
			size_t res = fwrite(row.data(), sizeof(uint8_t), 3 * dim.x, f);
			if(res != 3 * dim.x) throw std::runtime_error("cellgrid::write_ppm(): Error writing output!\n");
		}
	}
	
	/* save / load functionality */
	constexpr static uint64_t header = 0xc4cce680cc542db8UL;
	
	/* save the parameters of the grid to the given binary file;
	 * returns the number of bytes written, 0 on error
	 * (note that this does not save the carrying capacities) */
	size_t write(FILE* f) const {
		size_t ret = 0;
		size_t r1 = fwrite(&header, sizeof(uint64_t), 1, f);
		if(r1 != 1) return 0;
		ret += sizeof(uint64_t);
		r1 = fwrite(&dim.x, sizeof(uint32_t), 1, f);
		if(r1 != 1) return 0;
		ret += sizeof(uint32_t);
		r1 = fwrite(&dim.y, sizeof(uint32_t), 1, f);
		if(r1 != 1) return 0;
		ret += sizeof(uint32_t);
		return ret;
	}
	
	/* load the parameters from the given memory area (read from the file written by write())
	 * and initialize this class; note that this does not read carrying capacities (they are initialized to zero)
	 * returns the number of bytes read, or 0 on error */
	size_t read(const void* base, size_t base_size) {
		if(base_size < 2*sizeof(uint32_t) + sizeof(uint64_t)) return 0;
		const uint64_t* tmp = (const uint64_t*)base;
		if(*tmp != header) return 0;
		const uint8_t* tmp2 = (const uint8_t*)base;
		tmp2 += sizeof(uint64_t);
		const uint32_t* tmp3 = (const uint32_t*)tmp2;
		dim.x = tmp3[0];
		dim.y = tmp3[1];
		init();
		return 2*sizeof(uint32_t) + sizeof(uint64_t);
	}
	
protected:
	void init() {
		size_t size = dim.x;
		size *= dim.y;
		cells.resize(0); // erase any cells already stored
		cells.resize(size); // initialize all cells to default / zero state
	}
	
	/* helper for writing data as PGM */
	template<class T> void write_rows(FILE* f) {
		std::vector<T> row;
		row.resize(dim.x, 0);
		
		for(unsigned int y = 0; y < dim.y; y++) {
			for(unsigned int x = 0; x < dim.x; x++) {
				uint32_t N = at(x, y).N;
				if(N > std::numeric_limits<T>::max()) throw std::runtime_error("cellgrid::write_rows(): Too large population value!\n");
				row[x] = N;
			}
			size_t res = fwrite(row.data(), sizeof(T), dim.x, f);
			if(res != dim.x) throw std::runtime_error("cellgrid::write_rows(): Error writing output!\n");
		}
	}
};

#endif

