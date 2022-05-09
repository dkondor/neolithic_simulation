/*
 * neolithic_base.cpp -- base class for neolithic simulations with
 * 	shared functionality
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


#ifndef NEOLITHIC_BASE_HPP
#define NEOLITHIC_BASE_HPP

#include "cell.hpp"
#include "mmap.h"
#include "btprob.hpp"
#include "read_table_cpp.h"
static void rt_error(read_table2& rt) { rt.write_error(std::cerr); }
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <vector>
#include <utility>
#include <memory>
#include <mutex>
#include <thread>
#include <atomic>
#include <unordered_map>
#include <unordered_set>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

template<class cell_collection>
class nbase {
public:
	cell_collection g; /* simulation space -- should be cellgrid or cellnet */

protected:
	
	/* main parameters -- group migration */
	double cap_mig_exp = 4.0; /* exponent for split-off probabilities */
	double cap_mig_max = 0.25; /* probability for split-off if the population is equal to the carrying capacity */
	double cap_mig_min = 0.0; /* minimum share of people to trigger a split-off based on carrying capacity */
	double dunbar_limit = 150.0; /* Dunbar number for village populations */
	double dunbar_mig_exp = 2.0; /* exponent for Dunbar limit based split-off */
	double dunbar_mig_base = 0.1; /* probability for split-off if population is equal to the Dunbar-number */
	double dunbar_mig_min = 0.0; /* minimum share of people to trigger a migration (for Dunbar limit based split-off) */
	double group_mig_dist = 5.0; /* characteristic distance for group migration */
	double pempty = 1.0; /* extra factor for preference toward "emtpy" (unsettled) cells */
	bool target_only_free = false; /* if true, group migrations should only target cells where there is enough remaining capacity */
	bool group_pop_estimate = false; /* when estimating targets of group migrations, the remaining population is considered (if false, only the carrying capacity is considered) */
	double group_mig_pow = 0.0; /* if > 0.0, group migration distance weighting is done by a power-law function with this exponent,
		and group_mig_dist gives the limit on cells to consider */
	bool have_pars = false; /* set to true if any of pempty, group_mig_dist or group_mig_pow is given as a command line argument */
	unsigned int target_max_retries = 1000; /* try to select a target cell this many times for a split-off */
	bool scale_max_retries = false; /* if true, scale the above with the split-off probability */
	const bool group_mig_prob_pempty; /* whether the group_mig_prob function takes into account pempty */
	bool peaceful_migrations = false; /* if set to true, a farmer group can migrate to an occupied cell without conflict if there is
		enough remaining carrying capacity */
	
	/* demographic parameters */
	cell::cell_pars pars = cell::default_pars();
	
	/* base random number generator */
	std::mt19937_64 rng;
	std::vector<std::mt19937_64> par_rng; /* separate RNGs used by each thread if using OpenMP */
	int omp_threads = 0;
	unsigned int create_dist_threads = 1;
	
	bool use_global_neighbors = false; /* common parameter, but used by derived classes */
	
	/* more generic probability distribution interface */
	class prob_wrapper_base {
		public:
			virtual ~prob_wrapper_base() { }
			virtual double total_sum() const = 0; /* get the total sum -- used when sampling */
			virtual size_t upper_bound(double p1) const = 0; /* get the index of element corresponding to p1 in the CDF */
			virtual size_t size() const = 0; /* get the number of values actually stored -- used for saving */
			virtual double get_sum(size_t i) const = 0; /* get the ith stored value -- used for saving */
			virtual double check_sum() const = 0; /* check if the stored sums are consistent (return the difference, if any) */
			virtual prob_wrapper_base* clone() const = 0; /* create a copy of this instance */
	};
	
	template<class btprob_base>
	class btprob_wrapper : public prob_wrapper_base {
		protected:
			bool is_writable;
		public:
			btprob<btprob_base> btp;
			double total_sum() const override { return btp.total_sum(); }
			size_t upper_bound(double p1) const override { return btp.upper_bound(p1); }
			size_t size() const override { return btp.size() / 2UL; } /* TODO: this only works for btprob_cb_base */
			double get_sum(size_t i) const override { return btp.get_sum(i); }
			double check_sum() const override { return btp.check_sums(); }
			prob_wrapper_base* clone() const override {
				if(is_writable && !btp.can_clone_write())
					throw std::runtime_error("btprob_wrapper::clone(): attempt to copy an instance that is not allowed!\n");
				return new btprob_wrapper(*this);
			}
			explicit btprob_wrapper(bool should_be_writable = false) : is_writable(should_be_writable) { }
			~btprob_wrapper() { }
	};
	
	class cdf_wrapper : public prob_wrapper_base {
		protected:
			std::vector<double> cdfv;
			const double* cdf;
			const size_t size1;
			cdf_wrapper() = delete;
		public:
			cdf_wrapper(const double* cdf_, size_t size_) :
				cdf(cdf_), size1(size_) { }
			explicit cdf_wrapper(std::vector<double>&& cdfv_) :
				cdfv(std::move(cdfv_)), cdf(cdfv.data()), size1(cdfv.size()) { }
			double total_sum() const override { return cdf[size1-1]; }
			size_t upper_bound(double p1) const override {
				const double* it = std::upper_bound(cdf, cdf + size1, p1);
				size_t r = it - cdf;
				return r;
			}
			size_t size() const override { return size1; }
			double get_sum(size_t i) const override { return cdf[i]; }
			double check_sum() const override { return 0.0; } /* TODO: should there be a way to check the CDF here? */
			prob_wrapper_base* clone() const override { return new cdf_wrapper(*this); }
			~cdf_wrapper() { }
	};
	
	/* wrapper around the above for one source point */
	struct cell_prob_helper {
		/* IDs of points for which we have probabilities, ordered.
		 * Only used if use_global_neighbors == false */
		const typename cell_collection::point* neighbors = nullptr; /* pointer of neighbors, either loaded externally or stored in the vector below */
		size_t nb_size = 0; /* size of neighbors */
		std::vector<typename cell_collection::point> neighbors_vec;
		std::unique_ptr<prob_wrapper_base> btpw;
		std::vector<size_t> changed_ix;
		bool invalidate = false;
		
		cell_prob_helper() { }
		cell_prob_helper(const cell_prob_helper& h) : neighbors(h.neighbors), nb_size(h.nb_size),
			neighbors_vec(h.neighbors_vec), btpw(h.btpw->clone()) { }
		cell_prob_helper(cell_prob_helper&& h) = default;
	};
	
	std::unordered_map<typename cell_collection::point, cell_prob_helper, typename cell_collection::pointhash> cph;
	
	enum class cell_prob_helper_type { NONE, BTVEC, BTCB, CDF };
	cell_prob_helper_type use_prob_helper = cell_prob_helper_type::NONE;
	/* set to true by read_helper_matrix() -- if true, helpers are not erased and recreated when nodes are occupied / abandoned */
	bool external_helper = false;
	/* open file with an externally loaded helper matrix */
	std::shared_ptr<FileMapping> mp;
	
	explicit nbase(bool pe) : group_mig_prob_pempty(pe), mp(new FileMapping()) { rng.seed(time(0)); }
	
	auto get_cb_for_pt(const typename cell_collection::point& pt) {	
		return [this, pt] (size_t i) {
			const auto& pt2 = g.ptid(i);
			if(pt2 == pt) return 0.0; /* exclude self from targets */
			double dst = g.point_dist(pt, pt2);
			return group_mig_prob(pt2, dst);
		};
	}
	
	/* create a helper distribution for one point; it is assumed that use_prob_helper != NONE
	 * (note: if a non-NULL mutes is provided, it is locked when setting the btpw member of the helper) */
	void create_cell_prob_helper_one(const typename cell_collection::point& pt, bool should_exist = false, std::mutex* m = nullptr) {
		auto& helper = should_exist ? cph.at(pt) : cph[pt]; /* this will add the helper for this node */
		auto& nb = helper.neighbors_vec;
		std::unordered_map<typename cell_collection::point, double, typename cell_collection::pointhash> nb_dist;
		bool need_probs = (use_prob_helper == cell_prob_helper_type::BTVEC ||
			use_prob_helper == cell_prob_helper_type::CDF);
		
		if(!use_global_neighbors) {
			const double max_dist = group_mig_dist * (group_mig_pow > 0.0 ? 1.0 : 10.0);
			for(const auto& pt2 : g.neighbors(pt, max_dist)) {
				nb.push_back(pt2.first);
				if(need_probs) nb_dist.insert(pt2);
			}
			nb.shrink_to_fit();
			std::sort(nb.begin(), nb.end());
			helper.neighbors = nb.data();
			helper.nb_size = nb.size();
		}
		/* we need to generate the distances with a network search to support more complex cases */
		else if(need_probs) for(const auto& pt2 : g.neighbors(pt, 1e9)) nb_dist.insert(pt2);
		
		size_t nb_size = use_global_neighbors ? g.size() : nb.size();
		
		if(need_probs) {
			auto get_one_prob = [this, &nb, &nb_dist, &pt](size_t i) {
				const auto& pt1 = use_global_neighbors ? g.ptid(i) : nb[i];
				if(pt1 == pt) return 0.0;
				auto it = nb_dist.find(pt1);
				if(it == nb_dist.end()) return 0.0; /* unreachable nodes have zero migration probability */
				return group_mig_prob(pt1, it->second);
			};
			
			if(use_prob_helper == cell_prob_helper_type::CDF) {
				std::vector<double> pr;
				pr.resize(nb_size);
				pr[0] = get_one_prob(0);
				for(size_t i = 1; i < nb_size; i++) pr[i] = pr[i-1] + get_one_prob(i);
				if(m) m->lock();
				helper.btpw.reset(new typename nbase<cell_collection>::cdf_wrapper(std::move(pr)));
				if(m) m->unlock();
			}
			else { // BTVEC type
				auto btpw = new nbase<cell_collection>::btprob_wrapper<btprob_vec_base>(true);
				auto& btp = btpw->btp;
				btp.resize(nb_size);
				btp.set_all_probs_cb(get_one_prob);
				if(m) m->lock();
				helper.btpw.reset(btpw);
				if(m) m->unlock();
			}
		}
		else { // BTCB type
			auto btpw = new nbase<cell_collection>::btprob_wrapper<btprob_cb_base>(true);
			auto& btp = btpw->btp;
			btp.resize(nb_size);
			
			if(use_global_neighbors) btp.prob = get_cb_for_pt(pt);
			else btp.prob = [this, pt, &nb] (size_t i) {
				const auto& pt2 = nb[i]; /* note: here we know that pt != pt2 */
				double dst = g.point_dist(pt, pt2);
				return group_mig_prob(pt2, dst);
			};
			
			btp.create_sums();
			if(m) m->lock();
			helper.btpw.reset(btpw);
			if(m) m->unlock();
		}
	}
	
public:
	using point_type = typename cell_collection::point;
	using cell_collection_type = cell_collection;

	std::unordered_map<typename cell_collection::point, double, typename cell_collection::pointhash> orig_K; /* base carrying capacities of cells -- these are used if present and group_pop_estimate == false */
	
	/* parse command line options for parameters in this class;
	 * i is the current index in argv;
	 * return value:
	 * 	0 if the current option did not match
	 *  >0 if there was a match and i should be incremented by this amount
	 *  -1 if there is an error
	 * 
	 * derived classes can override this, but they should call this to
	 * parse to common options first */
	virtual int parse_option(int, char** argv, int i) {
		if(argv[i][0] != '-') return -1;
		int res = 0;
		switch(argv[i][1]) {
			case 's':
				rng.seed(atoi(argv[i+1]));
				res = 2;
				break;
			case 't':
				{
					int tmp1 = atoi(argv[i+1]);
					set_num_omp_threads(tmp1);
					if(tmp1 > 0) create_dist_threads = tmp1;
				}
				res = 2;
				break;
			case 'r':
				pars.r = atof(argv[i+1]);
				res = 2;
				break;
			case 'd':
				pars.delta = atof(argv[i+1]);
				res = 2;
				break;
			case 'G':
				group_mig_dist = atof(argv[i+1]);
				res = 2;
				have_pars = true;
				break;
			case 'P':
				group_mig_pow = atof(argv[i+1]);
				res = 2;
				have_pars = true;
				break;
			case 'E':
				pempty = atof(argv[i+1]);
				if(group_mig_prob_pempty) have_pars = true;
				res = 2;
				break;
/*			case 'e': -- not used
				empty_delay = atoi(argv[i+1]);
				res = 2;
				break; */
			case 'M':
				// max retries
				switch(argv[i][2]) {
					case 's':
						scale_max_retries = true;
						/* fallthrough */
					case 0:
						target_max_retries = atoi(argv[i+1]);
						res = 2;
						break;
					default:
						res = -1;
						break;
				}
				break;
			case 'D':
				switch(argv[i][2]) {
					case 'e':
						dunbar_mig_exp = atof(argv[i+1]);
						res = 2;
						break;
					case 'b':
						dunbar_mig_base = atof(argv[i+1]);
						res = 2;
						break;
					case 'l':
						dunbar_limit = atof(argv[i+1]);
						res = 2;
						break;
					case 'm':
						dunbar_mig_min = atof(argv[i+1]);
						res = 2;
						break;
					default:
						res = -1;
						break;
				}
				break;
			case 'C':
				switch(argv[i][2]) {
					case 'e':
						cap_mig_exp = atof(argv[i+1]);
						res = 2;
						break;
					case 'b':
						cap_mig_max = atof(argv[i+1]);
						res = 2;
						break;
					case 'm':
						cap_mig_min = atof(argv[i+1]);
						res = 2;
						break;
					default:
						res = -1;
						break;
				}
				break;
			case 'g':
				use_global_neighbors = true;
				res = 1;
				break;
		}
		return res;
	}
	
	/* estimate the base probability of migrating to cell pt2 from distance dst
	 * (not taking into account if the cell is occupied or not)
	 * this returns probabilities for empty cells -- for occupied cells, this
	 * needs to be multiplied by 1.0 / pempty */
	double group_mig_prob2(const typename cell_collection::point& pt2, const cell& c2, double dst) const {
		double p1 = 0.0;
		if(group_mig_pow > 0.0) p1 = pow(dst, -1.0 * group_mig_pow);
		else p1 = exp(-1.0 * dst / group_mig_dist);
		
		if(group_pop_estimate) p1 *= (c2.K - c2.N);
		else p1 *= orig_K.size() ? orig_K.at(pt2) : c2.K;
		return p1;
	}
	
	double group_mig_prob(const typename cell_collection::point& pt2, double dst, bool ignore_state = false) const {
		/* estiamte pt1 -> pt2 migration probability */
		const auto& c2 = g.at(pt2);
		if(group_pop_estimate && c2.N > c2.K) return 0.0;
		if(!ignore_state && target_only_free && c2.st != cell::state::EMPTY) return 0.0;
		
		double p1 = group_mig_prob2(pt2, c2, dst);
		if(!ignore_state && group_mig_prob_pempty && c2.st == cell::state::EMPTY) p1 *= pempty;
		return p1;
	}
	
	/* check group migration probabilities for the given cell;
	 * returns the size of the out-migrating group (0 means no migration)
	 * note: a random number generator is taken as a parameter, this way
	 * (1) this function can be declared const; (2) it can be called from
	 * multiple threads (e.g. when using OpenMP)
	 * returns the calculated (total) split-off probability in ps if a
	 * split-off actually occured */
	template<class RNG>
	unsigned int cell_group_mig(size_t pt_ix, RNG& rng2, double& ps) const {
		ps = 0.0;
		auto& c = g.at_ix(pt_ix);
		
		if(c.st != cell::state::FARMING || c.N == 0) return 0;
		if(cph.size()) {
			const auto& helper = cph.at(g.ptid(pt_ix));
			/* return 0 if there is nowhere to migrate to */
			if(helper.btpw && helper.btpw->total_sum() <= 0.0) return 0;
		}
				
		double pc = 0.0; // probability of split-off due to approaching carrying capacity
		if(cap_mig_exp >= 0.0) {
			double tmp = 1.0;
			if(cap_mig_exp > 0.0) {
				tmp = c.N / c.K;
				if(tmp >= cap_mig_min && cap_mig_min < 1.0)
					tmp = (tmp - cap_mig_min) / (1.0 - cap_mig_min);
				else tmp = 0.0;
				tmp = pow(tmp, cap_mig_exp);
			}
			pc = cap_mig_max * tmp;
		}
		double pd = 0.0; // probability of split-off due to the Dunbar limit
		if(dunbar_limit > 0.0 && dunbar_mig_exp >= 0.0) {
			double tmp = 1.0;
			if(dunbar_mig_exp > 0.0) {
				tmp = c.N / dunbar_limit;
				if(tmp >= dunbar_mig_min && dunbar_mig_min < 1.0)
					tmp = (tmp - dunbar_mig_min) / (1.0 - dunbar_mig_min);
				else tmp = 0.0;
				tmp = pow(tmp, dunbar_mig_exp);
			}
			pd = dunbar_mig_base * tmp;
		}
		
		if(pc == 0.0 && pd == 0.0) return 0;
		
		std::uniform_real_distribution<double> d1;
		double x1 = d1(rng2);
		double x2 = (x1 >= pc) ? d1(rng2) : 0.0;
		if(x1 >= pc && x2 >= pd) return 0;
		
		ps = pc + pd;
		/* choose a faction size */
		uint32_t N1 = c.N;
		double N2 = N1 / 2.0;
		std::poisson_distribution<unsigned int> pd1(N2);
		uint32_t N3 = pd1(rng2);
		if(N3 > N1) N3 = N1;
		return N3; /* note: N3 can be zero at this point */
	}
	
	
	struct group_mig_order {
		typename cell_collection::point src;
		typename cell_collection::point dst;
		uint32_t group_size;
		unsigned int seq; /* random sequence number used for sorting */
		double dist; /* distance -- used by raider attacks optionally */
		
		/* helpers to sort a vector of orders */
		static void sort_by_seq(std::vector<group_mig_order>& orders) {
			std::sort(orders.begin(), orders.end(),
				[](const group_mig_order& o1, const group_mig_order& o2) {
					return o1.seq < o2.seq;
				});
		}
		static void sort_by_ptid(std::vector<group_mig_order>& orders) {
			std::sort(orders.begin(), orders.end(),
				[](const group_mig_order& o1, const group_mig_order& o2) {
					return o1.dst < o2.dst || (o1.dst == o2.dst && o1.src < o2.src);
				});
		}
	};
	
	
	/* choose a destination for a migrating group using the internal helper
	 * that stores the target distribution;
	 * if N3 > 0, this checks that a group of this size fits the target */
	template<class RNG>
	bool choose_dst_helper(const typename cell_collection::point& pt1, group_mig_order& res,
			unsigned int max_retries, RNG& rng1, unsigned int N3 = 0) const {
		const auto& helper = cph.at(pt1);
		double max1 = helper.btpw->total_sum();
		if(max1 <= 0.0) return false; /* this can happen if all cells are occupied */
		size_t nb_size = use_global_neighbors ? g.size() : helper.nb_size;
		
		std::uniform_real_distribution<double> dd(0.0, max1);
		std::uniform_real_distribution<double> pf1;
		double occupied_prob = (!group_mig_prob_pempty && pempty > 0.0) ? 1.0 / pempty : 0.0;
		for(unsigned int j = 0; j < max_retries; j++) {
			size_t i = helper.btpw->upper_bound(dd(rng1));
			if(i >= nb_size) throw std::runtime_error("nbase::choose_dst_helper(): error selecting a target!\n");
			const auto& pt2 = use_global_neighbors ? g.ptid(i) : helper.neighbors[i];
			const auto& c2 = g.at(pt2);
			if(c2.st == cell::state::EMPTY && c2.delay_counter) continue; /* empty cells that cannot be resettled yet */
			if(c2.K < 1.0) continue; /* cells with very low carrying capacity are skipped as well (they should be excluded in the input as well) */
			
			if(!group_mig_prob_pempty) {
				bool nonempty = false;
				if(peaceful_migrations && N3) {
					if(c2.N + N3 > c2.K) nonempty = true;
				}
				else if(c2.st != cell::state::EMPTY) nonempty = true;
				if(nonempty) {
					double x = (occupied_prob > 0.0) ? pf1(rng1) : 1.0;
					if(x >= occupied_prob) continue;
				}
			}
			/* skip this cell if it cannot fit the migrating population
			 * (for simulation variant without conflict) */
			else if(N3 && c2.N + N3 > c2.K) continue;
			res.src = pt1;
			res.dst = pt2;
			return true;
		}
		return false;
	}
	
	/* set the state of a cell externally (e.g. for initialization);
	 * this also takes care of creating the helper distributions here, but not updating;
	 * it can be overriden by a derived class if there is a need to also update the helper distributions
	 * returns if the state actually changed
	 * (note: this needs to be a virtual function to work from set_start_pop() below) */
	virtual bool set_cell_state(const typename cell_collection::point& pt, cell::state st, unsigned int pop) {
		cell& c = g.at(pt);
		auto old_state = c.st;
		c.st = st;
		c.N = pop;
		
		/* cell became occupied, might need to create the helper distribution for it */
		if(st != old_state && old_state == cell::state::EMPTY)
			if(use_prob_helper != cell_prob_helper_type::NONE && !external_helper)
				create_cell_prob_helper_one(pt);
		return (old_state != st);
	}
	
	virtual void reset() {
		for(const auto& pt : g) set_cell_state(pt, cell::state::EMPTY, 0);
	}
	
	/* add a starting population; if reset == true, ensure that all cells are set to empty before
	 * (needed only if multiple realizations are run with the same instance) */
	void set_start_pop(const typename cell_collection::point& start_cell, unsigned int starting_pop,
			bool reset_pop = false, double fill_factor = 1.0, double pchoice = 1.0) {
		if(reset_pop) reset();
	
		auto set_state = [this] (const typename cell_collection::point& pt, unsigned int pop) {
			set_cell_state(pt, cell::state::FARMING, pop);
		};
		
		std::bernoulli_distribution bd(pchoice);
		
		if(starting_pop) {
			for(auto it = g.neighbors_start(start_cell, 1e9, true); !it.is_end() && starting_pop; ++it) {
				auto& c1 = g.at(it->first);
				if(pchoice < 1.0 && !bd(rng)) continue;
				unsigned int tmp = std::min((unsigned int)std::floor(c1.K * fill_factor), starting_pop);
				if(tmp) {
					set_state(it->first, tmp);
					starting_pop -= tmp;
				}
			}
		}
		else set_state(start_cell, 100);
	}
	
	/* get a random number (useful for seeding other RNGs) */
	uint64_t get_random() { return rng(); }
	
	void create_cell_prob_helper(bool create_all = false) {
		cph.clear();
		for(const auto& pt : g) if(create_all || g.at(pt).st != cell::state::EMPTY)
			create_cell_prob_helper_one(pt);
	}
	
	
protected:
	/* base binary file header */
	constexpr static uint64_t bin_file_header_base = 0xae50f3f45e30a440UL;
	/* additional flags that can be used */
	constexpr static uint64_t bin_file_cdf = 1UL; /* CDF values are stored directly */
	constexpr static uint64_t bin_file_nonbr = 2UL; /* true if global neighbors are not duplicated in the bin file */
	constexpr static uint64_t bin_file_use_pempty = 4UL; /* true if migration probabilities are multiplied by pempty in the helper */
	constexpr static uint64_t bin_file_local_neighbors = 8UL; /* true if global neighbors are NOT used (each helper stores a list of neighbors separately) */
	constexpr static uint64_t bin_file_has_helpers = 16UL; /* true if helper distributions are stored in the binary files */
	constexpr static uint64_t header_flags = 31UL; /* combination of all flags */
	


	void write_one_helper(const typename cell_collection::point& x, size_t sums_size, bool use_cdf, FILE* f) const {
		const auto& helper = cph.at(x);
		const auto& y = helper.btpw;
		
		if(!use_global_neighbors) {
			sums_size = helper.nb_size;
			/* write the neighbors preceeded by the size of the array */
			const auto* data = helper.neighbors;
			const size_t pt_size = sizeof(typename cell_collection::point);
			if(fwrite(&sums_size, sizeof(uint64_t), 1, f) != 1)
				throw std::runtime_error("nbase::write_combined_matrix(): error writing to the output file!\n");
			for(size_t i = 0; i < sums_size; i++) if(fwrite(data + i, pt_size, 1, f) != 1)
				throw std::runtime_error("nbase::write_combined_matrix(): error writing to the output file!\n");
			/* potential padding */
			size_t tmpsize = pt_size * sums_size;
			if(tmpsize % 8) {
				uint8_t padding[8] = {0};
				size_t pad = 8 - (tmpsize % 8);
				if(fwrite(padding, 1, pad, f) != pad)
					throw std::runtime_error("nbase::write_combined_matrix(): error writing to the output file!\n");
			}
			if(!use_cdf) sums_size /= 2;
		}
		if(y->size() != sums_size) throw std::runtime_error("nbase::write_combined_matrix(): unexpected size of the stored distribution!\n");
		
		
		for(size_t i = 0; i < sums_size; i++) {
			double z = y->get_sum(i);
			if(fwrite(&z, sizeof(double), 1, f) != 1)
				throw std::runtime_error("nbase::write_combined_matrix(): error writing to the output file!\n");
		}
	}

public:

	/* write a binary file with basic parameters, cells and probabilities */
	void write_combined_matrix(const char* fn, bool store_global_neighbors = false, bool use_global_neighbors_override = false,
			bool no_helpers = false) {
		if((store_global_neighbors || use_global_neighbors_override) && !no_helpers) {
			if(cph.size() && !use_global_neighbors) throw std::runtime_error("nbase::write_combined_matrix(): unsupported combination!\n");
			use_global_neighbors = true;
		}
		
		bool use_cdf = false;
		switch(use_prob_helper) {
			case cell_prob_helper_type::CDF:
				use_cdf = true;
				break;
			case cell_prob_helper_type::BTVEC:
				throw std::runtime_error("nbase::write_combined_matrix(): helper type not supported!\n");
			case cell_prob_helper_type::BTCB:
			case cell_prob_helper_type::NONE:
				break; /* default type: callback-based helper */
		}
		
		size_t s1 = g.size();
		if(s1 > std::numeric_limits<uint32_t>::max())
			throw std::runtime_error("nbase::write_combined_matrix(): too many cells!\n");
		uint32_t nnodes = (uint32_t)s1;
		
		FILE* f = fopen(fn, "w");
		if(!f) throw std::runtime_error("nbase::write_combined_matrix(): cannot open output file!\n");
		
		{
			uint64_t tmpheader = bin_file_header_base;
			if(use_cdf) tmpheader |= bin_file_cdf;
			if(!store_global_neighbors) tmpheader |= bin_file_nonbr;
			if(group_mig_prob_pempty) tmpheader |= bin_file_use_pempty;
			if(!use_global_neighbors) tmpheader |= bin_file_local_neighbors;
			if(!no_helpers) tmpheader |= bin_file_has_helpers;
			if(fwrite(&tmpheader, sizeof(uint64_t), 1, f) != 1)
				throw std::runtime_error("nbase::write_combined_matrix(): error writing to the output file!\n");
		}
		
		/* write the cell parameters / properties */
		{
			size_t tmp = g.write(f);
			if(!tmp) throw std::runtime_error("nbase::write_combined_matrix(): error writing to the output file!\n");
			
			/* pad to a 8 byte boundary */
			if(tmp % 8) {
				uint8_t padding[8] = {0};
				size_t pad = 8 - (tmp % 8);
				if(fwrite(padding, 1, pad, f) != pad)
					throw std::runtime_error("nbase::write_combined_matrix(): error writing to the output file!\n");
			}
		}
		
		/* parameters: group_mig_dist, group_mig_pow, potentially pempty */
		{
			double tmp1[2] = {group_mig_dist, group_mig_pow};
			if(fwrite(tmp1, sizeof(double), 2, f) != 2)
				throw std::runtime_error("nbase::write_combined_matrix(): error writing to the output file!\n");
			if(group_mig_prob_pempty) if(fwrite(&pempty, sizeof(double), 1, f) != 1)
				throw std::runtime_error("nbase::write_combined_matrix(): error writing to the output file!\n");
		}
		
		if(store_global_neighbors) {
			/* next field: number of cells (uint32_t) */
			if(fwrite(&nnodes, sizeof(uint32_t), 1, f) != 1)
				throw std::runtime_error("nbase::write_combined_matrix(): error writing to the output file!\n");
			
			/* write all cell IDs in the order used by global_neighbors -- TODO: this should be eliminated! */
			for(size_t i = 0; i < s1; i++) {
				const auto& pt1 = g.ptid(i);
				if(fwrite(&pt1, sizeof(typename cell_collection::point), 1, f) != 1)
					throw std::runtime_error("nbase::write_combined_matrix(): error writing to the output file!\n");
			}
			
			/* pad to 8 byte boundary */
			size_t s3 = sizeof(uint32_t) + sizeof(typename cell_collection::point) * s1;
			if(s3 % 8) {
				uint8_t tmp[8] = {0};
				size_t s4 = 8 - (s3 % 8);
				if(fwrite(tmp, 1, s4, f) != s4)
					throw std::runtime_error("nbase::write_combined_matrix(): error writing to the output file!\n");
			}
		}
		
		/* write all carrying capacities */
		for(size_t i = 0; i < s1; i++) {
			const auto& id = g.ptid(i);
			double K1 = this->orig_K.size() ? this->orig_K.at(id) : g.at(id).K;
			if(fwrite(&K1, sizeof(double), 1, f) != 1)
				throw std::runtime_error("nbase::write_combined_matrix(): error writing to the output file!\n");
		}
		
		if(!no_helpers) {
			/* write all partial sums in order */
			uint64_t sums_size = use_cdf ? s1 : s1 / 2;
			size_t j = 0;
			if(create_dist_threads > 1 && !cph.size()) {
				// use multiple threads to create helper distributions
				// 1. ensure that all points are present in the cph hashmap
				for(const auto& x : g) cph.emplace(std::piecewise_construct, std::forward_as_tuple(x), std::forward_as_tuple());
				
				// 2. counters and mutexes
				std::atomic<size_t> i_started{0}; // counter for which points we have started to create the helper distribution
				std::atomic<size_t> i_written{0}; // counter for which points we have written the helper distribution
				std::mutex mf; // mutex protecting setting the finished helpers (given to create_cell_prob_helper_one())
				size_t i_total = g.size();
				
				// 3. function executed in parallel
				auto tf = [this, &i_started, &i_written, i_total, s1, sums_size, use_cdf, f, &mf]() {
					while(true) {
						size_t i = i_started++;
						if(i >= i_total) break; // no more work to do
						auto x = g.ptid(i);
						create_cell_prob_helper_one(x, true, &mf);
						if(i == i_written) {
							/* it is our turn to write -- write out all contiguous finished results */
							while(true) {
								write_one_helper(x, sums_size, use_cdf, f);
								fprintf(stderr, "%lu / %lu cells processed\r", i + 1, s1);
								// erase the helper -- note: we cannot remove from the map to avoid invalidating references
								cph.at(x).btpw.reset(nullptr);
								mf.lock();
								i++;
								bool cont = false;
								if(i < i_total) {
									x = g.ptid(i);
									if(cph.at(x).btpw) cont = true;
								}
								if(!cont) i_written = i;
								mf.unlock();
								if(!cont) break;
							}
						}
					}
				};
				
				// 4. run the given number of threads
				std::vector<std::thread> threads;
				threads.reserve(create_dist_threads);
				for(unsigned int j = 0; j < create_dist_threads; j++) threads.emplace_back(tf);
				for(unsigned int j = 0; j < create_dist_threads; j++) threads[j].join();
				if(i_written != i_total) throw std::runtime_error("nbase::write_combined_matrix(): inconsistent state after writing helpers!\n");
				
				// 5. erase the helper distribution hashmap
				cph.clear();
			}
			else for(const auto& x : g) {
				bool erase_x = false;
				if(!cph.count(x)) {
					create_cell_prob_helper_one(x);
					erase_x = true;
				}
				
				write_one_helper(x, sums_size, use_cdf, f);
				
				if(erase_x) {
					cph.erase(x);
					fprintf(stderr, "%lu / %lu cells processed\r", j, s1);
				}
				j++;
			}
		}
		
		fclose(f);
	}
	
	/* read previously saved binary helper matrix, storing the carrying capacities in the cellgrid as well
	 * if no_helpers == true, do not read probability distributions even if they are present in the input file */
	void read_combined_matrix(const char* fn, bool no_helpers = false) {
		{
			bool tmp = mp->open_file(fn, FileMapping::Mode::ReadOnly);
			if(tmp) {
				/* depending on whether the distributions might change, we either want
				 * a private writeable or a shared readonly mapping */
				if(group_mig_prob_pempty && !no_helpers) tmp = mp->map_file(FileMapping::Mode::ReadWrite, false);
				else tmp = mp->map_file(FileMapping::Mode::ReadOnly, true);
			}
			if(!tmp) throw std::runtime_error("nbase::read_combined_matrix(): error opening the input file!\n");
		}
		
		cph.clear();
		
		uint8_t* data = (uint8_t*)mp->get_mapping();
		if(!data) throw std::runtime_error("nbase::read_combined_matrix(): error opening the input file!\n");
		
		bool store_global_neighbors = false;
		bool have_helpers = false;
		bool use_cdf = false; /* if helper distribution uses CDFs */
		bool use_pempty = false; /* if pempty is stored in the binary file */
		{
			uint64_t* tmp1 = (uint64_t*)data;
			uint64_t tmp2 = (*tmp1 & ~header_flags);
			if(tmp2 != bin_file_header_base)
				throw std::runtime_error("nbase::read_combined_matrix(): header does not match expected value!\n");
			have_helpers = (*tmp1 & bin_file_has_helpers) && !no_helpers;
			store_global_neighbors = !(*tmp1 & bin_file_nonbr);
			use_pempty = (*tmp1 & bin_file_use_pempty);
			
			if(have_helpers) {
				use_cdf = (*tmp1 & bin_file_cdf);
				use_global_neighbors = !(*tmp1 & bin_file_local_neighbors);
				if(store_global_neighbors && !use_global_neighbors)
					throw std::runtime_error("nbase::read_combined_matrix(): invalid parameter combination in header!\n");
				
				use_prob_helper = use_cdf ? cell_prob_helper_type::CDF : cell_prob_helper_type::BTCB;
			
				if(use_pempty != group_mig_prob_pempty)
					throw std::runtime_error("nbase::read_combined_matrix(): unsupported helper type!\n");
			}
			// slight optimization: if we are not using the helpers, we can always set the access mode to read-only
			else if(mp->mem_mode() != FileMapping::Mode::ReadOnly) mp->set_mem_mode(FileMapping::Mode::ReadOnly);
		}
		data += sizeof(uint64_t);
		
		if(have_helpers && have_pars)
			throw std::runtime_error("nbase::read_combined_matrix(): parameters (-P -G) cannot be set as command line arguments!\n");
		
		{
			/* read the cells */
			size_t size1 = g.read(data, mp->file_size() - sizeof(uint64_t));
			if(!size1) throw std::runtime_error("nbase::read_combined_matrix(): error reading the cell properties!\n");
			
			/* potential padding */
			if(size1 % 8) size1 += (8 - (size1 % 8));
			data += size1;
		}
		
		/* parameters: group_mig_dist, group_mig_pow
		 * only if we are reading the helpers or if these were not provided as command line arguments */
		if(have_helpers || !have_pars) {
			const double* tmp1 = (const double*)data;
			group_mig_dist = tmp1[0];
			group_mig_pow = tmp1[1];
			if(use_pempty) pempty = tmp1[2];
		}
		data += 2*sizeof(double); // need to advance the file position even if the parameters are not read / used
		if(use_pempty) data += sizeof(double);
		
		uint32_t s1 = (uint32_t)g.size();
		if(store_global_neighbors) {
			/* size of nodes */
			s1 = *(uint32_t*)data;
			data += sizeof(uint32_t);
			if(s1 != g.size()) throw std::runtime_error("nbase::read_combined_matrix(): number of nodes do not match!\n");
			
			/* check that the global_neighbors array in the matrix matches the nodes in g */
			typename cell_collection::point* tmp1 = (typename cell_collection::point*)data;
			for(uint32_t i = 0; i < s1; i++)
				if(tmp1[i] != g.ptid(i)) throw std::runtime_error("nbase::read_combined_matrix(): global_neighbors array does not match!\n");
			
			data += sizeof(typename cell_collection::point) * (uint64_t)s1;
			size_t s3 = sizeof(uint64_t) + sizeof(uint32_t) + sizeof(typename cell_collection::point) * (uint64_t)s1;
			if(s3 % 8) {
				size_t tmp1 = 8 - (s3 % 8);
				data += tmp1;
			}
		}
		
		{
			/* carrying capacities */
			const double* K1 = (const double*)data;
			for(uint32_t i = 0; i < s1; i++) g.at(g.ptid(i)).K = K1[i];
			data += sizeof(double) * (uint64_t)s1;
		}
		
		if(have_helpers) {
			size_t sums_size = use_cdf ? s1 : s1 / 2; /* size of sums array (if using global neighbors) */
			for(uint32_t i = 0; i < s1; i++) {
				const auto& pt = g.ptid(i);
				auto& x1 = cph[pt]; /* this will create the helper struct for this point */
				
				size_t nb_size = s1;
				if(!use_global_neighbors) {
					uint64_t* tmp1 = (uint64_t*)data;
					nb_size = *tmp1;
					sums_size = use_cdf ? nb_size : nb_size / 2;
					x1.nb_size = nb_size;
					data += sizeof(uint64_t);
					x1.neighbors = (const typename cell_collection::point*)data;
					size_t tmpsize = nb_size * sizeof(typename cell_collection::point);
					if(tmpsize % 8) tmpsize += 8 - (tmpsize % 8);
					data += tmpsize;
				}
				
				if(use_cdf) {
					const double* sums1 = (const double*)data;
					x1.btpw.reset(new typename nbase<cell_collection>::cdf_wrapper(sums1, sums_size));
					data += sums_size * sizeof(double);
				}
				else {
					/* note: if adjusting of probabilities is done inside group_mig_prob, then
					 * the helper needs to be writable (adjustable), we can set a flag for that here;
					 * this will make this class non-copyable (attempt to copy will throw an exception) */
					auto btpw = new nbase<cell_collection>::btprob_wrapper<btprob_cb_base>(group_mig_prob_pempty);
					double* sums1 = (double*)data;
					btpw->btp.set_sums_external(sums1, nb_size);
					btpw->btp.prob =  [this, pt, nb = x1.neighbors] (size_t j) {
						const auto& pt2 = use_global_neighbors ? g.ptid(j) : nb[j];
						if(pt2 == pt) return 0.0; /* exclude self from targets */
						double dst = g.point_dist(pt, pt2);
						return group_mig_prob(pt2, dst);
					};
					x1.btpw.reset(btpw);
					data += sums_size * sizeof(double);
				}
			}
			external_helper = true;
		}
	}
	
	double check_helper() const {
		double s1 = 0.0;
		for(const auto& pt : g) {
			const auto& x = cph.at(pt);
			s1 += x.btpw->check_sum();
		}
		return s1;
	}
	
	/* helper to write out a point as a set of unsigned integers */
	static void write_pt(const typename cell_collection::point& pt, FILE* fout, char delim = '\t') {
		auto t = cell_collection::get_point_as_tuple(pt);
		constexpr size_t s = std::tuple_size<decltype(t)>::value;
		// note: we assume that point IDs are made up of unsigned integers (this is true for cellgrid and cellnet)
		// also, we assume that the size of tuples is either 1 or 2
		static_assert(s == 1 || s == 2);
		unsigned int x = std::get<0>(t);
		fprintf(fout, "%u", x);
		if constexpr (s > 1) {
			x = std::get<1>(t);
			fprintf(fout, "%c%u", delim, x);
		}
	}

	/* diagnostic: write all neighbor relations in the helper maps to the given file */
	void write_all_neighbors(FILE* f) const {
		for(const auto& x : cph) {
			const auto& pt = x.first;
			const auto& helper = x.second;
			for(size_t i = 0; i < helper.nb_size; i++) {
				const auto& pt2 = helper.neighbors[i];
				write_pt(pt, f);
				fputc('\t', f);
				write_pt(pt2, f);
				fputc('\n', f);
			}
		}
	}

	/* get the total number of farming cells and total population numbers */
	struct totals {
		unsigned int farming_cells = 0;
		unsigned int raider_cells = 0;
		uint64_t farmers = 0;
		uint64_t raiders = 0;
	};
	totals get_totals() const {
		totals res;
		for(size_t i = 0; i < g.size(); i++) {
			const auto& c = g.at_ix(i);
			switch(c.st) {
				case cell::state::FARMING:
					res.farming_cells++;
					res.farmers += (uint64_t)c.N;
					break;
				case cell::state::RAIDERS:
					res.raider_cells++;
					res.raiders += (uint64_t)c.N;
					break;
				case cell::state::EMPTY:
					/* empty cells have no population, they are not counted */
					break;
			}
		}
		return res;
	}
	
protected:
	/* get the ID of the currently running OpenMP thread or 0 if compiled without OpenMP */
	static int get_thread_num() {
#ifdef _OPENMP
		return omp_get_thread_num();
#else
		return 0;
#endif
	}
	
	/* set the number of OpenMP threads in use -- noop if compiled without OpenMP */
	void set_num_omp_threads(int num_threads) {
#ifdef _OPENMP
		omp_threads = std::min(num_threads, omp_get_num_procs());
		if(omp_threads <= 0) omp_threads = 1;
		omp_set_dynamic(0);
		omp_set_num_threads(omp_threads);
#else
		omp_threads = 1;
#endif
	}
	
	/* ensure that the number of OMP threads is set to a sensible value
	 * and that each thread has a separate RNG */
	void ensure_omp_threads() {
		if(omp_threads <= 0) set_num_omp_threads(1);
		
		/* ensure that each thread has a separate RNG */
		if(par_rng.size() != (size_t)omp_threads) {
			par_rng.clear();
			for(int i = 0; i < omp_threads; i++)
				par_rng.emplace_back(rng());
		}
	}
	
	/* get a reference to a random number generator that is safe to use,
	 * depending on whether OMP is enabled and running on multiple threads */
	std::mt19937_64& prng() {
#ifdef _OPENMP
		if(omp_threads > 1) return par_rng.at(get_thread_num());
		else return rng;
#else
		return rng;
#endif
	}
	
public:
	/* re-seed the random number generator */
	void seed_rng(uint64_t s) {
		rng.seed(s);
		par_rng.clear();
		ensure_omp_threads();
	}
	
	virtual ~nbase() { }
};

/* helper derived class that is necessary to create new callback functions
 * for the helper objects */
template<class cell_collection>
class nbase2 : public nbase<cell_collection> {
protected:
	explicit nbase2(bool pe) : nbase<cell_collection>(pe) { }
	using cell_prob_helper_type = typename nbase<cell_collection>::cell_prob_helper_type;
public:
	nbase2(const nbase2& nb) : nbase<cell_collection>(nb) {
		if(this->use_prob_helper == cell_prob_helper_type::BTCB) {
			for(auto& x : this->cph) {
				auto& helper = x.second;
				auto& btpw = static_cast<typename nbase<cell_collection>::template btprob_wrapper<btprob_cb_base>&>(*helper.btpw);
				auto& btp = btpw.btp;
				btp.prob = this->get_cb_for_pt(x.first);
			}
		}
	}
};



/******************************************************************
 * Helper functionality for reading and processing weather data
 ******************************************************************/


/* shared information and functionality for updating weather data */
struct weather_update_base {
protected:
	std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, double> > > wmap;
	std::unordered_set<unsigned int> wcells; /* cells with weather data */
	
	template<class cell_collection>
	void prepare(nbase<cell_collection>& nn) const {
		if(weather_cut) {
			/* initialize carrying capacity to zero everywhere, so that
			 * only cells with valid climate data will be used in this step */
			for(const auto& pt : nn.g) nn.g.at(pt).K = 0.0;
		}
		else if(weather_factors) {
			/* only zero out cells which has associated weather data */
			for(const auto& pt : wcells) nn.g.at(pt).K = 0.0;
		}
	}
	
	template<class cell_collection>
	void update_one(nbase<cell_collection>& nn, unsigned int wid, double factor) const {
		{
			/* further scale the variations */
			double tmp = factor - 1.0;
			factor = std::max(1.0 + sd_factor * tmp, 0.0);
		}
		
		auto it = wmap.find(wid);
		if(it == wmap.end()) return;
		const auto& cells = it->second;
		for(const auto& x : cells) {
			if(weather_factors) nn.g.at(x.first).K += nn.orig_K.at(x.first) * factor * x.second;
			else nn.g.at(x.first).K = nn.orig_K.at(x.first) * factor;
		}
	}
	
public:
	bool weather_allow_missing = false; /* allow missing data in some years (last year's value is used) */
	bool weather_cut = false; /* cut study area to regions which have weather data (not compatible with the previous option) */
	bool weather_factors = false; /* true if match to weather cells is not one-to-one (a factor is read together with the matching) */
	int first_year = INT_MIN; /* "start year" of the simulation -- if given, ignore years before this */
	double sd_factor = 1.0; /* scale the standard deviation of weather variation by this factor -- it is assumed that factors have a mean of one */
	
	template<class cell_collection>
	void read_mapping(const char* fn, bool is_csv, nbase<cell_collection>& nn) {
		wmap.clear();
		wcells.clear();
		
		read_table2 rtwm(fn);
		if(is_csv) {
			rtwm.set_delim(',');
			rtwm.read_line();
		}
		while(rtwm.read_line()) {
			unsigned int cellid, wid;
			double factor = 1.0;
			if(!rtwm.read(cellid, wid)) break;
			if(weather_factors && !rtwm.read(factor)) break;
			/* note: only save mappings for cells that exist in the simulation space */
			auto& cells = wmap[wid]; /* note: ensure that the climate cell exists even if it has no associated cells in the simulation */
			if(nn.g.cell_exists(cellid)) {
				cells.push_back(std::make_pair(cellid, factor));
				wcells.insert(cellid);
			}
		}
		if(rtwm.get_last_error() != T_EOF) {
			rt_error(rtwm);
			throw std::runtime_error("Error reading the mapping between cells and weather data!\n");
		}
		
		/* save the original carrying capacity values to be used later
		 * (only if this was not created previously) */
		if(!nn.orig_K.size()) for(const auto& pt : nn.g) nn.orig_K[pt] = nn.g.at(pt).K;
	}
};

/* helper struct to update yields during the simulation by reading
 * weather data from a file */
struct weather_update_file : public weather_update_base {
	std::unique_ptr<read_table2> rtwp;
	int wf_last_year;
	
	void throw_error() {
		std::string err = rtwp->exception_string("Error reading weather data:\n");
		throw std::runtime_error(err);
	}
	
	void open_file(const char* weather_file) {
		rtwp.reset(new read_table2(weather_file));
		if(!(rtwp->read_line())) throw_error();
	}
	
	template<class cell_collection>
	void operator ()(nbase<cell_collection>& nn, unsigned int i) {
		if(!rtwp) return; /* allow running as no-op if no file was opened */
		auto& rtw = *rtwp; /* get a reference for convenience */
		if(rtw.get_last_error() != T_OK) throw_error();
		
		prepare(nn);
		
		bool first_line = true;
		size_t entries_read = 0;
		/* note: the "current" line is already read in rtw, except for the case when the file already ended */
		do {
			/* parse one line */
			int year;
			unsigned int wid;
			double factor;
			if(!rtw.read(year, wid, factor)) break;
			
			if(first_line) {
				if(!i && year < first_year) continue; /* skip ahead until the first year of the simulation is reached */
				if(i && wf_last_year + 1 != year) {
					fprintf(stderr, "Non-consecutive years in weather factor data (%d, %d; in line %" PRIu64 ")!\n",
						wf_last_year, year, rtw.get_line());
					throw std::runtime_error("Error processing weather data!\n");
				}
				wf_last_year = year;
				first_line = false;
			}
			else if(year != wf_last_year) {
				rtw.reset_pos();
				break;
			}
			
			/* update the carrying capacities based on the saved values */
			update_one(nn, wid, factor);
			
			entries_read++;
		} while(rtw.read_line());
		
		if( ! (rtw.get_last_error() == T_OK || rtw.get_last_error() == T_EOF) ) throw_error();
		
		if(!weather_allow_missing && entries_read != wmap.size()) {
			fprintf(stderr, "Missing entries in the climate variability data (line %" PRIu64 ")!\n", rtw.get_line());
			throw std::runtime_error("Error processing weather data!\n");
		}
	}
};

/* helper for the case when we don't have weather data */
struct weather_update_noop {
	template<class cell_collection>
	void operator ()(nbase<cell_collection>&, unsigned int) const { }
};

struct weather_update_all : public weather_update_base {
	int start_year = 0;
	std::unordered_map<int, std::pair<std::vector<unsigned int>, std::vector<double> > > data;
	bool have_data = false; /* simplification to allow using this case even in the case when there is no weather data */
	
	template<class cell_collection>
	void operator ()(nbase<cell_collection>& nn, unsigned int i) const {
		if(!have_data) return;
		int year = start_year + i;
		auto it = data.find(year);
		if(it == data.end()) throw std::runtime_error("weather_update_all: missing year in weather data!\n");
		const auto& tmp = it->second;
		const auto& wids = tmp.first;
		const auto& vals = tmp.second;
		
		prepare(nn);
		
		size_t entries_read = 0;
		for(size_t j = 0; j < wids.size(); j++) {
			unsigned int wid = wids[j];
			double factor = vals.at(j);
			
			/* update the carrying capacities based on the saved values */
			update_one(nn, wid, factor);
			
			entries_read++;
		}
		
		if(!weather_allow_missing && entries_read != wmap.size()) {
			fprintf(stderr, "Missing entries in the climate variability data (year %u)!\n", i);
			throw std::runtime_error("Error processing weather data!\n");
		}
	}
	
	void read_data(const char* fn, size_t size_hint = 0) {
		data.clear();
		read_table2 rtw(fn);
		
		while(rtw.read_line()) {
			int year;
			unsigned int wid;
			double factor;
			if(!rtw.read(year, wid, factor)) break;
			if(year < first_year) continue; // skip data before the expected first year of the simulation
			
			auto x = data.emplace(std::piecewise_construct, std::forward_as_tuple(year), std::forward_as_tuple());
			auto& tmp = x.first->second;
			if(size_hint && x.second) {
				tmp.first.reserve(size_hint);
				tmp.second.reserve(size_hint);
			}
			tmp.first.push_back(wid);
			tmp.second.push_back(factor);
			
			if(year < start_year) start_year = year;
		}
		
		if(rtw.get_last_error() != T_EOF) {
			std::string err = rtw.exception_string("Error reading weather data:\n");
			throw std::runtime_error(err);
		}
		
		have_data = true;
	}
};


/********************************************************************
 * Common options that are shared by different programs running
 * the simulation and facility to process them.
 ********************************************************************/
template<class cell_collection>
struct common_options {
	unsigned int steps = 2000; /* timesteps to run the simulation for */
	
	char* helper_matrix_fn = nullptr; /* previously created helper matrix to load */
	bool helper_only_net = false; /* if true, the above file is read, but the actual helper distributions are not used from it 
		(only the simulation space and the distances if included) */
	
	double K1 = 200.0;
	bool adjust_K = false; /* if have_K == true, adjust carrying capacities so that the average is the value given by K1 */
	bool adjust_K_only_nonzero = true; /* if adjust_K == true, only cells with nonzero carrying capacity are considered */
	double Kthres = 0.0; /* if > 0.0, leave out cells with capacitieis below this (after scaling) */
	
	char* out_base = nullptr; /* filename for (detailed) output */
	unsigned int out_period = 10; /* write detailed output this often (if running the simulation) */
	bool out_base_zip = false; /* set to true if the above file should be compressed with gzip */
	
	/* cell where to seed the simulation */
	typename cell_collection::point start_cell{};
	bool have_start_cell = false;
	unsigned int starting_pop = 0; /* if given, try to distribute this amount of population in the beginning */
	double start_pop_fill_factor = 1.0; /* starting population should fill ratio to carrying capacity */
	double start_pop_pchoice = 1.0; /* probability of each cell included in the start population */
	
	/* use a mapping from cell IDs to weather cell IDs */
	char* weather_mapping = nullptr;
	/* climate-based variability data (multiplicative factor for carrying capacities) */
	char* weather_file = nullptr;
	bool weather_allow_missing = false; /* allow missing data in some years (last year's value is used) */
	bool weather_cut = false; /* cut study area to regions which have weather data (not compatible with the previous option) */
	bool weather_factors = false; /* true if match to weather cells is not one-to-one (a factor is read together with the matching) */
	bool weather_mapping_csv = false; /* true if weather mapping file is in CSV format */
	int weather_first_year = INT_MIN; /* first year to consider in the weather file (i.e. "start date" of the simulation) */
	double weather_sd_factor = 1.0; /* further scale weather variations by this factor */
	
	bool track_avg_distance = false; /* track average distance of newly established cells */
	
	
	/* Try to parse a command line option, return the number of items
	 * consumed, or 0 if unknown, -1 if invalid. The flags parameter
	 * can be a bitwise OR of the following values: */
	struct option_flags {
		constexpr static uint32_t SIMULATION = 1; /* set if running a simulation; enables options related to that */
		constexpr static uint32_t WEATHER = 2; /* set if using the weather-related flags make sense */
		constexpr static uint32_t HELPER = 4; /* set if reading a previously generated helper matrix makes sense */
	};
	int parse_option(int argc, char** argv, int i, uint32_t flags = 0) {
		if(argv[i][0] != '-') return -1;
		int res = 0;
		switch(argv[i][1]) {
			case 'w':
				if(!(flags & option_flags::WEATHER)) return -1;
				res = 1;
				switch(argv[i][2]) {
					case 'm':
						weather_allow_missing = true;
						break;
					case 'c':
						weather_cut = true;
						break;
					case 'f':
						weather_factors = true;
						break;
					case 'C':
						weather_mapping_csv = true;
						break;
					case 0:						
						/* need two files */
						if(i + 2 >= argc) {
							fprintf(stderr, "-w needs two arguments!\n");
							return -1;
						}
						weather_mapping = argv[i+1];
						weather_file = argv[i+2];
						res += 2;
						break;
					case 'F':
						weather_first_year = atoi(argv[i+1]);
						res = 2;
						break;
					case 's':
						weather_sd_factor = atof(argv[i+1]);
						res = 2;
						break;
					default:
						return -1;
				}
				break;
			case 'i':
				if(!(flags & option_flags::SIMULATION)) return -1;
				res = cell_collection::try_parse_point(argc, argv, i+1, start_cell);
				if(res < 0) return res;
				have_start_cell = true;
				res++;
				break;
			case 'I':
				if(!(flags & option_flags::SIMULATION)) return -1;
				switch(argv[i][2]) {
					case 'f':
						start_pop_fill_factor = atof(argv[i+1]);
						break;
					case 'p':
						start_pop_pchoice = atof(argv[i+1]);
						break;
					case 0:
						starting_pop = atoi(argv[i+1]);
						break;
					default:
						return -1;
				}
				res = 2;
				break;
			case 'o':
				switch(argv[i][2]) {
					case 'p':
						if(!(flags & option_flags::SIMULATION)) return -1;
						out_period = atoi(argv[i+1]);
						if(!out_period) return -1;
						res = 2;
						break;
					case 'z':
						if(!(flags & option_flags::SIMULATION)) return -1;
						out_base_zip = true;
						res = 1;
						break;
					case 0:
						out_base = argv[i+1];
						res = 2;
						break;
				}
				break;
			case 'S':
				if(!(flags & option_flags::SIMULATION)) return -1;
				steps = atoi(argv[i+1]);
				res = 2;
				break;
			case 'K':
				switch(argv[i][2]) {
					case 0:
						K1 = atof(argv[i+1]);
						adjust_K = true;
						res = 2;
						break;
					case '0':
						adjust_K_only_nonzero = false;
						res = 1;
						break;
					case 't':
						Kthres = atof(argv[i+1]);
						res = 2;
						break;
					/* note: there could be other options handled by the caller */
				}
				break;
			case 'T':
				if(!(flags & option_flags::SIMULATION)) return -1;
				track_avg_distance = true;
				res = 1;
				break;
			case 'H':
				if(argv[i][2] == 'f') {
					helper_matrix_fn = argv[i+1];
					if(argv[i][3] == '0') helper_only_net = true;
					res = 2;
				}
				break;
		}
		return res;
	}
	
	
	/* some helpers */
	void read_weather_mapping(weather_update_base& wub, nbase<cell_collection>& nn) {
		wub.weather_allow_missing = weather_allow_missing;
		wub.weather_cut = weather_cut;
		wub.weather_factors = weather_factors;
		wub.first_year = weather_first_year;
		wub.sd_factor = weather_sd_factor;
		if(weather_mapping) wub.read_mapping(weather_mapping, weather_mapping_csv, nn);
	}
	
	void ensure_start_cell(nbase<cell_collection>& nn) {
		if(!have_start_cell) {
			auto it = nn.g.begin();
			start_cell = *it;
			have_start_cell = true;
		}
	}
	
	void do_adjust_K(nbase<cell_collection>& nn, bool save_orig_K = false) {
		if(Kthres > 0.0 && helper_matrix_fn)
			throw std::runtime_error("common_options::do_adjust_K(): cannot use threshold if probabilities were loaded from a binary file!\n");
		if(adjust_K) for(unsigned int i = 0; ; i++) {
			double s1 = 0.0;
			double n1 = 0.0;
			unsigned int nthres = 0;
			for(const auto& pt : nn.g) {
				double tmpK = nn.g.at(pt).K;
				if(save_orig_K && i == 0) nn.orig_K[pt] = tmpK;
				if(tmpK > 0.0 || !adjust_K_only_nonzero) {
					if(i && tmpK < Kthres) {
						nn.g.at(pt).K = 0.0;
						nthres++;
					}
					else {
						s1 += tmpK;
						n1 += 1.0;
					}
				}
			}
			s1 /= n1;
			
			double factor = K1 / s1;
			for(const auto& pt : nn.g) nn.g.at(pt).K *= factor;
			if(Kthres <= 0.0 || (i && !nthres)) break;
		}
	}
};




/********************************************************************
 * Common facilities for running the simulation and outputting or
 * aggregating the results
 ********************************************************************/

/* the "results" of one step of the simulation */
template<class nclass>
struct step_res {
	unsigned int farmer_targets[3] = {0, 0, 0}; /* empty, farmer, raider */
	unsigned int raider_targets[3] = {0, 0, 0};
	
	double sum_dist = 0.0;
	unsigned int cnt_dist = 0;
	double sum_dist_new = 0.0;
	unsigned int cnt_dist_new = 0;
	unsigned int total_mig = 0;
	
	unsigned int group_mig_rejected = 0;
	
	typename nclass::totals t;
};


/* run one realization of the simulation (of the given steps), with class nn
 * (that should be a derived class of nbase, i.e. neolithic or neolithic_w)
 * 
 * wfun shuold be one of the weather update helper classes and ofun should
 * be a callable that takes the result of one simulation step and writes
 * or aggregates the relevant statistics (see below)
 * track_new_distance: track the average distance of cells settled for the first time separately */
template<class nclass, class wfun, class ofun>
void do_one_run(nclass& nn, unsigned int steps, typename nclass::point_type start_cell, wfun&& weather_update, ofun&& output,
		bool track_new_distance = false) {
	
	/* TODO: hash type for the case of points with no std::hash implementation! */
	std::unordered_set<typename nclass::point_type> settled;
	
	for(unsigned int i = 0; i < steps; i++) {
		weather_update(nn, i);
		
		/* step the demographic simulation */
		nn.cells_step(1.0);
		
		/* evaluate migrations and raider attacks */
		nn.create_group_orders();
		
		/* count the types of targets for all orders */
		step_res<nclass> sr;
		for(const auto& o : nn.orders) {
			const auto& cs = nn.g.at(o.src);
			const auto& cd = nn.g.at(o.dst);
			
			unsigned int t1 = 0;
			switch(cd.st) {
				case cell::state::EMPTY:
					t1 = 0;
					break;
				case cell::state::FARMING:
					t1 = 1;
					break;
				case cell::state::RAIDERS:
					t1 = 2;
					break;
			}
			if(cs.st == cell::state::FARMING) sr.farmer_targets[t1]++;
			else if(cs.st == cell::state::RAIDERS) sr.raider_targets[t1]++;
		}
		sr.group_mig_rejected = nn.orders_rejected;
		
		nn.apply_group_mig([&sr, &nn, start_cell, &settled, track_new_distance]
		(const auto&, const auto& pt2, unsigned int N) {
			const auto& c2 = nn.g.at(pt2);
			if(c2.st != cell::state::EMPTY) return;
			double dist1 = nn.g.coords_dist(start_cell, pt2);
			sr.sum_dist += dist1;
			sr.cnt_dist++;
			if(track_new_distance && !settled.count(pt2)) {
				settled.insert(pt2);
				sr.sum_dist_new += dist1;
				sr.cnt_dist_new++;
			}
			sr.total_mig += N;
		});
		
		sr.t = nn.get_totals();
		
		output(nn, i, sr);
	}
}

/* simulation output aggregating the population of cells */
template<class nclass>
struct output_avg {
	struct record {
		double nfarmers = 0.0; /* number of farmer cells */
		double pfarmers = 0.0; /* farmer population */
		double nraiders = 0.0; /* number of raider cells */
		double praiders = 0.0; /* raider population */
		double avg_dist = 0.0; /* average distance of newly settled cells from the start */
	};
	
	std::vector<record> sums; /* sum of records per run */
	std::vector<record> sum2; /* sum of squares */
	
	struct cell_record {
		double pop = 0.0; /* population */
		unsigned int n = 0; /* number of times this cell was occupied */
		bool nonzero() const {
			return (pop > 0.0) || (n > 0);
		}
	};
	
	std::vector<std::vector<cell_record> > cells_farmers; /* sum of farmer population in each cell */
	std::vector<std::vector<cell_record> > cells_raiders; /* sum of raider population in each cell */
	unsigned int out_period = 1; /* output period for the above */
	
	std::mutex* m = nullptr; /* mutex to protect access to the above arrays -- separate mutex for each time step */
	
	/* mapping from cell IDs to regions */
	std::unordered_map<unsigned int, unsigned int> cells_regions;
	
	void operator () (const nclass& nn, unsigned int i, const step_res<nclass>& sr) {
		if(i >= sums.size() || i >= sum2.size()) throw std::runtime_error("output_avg: no space allocated for output!\n");
		unsigned int ip = i / out_period;
		if((cells_farmers.size() || cells_raiders.size()) && (ip >= cells_farmers.size() || ip >= cells_farmers.size()))
			throw std::runtime_error("output_avg: no space allocated for output!\n");
		
		std::unique_lock<std::mutex> lock;
		if(m) lock = std::unique_lock<std::mutex>(m[ip]);
		
		record& r1 = sums[i];
		record& r2 = sum2[i];
		
		r1.nfarmers += sr.t.farming_cells;
		r1.pfarmers += sr.t.farmers;
		r1.nraiders += sr.t.raider_cells;
		r1.praiders += sr.t.raiders;
		
		r2.nfarmers += sr.t.farming_cells * sr.t.farming_cells;
		r2.pfarmers += sr.t.farmers * sr.t.farmers;
		r2.nraiders += sr.t.raider_cells * sr.t.raider_cells;
		r2.praiders += sr.t.raiders * sr.t.raiders;
		
		if(cells_farmers.size() && (i % out_period == 0)) {
			auto& c1 = cells_farmers[ip];
			auto& c2 = cells_raiders[ip];
			for(size_t j = 0; j < nn.g.size(); j++) {
				const auto& pt = nn.g.ptid(j);
				const cell& c = nn.g.at(pt);
				if(c.st == cell::state::EMPTY) continue;
				size_t j2 = j;
				if(cells_regions.size()) {
					auto it = cells_regions.find(pt);
					if(it == cells_regions.end()) continue;
					j2 = it->second;
				}
				if(c.st == cell::state::FARMING) {
					c1[j2].n++;
					c1[j2].pop += c.N;
				}
				else { /* raiders */
					c2[j2].n++;
					c2[j2].pop += c.N;
				}
			}
		}
	}
	
	std::vector<unsigned int> read_region_ids(read_table2& rtr) {
		std::unordered_map<unsigned int, unsigned int> region_ids;
		std::vector<unsigned int> region_ids_r;
		unsigned int nregions = 0;
		
		while(rtr.read_line()) {
			unsigned int cell_id, region_id;
			if(!rtr.read(cell_id, region_id)) break;
			
			auto r = region_ids.insert(std::make_pair(region_id, nregions));
			unsigned int rix = r.first->second;
			if(r.second) nregions++;
			
			cells_regions[cell_id] = rix;
		}
		
		if(rtr.get_last_error() != T_EOF)
			throw std::runtime_error(rtr.exception_string("Error reading cell -- region matching:\n"));
		
		region_ids_r.resize(nregions);
		for(const auto& x : region_ids) region_ids_r[x.second] = x.first;
		return region_ids_r;
	}
	std::vector<unsigned int> read_region_ids(read_table2&& rtr) { return read_region_ids(rtr); }
	
	void allocate(unsigned int steps, unsigned int ncells, bool have_raiders) {
		sums.clear();
		sum2.clear();
		cells_farmers.clear();
		cells_raiders.clear();
		
		sums.resize(steps);
		sum2.resize(steps);
		if(ncells) {
			unsigned int ns2 = steps / out_period + 1;
			cells_farmers.resize(ns2);
			cells_raiders.resize(ns2);
		
			for(unsigned int i = 0; i < ns2; i++) {
				cells_farmers[i].resize(ncells);
				if(have_raiders) cells_raiders[i].resize(ncells);
			}
		}
	}
	
	
	template<class wfun>
	void runmultiple(nclass& nn, const common_options<typename nclass::cell_collection_type>& opts, wfun&& wd,
			unsigned int nrepeat, unsigned int nthreads) {
		/* if running multi-threaded, we need to use a mutex */
		m = new std::mutex[opts.steps / out_period + 1];
		
		std::vector<uint64_t> seeds;
		std::vector<std::thread> threads;
		for(unsigned int j = 0; j < nthreads; j++) seeds.push_back(nn.get_random()); /* TODO: rng should be public member! */
		std::atomic<unsigned int> ix{0};
		auto tf = [&nn, &opts, nrepeat, &wd, &ix, this] (uint64_t seed) {
			nclass nn2(nn);
			nn2.seed_rng(seed);
			unsigned int k = 0;
			while(true) {
				unsigned int j = ix++;
				if(j >= nrepeat) break;
				if(k) nn2.set_start_pop(opts.start_cell, opts.starting_pop, true, opts.start_pop_fill_factor, opts.start_pop_pchoice);
				do_one_run(nn2, opts.steps, opts.start_cell, wd, *this);
				k++;
			}
		};
		threads.reserve(nthreads);
		for(unsigned int j = 0; j < nthreads; j++) threads.emplace_back(tf, seeds[j]);
		for(unsigned int j = 0; j < nthreads; j++) threads[j].join();
		
		delete[]m;
		m = nullptr;
	}
};


#endif

