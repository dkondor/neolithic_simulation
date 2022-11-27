/*
 * neolithic_w.hpp -- simulation of neolithic farming communities with
 * 	simple model of warfare and conflicts
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



#ifndef NEOLITHIC_W_HPP
#define NEOLITHIC_W_HPP

#include "cell.hpp"
#include "neolithic_base.hpp"
#include <vector>
#include <unordered_set>
#include <utility>
#include <memory>
#include <time.h>

template<class cell_collection>
class neolithic_w : public nbase2<cell_collection> {
protected:
	/* main parameters -- group migration */
	bool mobile_raiders = false; /* if true, "raiders" move to a new cell after each attack */
	bool stationary_raiders_split = true; /* if true, even for stationary raiders, the origin cell is split for the first conflict */
	bool defenders_create_raiders = false; /* if true, raiders are created if a cell is successfully defended from an attack as well */
	bool raiders_can_revert = true; /* if true, raiders can revert to farming if they cannot find a suitable target, or if their attack fails */
	bool raider_dist_real_pop = false; /* if true, raiders consider the real population of target cells for attack */
	bool new_stationary_model = false; /* new model with stationary warriors */
	double raider_dist = 5.0; /* characteristic distance for raiders (if mobile_raiders == false) */
	double raider_max_dist = -1.0; /* separate variable controlling the maximum distance raiders can look for targets */
	double attack_success_prob = 0.5; /* probability of an attacker being successfull */
	double raider_success_dist = -1.0; /* if > 0.0, probability of successful attack depends on the distance for raiders */
	double survivor_ratio = 0.2; /* rate of survivors after (1) a successfull attack by non-mobile raiders; (2) an unsuccessful attack of raiders */
	double praiders = 0.5; /* chance of attackers converting to raiders */
	double raiders_revert_rate = 0.0; /* rate of raiders reverting (only for the new stationary warrior model) */
	double pattack = 1.0; /* rate at which raiders attack in the new stationary model */
	unsigned int empty_delay = 0; /* if > 0, cells that became abandoned cannot be settled for this many years */
	unsigned int raider_attack_threshold = 0; /* raiders don't attack cells below this threshold */
	unsigned int raider_revert_threshold = 0; /* raiders revert back to farming if their population falls below this threshold (after an attack by other raiders) */
	bool peaceful_migrations2 = false; /* if set to true, a farmer group can migrate to an occupied cell without conflict if there is
		enough remaining carrying capacity (but selection can be still preferred for empty cells) */
	
	using nbase<cell_collection>::rng;
	using nbase<cell_collection>::pempty;
	using nbase<cell_collection>::use_global_neighbors;
	using nbase<cell_collection>::group_mig_dist;
	using nbase<cell_collection>::group_mig_pow;
	using nbase<cell_collection>::cph;
	using nbase<cell_collection>::group_mig_prob;
	using nbase<cell_collection>::use_prob_helper;
	using nbase<cell_collection>::peaceful_migrations;
	using group_mig_order = typename nbase<cell_collection>::group_mig_order;
	using cell_prob_helper_type = typename nbase<cell_collection>::cell_prob_helper_type;
	
public:
	using nbase<cell_collection>::g;
	// using nbase<cell_collection>::orig_K; -- not sure if this is required
	
	neolithic_w() : nbase2<cell_collection>(false) {
		use_prob_helper = cell_prob_helper_type::BTCB;
	}

	unsigned int attack_success = 0; /* count the number of successful and failed attacks in each time step */
	unsigned int attack_fail = 0;
	unsigned int raiders_created = 0; /* count the number of raider cells created and reverted */
	unsigned int raiders_reverted = 0;
	
	/* parse command line options for parameters in this class;
	 * i is the current index in argv;
	 * return value:
	 * 	0 if the current option did not match
	 *  >0 if there was a match and i should be incremented by this amount
	 *  -1 if there is an error */
	int parse_option(int argc, char** argv, int i) override {
		int res = nbase<cell_collection>::parse_option(argc, argv, i);
		if(res != 0) return res;
		
		/* at this point, argv[i][0] == '-' */
		switch(argv[i][1]) {
			case 'a':
				attack_success_prob = atof(argv[i+1]);
				res = 2;
				break;
			case 'R':
				switch(argv[i][2]) {
					case 'c':
						defenders_create_raiders = true;
						res = 1;
						break;
					case 'r':
						raiders_can_revert = false;
						res = 1;
						break;
					case 'R':
						raiders_revert_rate = atof(argv[i+1]);
						res = 2;
						break;
					case 'p':
						praiders = atof(argv[i+1]);
						res = 2;
						break;
					case 'P':
						raider_dist_real_pop = true;
						res = 1;
						break;
					case 'A':
						pattack = atof(argv[i+1]);
						res = 2;
						break;
					case 's':
						survivor_ratio = atof(argv[i+1]);
						res = 2;
						break;
					case 'm':
						raider_max_dist = atof(argv[i+1]);
						res = 2;
						break;
					case 'a':
						raider_success_dist = atof(argv[i+1]);
						res = 2;
						break;
					case 't':
						raider_attack_threshold = atoi(argv[i+1]);
						res = 2;
						break;
					case 'T':
						raider_revert_threshold = atoi(argv[i+1]);
						res = 2;
						break;
					case '2':
						new_stationary_model = true;
						res = 1;
						break;
					case 'S':
						stationary_raiders_split = false;
						res = 1;
						break;
					case 0:
						raider_dist = atof(argv[i+1]);
						res = 2;
						break;
					default:
						res = -1;
						break;
				}
				break;
			case 'p':
				peaceful_migrations2 = true;
				if(argv[i][2] == '2') peaceful_migrations = true;
				res = 1;
				break;
			case 'm':
				mobile_raiders = true;
				res = 1;
				break;
			case 'e':
				empty_delay = atoi(argv[i+1]);
				res = 2;
				break;
			case 'b':
				use_prob_helper = cell_prob_helper_type::CDF;
				res = 1;
				break;
		}
		if(mobile_raiders && new_stationary_model) throw std::runtime_error("neolithic_w::parse_options(): Invalid combination (-n and -R2)!\n");
		return res;
	}
	
	std::vector<group_mig_order> orders;
	unsigned int orders_rejected = 0; // number of cells where a split-off was rejected due to not finding a target
	std::vector<typename cell_collection::point> raiders_revert;
	
	/* Create a set of "orders", i.e. migration choices for split-off groups,
	 * to be carried out later. */
	void create_group_orders() {
		orders.clear();
		raiders_revert.clear();
		
		this->ensure_omp_threads();
		
		/* init orders and reverted raiders to the maximum possible size */
		size_t tmpsize = g.size();
		orders.resize(tmpsize);
		raiders_revert.resize(tmpsize);
		size_t norders = 0;
		size_t nrr = 0;
		orders_rejected = 0;

#ifdef _OPENMP		
#pragma omp parallel for
#endif
		for(size_t j = 0; j < tmpsize; j++) {
			auto& rng2 = this->prng();
			const auto& pt1 = g.ptid(j);
			auto& c = g.at_ix(j);
			if(c.st == cell::state::EMPTY) continue;
			unsigned int N3 = 0; // faction size (in case of farming cell)
			
			/* estimate split-off groups and their targets, apply these changes */
			std::uniform_int_distribution<unsigned int> osd;
			
			std::vector<typename cell_collection::point> targets;
			std::vector<double> w;
			std::vector<double> target_dst;
			
			if(c.st == cell::state::FARMING) {
				double ps;
				N3 = this->cell_group_mig(j, rng2, ps);
				/* no migration (very unlikely) -- or should we set a lower limit on the split-off faction? */
				if(!N3) continue;
				/* note: only farmers choose from the long-range distribution,
				 * raiders make separate choices among only farming cells (see below) */
				group_mig_order o;
				unsigned int max_retries = this->target_max_retries;
				if(this->scale_max_retries) {
					max_retries = (unsigned int)std::round(max_retries * ps);
					if(!max_retries) max_retries = 1;
				}
				if(this->choose_dst_helper(pt1, o, max_retries, rng2, N3)) {
					o.group_size = N3;
					o.seq = osd(rng2);
					size_t tmpo;
#ifdef _OPENMP		
#pragma omp atomic capture
#endif
					tmpo = norders++;
					orders[tmpo] = o;
				}
				else {
#ifdef _OPENMP
#pragma omp atomic update
#endif
					orders_rejected++;
				}
			}
			else {
				/* possibility for a raider cell reverting */
				bool revert = false;
				if(raiders_revert_rate > 0.0) {
					std::uniform_real_distribution<double> rrd(0.0, 1.0);
					double x = rrd(rng2);
					if(x < raiders_revert_rate) revert = true;
				}
				
				if(!revert && pattack < 1.0) {
					/* in this case, raiders not necessarily attack in each year */
					std::uniform_real_distribution<double> da(0.0, 1.0);
					double x = da(rng2);
					if(x >= pattack) continue;
				}
				
				/* here, this is a raider cell, so we have to choose
				 * from a "local" distribution */
				double max_dist = raider_max_dist > 0.0 ? raider_max_dist : raider_dist * 10.0;
				if(!revert) for(const auto& tmp1 : g.neighbors(pt1, max_dist)) {
					const auto& pt2 = tmp1.first;
					const auto& c2 = g.at(pt2);
					if(c2.st != cell::state::FARMING && !new_stationary_model) continue; /* only farming cells can be targets */
					if(c2.N < raider_attack_threshold) continue; /* too low population for an attack target */
					double p1 = exp(-1.0 * tmp1.second / raider_dist);
					if(raider_dist_real_pop) p1 *= c2.N;
					else p1 *= c2.K;
					
					if(p1 > 0.0) {
						targets.push_back(pt2);
						w.push_back(p1);
						if(raider_success_dist > 0.0) target_dst.push_back(tmp1.second);
					}
				}
				
				if(w.size()) {
					/* raiding this way is only possible if there are possible target cells */
					std::discrete_distribution<unsigned int> td(w.cbegin(), w.cend());
					group_mig_order o;
					o.src = pt1;
					unsigned int tix = td(rng2);
					o.dst = targets[tix];
					if(raider_success_dist > 0.0) o.dist = target_dst[tix];
					o.group_size = 0;
					o.seq = osd(rng2);
					size_t tmpo;
#ifdef _OPENMP		
#pragma omp atomic capture
#endif
					tmpo = norders++;
					orders[tmpo] = o;
					targets.clear();
					w.clear();
					target_dst.clear();
				}
				else {
					size_t tmpr;
#ifdef _OPENMP		
#pragma omp atomic capture
#endif
					tmpr = nrr++;
					raiders_revert[tmpr] = pt1;
				}
			}
		}
		
		orders.resize(norders);
		raiders_revert.resize(nrr);
	}
	
	template<class CB>
	void apply_group_mig(CB&& cb) {
		std::unordered_set<typename cell_collection::point, typename cell_collection::pointhash> already_targeted;
		
		attack_success = 0;
		attack_fail = 0;
		raiders_created = 0;
		raiders_reverted = 0;
		
		/* process reverted raiders */
		for(const auto& pt1 : raiders_revert) {
			auto& c = g.at(pt1);
			if(raiders_can_revert) c.st = cell::state::FARMING;
			else {
				c.st = cell::state::EMPTY;
				c.N = 0;
				if(empty_delay) c.delay_counter = empty_delay;
			}
			raiders_reverted++;
		}
		
		/* sort the orders randomly */
		group_mig_order::sort_by_seq(orders);
		
		/* process migration events */
		std::uniform_real_distribution<double> pf1;
		for(size_t i = 0; i < orders.size(); i++) {
			const auto& src = orders[i].src;
			const auto& dst = orders[i].dst;
			if(already_targeted.count(src)) {
				/* source cell was attacked already, simply skip this order */
				continue;
			}
			auto& c1 = g.at(src);
			auto& c2 = g.at(dst);
			if(already_targeted.count(dst) || 
				(c2.st == cell::state::EMPTY && c1.st == cell::state::RAIDERS) )  {
				/* target was already attacked, or it was originally a farming cell, but
				 * everyone migrated away -- cancel this order */
				if(c1.st == cell::state::RAIDERS) {
					raiders_reverted++;
					if(raiders_can_revert) c1.st = cell::state::FARMING;
					else {
						c1.st = cell::state::EMPTY;
						c1.N = 0;
						if(empty_delay) c1.delay_counter = empty_delay;
					}
				}
				continue;
			}
			already_targeted.insert(dst);
			
			uint32_t N3 = orders[i].group_size;
			cb(orders[i].src, orders[i].dst, N3);
			
			cell::state dsts = c2.st;
			
			bool conflict = false;
			bool migration = false;
			
			if(c2.st != cell::state::EMPTY) {
				if(c1.st == cell::state::FARMING &&
					peaceful_migrations2 && c2.N + N3 < c2.K)
						migration = true;
				else conflict = true;
			}
			else migration = true; /* note: the case when raiders selected an empty cell was already handled above */
			
			if(conflict) {
				double prob1 = attack_success_prob;
				if(raider_success_dist > 0.0) prob1 *= exp(-1.0 * orders[i].dist / raider_success_dist);
				double x1 = (prob1 < 1.0) ? pf1(rng) : 0.0;
				if(x1 < prob1) {
					/* successful attack */
					attack_success++;
					if(c1.st == cell::state::RAIDERS) migration = mobile_raiders;
					else migration = stationary_raiders_split;
				}
				else {
					conflict = false;
					attack_fail++;
					if(!new_stationary_model && c1.st == cell::state::RAIDERS) {
						/* unsuccessful attack, raiders revert to farming or starve */
						if(raiders_can_revert) {
							c1.N = (unsigned int)std::max(1.0, std::round(survivor_ratio * c1.N));
							c1.st = cell::state::FARMING;
						}
						else {
							c1.N = 0;
							c1.st = cell::state::EMPTY;
							if(empty_delay) c1.delay_counter = empty_delay;
						}
						raiders_reverted++;
					}
					if(defenders_create_raiders) {
						x1 = pf1(rng);
						if(x1 < praiders) {
							c2.st = cell::state::RAIDERS;
							raiders_created++;
						}
					}
				}
			}
			
			if(migration) {
				if(c1.st == cell::state::FARMING) {
					/* this is the case of mobile split-off from a farming cell */
					c2.N = N3;
					c1.N -= N3;
					c2.st = cell::state::FARMING;
					if(conflict) {
						double x1 = pf1(rng);
						if(x1 < praiders) {
							c2.st = cell::state::RAIDERS;
							raiders_created++;
						}
					}
				}
				else { /* c1.st == cell::state::RAIDERS, in this case conflict == true */
					/* this is the case of mobile raiders */
					c2.N = c1.N;
					c2.st = cell::state::RAIDERS;
					c1.st = cell::state::EMPTY;
					c1.N = 0;
					if(empty_delay) c1.delay_counter = empty_delay;
				}
			}
			else if(conflict) {
				/* conflict with stationary raiders or farmers */
				/* target cell is depopulated */
				c2.N = (unsigned int)std::round(survivor_ratio * c2.N);
				if(!c2.N) {
					c2.st = cell::state::EMPTY;
					if(empty_delay) c2.delay_counter = empty_delay;
				}
				else if(c2.st == cell::state::RAIDERS && c2.N < raider_revert_threshold) {
					c2.st = cell::state::FARMING;
				}
				if(c1.st == cell::state::FARMING) {
					/* chance to convert to raiders in this case */
					double x1 = pf1(rng);
					if(x1 < praiders) {
						c1.st = cell::state::RAIDERS;
						raiders_created++;
					}
				}
			}
			
			if(c1.N == 0) {
				/* in (the very unlikely) case that everyone migrated away, mark the source cell as empty */
				c1.st = cell::state::EMPTY;
				if(empty_delay) c1.delay_counter = empty_delay;
			}
			
			if(c2.st != dsts && dsts == cell::state::EMPTY && !this->external_helper) {
				/* cell became occupied might need to create the helper distribution for it */
				if(cph.count(dst) == 0) this->create_cell_prob_helper_one(dst);
			}
			
		}
	}
	
	template<class point>
	struct noop_cb {
		void operator ()(const point& pt1, const point& pt2, unsigned int N) const { }
	};
	
	void apply_group_mig() { apply_group_mig(noop_cb<typename cell_collection::point>()); }
	
	void do_group_mig() {
		create_group_orders();
		apply_group_mig(noop_cb<typename cell_collection::point>());
	}
	
	template<class CB>
	void do_group_mig(CB&& cb) {
		create_group_orders();
		apply_group_mig(cb);
	}
	
	/* step the demographic simulation in each cell */
	void cells_step(double t) {
		this->ensure_omp_threads();
		size_t tmpsize = g.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for(size_t i = 0; i < tmpsize; i++) {
			auto& c = g.at_ix(i); /* cell to update */
			if(c.st == cell::state::EMPTY) continue;
			c.step(t, this->prng(), this->pars);
			/* in the case this cell became empty, we need to update the state */
			if(c.N == 0) c.st = cell::state::EMPTY;
		}
	}
};

#endif

