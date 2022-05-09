/*
 * neolithic.hpp -- implementation class for agent-based simulation of
 * 	neolithic settlements; version without warfare
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

#ifndef NEOLITHIC_HPP
#define NEOLITHIC_HPP

#include "neolithic_base.hpp"

template<class cell_collection>
class neolithic : public nbase2<cell_collection> {
protected:
	bool helper_batch_update = false;
	
	using nbase<cell_collection>::group_pop_estimate;
	using nbase<cell_collection>::target_only_free;
	using nbase<cell_collection>::rng;
	using nbase<cell_collection>::pempty;
	using nbase<cell_collection>::use_global_neighbors;
	using nbase<cell_collection>::group_mig_dist;
	using nbase<cell_collection>::group_mig_pow;
	using nbase<cell_collection>::group_mig_prob;
	using nbase<cell_collection>::cph;
	using nbase<cell_collection>::use_prob_helper;
	using group_mig_order = typename nbase<cell_collection>::group_mig_order;
	using cell_prob_helper_type = typename nbase<cell_collection>::cell_prob_helper_type;
	
public:
	using nbase<cell_collection>::g;
	
	neolithic() : nbase2<cell_collection>(true) {
		group_pop_estimate = true;
		this->target_max_retries = 10;
		this->target_only_free = true;
	}

	int parse_option(int argc, char** argv, int i) override {
		int res = nbase<cell_collection>::parse_option(argc, argv, i);
		if(res != 0) return res;
		
		/* at this point, argv[i][0] == '-' */
		switch(argv[i][1]) {
			case 'p':
				group_pop_estimate = false;
				res = 1;
				break;
			case 'H':
				if(argv[i][2] == 0 || argv[i][2] == '2') {
					cell_prob_helper_type tmp = (argv[i][2] == '2') ?
						cell_prob_helper_type::BTCB : cell_prob_helper_type::BTVEC;
					if(use_prob_helper != cell_prob_helper_type::NONE && tmp != use_prob_helper)
						throw std::runtime_error("neolithic::parse_option(): conflicting options!\n");
					use_prob_helper = tmp;
					res = 1;
				}
				else if(argv[i][2] == 'b') {
					helper_batch_update = true;
					res = 1;
				}
				break;
		}
		return res;
	}
	
	std::vector<group_mig_order> orders;
	unsigned int orders_rejected = 0; // number of cells where a split-off was rejected due to not finding a target
	
	/* Create a set of "orders", i.e. migration choices for split-off groups,
	 * to be carried out later. */
	void create_group_orders() {
		orders.clear();
		
		this->ensure_omp_threads();
		
		size_t tmpsize = g.size();
		size_t norders = 0;
		orders.resize(tmpsize);
		unsigned int nrejected = 0;

#ifdef _OPENMP		
#pragma omp parallel
#endif
		{ /* parallel scope to have targets and w as thread-local variables */
			std::vector<typename cell_collection::point> targets;
			std::vector<double> w;
		
#ifdef _OPENMP		
#pragma omp for
#endif
			for(size_t pt_ix = 0; pt_ix < tmpsize; pt_ix++) {
				const auto& pt1 = g.ptid(pt_ix);
				double ps; /* unused */
				unsigned int N3 = this->cell_group_mig(pt_ix, this->prng(), ps);
				if(!N3) continue;
				
				if(cph.size()) {
					/* use the helper probability distributions -- ensure that it exists */
					auto& helper = cph.at(pt1);
					if(!helper.btpw) this->create_cell_prob_helper_one(pt1);
					else if(helper_batch_update) {
						if(helper.invalidate) {
							this->create_cell_prob_helper_one(pt1);
							helper.invalidate = false;
							helper.changed_ix.clear();
						}
						else if(helper.changed_ix.size() > 0) adjust_helper_batched(helper);
					}
					group_mig_order o;
					if(this->choose_dst_helper(pt1, o, this->target_max_retries, this->prng(), target_only_free ? N3 : 0)) {
						o.group_size = N3;
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
						nrejected++;
					}
				}
				else {
					/* this cell has a split-off event, we need to choose a target
					 * for the migrating group */
					double max_dist = group_mig_dist * (group_mig_pow > 0.0 ? 1.0 : 10.0);
					for(const auto& tmp1 : g.neighbors(pt1, max_dist)) {
						const auto& pt2 = tmp1.first;
						/* estiamte pt1 -> pt2 migration probability */
						double p1 = group_mig_prob(pt2, tmp1.second);
						if(p1 > 0.0) {
							const auto& c2 = g.at(pt2);
							if(target_only_free && c2.N + N3 > c2.K) continue;
							targets.push_back(pt2);
							w.push_back(p1);
						}
					}
					
					if(w.size()) {
						/* migration is only possible if there are neighbor cells with free capacities (this should always be the case) */
						std::discrete_distribution<unsigned int> td(w.cbegin(), w.cend());
						group_mig_order o;
						o.src = pt1;
						o.dst = targets[td(this->prng())];
						o.group_size = N3;
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
						nrejected++;
					}
				
					targets.clear();
					w.clear();
				}
			}
		} /* end of pragma omp parallel -- scope of thread-local targets and w */
		
		orders.resize(norders);
		orders_rejected = nrejected;
	}
	
	/* estimate split-off groups and their targets, apply these changes
	 * return the number of successful migrations and attempted, but rejected ones (due to no target) */
	template<class CB>
	void apply_group_mig(CB&& cb) {
		/* sort the orders, choose one randomly for all cases where multiple groups have targeted the same cell */
		group_mig_order::sort_by_ptid(orders);
		size_t k = 0;
		for(size_t i = 0; i < orders.size();) {
			size_t j = i + 1;
			for(; j < orders.size(); j++) if(orders[j].dst != orders[i].dst) break;
			if(j == i + 1) {
				/* no duplicates */
				if(i > k) orders[k] = orders[i];
				i++;
				k++;
			}
			else {
				/* we have a set of duplicates, choose only one */
				std::uniform_int_distribution<size_t> d1(i, j - 1);
				size_t x = d1(rng);
				if(x > k) orders[k] = orders[x];
				k++;
				i = j;
			}
		}
		if(k < orders.size()) orders.erase(orders.begin() + k, orders.end());
		
		/* process migration events */
		for(size_t i = 0; i < orders.size(); i++) {
			auto& c1 = g.at(orders[i].src);
			auto& c2 = g.at(orders[i].dst);
			uint32_t N3 = orders[i].group_size;
			
#ifdef MIG_DEBUG
			{
				auto coords1 = g.get_coords(orders[i].dst);
				if(coords1.second > 52.0) {
					fprintf(stderr, "Settled cell at %f, %f\n", coords1.first, coords1.second);
				}
			}
#endif
			
			cb(orders[i].src, orders[i].dst, N3);
			
			if(c2.st == cell::state::EMPTY) {
				if(use_prob_helper != cell_prob_helper_type::NONE)
					if(cph.count(orders[i].dst) == 0) cph.emplace(std::piecewise_construct,
						std::forward_as_tuple(orders[i].dst), std::forward_as_tuple());
				
				c2.st = cell::state::FARMING;
				if(use_prob_helper != cell_prob_helper_type::NONE && !group_pop_estimate)
					update_cell_prob_helper(orders[i].dst);
			}
			else if(c2.N + N3 > c2.K && target_only_free) continue;
			
			c2.N += N3; /* add to the new cell's population */
			c1.N -= N3; /* apply population decrease */
			if(c1.N == 0) c1.st = cell::state::EMPTY;
			
			if(use_prob_helper != cell_prob_helper_type::NONE) {
				if(group_pop_estimate) {
					/* we need to update the helper probabilities */
					update_cell_prob_helper(orders[i].src);
					update_cell_prob_helper(orders[i].dst);
				}
				else if(c1.N == 0) update_cell_prob_helper(orders[i].src);
			}
		}
	}
	
	void apply_group_mig() {
		apply_group_mig([](const auto&, const auto&, unsigned int) { });
	}
	
	/* step the demographic simulation in each cell, updating the
	 * migration base probabilities if needed */
	void cells_step(double t) {
		size_t tmpsize = g.size();
		this->ensure_omp_threads();
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for(size_t i = 0; i < tmpsize; i++) {
			auto& c = g.at_ix(i); /* cell to update */
			if(c.st == cell::state::EMPTY) continue;
			unsigned int N1 = c.N;
			c.step(t, this->prng(), this->pars);
			/* if the population changed, we need to update the migration probabilities */
			if(c.N == 0) c.st = cell::state::EMPTY;
			if(use_prob_helper != cell_prob_helper_type::NONE &&
					((N1 != c.N && group_pop_estimate) || c.N == 0)) {
#ifdef _OPENMP		
#pragma omp critical
#endif
				update_cell_prob_helper(g.ptid(i));
			}
		}
	}
	
	/* override that updates the helpers if necessary */
	bool set_cell_state(const typename cell_collection::point& pt, cell::state st, unsigned int pop) override {
		bool tmp = nbase<cell_collection>::set_cell_state(pt, st, pop);
		if(tmp || group_pop_estimate) update_cell_prob_helper(pt);
		return tmp;
	}
	
	/* reset all cells to zero population -- this is an optimization to update all helpers at once */
	void reset() override {
		/* in this case, it is easier to clear everything instead of updating the
		 * distributions, which would be an O(N*N*log(n)) operation */
		cph.clear();
		for(size_t i = 0; i < g.size(); i++) {
			auto& c = g.at_ix(i);
			c.st = cell::state::EMPTY;
			c.N = 0;
		}
	}
	
	uint64_t helper_total_updates = 0;
	
protected:
	
	void adjust_one_helper(const typename cell_collection::point& pt, const typename cell_collection::point& pt2,
			size_t ix, typename nbase<cell_collection>::cell_prob_helper& helper) {
		switch(use_prob_helper) {
			case cell_prob_helper_type::BTCB:
				{
					auto& btpw = static_cast<typename nbase<cell_collection>::template btprob_wrapper<btprob_cb_base>&>(*helper.btpw);
					auto& btp = btpw.btp;
/*					if(!group_pop_estimate) {
						const cell& c = g.at(pt);
						double dst = g.point_dist(pt, pt2);
						/ * note: w is the probability without caring whether the cell is empty;
						 * i.e. it is not zeroed or multiplied by pempty *
						double w = group_mig_prob(pt, dst, true);
						
						if(target_only_free) {
							if(c.st != cell::state::EMPTY) w *= -1;
						}
						else {
							if(c.st != cell::state::EMPTY) w *= (1.0 - pempty);
							else w *= (pempty - 1.0);
						}
						
						btp.adjust_sums(ix, w);
					}
					else */ btp.recalculate_sum(ix);
				}
				break;
			case cell_prob_helper_type::BTVEC:
				{
					auto& btpw = static_cast<typename nbase<cell_collection>::template btprob_wrapper<btprob_vec_base>&>(*helper.btpw);
					auto& btp = btpw.btp;
					double dst = g.point_dist(pt, pt2);
					btp.set_prob(ix, group_mig_prob(pt, dst));
				}
				break;
			default:
				throw std::runtime_error("neolithic::update_cell_prob_helper(): unsupported helper type!\n");
		}
		helper_total_updates++;
	}
	
	void adjust_helper_batched(typename nbase<cell_collection>::cell_prob_helper& helper) {
		switch(use_prob_helper) {
			case cell_prob_helper_type::BTCB:
				{
					auto& btpw = static_cast<typename nbase<cell_collection>::template btprob_wrapper<btprob_cb_base>&>(*helper.btpw);
					auto& btp = btpw.btp;
					btp.recalculate_sums(helper.changed_ix);
					helper.changed_ix.clear();
				}
				break;
			case cell_prob_helper_type::BTVEC:
				{
					throw std::runtime_error("neolithic::adjust_helper_batched(): not implemented yet for probability vectors!\n");
					/* TODO: batch update of probabilities *
					auto& btpw = static_cast<typename nbase<cell_collection>::template btprob_wrapper<btprob_vec_base>&>(*helper.btpw);
					auto& btp = btpw.btp;
					double dst = g.point_dist(pt, pt2);
					btp.set_prob(ix, group_mig_prob(pt, dst)); */
				}
				break;
			default:
				throw std::runtime_error("neolithic::update_cell_prob_helper(): unsupported helper type!\n");
		}
		helper_total_updates++;
	}
	
	/* update the probability helper for all neighbors of pt
	 * (after a change in pt) */
	void update_cell_prob_helper(const typename cell_collection::point& pt) {
		double max_dist = group_mig_dist * (group_mig_pow > 0.0 ? 1.0 : 10.0);
		if(use_global_neighbors) {
			size_t global_ix = g.ptix(pt);
			
			for(const auto& pt2 : g) {
				auto it = cph.find(pt2);
				if(it == cph.end()) continue;
				auto& helper = it->second;
				if(!helper.btpw) continue;
				adjust_one_helper(pt, pt2, global_ix, helper);
			}
		}
		else for(const auto& tmp : g.neighbors(pt, 1.2 * max_dist)) {
			const auto& pt2 = tmp.first;
			auto it2 = cph.find(pt2);
			if(it2 == cph.end()) continue;
			auto& helper = it2->second;
			if(!helper.btpw) continue;
			
			auto it = std::lower_bound(helper.neighbors, helper.neighbors + helper.nb_size, pt);
			if(it != helper.neighbors + helper.nb_size && *it == pt) {
				/* we really need to make an update */
				size_t i = it - helper.neighbors;
				if(helper_batch_update) {
					size_t max_size;
					if(helper.nb_size) max_size = std::max(1UL, helper.nb_size / 10UL);
					else max_size = std::max(1UL, g.size() / 10UL);
					if(helper.changed_ix.size() >= max_size) {
						helper.invalidate = true;
						helper.changed_ix.clear();
					}
					else helper.changed_ix.push_back(i);
				}
				else adjust_one_helper(pt, pt2, i, helper);
			}
		}
	}
};

#endif

