/*
 * btprob.hpp -- maintain a discrete probability distribution with a
 * 	possibility of individual probabilities changing
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
 * This header contains implementations for efficiently sampling from a
 * discrete probability distribution where individual probabilities can
 * change in time. Both changing one of the probabilities and sampling
 * the distribution is achieved in O(log(N)) time (where N is the number
 * of discrete choices). Initial creation of the distribution requires
 * O(N) time. Changing the number of discrete choices is not supported
 * (but elements can be "disabled" by setting their weights to zero).
 * 
 * Note that the distribution stores raw weights / frequencies, i.e.
 * these do not need to be normalized. This is essential, since this
 * way, re-normalization need not be performed when changing one of the
 * weights.
 * 
 * Weights are assumed to be floating point (i.e. double type) numbers.
 * Inaccuracies could arise from rounding errors when updating the
 * distribution or when sampling.
 * 
 * The implementation internally arranges the weights in a binary tree,
 * with each node storing the sum of weights in its induced subtree.
 * This allows searching for any elements based on the cumulative
 * distribution function (CDF) values (order of elements in this case is
 * implicit based on the usual ordering of binary trees). Since the
 * total number of elements is assumed to be known at creation and does
 * not change, the binary tree itself is stored in flat arrays instead
 * of an actual tree structure to save memory.
 * 
 * Sampling from the distribution then requires choosing a value
 * randomly between 0 and the total sum of probabilities (i.e. the value
 * stored in the root of the tree). The corresponding element can then
 * be found efficiently using a binary search on the tree.
 * 
 */

#ifndef BTPROB_HPP
#define BTPROB_HPP

#include <vector>
#include <stdexcept>
#include <functional>

/* Base class that stores both the probabilities and the partial sums in
 * std::vectors; it is assumed that element i has children at index
 * 2*i + 1 and 2*i + 2 */
class btprob_vec_base {
	protected:
		std::vector<double> prob; /* probabilities */
		std::vector<double> sums; /* partial sums (not stored for lead nodes) */
		
		void set_prob1(size_t i, double p) { prob[i] = p; }
		void set_sum1(size_t i, double s) { sums[i] = s; }
	
	public:
		/* access to array elements is encapsulated for derived classes */
		double get_prob(size_t i) const { return prob[i]; }
		double get_sum(size_t i) const { return sums[i]; }
		
		size_t size() const { return prob.size(); }
		/* resize the containers and reset all probabilities to zero */
		void resize(size_t new_size) {
			prob.clear();
			sums.clear();
			
			prob.resize(new_size, 0.0);
			sums.resize(new_size / 2, 0.0);
			/* TODO: should we try to call shrink_to_fit() here? */
		}
		
		size_t prob_size() const { return prob.size(); }
		size_t sums_size() const { return sums.size(); }
		
		explicit btprob_vec_base(size_t new_size = 0) { resize(new_size); }
		
		/* return if it is OK to modify this class after copying
		 * (in this case, this is always true, since the vectors will be copied) */
		bool can_clone_write() const { return true; }
};

/* class to store a possibly changing set of probabilities in a binary
 * tree and sample them
 * note: base should be btprob_vec_base, btprob_cb_base or a compatible
 * base class (e.g. with array based storage)
 */
template<class base>
class btprob : public base {
	protected:
		
		/* element access and getting the size should be implemented
		 * by the base class */
		using base::set_prob1;
		using base::set_sum1;
		using base::prob_size;
		using base::sums_size;
		
		/* convenience wrapper to get the sum of a node -- for leaf
		 * nodes it returns the weights */
		double get_sum_or_prob(size_t i) const {
			if(i < sums_size()) return get_sum(i);
			else return get_prob(i);
		}
		
		constexpr static double adjust_sums_eps = 1e-10; /* numerical precision when adjusting sums */
		
	public:
		
		using base::get_prob;
		using base::get_sum;
		using base::size;
		
		/* get the sum of all elements -- need to be used when sampling
		 * the distribution */
		double total_sum() const {
			if(prob_size()) return get_sum_or_prob(0);
			else return 0.0;
		}
		
		/* get current the weight / frequency of element i */
		double at(size_t i) const {
			if(i >= prob_size()) throw std::out_of_range("Element index out of range!\n");
			return get_prob(i);
		}
		
		/* get current the weight / frequency of element i */
		double operator [](size_t i) const {
			return get_prob(i);
		}
		
		/* return the index of the first element where the sum is
		 * larger than p1, or size() if no such element is found
		 * 
		 * sampling from the distribution is achieved by choosing a
		 * value uniformly at random in the range [0, total_sum() ),
		 * and then calling upper_bound() with this value */
		size_t upper_bound(double p1) const {
			const size_t psize = prob_size();
			if(p1 >= total_sum()) return psize;
			
			double sum1 = 0.0;
			size_t i = 0;
			size_t last_i = 0;
			while(2*i + 1 < psize) { /* while we have children */
				double tmp1 = get_sum_or_prob(2*i + 1); /* sum in the left child */
				double tmp2 = get_prob(i); /* value of i */
				double sum2 = sum1 + tmp1 + tmp2; /* value of i in the sorted list of sums */
				if(p1 < sum2) {
					/* i might be a good candidate or we have to descend to the left */
					if(p1 >= sum1 + tmp1) break; /* i is the good candidate */
					last_i = i;
					i = 2*i + 1;
				}
				else {
					/* in this case, we always have to go to the right */
					if(2*i + 2 < psize) {
						sum1 = sum2;
						i = 2*i + 2;
					}
					else {
						/* this probably should not happen, but is possible due to rounding errors? */
						i = last_i;
						break;
					}
				}
			}
			return i;
		}
		
		/* set the ith probability to the given value
		 * (and update the sums accordingly) */
		void set_prob(size_t i, double p1) {
			set_prob1(i, p1);
			recalculate_sum(i);
		}
		
		/* adjust the stored partial sums after a probability was updated
		 * (only needs to be called if weights are adjusted outside, i.e.
		 * not with the set_prob() function) */
		void recalculate_sum(size_t i) {
			while(true) {
				if(2*i + 1 < prob_size()) {
					double tmp = get_prob(i) + get_sum_or_prob(2*i + 1);
					if(2*i + 2 < prob_size()) tmp += get_sum_or_prob(2*i + 2);
					set_sum1(i, tmp);
				}
				if(!i) break;
				i = (i - 1) / 2;
			}
		}
		
		/* adjust the stored partial sums after multiple updates in the
		 * probabilities for elements in the vector given as argument;
		 * note that the contents of the vector are altered (destroyed) */
		void recalculate_sums(std::vector<size_t>& changed) {
			if(prob_size() <= 1 || changed.size() == 0) return; /* protect against corner cases */
			/* note: we only have to do something for cases where the
			 * index is below prob_size() / 2 -- we alter indices first
			 * to be able to filter out duplicates */
			for(size_t& x : changed)
				if(2*x + 1 >= prob_size()) x = (x - 1) / 2;
			/* sort the indices */
			std::sort(changed.begin(), changed.end());
			/* remove duplicates */
			{
				auto it = std::unique(changed.begin(), changed.end());
				changed.erase(it, changed.end());
			}
			
			/* we process changes per level, going from bottom up */
			
			/* identify the last level */
			size_t cur_level_min = 0;
			while(true) {
				size_t next_level_min = cur_level_min * 2 + 1;
				if(next_level_min > changed.back()) break;
				cur_level_min = next_level_min;
			}
			
			
			while(true) {
				size_t j = changed.size() - 1;
				
				/* 1. process nodes on the current level */
				while(true) {
					size_t i = changed[j];
					double tmp = get_prob(i) + get_sum_or_prob(2*i + 1);
					if(2*i + 2 < prob_size()) tmp += get_sum_or_prob(2*i + 2);
					set_sum1(i, tmp);
					if(!i) break;
					changed[j] = (i - 1) / 2; /* update the index with the parent node which need to be updated next */
					if(!j) break; /* we are at the end of the array */
					if(changed[j-1] < cur_level_min) break; /* next node is already on the next level */
					j--;
				}
				
				if(!cur_level_min) break; /* we finished with the last level */
				
				/* 2. continue to find nodes on the previous level */
				cur_level_min = (cur_level_min - 1) / 2;
				size_t mid_ix = j; /* save the index of the start of the level just processed */
				while(j && changed[j-1] >= cur_level_min) j--;
				
				/* 3. merge the two levels, remove duplicates */
				if(j < mid_ix) std::inplace_merge(changed.begin() + j, changed.begin() + mid_ix, changed.end());
				auto it = std::unique(changed.begin() + j, changed.end());
				changed.erase(it, changed.end());
			}
		}
		
		
		/* Adjust weights after adding to or subtracting from the weight
		 * of a given element. This function only makes sense if the
		 * probabilities are not stored internally. It is the caller's
		 * responsibility to ensure that weights stay positive (negative
		 * weights will result in invalid results when sampling the
		 * distribution). */
		void adjust_sums(size_t i, double diff) {
			/* skip if we are at the last level where sums are not stored */
			if(2*i + 1 >= prob_size()) {
				if(i) i = (i - 1)/2;
				else return; /* corner case: only one element */
			}
			double diff_eps = std::max(std::abs(diff * adjust_sums_eps),
				std::abs(get_sum(0) * adjust_sums_eps));
			while(true) {
				/* here, 2*i + 1 < prob_size(), so sums[i] exists */
				double tmp = get_sum(i);
				tmp += diff;
				if(std::abs(tmp) < diff_eps) tmp = 0.0;
				if(tmp < 0.0) throw std::runtime_error("btprob::adjust_sums(): negative sum detected!\n");
				set_sum1(i, tmp);
				if(!i) break;
				i = (i - 1) / 2;
			}
		}
		
		/* add to or subtract from the weight of element i,
		 * also adjusting the partial sums */
		void adjust_prob(size_t i, double diff) {
			double tmp = get_prob(i);
			tmp += diff;
			if(tmp < 0.0) throw std::runtime_error("btprob::adjust_sums(): negative sum detected!\n");
			set_prob1(i, tmp);
			adjust_sums(i);
		}
		
		/* update all probabilities together and recalculate all sums
		 * (need to provide an iterator that returns the weight of
		 * elements) */
		template<class It, class Sent>
		void set_all_probs(It it, Sent s) {
			size_t i = 0;
			for(; it != s && i < size(); ++i, ++it) set_prob1(i, *it);
			if(i < size()) throw std::runtime_error("btprob::set_all_probs(): not all probabilities given!\n");
			if(it != s) throw std::runtime_error("btprob::set_all_probs(): too many probabilities given!\n");
			create_sums();
		}
		
		/* same, but no upper limit given */
		template<class It>
		void set_all_probs(It it) {
			for(size_t i = 0; i < size(); ++i, ++it) set_prob1(i, *it);
			create_sums();
		}
		
		/* same, but the given function object is called for each possible value of i */
		template<class CB>
		void set_all_probs_cb(CB&& cb) {
			for(size_t i = 0; i < size(); i++) set_prob1(i, cb(i));
			create_sums();
		}
		
		/* re-create the sums array, based on the probabilities;
		 * note: caller has to ensure that probabilities were already
		 * set with some external method */
		void create_sums() {
			size_t size1 = prob_size() / 2; /* size of sums array */
			if(!size1) return;
			for(size_t j = 0; j < size1; j++) {
				size_t i = size1 - j - 1;
				double tmp = get_prob(i) + get_sum_or_prob(2*i + 1);
				if(2*i + 2 < prob_size()) tmp += get_sum_or_prob(2*i + 2);
				set_sum1(i, tmp);
			}
		}
		
		/* checks that the sums are in a consistent state;
		 * returns the absolute difference over the whole tree */
		double check_sums() const {
			double r = 0.0;
			size_t size1 = prob_size() / 2;
			for(size_t i = 0; i < size1; i++) {
				double tmp = get_prob(i) + get_sum_or_prob(2*i + 1);
				if(2*i + 2 < prob_size()) tmp += get_sum_or_prob(2*i + 2);
				r += std::abs(tmp - get_sum(i));
			}
			return r;
		}
};


/* Version of the base class where probabilities are not stored explicitly,
 * but calculated with a callback function when needed. This reduces memory
 * usage to 1/3, but at the cost of recalculating some of the probabilities
 * when sampling the distribution. */
class btprob_cb_base {
	protected:
		size_t size1; /* size of the underlying array (number of discrete choices) */
		std::vector<double> sums_v;
		double* sums = nullptr;
		
		void set_prob1(size_t, double) { throw std::runtime_error("btprob_cb_base::set_prob1(): Setting probabilities directly is not supported!\n"); }
		void set_sum1(size_t i, double s) { sums[i] = s; }
	
	public:
		
		/* access to array elements is encapsulated for derived classes */
		double get_prob(size_t i) const { return prob(i); }
		double get_sum(size_t i) const { return sums[i]; }
		
		/* function that calculates the probabilities, should be set be the user;
		 * note: need to call create_sums() after setting this! */
		std::function<double(size_t)> prob; 
		size_t size() const { return size1; }
		
		/* resize the containers and reset all sums to zero */
		void resize(size_t new_size) {
			sums_v.clear();
			sums_v.resize(new_size / 2, 0.0);
			size1 = new_size;
			sums = sums_v.data();
		}
		
		size_t prob_size() const { return size1; }
		size_t sums_size() const { return size1 / 2; }
		
		/* set the sums based on an externally supplied array;
		 * s must be a new_size / 2 size array that will be used subsequently;
		 * note that it is not freed by this class (that is the caller's responsibility) */
		void set_sums_external(double* s, size_t new_size) {
			sums_v.clear();
			sums_v.shrink_to_fit();
			sums = s;
			size1 = new_size;
		}
		
		explicit btprob_cb_base(size_t new_size = 0) { resize(new_size); }
		btprob_cb_base(const btprob_cb_base& bt) : size1(bt.size1), sums_v(bt.sums_v), sums(bt.sums) {
			if(sums_v.size()) sums = sums_v.data(); /* in this case, we made a copy of the values */
		}
		btprob_cb_base(btprob_cb_base&& bt) = default;
		
		/* return if it is OK to modify this class after copying
		 * (in this case, this depends on whether the sums array was supplied externally) */
		bool can_clone_write() const {
			if(!size1) return true; /* if nothing is allocated, it is OK */
			return (sums_v.size() > 0); /* otherwise, it is only OK, if we are using a vector that is copied */
		}
};


#endif

