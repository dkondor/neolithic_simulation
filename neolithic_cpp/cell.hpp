/*
 * cell.hpp -- a cell represents a local area where demographic simulations take place
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

#ifndef CELL_HPP
#define CELL_HPP

#include <stdint.h>
#include <cmath>
#include <algorithm>
#include <random>
#include <stdexcept>

/**
 * A cell is the basic unit of the simulation. Demographics in cells are
 * simulated using a logistic equation (individual people are not considered).
 * This class contains the variables defining the current state of a cell
 * and basic functionality to update the local population.
 * 
 * Cells are organized in a landscape in the container classes:
 * cellnet and cellgrid (in their respective header)
 * Typically they are not used directly.
 */
struct cell {
	/* possible states of a cell */
	enum class state {
		EMPTY = 0, /* cell is empty, no agricultural village present */
		FARMING, /* farming village is present */
		RAIDERS /* group of warriors who live by raiding nearby farming villages */
		/* additional states could represent foragers or warrior groups */
	};
	
	uint32_t N = 0; // current population
	double K = 0.0; // carrying capacity (number of people) -- this can change in time generally, or multiplied by factors dependant on local conditions
	state st = state::EMPTY; // current state (note: st == EMPTY implies N == 0 and vice versa)
	uint32_t delay_counter = 0; // counter that can be used for certain types of delays -- not used in the main simulation
	
	uint32_t out_mig = 0; // current estimate of out-migration -- not used in the main simulation
	uint32_t in_mig = 0; // current estimate of in-migration -- not used in the main simulation
	
	uint32_t settled_year = 0; // date when this cell was first settled
	uint32_t settled_pop = 0; // initial population -- these can be used to estimate an "expected" population number at any time
	
	// some simulation parameters (note: these are not stored for each cell, but defined here since they belong here logically)
	struct cell_pars {
		double r = 0.0135; // birth rate (per year)
		double delta = 0.8; // population collapse if larger than carrying capacity
		double cr = 0.1; // probability of a cultural change happening (one bit flipping) in a year
	};
	
	// advance the demographic model by t time (years); above parameters are given as a variable
	template<class RNG>
	void step(double t, RNG& rng, const cell_pars& pars = default_pars());
	
	// default cell parameters
	static constexpr cell_pars default_pars() {
		return cell_pars();
	}
	
	/* cells can be created by specifying their carrying capacity */
	explicit cell(double K_ = 0.0) : K(K_) { }
};


template<class RNG>
void cell::step(double t, RNG& rng, const cell_pars& pars) {
	double target = K; /* note: actual carrying capacity could be a function of K */
	double r = pars.r;
	double delta = pars.delta;
	
	if(delay_counter) delay_counter--;
	
	if(N == 0) throw std::runtime_error("cell::step(): empty cell found!\n");
	
	if(N > target) {
		/* population collapses to delta * target */
		N = std::round(delta * target);
		if(N > target) N = std::floor(target);
	}
	else {
		/* logistic growth */
		double tmp = N;
		tmp = (target - tmp) / tmp;
		tmp *= exp(-1.0 * r * t);
		tmp = target / (1 + tmp);
		/* tmp is the expected population based on the solution of the 
		 * logistic differential equation
		 * we calculate actual growth based on a Poisson distribution */
		double gr = tmp - N;
		if(gr > 0.0) {
			std::poisson_distribution<unsigned int> pd(gr);
			N += pd(rng);
			if(N > target) N = std::floor(target);
		}
	}
}


#endif

