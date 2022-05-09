/*
 * cellnet_img.hpp -- helper function to write out the state of the simulation
 * (old style)
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


#ifndef CELLNET_IMG_HPP
#define CELLNET_IMG_HPP
#include "cellnet.hpp"
#include "cellnet_poly.hpp"

/* output parameters */
struct write_image_pars {
	double max_pop;
	double max_pix_val = 1.0; // scale maximum values
	double gamma = 1.0;
	bool raider_pop = false;
};

	
void cellnet_write_image(const cellnet& cn, FILE* f, unsigned int w, unsigned int h, poly_output::image_format fmt,
	const write_image_pars& pars, const std::vector<poly_output::image_label>& labels, poly_output& pout);

#endif

