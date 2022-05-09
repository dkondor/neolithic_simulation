/*
 * cellnet_img.cpp -- write current simulation state as a png image
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


#include "cellnet_img.hpp"


void cellnet_write_image(const cellnet& cn, FILE* f, unsigned int w, unsigned int h, poly_output::image_format fmt,
		const write_image_pars& pars, const std::vector<poly_output::image_label>& labels, poly_output& pout) {
	pout.start_img(w, h);
	
	/* go through all cells, draw them as filled polygons */
	for(size_t j = 0; j < cn.size(); j++) {
		unsigned int pt = cn.ptid(j);
		const cell& c = cn.at_ix(j);
		
		if(c.st == cell::state::EMPTY || c.N == 0) {
			if(c.N) throw std::runtime_error("cellnet::write_image(): empty cell found with nonzero population!\n");
			continue;
		}
		
		double val = std::min((double)c.N, pars.max_pop) / pars.max_pop;
		if(pars.gamma != 1.0) val = std::pow(val, pars.gamma);
		val *= pars.max_pix_val;
		
		pout.draw_cell(pt, val, c.st == cell::state::FARMING ? poly_output::color_scale_t::BLUE :
			poly_output::color_scale_t::RED);
	}
	
	/* write any labels */
	pout.add_labels(labels);
	
	/* write out the image */
	pout.write_image(f, fmt);
}

