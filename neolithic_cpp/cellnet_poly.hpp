/*
 * cellnet_poly.hpp -- helper class to create a series of images from
 * 	the simulation state
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


#ifndef CELLNET_POLY_HPP
#define CELLNET_POLY_HPP
#include "read_table_cpp.h"
#include <iostream>
#include <vector>
#include <unordered_map>
#include <utility>
#include <math.h>
#include <stdexcept>
#include <stdint.h>
#include <stdio.h>
#include <cairo.h>


static cairo_status_t write_cb(void* x, const unsigned char* data, unsigned int length) {
	FILE* f = (FILE*)x;
	if(fwrite(data, 1, length, f) != length) return CAIRO_STATUS_WRITE_ERROR;
	else return CAIRO_STATUS_SUCCESS;
}


class poly_output {
public:
	/* possible image output types */
	enum class image_format {
		RAW, /* raw rgba data (equivalent to CAIRO_FORMAT_ARGB32) */
		PNG /* png, using Cairo's default options */
		/* PPM -- TODO: PPM with header */
	};
	struct image_label {
		std::string label;
		unsigned int x;
		unsigned int y;
		unsigned int size = 12;
	};
	
	enum class color_scale_t {
		BLUE, /* dark to light blue color scale */
		RED, /* dark to light red color scale */
		DIVERGING /* scale going blue->white->red */
	};
	
protected:
	// separate vector for polygons (used when creating images); this is not required to run the simulation
	std::unordered_map<unsigned int, std::vector<std::pair<double, double> > > poly;
	cairo_surface_t* s = nullptr;
	cairo_t* cr = nullptr;
	
	void deallocate() {
		if(cr) { cairo_destroy(cr); cr = nullptr; }
		if(s) { cairo_surface_destroy(s); s = nullptr; }
	}
	
	static void hsv2rgb(double h, double s, double v, double& r, double& g, double& b) {
		double c = s * v;
		h /= 60.0;
		double hmod = h - 2.0 * floor(h / 2.0);
		double x = c * (1 - fabs(hmod - 1));
		if(h < 0) throw std::runtime_error("poly_output::hsv2rgb(): invalid hue value!\n");
		if(h < 1)       { r = c; g = x; b = 0; }
		else if(h < 2)  { r = x; g = c; b = 0; }
		else if(h < 3)  { r = 0; g = c; b = x; }
		else if(h < 4)  { r = 0; g = x; b = c; }
		else if(h < 5)  { r = x; g = 0; b = c; }
		else if(h <= 6) { r = c; g = 0; b = x; }
		else throw std::runtime_error("poly_output::hsv2rgb(): invalid hue value!\n");
		double m = v - c;
		r += m;
		g += m;
		b += m;
	}
	
	static void set_color_from_val(cairo_t* cr, double val, color_scale_t cs) {
		if(cs == color_scale_t::BLUE || cs == color_scale_t::RED) {
			double r, g, b;
			if(val < 0.5) {
				double val1 = 2.0 * val;
				r = 0.196 * val1 + 0.086 * (1.0 - val1);
				g = 0.408 * val1 + 0.192 * (1.0 - val1);
				b = 0.584 * val1 + 0.294 * (1.0 - val1);
			}
			else {
				double val1 = 2.0 * (val - 0.5);
				r = 0.333 * val1 + 0.196 * (1.0 - val1);
				g = 0.686 * val1 + 0.408 * (1.0 - val1);
				b = 0.957 * val1 + 0.584 * (1.0 - val1);
			}
			if(cs == color_scale_t::RED) cairo_set_source_rgba(cr, b, r, g, 1.0);
			else cairo_set_source_rgba(cr, r, g, b, 1.0);
		}
		else {
			/* diverging color scale */
			double v = 0.8;
			double s = 2.0 * fabs(val - 0.5);
			double h = (val < 0.5) ? 0.0 : 240.0;
			double r, g, b;
			hsv2rgb(h, s, v, r, g, b);
			cairo_set_source_rgba(cr, r, g, b, 1.0);
		}
	}
	
	/* current image width and height */
	unsigned int w = 0;
	unsigned int h = 0;
	
public:
	
	~poly_output() { deallocate(); }
	
	/* parameters for image scaling -- these should be set to the extents of the image
	 * and not modified during rendering */
	double lonmin = -180.0;
	double lonmax = 180.0;
	double latmin = -90.0;
	double latmax = 90.0;
	
	/* read polygon coordinates from a given CSV file -- note: these should be ordered */
	void read_poly(const char* fn, bool header = false) {
		poly.clear();
		read_table2 rt(fn);
		rt.set_delim(',');
		if(header) rt.read_line();
		while(rt.read_line()) {
			uint32_t id;
			double lon, lat;
			if(!rt.read(id, lon, lat)) break;
			poly[id].push_back(std::make_pair(lon, lat));
		}
		if(rt.get_last_error() != T_EOF) {
			std::string err = rt.exception_string("poly_output::read_poly: error reading coordinates:\n");
			throw std::runtime_error(err);
		}
	}
	
	
	/* helper to figure out the extents of coordinates */
	template<class F>
	void get_poly_extents(F&& filter) {
		lonmin = 360.0;
		lonmax = -180.0;
		latmin = 90.0;
		latmax = -90.0;
		for(const auto& p1 : poly) if(filter(p1.first)) for(const auto& c : p1.second) {
			if(c.first  < lonmin) lonmin = c.first;
			if(c.first  > lonmax) lonmax = c.first;
			if(c.second < latmin) latmin = c.second;
			if(c.second > latmax) latmax = c.second;
		}
	}
	
	void get_poly_extents() {
		get_poly_extents([] (unsigned int) { return true; });
	}
	
	void start_img(unsigned int w_, unsigned int h_) {
		if(!(w_ && h_)) throw std::runtime_error("poly_output::start_img(): invalid size!\n");
		deallocate();
		w = w_;
		h = h_;
		s = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, w, h);
		if(cairo_surface_status(s) != CAIRO_STATUS_SUCCESS) throw std::runtime_error("poly_output::start_img(): error creating image!\n");
		
		cr = cairo_create(s);
		cairo_set_source_rgba(cr, 0.0, 0.0, 0.0, 1.0);
		cairo_paint(cr); /* background */
	}
	
	void draw_cell(unsigned int id, double val, color_scale_t cs) {
		if(!(cr && s)) throw std::runtime_error("poly_output::draw_cell(): no image (call start_img() first)!\n");
		set_color_from_val(cr, val, cs);
		
		const auto& coords = poly.at(id);
		for(size_t i = 0; i < coords.size(); i++) {
			const auto& y = coords[i];
			double lon = y.first;
			double lat = y.second;
			double x1 = (lon - lonmin) / (lonmax - lonmin);
			double y1 = 1.0 - (lat - latmin) / (latmax - latmin);
			x1 *= (double)w;
			y1 *= (double)h;
			if(i == 0) cairo_move_to(cr, x1, y1);
			else cairo_line_to(cr, x1, y1);
		}
		cairo_fill(cr);
		
		if(cairo_status(cr) != CAIRO_STATUS_SUCCESS)
			throw std::runtime_error("poly_output::draw_cell(): error creating image!\n");
	}
	
	void add_labels(const std::vector<image_label>& labels) {
		if(!(cr && s)) throw std::runtime_error("poly_output::add_labels(): no image (call start_img() first)!\n");
		cairo_select_font_face(cr, "sans-serif", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
		for(const auto& l : labels) {
			cairo_set_font_size(cr, l.size);
			cairo_move_to(cr, l.x, l.y);
			cairo_set_source_rgba(cr, 1.0, 1.0, 1.0, 1.0);
			cairo_show_text(cr, l.label.c_str());
		}
		
		if(cairo_status(cr) != CAIRO_STATUS_SUCCESS)
			throw std::runtime_error("poly_output::add_labels(): error creating image!\n");
	}
	
	void write_image(FILE* f, image_format fmt) {
		if(!(cr && s)) throw std::runtime_error("poly_output::write_image(): no image (call start_img() first)!\n");
		cairo_surface_flush(s);
		
		bool res = true;
		switch(fmt) {
			case image_format::PNG:
				if(cairo_surface_write_to_png_stream(s, write_cb, f) != CAIRO_STATUS_SUCCESS)
					res = false;
				break;
			case image_format::RAW:
				{
					unsigned char* data = cairo_image_surface_get_data(s);
					int stride = cairo_image_surface_get_stride(s);
					// size_t rowsize = 4UL*w;
					for(unsigned int i = 0; i < h; i++) {
						ssize_t off = stride;
						off *= i;
						unsigned char* tmp = data + off;
						if(fwrite(tmp, 4U, w, f) != w) {
							res = false;
							break;
						}
					}
				}
				break;
			default:
				throw std::runtime_error("poly_output::write_image(): unsupported image format!\n");
		}
		if(!res) throw std::runtime_error("poly_output::write_image(): error writing output!\n");
	}
};

#endif

