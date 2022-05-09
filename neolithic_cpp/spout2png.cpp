/*
 * spout2png.cpp -- create a set of png images of spatial distribution outputs
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
#include <stdio.h>
#include <iostream>

static std::string get_label(int year, bool use_bce) {
	if(!use_bce) return std::to_string(year);
	else if(year < 0) return std::to_string(-1*year) + " BCE";
	else return std::to_string(year) + " CE";
}

int main(int argc, char **argv)
{
	int year_step = 10;
	write_image_pars pars;
	poly_output pout;
	pars.raider_pop = true;
	unsigned int width = 1200;
	unsigned int height = 800;
	char* coords_file = nullptr;
	char* net_file = nullptr;
	char* polygon_file = nullptr;
	bool coords_file_csv = false;
	bool net_file_edges = false;
	char* out_base = nullptr;
	bool use_ffmpeg = false;
	const char* codec = "vp9";
	unsigned int video_quality = 15;
	bool have_raiders = true;
	int start_year = 0;
	
	for(int i = 1; i < argc; i++) if(argv[i][0] == '-') switch(argv[i][1]) {
		case 'o':
			switch(argv[i][2]) {
				case 'y':
					start_year = atoi(argv[i+1]);
					i++;
					break;
				case 'p':
					year_step = atoi(argv[i+1]);
					if(!year_step) year_step = 1;
					i++;
					break; /*
				case 's':
					ffmpeg_scale = atoi(argv[i+1]);
					i++;
					break; */
				case 'g':
					pars.gamma = atof(argv[i+1]);
					i++;
					break;
				case 'm':
					pars.max_pop = atoi(argv[i+1]);
					i++;
					break;
				case 'w':
					width = atoi(argv[i+1]);
					i++;
					break;
				case 'h':
					height = atoi(argv[i+1]);
					i++;
					break;
				case 'r':
					pars.raider_pop = false;
					break;
				case 'c':
					/* video codec to use with ffmpeg */
					codec = argv[i+1];
					i++;
					break;
				case 'q':
					video_quality = atoi(argv[i+1]);
					i++;
					break;
				case 'f':
					use_ffmpeg = true;
					/* fallthrough */
				case 0:
					out_base = argv[i+1];
					i++;
					break;
				default:
					fprintf(stderr, "Unknown parameter: %s!\n", argv[i]);
					return 1;
			}
			break;
		case 'c':
			if(argv[i][2] == 'p') polygon_file = argv[i+1];
			else {
				coords_file = argv[i+1];
				if(argv[i][2] == 'c') coords_file_csv = true;
			}
			i++;
			break;
		case 'n':
			net_file = argv[i+1];
			if(argv[i][2] == 'e') net_file_edges = true;
			i++;
			break;
		case 'R':
			have_raiders = false;
			break;
		default:
			fprintf(stderr, "Unknown parameter: %s!\n", argv[i]);
			break;
	}
	else fprintf(stderr, "Unknown parameter: %s!\n", argv[i]);
	
	if(!(polygon_file && out_base)) {
		fprintf(stderr, "Missing parameters!\n");
		return 1;
	}
	
	cellnet n;
	bool use_cellnet = false;
	if(coords_file && net_file) {
		n.load_net(net_file, coords_file, net_file_edges ? cellnet::load_net_flags::read_edges : 0U, coords_file_csv);
		pout.read_poly(polygon_file, false);
		/* figure out the extents of the simulation */
		pout.get_poly_extents([&n](unsigned int id) { return n.cell_exists(id); });
		use_cellnet = true;
	}
	else {
		pout.read_poly(polygon_file, false);
		pout.get_poly_extents();
	}
	pout.lonmin -= 0.5;
	pout.lonmax += 0.5;
	pout.latmin -= 0.5;
	pout.latmax += 0.5;
	
	char* out_fn = nullptr;
	
	FILE* outp = nullptr;
	if(use_ffmpeg) {
		const char* base_cmd = "ffmpeg -y -f rawvideo -pix_fmt bgra -s %ux%u -r 24 -i - -codec:v %s -crf %u %s";
		size_t len1 = strlen(base_cmd) + strlen(out_base) + strlen(codec) + 50;
		out_fn = new char[len1];
		for(int j = 0; j < 2; j++) {
			int tmps = snprintf(out_fn, len1, base_cmd, width, height, codec, video_quality, out_base);
			
			if(tmps > 0 && (size_t)tmps < len1) break;
			delete[]out_fn;
			if(tmps <= 0 || j == 1) {
				fprintf(stderr, "Error creating the output command!\n");
				return 1;
			}
			len1 = tmps + 1;
			out_fn = new char[len1];
		}
		outp = popen(out_fn, "w");
		if(!outp) {
			fprintf(stderr, "Error starting ffmpeg!\n");
			return 1;
		}
		delete[]out_fn;
		out_fn = nullptr;
	}
	else out_fn = new char[strlen(out_base) + 16];
	
	read_table2 rt(std::cin);
	int last_year = start_year;
	std::vector<poly_output::image_label> labels;
	labels.resize(1);
	labels[0].x = 0.85 * width;
	labels[0].y = 0.15 * height;
	labels[0].size = 36;
	
	bool pout_have_img = false;
	
	while(true) {
		unsigned int year0 = 0, id, pop, raiders = 0;
		bool have_line = rt.read_line();
		if(have_line) {
			if(!rt.read(year0, id, pop)) break;
			if(have_raiders && !rt.read(raiders)) break;
		}
		int year = start_year + (int)year0;
		if(!have_line || year != last_year) {
			if(last_year % year_step == 0) {
				/* write out last year */
				labels[0].label = get_label(last_year, start_year < 0);
				if(!use_cellnet) pout.add_labels(labels);
				if(use_ffmpeg) {
					if(use_cellnet) cellnet_write_image(n, outp, width, height, poly_output::image_format::RAW, pars, labels, pout);
					else pout.write_image(outp, poly_output::image_format::RAW);
				}
				else {
					sprintf(out_fn, "%s_%06d.png", out_base, last_year);
					FILE* f = fopen(out_fn, "w");
					if(!f) {
						fprintf(stderr, "Error opening output file %s\n", out_fn);
						delete[]out_fn;
						return 1;
					}
					if(use_cellnet) cellnet_write_image(n, f, width, height, poly_output::image_format::PNG, pars, labels, pout);
					else pout.write_image(f, poly_output::image_format::PNG);
					fclose(f);
				}
				
				/* reset all cells to empty (unfortunately this needs to be done by hand) */
				if(use_cellnet) for(size_t i = 0; i < n.size(); i++) {
					cell& c = n.at_ix(i);
					c.st = cell::state::EMPTY;
					c.N = 0;
				}
				else pout_have_img = false;
			}
			last_year = year;
			if(!have_line) break;
		}
		
		if(year % year_step == 0) {
			if(use_cellnet) {
				cell& c = n.at(id);
				c.st = raiders ? cell::state::RAIDERS : cell::state::FARMING;
				c.N = pop;
			}
			else {
				if(!pout_have_img) {
					pout.start_img(width, height);
					pout_have_img = true;
				}
				double val = std::min((double)pop, pars.max_pop) / pars.max_pop;
				if(!pars.raider_pop && raiders) val = 1.0;
				pout.draw_cell(id, val, raiders ? poly_output::color_scale_t::RED : poly_output::color_scale_t::BLUE);
			}
		}
	}
	if(out_fn) delete[]out_fn;
	if(outp) pclose(outp);
	
	if(rt.get_last_error() != T_EOF) {
		std::cerr<<"Error reading input:\n";
		rt.write_error(std::cerr);
	}
	
	return 0;
}

