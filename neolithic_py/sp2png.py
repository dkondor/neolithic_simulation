#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  sp2png.py -- read the output of a simulation, convert to a series of
# 	PNG images
#  
#  Copyright 2021 Daniel Kondor <kondor@csh.ac.at>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

import os, math
import cairo as cr

def set_color_from_val(val, red = False):
	r = 0.0
	g = 0.0
	b = 0.0
	if val < 0.5:
		val1 = 2.0 * val
		r = 0.196 * val1 + 0.086 * (1.0 - val1)
		g = 0.408 * val1 + 0.192 * (1.0 - val1)
		b = 0.584 * val1 + 0.294 * (1.0 - val1)
	else:
		val1 = 2.0 * (val - 0.5)
		r = 0.333 * val1 + 0.196 * (1.0 - val1)
		g = 0.686 * val1 + 0.408 * (1.0 - val1)
		b = 0.957 * val1 + 0.584 * (1.0 - val1)
	if red:
		tmp = r
		r = b
		b = g
		g = tmp
	return (r, g, b)



def main(args):
	
	polygons = dict()
	poly_filter = None
	
	poly_fn = None
	poly_filter_fn = None
	poly_filter_csv = False
	
	out_base = None
	
	width = 1200
	height = 800
	use_geo_ratio = False
	
	max_pop = 1
	gamma = 1.0
	
	i = 1
	while i < len(args):
		if args[i][0] != '-':
			raise BaseException('Invalid argument: {}!\n'.format(args[i]))
		if args[i][1] == 'p':
			poly_fn = args[i+1]
			i += 1
		elif args[i][1] == 'f':
			if len(args[i]) > 2 and args[i][2] == 'c':
				poly_filter_csv = True
			poly_filter_fn = args[i+1]
			i += 1
		elif args[i][1] == 'o':
			out_base = args[i+1]
			i += 1
		elif args[i][1] == 'w':
			width = int(args[i+1])
			i += 1
		elif args[i][1] == 'h':
			height = int(args[i+1])
			i += 1
		elif args[i][1] == 'g':
			use_geo_ratio = True
		elif args[i][1] == 'm':
			max_pop = float(args[i+1])
			i += 1
		elif args[i][1] == 'G':
			gamma = float(args[i+1])
			i += 1
		else:
			raise BaseException('Invalid argument: {}!\n'.format(args[i]))
		i += 1
	
	if poly_fn is None or out_base is None or width <= 0 or height <= 0:
		raise BaseException('Invalid or missing arguments!\n')
	
	# read a list of "filter" (if not all polygons should be used)
	if poly_filter_fn is not None:
		poly_filter = set()
		with open(poly_filter_fn, 'r') as f:
			if poly_filter_csv:
				f.readline()
			for l in f.readlines():
				l2 = l.split(',' if poly_filter_csv else None)
				poly_filter.add(int(l2[0]))
	
	lonmin = 180.0
	lonmax = -180.0
	latmin = 90.0
	latmax = -90.0
	with open(poly_fn, 'r') as f:
		for l in f.readlines():
			l2 = l.split(',')
			id1 = int(l2[0])
			if poly_filter and id1 not in poly_filter:
				continue
			c1 = (float(l2[1]), float(l2[2]))
			lonmin = min(c1[0], lonmin)
			lonmax = max(c1[0], lonmax)
			latmin = min(c1[1], latmin)
			latmax = max(c1[1], latmax)
			if id1 not in polygons:
				polygons[id1] = list()
			polygons[id1].append(c1)
	
	# calculate scaling factors
	dlon = lonmax - lonmin
	dlat = latmax - latmin
	rlon = width / (dlon + 0.2)
	rlat = height / (dlat + 0.2)
	xpad = 0
	ypad = 0
	
	if use_geo_ratio:
		clat = (latmin + latmax) / 2.0
		gr = math.cos(math.pi * clat / 180.0)
		rlon2 = rlat * gr
		if rlon2 < rlon:
			rlon = rlon2
			w2 = rlon * (dlon + 0.2)
			# note: w2 < width
			xpad = (width - w2) / 2.0
		elif rlon2 > rlon:
			rlat = rlon / gr
			h2 = rlat * (dlat + 0.2)
			ypad = (height - h2) / 2.0
	
	def write_image(out_base, year, s):
		if out_base == '-':
			# write raw image to stdout
			sys.stdout.buffer.write(s.get_data())
		else:
			# write to png
			fn = '{}_{:06d}.png'.format(out_base, year)
			s.write_to_png(fn)
	
	s = None
	ctx = None
	last_year = None
	
	for l in sys.stdin.readlines():
		l2 = l.split()
		year = int(l2[0])
		id1 = int(l2[1])
		pop = int(l2[2])
		raiders = bool(int(l2[3]))
		
		if not last_year or year > last_year:
			if last_year:
				ctx = None
				s.flush()
				write_image(out_base, last_year, s)
			s = cr.ImageSurface(cr.FORMAT_ARGB32, width, height)
			ctx = cr.Context(s)
			ctx.set_source_rgba(0.0, 0.0, 0.0, 1.0)
			ctx.paint()
			last_year = year
		elif year < last_year:
			raise BaseException('Input not ordered!\n')
		
		poly1 = polygons[id1]
		for i in range(len(poly1)):
			lon = poly1[i][0]
			lat = poly1[i][1]
			x = (lon - lonmin + 0.1) * rlon + xpad
			y = height - ypad - (lat - latmin + 0.1) * rlat
			if i == 0:
				ctx.move_to(x, y)
			else:
				ctx.line_to(x, y)
		
		val1 = min(pop, max_pop) / max_pop
		if gamma != 1.0:
			val1 = math.pow(val1, gamma)
		c1 = set_color_from_val(val1, raiders)
		
		ctx.set_source_rgba(*c1, 1.0)
		ctx.fill()
	
	if ctx:
		ctx = None
		s.flush()
		write_image(out_base, last_year, s)
	
	return 0


if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))


