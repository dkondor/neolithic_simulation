#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  cellnet.py -- classes implementing the basic units (cells) of the simulation
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

from enum import Enum
import math, struct
import sortedcontainers as sc
import numpy as np


class cell:
	"""
	Base unit of the simulation. Represents a geographic area roughly on the scale of one village.
	"""
	def __init__(self, K = 0.0):
		self.N = 0 # population
		self.K = K # carrying capacity
		self.st = cell.state.EMPTY # state
		

	class state(Enum):
		EMPTY = 0
		FARMING = 1
		RAIDERS = 2
	
	class pars:
		def __init__(self):
			self.r = 0.0135
			self.delta = 1.0
	
	def step(self, t, rng, pars):
		"""
		Advance the demographic simulation in this cell by t years.
		rng should be a random number generator to use
		pars should be an object of cell.pars with the demographic simulator parameters
		"""
		target = self.K
		if self.N > target:
			# population collapse
			tmp = int(pars.delta * target)
			if tmp > target:
				tmp = int(target)
			self.N = tmp
		else:
			# population growth
			tmp = (target - self.N) / self.N * math.exp(-1.0 * t * pars.r)
			tmp = target / (1 + tmp) # tmp is the expected value of the population after this step
			d1 = tmp - self.N
			if d1 <= 0:
				return
			self.N += rng.poisson(d1)
			if self.N > target:
				self.N = int(target)
		if self.N == 0:
			# if the population got wiped out (e.g. since the carrying capacity was reset to < 1),
			# change the cell's state to empty as well
			self.st = cell.state.EMPTY
	

class spdist:
	"""
	Class storing a symmetric distance matrix, optionally loaded from a file
	"""
	def __init__(self):
		self.mp1 = None # file mapping storing the distance matrix
		self.base_offset = 0 # offset into mp1 where the matrix starts
		self.reader = None # struct.Struct to read from mp1
		self.size = 0 # number of nodes
		self.dists = None # array with stored elements (if allocated by us)
	
	def get_ix(self, i, j):
		"""
		Get the index of the distance between two nodes, given by their
		0-based integer indices (used internally).
		"""
		if i >= self.size or j >= self.size:
			raise IndexError('spdist.get_ix(): node indices are out of range!\n')
		if i < j:
			tmp = i
			i = j
			j = tmp
		# here i >= j
		return (i*(i+1)) // 2 + j;
		
	
	def get_dist(self, i, j):
		"""
		Get the distance between two nodes, given by their 0-based integer indices
		"""
		ix1 = self.get_ix(i, j)
		if self.mp1 is not None:
			tmp1 = self.reader.unpack_from(self.mp1, self.base_offset + ix1 * self.reader.size)
			return tmp1[0]
		else:
			return self.dists[ix1]
	
	def set_dist(self, i, j, d1):
		"""
		Set the distance between two nodes, given by their 0-based integer indices.
		Only works if the storage was reserved by allocate()
		"""
		if self.dists is None:
			raise BaseException('spdist.set_dist(): need to allocate storage first!\n')
		ix1 = self.get_ix(i, j)
		self.dists[ix1] = d1
	
	def allocate(self, size, dtype = np.double):
		"""
		Allocate space for storing a distance matrix among size nodes,
		i.e. a matrix of size*(size+1)/2, with the given datatype
		(default is double).
		"""
		self.mp1 = None
		self.base_offset = 0
		self.reader = None
		self.size = size
		self.dists = np.zeros((size * (size+1)) // 2, dtype = dtype)
	
	# binary file header base
	HEADER_BASE = 0x635b054329591700
	
	# binary file header bits identifying possible data types
	HEADER_INT16 = 0x1
	HEADER_UINT16 = 0x2
	HEADER_INT32 = 0x3
	HEADER_UINT32 = 0x4
	HEADER_INT64 = 0x5
	HEADER_UINT64 = 0x6
	HEADER_FLOAT = 0x7 # np.float32
	HEADER_DOUBLE = 0x8 # np.double == np.float64
	HEADER_FLAGS = 0xf # combination of all above flags
	

	
	def write(self, f):
		"""
		Write this matrix to a binary file opened by the caller as f.
		Returns the number of bytes written or 0 on error.
		"""
		
		if self.size == 0 or self.dists is None:
			raise BaseException('spdist.write(): nothing to write!\n')
		
		ret = 0
		header1 = self.HEADER_BASE
		tmp1 = self.dists.dtype
		tmp2 = None
		if tmp1 == np.int16:
			header1 += self.HEADER_INT16
			tmp2 = '=h'
		elif tmp1 == np.uint16:
			header1 += self.HEADER_UINT16
			tmp2 = '=H'
		elif tmp1 == np.int32:
			header1 += self.HEADER_INT32
			tmp2 = '=i'
		elif tmp1 == np.uint32:
			header1 += self.HEADER_UINT32
			tmp2 = '=I'
		elif tmp1 == np.int64:
			header1 += self.HEADER_INT64
			tmp2 = '=q'
		elif tmp1 == np.uint64:
			header1 += self.HEADER_UINT64
			tmp2 = '=Q'
		elif tmp1 == np.float32:
			header1 += self.HEADER_FLOAT
			tmp2 = '=f'
		elif tmp1 == np.double:
			header1 += self.HEADER_DOUBLE
			tmp2 = '=d'
		else:
			raise BaseException('spdist.write(): unsupported data type!\n')
	
		# 1. write the header
		ret2 = f.write(struct.pack('=Q', header1))
		if ret2 == 0:
			return 0
		ret += ret2
		# 2. write the size
		ret2 = f.write(struct.pack('=Q', self.size))
		if ret2 == 0:
			return 0
		ret += ret2
		
		# TODO: padding? (should not be required, but the C++ code has it)
		tmp1 = (self.size * (self.size + 1)) // 2
		helper1 = struct.Struct(tmp2)
		for i in range(tmp1):
			ret2 = f.write(helper1.pack(self.dists[i]))
			if ret2 == 0:
				return 0
			ret += ret2
		
		return ret
	
	def read(self, base, base_off, base_size):
		"""
		set up the distance matrix from the given memory area
		(corresponding to a file written by write() and opened by mmap());
		returns the number of bytes consumed
		note: it is the caller's responsibility to ensure that the
		given memory stays valid while this class is used
		"""
		self.mp1 = None
		self.base_offset = 0
		self.reader = None
		self.size = 0
		self.dists = None
		
		if base_size < base_off + 16:
			raise BaseException('spdist.read(): input is too small!\n')
		header1 = struct.unpack_from('=Q', base, base_off)
		if header1[0] & (~self.HEADER_FLAGS) != self.HEADER_BASE:
			raise BaseException('spdist.read(): invalid header!\n')
		
		header1 = header1[0] - self.HEADER_BASE
		tmp1 = None
		data_size = 0 # size of one element, in bytes
		if header1 == self.HEADER_INT16:
			tmp1 = '=h'
			data_size = 2
		elif header1 == self.HEADER_UINT16:
			tmp1 = '=H'
			data_size = 2
		elif header1 == self.HEADER_INT32:
			tmp1 = '=i'
			data_size = 4
		elif header1 == self.HEADER_UINT32:
			tmp1 = '=I'
			data_size = 4
		elif header1 == self.HEADER_INT64:
			tmp1 = '=q'
			data_size = 8
		elif header1 == self.HEADER_UINT64:
			tmp1 = '=q'
			data_size = 8
		elif header1 == self.HEADER_FLOAT:
			tmp1 = '=f'
			data_size = 4
		elif header1 == self.HEADER_DOUBLE:
			tmp1 = '=d'
			data_size = 8
		else:
			raise BaseException('spdist.read(): invalid header or unsupported data type!\n')
		
		size1 = struct.unpack_from('=Q', base, base_off + 8)
		size1 = size1[0]
		total_size = 16 + data_size * (size1 * (size1 + 1)) // 2
		if base_off + total_size > base_size:
			raise BaseException('spdist.read(): input is too small!\n')
		
		self.mp1 = base
		self.base_offset = base_off + 16
		self.size = size1
		self.reader = struct.Struct(tmp1)
		
		return total_size


	


class cellnet:
	"""
	A collection of cells organized by neighbor connections.
	"""
	
	def __init__(self):
		self.cells = dict() # main storage of cells
		self.n = dict() # neighbors for each cell
		self.coords = dict() # center coordinates of cells
		self.coastal_cells = set() # set of coastal cells -- travel among these is assumed to be faster by the below factor
		self.coastal_travel_factor = 1.0 # (by default it is the same as between land cells)
		self.land_travel_factors = dict() # extra factors affecting travel
		self.use_net_dist = False # set to true if distances should be measured along the edges of the cell network
		self.spd = None # pre-calculated network distances
		self.spd_ids = None # dictionary of indices to use for the above
		self.edge_dists = None # distances calculated or read for each edge
		
	def __getitem__(self, key):
		"""
		Access cells based on their IDs
		"""
		return self.cells[key]
	
	def __len__(self):
		return len(self.cells)
	
	def __iter__(self):
		return self.cells.__iter__()
	
	def __contains__(self, key):
		return key in self.cells
	
	def step(self, t, rng, pars):
		"""
		Step the demographic simulation in each cell with t years.
		"""
		for x in self.cells:
			c1 = self.cells[x]
			if c1.st != cell.state.EMPTY:
				if c1.N <= 0:
					raise BaseException('Empty cell found: {}, {}, {}, {}!\n'.format(x,
						'farmers' if c1.st == cell.state.FARMING else 'raiders', c1.N, c1.K))
				c1.step(t, rng, pars)
			
			
	###############################################################
	# functions for travel distances and iterating over neighbors
	
	def edge_dist_helper(self, id1, id2):
		"""
		Calculate the appropriately scaled travel distance along a network edge
		"""
		dist = self.coords_dist(id1, id2)
		if id1 in self.coastal_cells and id2 in self.coastal_cells:
			dist *= self.coastal_travel_factor;
		else:
			f1 = 1.0
			f2 = 1.0
			if id1 in self.land_travel_factors:
				f1 = self.land_travel_factors[id1]
			if id2 in self.land_travel_factors:
				f2 = self.land_travel_factors[id2]
			dist *= (f1 + f2) / 2.0;
		return dist
	
	def get_edge_dist(self, id1, id2):
		if self.edge_dists is not None:
			if id1 < id2:
				return self.edge_dists[(id1, id2)]
			else:
				return self.edge_dists[(id2, id1)]
		else:
			return self.edge_dist_helper(id1, id2)
	
	def calculate_edge_dists(self):
		"""
		Calculate scaled distances along each edge to be used later
		"""
		self.edge_dists = dict()
		for id1 in self.n:
			for id2 in self.n[id1]:
				d1 = self.edge_dist_helper(id1, id2)
				if id1 < id2:
					self.edge_dists[(id1, id2)] = d1
				else:
					self.edge_dists[(id2, id1)] = d1

			
	def coords_dist(self, id1, id2):
		"""
		Return the great circle distance between the center of the two cells given by their IDs.
		"""
		c1 = self.coords[id1]
		c2 = self.coords[id2]
		RADIUS = 6371.0
		lon1 = math.pi * c1[0] / 180.0
		lat1 = math.pi * c1[1] / 180.0
		lon2 = math.pi * c2[0] / 180.0
		lat2 = math.pi * c2[1] / 180.0

		s1 = math.sin((lat1 - lat2) / 2.0)
		s2 = math.sin((lon1 - lon2) / 2.0)
		r1 = s1 * s1 + s2 * s2 * math.cos(lat1) * math.cos(lat2)
		return RADIUS * 2.0 * math.asin(math.sqrt(r1))
		
	def point_dist(self, id1, id2):
		"""
		Return the distance among two cells.
		"""
		if self.spd is not None:
			return self.spd.get_dist(self.spd_ids[id1], self.spd_ids[id2])
		else:
			return self.coords_dist(id1, id2)
	
	# neighbor iterator for finding nodes by a BFS starting from a given node
	def neighbors(self, start_node, max_dist, include_self = False):
		current = start_node
		cd = 0.0
		if include_self:
			yield (current, cd)
		
		node_dist = dict()
		node_dist[current] = cd
		nodes_visited = set()
		nodes_visited.add(current)
		queue = sc.SortedList()
		
		while True:
			# add all the neighbors of current to the queue
			for x in self.n[current]:
				if x not in nodes_visited:
					cd2 = cd + self.get_edge_dist(current, x) if self.use_net_dist else self.point_dist(start_node, x)
					if x not in node_dist:
						node_dist[x] = cd2
						queue.add((cd2, x))
					else:
						cdtmp = node_dist[x]
						if cd2 < cdtmp:
							queue.remove((cdtmp, x))
							queue.add((cd2, x))
							node_dist[x] = cd2
			# get the next node in the queue unless we are at the end
			if len(queue) == 0:
				break
			cd, current = queue.pop(0)
			if cd > max_dist:
				break
			nodes_visited.add(current)
			yield (current, cd)
				
	
	def load_net(self, net_fn, coords_fn, coords_csv = False, net_edges = False, read_edge_dists = False,
			read_types = False, read_travel_factors = False):
		"""
		Read a set of cells from the given input files.
		"""
		self.cells.clear()
		self.n.clear()
		self.coords.clear()
		self.coastal_cells.clear()
		self.land_travel_factors.clear()
		self.spd = None
		self.edge_dists = None
		if read_edge_dists:
			self.edge_dists = dict()
		
		# read the coordinates first and create all cell objects
		with open(coords_fn, 'r') as cf:
			if coords_csv:
				cf.readline()
			for l in cf.readlines():
				l2 = l.split(',' if coords_csv else None)
				if len(l2) < 4:
					raise BaseException('Invalid data in input file {}!\n'.format(coords_fn))
				id1 = int(l2[0])
				lon = float(l2[1])
				lat = float(l2[2])
				K = float(l2[3])
				c1 = cell(K)
				if id1 in self.cells or id1 in self.coords:
					raise BaseException('Duplicate cell ID: {}!\n'.format(id1))
				self.cells[id1] = c1
				self.coords[id1] = (lon, lat)
				self.n[id1] = set() # create the neighbor list so it can be used later
				
				i = 4
				try:
					if read_types:
						type1 = int(l2[i])
						if type1 == 1:
							self.coastal_cells.add(id1)
						i += 1
					if read_travel_factors:
						factor1 = float(l2[i])
						self.land_travel_factors[id1] = factor1
						i += 1
				except:
					raise BaseException('Invalid data in input file {}!\n'.format(coords_fn))
		
		# read the neighbor connections
		with open(net_fn, 'r') as nf:
			for l in nf.readlines():
				l2 = l.split()
				id1 = int(l2[0])
				if id1 not in self.cells:
					continue
				if net_edges:
					id2 = int(l2[1])
					if id2 in self.cells:
						self.n[id1].add(id2)
						if read_edge_dists:
							d1 = float(l2[2])
							if id1 < id2:
								self.edge_dists[(id1, id2)] = d1
							else:
								self.edge_dists[(id2, id1)] = d1
				else:
					for i in range(1, len(l2)):
						id2 = int(l2[i])
						if id2 in self.cells:
							self.n[id1].add(id2)
		
		# TODO: check that all nodes have neighbors and the network is symmetric!
	
	def read_edges(self, fn, read_dist = False, make_symmetric = True):
		"""
		Read a set of additional edges from a given file, supplied as an edgelist.
		"""
		with open(fn, 'r') as f:
			for l in f.readlines():
				l2 = l.split()
				if len(l2) == 0:
					continue
				if len(l2) < (3 if read_dist else 2):
					raise BaseException('Invalid data in input file {}!\n'.format(fn))
				id1 = int(l2[0])
				id2 = int(l2[1])
				if id1 in self.cells and id2 in self.cells:
					self.n[id1].add(id2)
					if make_symmetric:
						self.n[id2].add(id1)
				
				if read_dist:
					d1 = float(l2[2])
					if id1 < id2:
						self.edge_dists[(id1, id2)] = d1
					else:
						self.edge_dists[(id2, id1)] = d1
		
	
	
	def create_distance_matrix(self):
		"""
		Create a matrix of distances among all cells.
		"""
		self.spd = None
		self.spd_ids = dict()
		
		tmpids = sorted(list(self.cells)) # list of cell IDs
		ncells = len(tmpids)
		for i in range(ncells):
			self.spd_ids[tmpids[i]] = i
		
		spd1 = spdist()
		spd1.allocate(ncells, np.uint16)
		
		for i in range(ncells):
			for id2, dst in self.neighbors(tmpids[i], 1e9, False):
				d2 = round(dst)
				if d2 > 65535 or d2 < 0:
					raise BaseException('cellnet.create_distance_matrix(): invalid distance!\n')
				spd1.set_dist(i, self.spd_ids[id2], d2)
		self.spd = spd1
		self.use_net_dist = False
	
	
	def scale_K_lat(self):
		"""
		Scale carrying capacities by a factor dependent on the latitude. This is intended for the
		case of rectangular grids in geographical coordinates where the size of cells depends on
		the latitude.
		"""
		for pt in self.cells:
			c1 = self.cells[pt]
			lat = self.coords[pt][1]
			factor = math.cos(math.pi * lat / 180.0)
			c1.K *= factor
	
	def scale_K(self, fn, zero_missing = False, is_csv = True, header = True):
		"""
		Scale carrying capacities based on the factors read from the given file.
		If zero_missing == True, capacities for all cells that do not appear in the
		file are set to zero.
		"""
		with open(fn, 'r') as f:
			all_ids = set(self.cells) if zero_missing else None
			if header:
				f.readline()
			for l in f.readlines():
				l2 = l.split(',' if is_csv else None)
				id1 = int(l2[0])
				factor = float(l2[1])
				if id1 in self.cells:
					c1 = self.cells[id1]
					c1.K *= factor
					if zero_missing:
						all_ids.remove(id1)
			if zero_missing:
				for id2 in all_ids:
					self.cells[id2].K = 0.0
	
	def adjust_K(self, K1, Kthresh = 0.0):
		"""
		Adjust carrying capacities so that the average is K1. Note that
		empty cells (with zero capacity) are not included in the average
		calculation.
		If Kthresh > 0.0, cells with carrying capacity below this value
		are iteratively set to zero.
		"""
		i = 0
		while True:
			sumK = 0.0
			sumN = 0.0
			Kt = 0
			for pt in self.cells:
				c1 = self.cells[pt]
				if c1.K > 0.0:
					if i > 0 and c1.K < Kthresh:
						c1.K = 0.0
						Kt += 1
					elif K1 is not None:
						sumK += c1.K
						sumN += 1
			if (i > 0 or K1 is None) and Kt == 0:
				break
			if sumN == 0.0:
				break
			factor = K1 / (sumK / sumN)
			for pt in self.cells:
				c1 = self.cells[pt]
				if c1.K > 0.0:
					c1.K *= factor
			if Kthresh == 0.0:
				break
			i += 1
	
	HEADER = 0x0b5c22f6d3deb068 # expected header value
	HEADER_HAVE_SPD = 0x1 # full distance matrix is stored
	HEADER_HAVE_EDGE_DISTS = 0x2 # edge distances are stored along with the edges
	HEADER_FLAGS = 0x3 # combination of all possible flags
	
	def read(self, base, base_offset, base_size):
		"""
		Read the cell IDs, coordinates and carrying capacities from the (memory-mapped) binary file
		opened at base, starting from base_offset, up to maximum base_size.
		Return the number of bytes read or zero on error.
		"""
		self.cells.clear()
		self.coords.clear()
		self.n.clear()
		self.coastal_cells.clear()
		self.land_travel_factors.clear()
		self.spd = None
		self.spd_ids = None
		self.edge_dists = None
		
		ret = 0
		if base_offset + 24 > base_size:
			raise BaseException('cellnet.read(): input file too small!\n')
		
		
		tmp1 = struct.unpack_from('=Q', base, base_offset)
		if tmp1[0] & (~self.HEADER_FLAGS) != self.HEADER:
			raise BaseException('cellnet.read(): Unexpected file header!\n')
		have_spd = bool(tmp1[0] & self.HEADER_HAVE_SPD)
		have_edge_dists = bool(tmp1[0] & self.HEADER_HAVE_EDGE_DISTS)
		tmp1 = struct.unpack_from('=Q', base, base_offset + 8)
		nnodes = tmp1[0]
		tmp1 = struct.unpack_from('=Q', base, base_offset + 16)
		nedges = tmp1[0]
		ret = 24
		
		size1 = ret + nnodes*20 + nedges*8
		if have_edge_dists:
			size1 += 8*nedges
		if nnodes % 2:
			size1 += 4 # padding
		
		if base_offset + size1 > base_size:
			raise BaseException('cellnet.read(): input file too small!\n')
		
		# 1. read the node IDs
		node_ids = list()
		helper1 = struct.Struct('=I')
		for i in range(nnodes):
			tmp1 = helper1.unpack_from(base, base_offset + ret + i*4)
			id1 = tmp1[0]
			node_ids.append(id1)
		ret += nnodes * 4
		
		if nnodes % 2:
			ret += 4
		
		helper2 = struct.Struct('=dd')
		for i in range(nnodes):
			tmp1 = helper2.unpack_from(base, base_offset + ret + i*16)
			id1 = node_ids[i]
			self.cells[id1] = cell()
			self.coords[id1] = tmp1 # note: tmp1 is a tuple with two coordinates
			self.n[id1] = set()
		ret += nnodes * 16
		
		helper3 = struct.Struct('=II')
		helper4 = struct.Struct('=d')
		tmpoff = base_offset + ret
		if have_edge_dists:
			self.edge_dists = dict()
		for i in range(nedges):
			tmp1 = helper3.unpack_from(base, tmpoff)
			id1 = tmp1[0]
			id2 = tmp1[1]
			self.n[id1].add(id2)
			tmpoff += 8
			
			if have_edge_dists:
				d1 = helper4.unpack_from(base, tmpoff)
				tmpoff += 8
				if id1 < id2:
					self.edge_dists[(id1, id2)] = d1
				else:
					self.edge_dists[(id2, id1)] = d1
			
		ret += nedges * 16 if have_edge_dists else nedges * 8
		
		if have_spd:
			self.spd = spdist()
			ret2 = self.spd.read(base, tmpoff, base_size)
			if ret2 == 0:
				raise BaseException('cellnet.read(): Error reading the distance matrix!\n')
			ret += ret2
			# need to create a mapping from node IDs to indices
			# which are based on the position of sorted node IDs
			if self.spd.size != len(node_ids):
				raise BaseException('cellnet.read(): Size of the distance matrix does not match!\n')
			node_ids.sort()
			self.spd_ids = dict()
			for i in range(len(node_ids)):
				self.spd_ids[node_ids[i]] = i
		
		return ret
		
	def write(self, f):
		"""
		Write this class (coordinates and network) to a binary file opened as f).
		Returns the number of bytes written
		"""
		ret = 0
		header1 = self.HEADER
		have_edge_dists = self.edge_dists is not None
		if have_edge_dists:
			header1 |= self.HEADER_HAVE_EDGE_DISTS
		if self.spd is not None:
			header1 |= self.HEADER_HAVE_SPD
		ret += f.write(struct.pack('=Q', header1))
		nnodes = len(self.cells)
		nedges = sum(len(self.n[pt]) for pt in self.cells)
		ret += f.write(struct.pack('=Q', nnodes))
		ret += f.write(struct.pack('=Q', nedges))
		
		# write node IDs
		helper1 = struct.Struct('=I')
		for x in self.n:
			ret += f.write(helper1.pack(x))
		if ret % 8:
			pad = 8 - (ret % 8)
			for i in range(pad):
				ret += f.write(struct.pack('x'))
		
		# write the coordinates (note: we need to use the same order as before)
		helper2 = struct.Struct('=dd')
		for x in self.n:
			c1 = self.coords[x]
			ret += f.write(helper2.pack(*c1))
		
		# write the network (as a sequence of edges)
		helper3 = struct.Struct('=II')
		helper4 = struct.Struct('=d')
		for x in self.n:
			for y in self.n[x]:
				ret += f.write(helper3.pack(x, y))
				if have_edge_dists:
					ret += f.write(helper4.pack(self.get_edge_dist(x, y)))
		
		if self.spd is not None:
			ret2 = self.spd.write(f)
			if ret2 == 0:
				# TODO: error handling?
				return 0
			ret += ret2
		
		return ret
		
		
