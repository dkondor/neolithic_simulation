#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  neolithic.py
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

import os, sys, math, struct, bisect, mmap
import numpy as np
import btprob as btp
from cellnet import cell, cellnet


class mmap_array:
	"""
	Helper class for a (read-only) array in a memory-mapped file.
	"""
	def __init__(self):
		self.mp1 = None
		self.off1 = 0
		self.r = struct.Struct('=I')
		self.size1 = 0
	
	def set_array(self, mp1, off, size):
		self.mp1 = mp1
		self.off1 = off
		self.size1 = size
	
	def __getitem__(self, i):
		if i >= self.size1:
			raise BaseException('mmap_array.__getitem__: index out of range!\n')
		if self.mp1 is None:
			raise BaseException('mmap_array.__getitem__: no data!\n')
		x = self.r.unpack_from(self.mp1, self.off1 + 4*i)
		return x[0]


class neolithic:
	"""
	Main class holding the parameters used in the simulation and calculating the interactions among cells.
	"""
	def __init__(self):
		# network of cells that the simulation is run on
		self.g = cellnet()
		# main parameters
		self.cap_mig_exp = 4.0 # exponent for split-off probabilities 
		self.cap_mig_max = 0.25 # probability for split-off if the population is equal to the carrying capacity 
		self.cap_mig_min = 0.0 # minimum share of people to trigger a split-off based on carrying capacity 
		self.dunbar_limit = 150.0 # Dunbar number for village populations 
		self.dunbar_mig_exp = 2.0 # exponent for Dunbar limit based split-off 
		self.dunbar_mig_base = 0.1 # probability for split-off if population is equal to the Dunbar-number 
		self.dunbar_mig_min = 0.0 # minimum share of people to trigger a migration (for Dunbar limit based split-off) 
		self.group_mig_dist = None # characteristic distance for group migration 
		self.target_only_free = False # if true, migrations only target empty cells -- in this case, the probability distributions should be regularly updated
		self.pempty = 1.0 # extra factor for preference toward "emtpy" (unsettled) cells 
		self.group_mig_pow = None # if > 0.0, group migration distance weighting is done by a power-law function with this exponent,
									# and group_mig_dist gives the limit on cells to consider 
		self.have_pars = False # set to true if any of pempty, group_mig_dist or group_mig_pow is given as a command line argument 
		self.mobile_raiders = False # if true, "raiders" move to a new cell each year 
		self.defenders_create_raiders = False # if true, raiders are created if a cell is successfully defended from an attack as well 
		self.raiders_can_revert = True # if true, raiders can revert to farming if they cannot find a suitable target, or if their attack fails 
		self.raider_dist = 5.0 # characteristic distance for raiders (if mobile_raiders == false) 
		self.attack_success_prob = 0.5 # probability of an attacker being successfull 
		self.survivor_ratio = 0.2 # rate of survivors after (1) a successfull attack by non-mobile raiders; (2) an unsuccessful attack of raiders 
		self.praiders = 0.5 # chance of attackers converting to raiders 
		# self.empty_delay = 0 # if > 0, cells that became abandoned cannot be settled for this many years
		self.pars = cell.pars()
		
		# additional parameters
		self.raiders_revert_rate = 0.0 # if > 0 every raider cell has a yearly chance of reverting back to farming by itself
		self.pattack = 1.0 # chance of a raider cell doing an attack every year (i.e. if < 1, raider do not attack every year)
		self.raiders_dist_real_pop = False # if True, raiders consider the total population when selecting an attack target
		self.new_stationary_model = False # new model in which raiders are stationary, but they are not eliminated on a failed attack
		
		# additional parameters not yet implemented
		# raider_max_dist -- '-Rm' custom setting for the maximum distance that the raiders look for targets
		# raider_success_dist -- '-Ra' attack success depends on the distance (exponentially), this is the characteristic distance

		
		self.max_retries = 1000 # number of times choosing a destination is attempted
		
		self.orig_K = dict() # original yield values (saved if needed)
		
		# random number generator
		self.rng = np.random.default_rng()
		
		# memory mapped file with the pre-calculated probability distribution sums -- need to be opened later
		self.matrix_f = None
		self.mp1 = None
		# list of probability helpers
		self.helpers = dict()
		self.helpers_cdf = False # set to true if helpers should use CDF instead of binary trees
		
		# list of node IDs in the order used in the above distributions
		self.global_ids = list()
		
		# per-node lists of neighbors (a dict if used, None otherwise)
		self.neighbors = None
		
		# orders -- group migration events in one step
		self.orders = list()
		self.raiders_revert = list()
	
	def parse_arg(self, args, i):
		"""
		Try to parse one command line argument in args, starting from position i.
		Returns the number of elements taken (including the option argument) or
		0 if the argument is not recognized, or -1 if it is invalid.
		"""
		if len(args[i]) < 2 or args[i][0] != '-':
			return -1
		if args[i][1] == 's':
			self.rng = np.random.default_rng(int(args[i+1]))
			return 2
		elif args[i][1] == 'r':
			self.pars.r = float(args[i+1])
			return 2
		elif args[i][1] == 'd':
			self.pars.delta = float(args[i+1])
			return 2
		elif args[i][1] == 'G':
			self.group_mig_dist = float(args[i+1])
			return 2
		elif args[i][1] == 'P':
			self.group_mig_pow = float(args[i+1])
			return 2
		elif args[i][1] == 'R':
			if len(args[i]) == 2:
				self.raider_dist = float(args[i+1])
				return 2
			elif args[i][2] == 'c':
				self.defenders_create_raiders = True
				return 1
			elif args[i][2] == 'r':
				self.raiders_can_revert = False
				return 1
			elif args[i][2] == 'p':
				self.praiders = float(args[i+1])
				return 2
			elif args[i][2] == 'R':
				self.raiders_revert_rate = float(args[i+1])
				return 2
			elif args[i][2] == 'A':
				self.pattack = float(args[i+1])
				return 2
			elif args[i][2] == 'P':
				self.raiders_dist_real_pop = True
				return 1
			elif args[i][2] == '2':
				self.new_stationary_model = True
				if self.mobile_raiders:
					raise BaseException('Incompatible parameter combination (only one of -m or -R2 can be given)!\n')
				return 1
			elif args[i][2] == 's':
				self.survivor_ratio = float(args[i+1])
				return 2
			else:
				return -1
		elif args[i][1] == 'm':
			self.mobile_raiders = True
			if self.new_stationary_model:
				raise BaseException('Incompatible parameter combination (only one of -m or -R2 can be given)!\n')
			return 1
		elif args[i][1] == 'E':
			self.pempty = float(args[i+1])
			return 2
		elif args[i][1] == 'e':
			raise BaseException('neolithic.parse_arg(): delay for re-settling abandoned cells is not supported!\n')
		elif args[i][1] == 'D':
			if len(args[i]) < 3:
				return -1
			elif args[i][2] == 'e':
				self.dunbar_mig_exp = float(args[i+1])
				return 2
			elif args[i][2] == 'b':
				self.dunbar_mig_base = float(args[i+1])
				return 2
			elif args[i][2] == 'l':
				self.dunbar_limit = float(args[i+1])
				return 2
			elif args[i][2] == 'm':
				self.dunbar_mig_min = float(args[i+1])
				return 2
			else:
				return -1
		elif args[i][1] == 'C':
			if len(args[i]) < 3:
				return -1
			elif args[i][2] == 'e':
				self.cap_mig_exp = float(args[i+1])
				return 2
			elif args[i][2] == 'b':
				self.cap_mig_max = float(args[i+1])
				return 2
			elif args[i][2] == 'm':
				self.cap_mig_min = float(args[i+1])
				return 2
			else:
				return -1
		elif args[i][1] == 'a':
			self.attack_success_prob = float(args[i+1])
			return 2
		elif args[i][1] == 'F':
			return 1 # -F is accepted, but has no effect (we always use fixed probabilities)
		elif args[i][1] == 'H':
			# settings for the probability helper distributions
			if len(args[i]) > 2:
				if args[i][2] == 'C':
					self.helpers_cdf = True
					return 1
				elif args[i][2] == '0':
					self.helpers = None # disable helpers
					return 1
			return 0 # the caller might want to interpret this option
		elif args[i][1] == 'l':
			self.target_only_free = True
			return 1
		elif args[i][1] == 'L':
			self.neighbors = dict()
			return 1
		else:
			return 0
	
	def __del__(self):
		self._close_file()
	
	def _close_file(self):
		if self.mp1 is not None:
			self.mp1.close()
			self.mp1 = None
		if self.matrix_f is not None:
			os.close(self.matrix_f)
			self.matrix_f = None
		self.helpers = dict()
	
	
	def group_mig_prob(self, id1, id2):
		"""
		Get the probability of migrating from cell id1 to id2
		"""
		if self.target_only_free and self.g[id2].st != cell.state.EMPTY:
			return 0.0
		p1 = 0.0
		dst = self.g.point_dist(id1, id2)
		if self.group_mig_pow is not None and self.group_mig_pow > 0.0:
			p1 = pow(dst, -1.0 * self.group_mig_pow)
		else:
			p1 = np.exp(-1.0 * dst / self.group_mig_dist)
		
		if len(self.orig_K) > 0:
			p1 *= self.orig_K[id2]
		else:
			p1 *= self.g[id2].K
		return p1
	
	# facilities for calculating the probabilities and choosing targets
	def _get_ix(self, id1):
		# get the index of the given id in the probability matrix
		ix1 = bisect.bisect_left(self.global_ids, id1, 0, len(self.global_ids))
		if ix1 >= len(self.global_ids) or self.global_ids[ix1] != id1:
			raise BaseException('Node ID not found: {}!\n'.format(id1))
		return ix1
	
	
	def create_one_helper(self, src, use_global_neighbors):
		"""
		Create the migration probability distribution for one source node.
		Returns a tuple with the helper and (optionally) the list of neighbors used.
		Note: this does not store the created helper in self.helpers (and self.neighbors),
		the caller should do this if required!
		"""
		helper = None
		neighbors = None
		local_neighbors = None
		if use_global_neighbors:
			neighbors = self.global_ids
		else:
			local_neighbors = list()
			max_dist = self.group_mig_dist 
			if self.group_mig_pow is None or self.group_mig_pow <= 0.0:
				max_dist *= 10.0
			for pt1, dst in self.g.neighbors(src, max_dist):
				local_neighbors.append(pt1)
			local_neighbors.sort()
			neighbors = local_neighbors
			
			
		if self.helpers_cdf:
			helper = btp.btprob_cdf()
			helper.set_probs(lambda ix: self.group_mig_prob(src, neighbors[ix]), len(neighbors))
		else:
			helper = btp.btprob_tree()
			helper.set_all_probs(len(neighbors), lambda ix, src=src, neighbors=neighbors: self.group_mig_prob(src, neighbors[ix]))
		return (helper, local_neighbors)
	
	def update_helper_neighbors(self, src):
		"""
		Update migration probabilities after a change in the state of the given cell.
		"""
		if self.helpers_cdf:
			raise BaseException('neolithic.update_helper_neighbors(): Cannot update CDF distribution!\n')
		if self.neighbors is not None:
			# we are using local neighbors, we only have to go through the neighbors of the current node
			# (note: we assume that distance and such neighbor relations are symmetric)
			tmpn = None
			if src in self.neighbors:
				tmpn = ((x, 0) for x in self.neighbors[src])
			else:
				max_dist = self.group_mig_dist 
				if self.group_mig_pow is None or self.group_mig_pow <= 0.0:
					max_dist *= 10.0
				tmpn = self.g.neighbors(src, max_dist)
			for x,_ in tmpn:
				if x in self.helpers:
					h = self.helpers[x]
					n2 = self.neighbors[x]
					ix = bisect.bisect_left(n2, src, 0, len(n2))
					if ix >= len(n2) or n2[ix] != src:
						raise BaseException('neolithic.update_helper_neighbors(): neighbor relation ID not found: {} -- {}!\n'.format(x, src))
					h.recalculate_sum(ix)
					
		else:
			# global neighbors, we have to go through all existing helpers
			# find the index of the current node
			ix = self._get_ix(src)
			for y, h in self.helpers:
				if y != src:
					h.recalculate_sum(ix)
				

	def choose_dst_nohelper(self, src, N3 = 0):
		"""
		Choose a destination for a group migrating out of src and return it or None if none found.
		Version without using pre-computed probability distributions.
		"""
		targets = list()
		w = list()
		sumw = 0.0
		for pt1, dst in self.g.neighbors(src, 1e9):
			# go through all possible points
			c2 = self.g[pt1]
			if N3 > 0 and c2.K < N3:
				continue
			p1 = 0.0
			if self.group_mig_pow is not None and self.group_mig_pow > 0.0:
				p1 = pow(dst, -1.0 * self.group_mig_pow)
			else:
				p1 = np.exp(-1.0 * dst / self.group_mig_dist)
			
			if len(self.orig_K) > 0:
				p1 *= self.orig_K[pt1]
			else:
				p1 *= c2.K
			
			if c2.st == cell.state.EMPTY:
				p1 *= self.pempty
			
			targets.append(pt1)
			w.append(p1)
			sumw += p1
		
		if len(targets) > 0:
			for i in range(len(w)):
				w[i] /= sumw
			return self.rng.choice(targets, p = w)
		else:
			return None

	
	def choose_dst(self, src, N3 = 0):
		"""
		Choose a destination for a group migrating out of src and return it or None if none found.
		If N3 > 0, an empty cell is only chosen if it has capacity for at least N3 people.
		"""
		if self.helpers is None:
			return self.choose_dst_nohelper(self, src, N3)
		
		tmp1 = 1.0 / self.pempty
		helper = None
		neighbors = None # if local neighbors are used
		
		if src in self.helpers:
			helper = self.helpers[src]
			if self.neighbors is not None:
				if src in self.neighbors:
					neighbors = self.neighbors[src]
				else:
					raise BaseException('Inconsistent state in the migration probability distributions (src: {})!'.format(src))
		else:
			(helper, neighbors) = self.create_one_helper(src, self.neighbors is None)
			self.helpers[src] = helper
			if self.neighbors is not None:
				self.neighbors[src] = neighbors
		
		s1 = helper.total_sum()
		if s1 <= 0.0:
			return None
		for _ in range(self.max_retries):
			pr1 = self.rng.random() * helper.total_sum()
			ix2 = helper.upper_bound(pr1)
			id2 = neighbors[ix2] if neighbors is not None else self.global_ids[ix2]
			c2 = self.g[id2]
			if c2.K < 1.0:
				continue # very low carrying capacity cells are skipped
			if c2.st != cell.state.EMPTY:
				p1 = self.rng.random()
				if p1 >= tmp1:
					continue
			elif N3 > 0 and c2.K < N3:
				continue
			return id2
		return None
	
	def create_group_orders(self):
		self.orders.clear()
		self.raiders_revert.clear()
		
		for pt in self.g:
			c = self.g[pt]
			if c.st == cell.state.EMPTY:
				continue
			N3 = 0 # split-off faction size (in case of farming cell)
			
			if c.st == cell.state.FARMING:
				pc = 0.0 # probability of split-off due to approaching carrying capacity
				if self.cap_mig_exp >= 0.0:
					tmp = 1.0
					if self.cap_mig_exp > 0.0:
						tmp = c.N / c.K
						if tmp >= self.cap_mig_min and self.cap_mig_min < 1.0:
							tmp = (tmp - self.cap_mig_min) / (1.0 - self.cap_mig_min)
						else:
							tmp = 0.0
						tmp = pow(tmp, self.cap_mig_exp)
					pc = self.cap_mig_max * tmp
				
				pd = 0.0 # probability of split-off due to the Dunbar limit
				if self.dunbar_limit > 0.0 and self.dunbar_mig_exp >= 0.0:
					tmp = 1.0
					if self.dunbar_mig_exp > 0.0:
						tmp = c.N / self.dunbar_limit
						if tmp >= self.dunbar_mig_min and self.dunbar_mig_min < 1.0:
							tmp = (tmp - self.dunbar_mig_min) / (1.0 - self.dunbar_mig_min)
						else:
							tmp = 0.0
						tmp = pow(tmp, self.dunbar_mig_exp)
					pd = self.dunbar_mig_base * tmp
				
				if pc == 0.0 and pd == 0.0:
					continue
				
				x1 = self.rng.random()
				if x1 >= pc:
					x1 = self.rng.random()
					if x1 >= pd:
						continue
				
				# choose a faction size
				N1 = c.N
				# very unlikely, but possible if previous migrations resulted in everyone leaving
				if N1 == 0:
					continue
				N2 = N1 / 2.0
				N3 = self.rng.poisson(N2)
				if N3 > N1:
					N3 = N1
				if N3 == 0:
					continue
				
				# choose a target from the long-range distribution
				dst = self.choose_dst(pt, N3)
				if dst is not None:
					self.orders.append((pt, dst, N3))
			else:
				revert = False
				if self.raiders_revert_rate > 0.0:
					x = self.rng.random()
					if x < self.raiders_revert_rate:
						revert = True
				
				if not revert and self.pattack < 1.0:
					x = self.rng.random()
					if x >= self.pattack:
						continue # this raider cell is skipped for this year
					
				
				# here, this is a raider cell, so we have to choose from a "local" distribution for the next target
				targets = list()
				w = list()
				
				if not revert:
					max_dist = self.raider_dist * 10.0
					sumw = 0.0
					for pt2, d1 in self.g.neighbors(pt, max_dist):
						c2 = self.g[pt2]
						if not self.new_stationary_model and c2.st != cell.state.FARMING:
							continue # only farming cells can be targets
						tmp1 = c2.N if self.raiders_dist_real_pop else c2.K
						p1 = np.exp(-1.0 * d1 / self.raider_dist) * tmp1
						
						if p1 > 0.0:
							targets.append(pt2)
							w.append(p1)
							sumw += p1
				
				if len(w) > 0:
					for i in range(len(w)):
						w[i] /= sumw
					dst = self.rng.choice(targets, p = w)
					self.orders.append((pt, dst, 0))
				else:
					self.raiders_revert.append(pt)
	
	def apply_group_mig(self, nsteps):
		already_targeted = set()
		
		attack_success = 0
		attack_fail = 0
		raiders_created = 0
		raiders_reverted = 0
		
		# process reverted raiders
		for pt1 in self.raiders_revert:
			c = self.g[pt1]
			if self.raiders_can_revert:
				c.st = cell.state.FARMING
			else:
				c.st = cell.state.EMPTY
			raiders_reverted += 1
		
		# shuffle the orders randomly
		self.rng.shuffle(self.orders)
		
		# process migration events
		for src, dst, N3 in self.orders:
			if src in already_targeted:
				# source cell was attacked already, simply skip this order
				continue
			
			c1 = self.g[src]
			c2 = self.g[dst]
			if dst in already_targeted or (c2.st == cell.state.EMPTY and c1.st == cell.state.RAIDERS):
				# target was already attacked, cancel this order
				# other, unlikely case: target was a farmer cell, but everyone migrated away
				if c1.st == cell.state.RAIDERS:
					raiders_reverted += 1
					if self.raiders_can_revert:
						c1.st = cell.state.FARMING
					else:
						c1.st = cell.state.EMPTY
						c1.N = 0 # raiders are assumed to starve in this case
				continue
			
			already_targeted.add(dst)
			
			
			
			if c1.st == cell.state.RAIDERS:
				# source cell is raiders
				x1 = self.rng.random()
				if x1 < self.attack_success_prob:
					# successful attack
					attack_success += 1
					if self.mobile_raiders:
						# cell is abandoned, raiders move to the new cell
						c2.st = cell.state.RAIDERS
						c2.N = c1.N
						c1.st = cell.state.EMPTY
						c1.N = 0
					else:
						# target cell is depopulated
						c2.N = int(round(self.survivor_ratio * c2.N))
						if c2.N == 0:
							c2.st = cell.state.EMPTY
				else:
					# unsuccessful attack
					attack_fail += 1
					if not self.new_stationary_model:
						# raiders revert to farming or starve
						if self.raiders_can_revert:
							c1.N = max(1, round(self.survivor_ratio * c1.N))
							c1.st = cell.state.FARMING
						else:
							c1.N = 0
							c1.st = cell.state.EMPTY
						raiders_reverted += 1
					if self.defenders_create_raiders:
						x1 = self.rng.random()
						if x1 < self.praiders:
							c2.st = cell.state.RAIDERS
							raiders_created += 1
			else:
				# split-off from a farmer cell
				if c2.st == cell.state.EMPTY:
					c1.N -= N3
					c2.N = N3
					c2.st = cell.state.FARMING
					if self.target_only_free:
						# if only empty cells can be targeted, we need to update the migration probability distributions
						self.update_helper_neighbors(dst)
				else:
					# non-empty cell, in this case, there is a conflict
					x1 = self.rng.random()
					if x1 < self.attack_success_prob:
						# attackers succeed
						c1.N -= N3
						c2.N = N3
						x1 = self.rng.random()
						if x1 < self.praiders:
							c2.st = cell.state.RAIDERS
							raiders_created += 1
							sys.stderr.write('Raiders created: {} -> {} (step: {})\n'.format(src, dst, nsteps))
						else:
							c2.st = cell.state.FARMING
						attack_success += 1
					else:
						# attack failed
						c1.N -= N3
						c1.N += round(N3 * self.survivor_ratio)
						attack_fail += 1
						if self.defenders_create_raiders:
							x1 = self.rng.random()
							if x1 < self.praiders:
								c2.st = cell.state.RAIDERS
								raiders_created += 1
				if c1.N == 0:
					# if everyone migrated away (very unlikely), mark the cell as empty
					c1.st = cell.state.EMPTY
					if self.target_only_free:
						self.update_helper_neighbors(src)
		
		return (attack_success, attack_fail, raiders_created, raiders_reverted)
		
	def step(self, t):
		"""
		Step the demographic simulation in each cell.
		"""
		for x in self.g:
			c1 = self.g[x]
			if c1.st != cell.state.EMPTY:
				if c1.N <= 0:
					raise BaseException('Empty cell found: {}, {}, {}, {}!\n'.format(x,
						'farmers' if c1.st == cell.state.FARMING else 'raiders', c1.N, c1.K))
				c1.step(t, self.rng, self.pars)
				if self.target_only_free and c1.st == cell.state.EMPTY:
					# if this cell became empty, we need to update migration probabilities
					self.update_helper_neighbors(x)
	
	
	def create_global_ids(self):
		self.global_ids = sorted(list(self.g))
	
	
	# binary file header
	HEADER = 0xae50f3f45e30a440
	# supported flags
	BIN_FILE_CDF = 0x1 # CDF values are stored directly
	BIN_FILE_NONBR = 0x2 # true if global neighbors are not duplicated in the bin file
	BIN_FILE_USE_PEMPTY = 0x4 # not supported (true if migration probabilities are multiplied by pempty in the helper)
	BIN_FILE_LOCAL_NEIGHBORS = 0x8 # true if global neighbors are NOT used (each helper stores a list of neighbors separately)
	BIN_FILE_HAS_HELPERS = 0x10 # true if the binary file has migration probability distributions (this is the default)
	
	HEADER_FLAGS = 0x1a # combination of all flags
		
	def read(self, fn, no_helpers = False):
		"""
		Load the saved data about cells and the probability matrix from the given file.
		"""
		if self.group_mig_dist is not None or self.group_mig_pow is not None:
			raise BaseException('neolithic.read(): setting the migration distance parameters from the command line is not supported!\n')
		
		self._close_file()
		self.global_ids.clear()
		self.orig_K.clear()
		
		f1 = os.open(fn, os.O_RDONLY) # note: we might need to add os.O_CLOEXEC here as well?
		if f1 < 0:
			raise BaseException('neolithic.read(): error opening the input file!\n')
		self.matrix_f = f1
		self.mp1 = mmap.mmap(f1, 0, access = mmap.ACCESS_READ)
		if len(self.mp1) < 8:
			raise BaseException('neolithic.read(): input file too short!\n')
		
		tmp1 = struct.unpack_from('=Q', self.mp1, 0)
		if tmp1[0] & (~self.HEADER_FLAGS) != self.HEADER:
			raise BaseException('cellnet.read(): Unexpected file header!\n')
		use_cdf = bool(tmp1[0] & self.BIN_FILE_CDF)
		nonbr = bool(tmp1[0] & self.BIN_FILE_NONBR)
		if bool(tmp1[0] & self.BIN_FILE_USE_PEMPTY):
			raise BaseException('cellnet.read(): Unsupported options!\n')
		use_local_neighbors = bool(tmp1[0] & self.BIN_FILE_LOCAL_NEIGHBORS)
		
		have_helpers = False
		if not no_helpers:
			have_helpers = bool(tmp1[0] & self.BIN_FILE_HAS_HELPERS)
		
		if have_helpers:
			if use_local_neighbors:
				self.neighbors = dict()
			self.helpers_cdf = use_cdf
		
		# read the cell network
		off1 = 8 + self.g.read(self.mp1, 8, len(self.mp1))
		if off1 % 8:
			off1 += (8 - off1 % 8)
		
		# check the file size
		if off1 + 20 > len(self.mp1):
			raise BaseException('neolithic.read(): input file too short!\n')
		
		# parameters: group_mig_dist, group_mig_pow
		tmp1 = struct.unpack_from('=d', self.mp1, off1)
		self.group_mig_dist = tmp1[0]
		off1 += 8
		tmp1 = struct.unpack_from('=d', self.mp1, off1)
		self.group_mig_pow = tmp1[0]
		off1 += 8
		
		nnodes1 = len(self.g)
		self.create_global_ids()
		
		if not nonbr:
			# node IDs in the global neighbors list
			helper1 = struct.Struct('=I')
			tmp1 = helper1.unpack_from(self.mp1, off1)
			N1 = tmp1[0]
			if N1 != nnodes1:
				raise BaseException('neolithic.read(): unexpected number of nodes!\n')
			off1 += 4
			# check the remaining size
			tmp1 = off1 + nnodes1 * 12
			if tmp1 % 8:
				tmp1 += (8 - tmp1 % 8)
			tmp2 = nnodes1 if use_cdf else int(nnodes1 / 2)
			tmp1 += tmp2 * nnodes1 * 8
			if tmp1 > len(self.mp1):
				raise BaseException('neolithic.read(): input file too short!\n')
			# actually read the node IDs and compare them with the ones from the cellnet
			for i in range(nnodes1):
				tmp1 = helper1.unpack_from(self.mp1, off1)
				if self.global_ids[i] != tmp1[0]:
					raise BaseException('neolithic.read(): global_neighbors array does not match!\n')
				off1 += 4
			if off1 % 8:
				off1 += (8 - off1 % 8)
		else:
			# just check that the size is OK
			tmp1 = off1 + nnodes1 * 8
			tmp2 = nnodes1 if use_cdf else int(nnodes1 / 2)
			tmp1 += tmp2 * nnodes1 * 8
			if tmp1 > len(self.mp1):
				raise BaseException('neolithic.read(): input file too short!\n')
			
			
		# carrying capacities
		helper2 = struct.Struct('=d')
		for i in range(nnodes1):
			tmp1 = helper2.unpack_from(self.mp1, off1)
			self.orig_K[self.global_ids[i]] = tmp1[0]
			off1 += 8
		
		if have_helpers:
			self.helpers = dict()
			# create the probability helpers
			helper3 = struct.Struct('=Q')
			for i in range(nnodes1):
				src = self.global_ids[i]
				size2 = nnodes1
				if use_local_neighbors:
					size2 = helper3.unpack_from(self.mp1, off1)
					off1 += 8
					arr = mmap_array()
					arr.set_array(self.mp1, off1, size2)
					self.neighbors[src] = arr
					off1 += 4 * size2
					if size2 % 2 == 1:
						off1 += 4 # padding so that the following distribution starts at an 8-byte boundary
					
				if use_cdf:
					tmp = btp.btprob_cdf()
					tmp.set_cdf(self.mp1, off1, size2)
					self.helpers[src] = tmp
					off1 += 8 * size2
				else:
					neighbors = arr if use_local_neighbors else self.global_ids
					tmp = btp.btprob_tree(lambda ix, src=src, neighbors=neighbors : self.group_mig_prob(src, neighbors[ix]))
					tmp.set_external(self.mp1, off1, size2)
					self.helpers[src] = tmp
					off1 += 4 * size2
				
				
	def write(self, fn):
		"""
		Write the cell network and the migration distributions as a binary file.
		"""
		with open(fn, 'wb') as f:
			header1 = self.HEADER | self.BIN_FILE_NONBR | self.BIN_FILE_HAS_HELPERS
			if self.helpers_cdf:
				header1 |= self.BIN_FILE_CDF
			use_local_neighbors = False
			if self.neighbors is not None:
				header1 |= self.BIN_FILE_LOCAL_NEIGHBORS
				use_local_neighbors = True
			f.write(struct.pack('=Q', header1))
			
			tmp = self.g.write(f)
			if tmp % 8:
				pad = 8 - (tmp % 8)
				for i in range(pad):
					f.write(struct.pack('x'))
			
			# parameters: group_mig_dist, group_mig_pow
			f.write(struct.pack('=d', self.group_mig_dist))
			f.write(struct.pack('=d', self.group_mig_pow or 0.0)) # note: group_mig_pow can be None, in this case we write 0.0
			
			# note: we do not store the neighbor list
			# we write the carrying capacities, in the same order as we have the neighbors
			helper1 = struct.Struct('=d')
			if len(self.global_ids) == 0:
				self.create_global_ids()
			for x in self.global_ids:
				f.write(helper1.pack(self.g[x].K))
			
			# write the probability helper distributions
			size1 = len(self.global_ids)
			helper2 = struct.Struct('=Q')
			helper3 = struct.Struct('=I')
			for x in self.global_ids:
				helper = None
				size2 = size1
				nbr = None
				
				if x in self.helpers:
					helper = self.helpers[x]
					if use_local_neighbors:
						nbr = self.neighbors[x]
				else:
					(helper, nbr) = self.create_one_helper(x, not use_local_neighbors)
				
				if use_local_neighbors:
					size2 = len(nbr)
					f.write(helper2.pack(size2))
					for i in range(size2):
						f.write(helper3.pack(nbr[i]))
					if size2 % 2 == 1:
						f.write(helper3.pack(0)) # padding to 8 byte boundary
				
				if self.helpers_cdf:
					size2 //= 2
				
				tmp = helper.get_array()
				for i in range(size2):
					f.write(helper1.pack(tmp[i]))
						
			
class weather_update:
	"""
	Helper class for updating carrying capacities during the simulation
	run based on climatic variability.
	Weather data is given in spatial units ("weather cells") that can be
	different from the cells used in the simulation. This is handled by
	reading and storing a "dictionary" that establishes the matching
	among these. This matching need not be 1-1: the scaling factor for a
	(simulation) cell can be the weighted average of the current climatic
	factor in multiple weather cells.
	"""
	
	def __init__(self):
		self.wmap = dict() # matching from weather cells to a list of simulation cells
		self.wcells = set() # simulation cells that have an associated weather cell
		
		# main parameters -- these should be set before starting the simulation
		self.weather_allow_missing = False # allow missing data in some years (last year's value is used)
		self.weather_cut = False # cut study area to regions which have weather data (not compatible with the previous option)
		self.weather_factors = False # true if match to weather cells is not one-to-one (a factor is read together with the matching)
		self.first_year = None # "start year" of the simulation -- if given, ignore years before this
		self.sd_factor = 1.0 # scale the standard deviation of weather variation by this factor -- it is assumed that factors have a mean of one
		
		# file with the weather data
		self.wf = None
		self.wf_last_year = None # current year in the above file
		self.last_line = None # line read last
		
	def read_mapping(self, fn, is_csv, nn):
		"""
		Read a mapping between weather cells and simulation cells to be
		used from the given file (fn).
		The parameter nn should be a class of type neolithic (above) and
		it is used to determine which simulation cells exist.
		"""
		
		self.wmap = dict()
		self.wcells = set()
		
		with open(fn, 'r') as f:
			if is_csv:
				f.readline() # skip header
			for l in f.readlines():
				factor = 1.0
				l2 = l.split(',' if is_csv else None)
				if len(l2) < (3 if self.weather_factors else 2):
					raise BaseException('weather_update.read_mappings(): Invalid data in input file!\n')
				cellid = int(l2[0])
				wid = int(l2[1])
				if self.weather_factors:
					factor = float(l2[2])
				# note: only save mappings for cells that exist in the simulation space
				if wid not in self.wmap:
					# note: ensure that the climate cell exists even if it has no associated cells in the simulation
					self.wmap[wid] = list()
				
				if cellid in nn.g:
					self.wmap[wid].append((cellid, factor))
					self.wcells.add(cellid)
		
		# save the original carrying capacity values to be used later
		# (only if this was not created previously)
		if len(nn.orig_K) == 0:
			for pt in nn.g:
				nn.orig_K[pt] = nn.g[pt].K
		
	
	def prepare(self, nn):
		"""
		Function called once per simulation step to clear the carrying
		capacities before reading the new values
		"""
		if self.weather_cut:
			# initialize carrying capacity to zero everywhere, so that
			# only cells with valid climate data will be used in this step
			for pt in nn.g:
				nn.g[pt].K = 0.0
		elif self.weather_factors:
			# only zero out cells which has associated weather data
			# (this is necessary in this case, since the results will be
			# the sum from multiple factors)
			for pt in self.wcells:
				nn.g[pt].K = 0.0
	
	def update_one(self, nn, wid, factor):
		"""
		Update the carrying capacity value in one cell, based on one line
		in the input data
		"""
		
		if wid not in self.wmap:
			return
		
		# further scale the variations
		tmp = factor - 1.0
		factor = max(1.0 + self.sd_factor * tmp, 0.0)
		
		for x in self.wmap[wid]:
			cellid = x[0]
			if self.weather_factors:
				nn.g[cellid].K += nn.orig_K[cellid] * factor * x[1]
			else:
				nn.g[cellid].K = nn.orig_K[cellid] * factor
	
	
	def step_one_year(self, nn, i):
		"""
		Update carrying capacities in all cells by stepping ahead one
		year in the weather data.
		The parameter i is the current simulation step number; this
		should be advanced by 1 each time this function is called.
		"""
		if self.wf is None:
			return # allow running as no-op if no file was opened
		
		self.prepare(nn)
		
		first_line = True
		entries_read = 0
		# note: the "current" line is already read in last_line, except
		# for the case when calling this function for the first time or
		# when the file has already ended
		if self.last_line is None:
			self.last_line = self.wf.readline()
			if self.last_line == "":
				self.last_line = None
				raise BaseException('weather_update.step_one_year(): no data in input file!\n')
		while True:
			# parse one line 
			l2 = self.last_line.split()
			if len(l2) < 3:
				raise BaseException('weather_update.step_one_year(): invalid data in input file!\n')
			year = int(l2[0])
			wid = int(l2[1])
			factor = float(l2[2])
			
			skip = False
			if first_line: # first line read in the current step
				if i == 0:
					if self.first_year and year < self.first_year: # if this is the start year
						skip = True # skip ahead until the first year of the simulation is reached
				elif self.wf_last_year + 1 != year:
					raise BaseException('weather_update.step_one_year(): Non-consecutive years in weather data!\n')
				self.wf_last_year = year;
				first_line = False;
			elif year != self.wf_last_year:
				break
			
			if not skip:
				# update the carrying capacities based on the values read from the current line
				self.update_one(nn, wid, factor);
				entries_read += 1
			
			# read the next line
			self.last_line = self.wf.readline()
			if self.last_line == "":
				# end of file
				self.last_line = None
				break
		
		if not self.weather_allow_missing and entries_read != len(self.wmap):
			raise BaseException('weather_update.step_one_year(): missing entries in the climate variability data!\n')
	
		
		
