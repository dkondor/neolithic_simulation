#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  btprob.py -- store a discrete probability distribution in a data
#	structure that potentially allows efficient change of individual
#	members and also storage on disk in a binary format
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


import os, math, struct, bisect, mmap

class btprob_base:
	"""
		Base class providing the interface for storing a probability
		distributions. This class should not be used directly; instead,
		use one of the derived classes that provide actual implementations.
	"""
	
	def __init__(self):
		self._N = 0 # size of the distribution (number of nodes)
	
	def _get_sum(self, ix):
		raise NotImplementedError()
	
	def total_sum(self):
		raise NotImplementedError()
	
	def upper_bound(self, p1):
		"""
			Find the largest index such that the sum of probabilities up
			to it is below p1.
		"""
		raise NotImplementedError()


class btprob_cdf(btprob_base):
	"""
		Version that stores a static probability distribution that cannot be changed.
		Supports using an external read-only memory area as the prepared CDF
		(e.g. opened with mmap).
	"""
	def __init__(self):
		btprob_base.__init__(self)
		self._cdf = []
		self._mp = None
		self._N = 0
		self._mp_off = 0
		self._mps = struct.Struct('=d')
	
	def _get_sum(self, ix):
		if self._mp is not None:
			off1 = self._mp_off + 8 * ix
			tmp1 = self.mp1s.unpack_from(self._mp, off1)
			return tmp1[0]
		else:
			return self._cdf[ix] # this will raise an exception if the index is out of range
	
	def total_sum(self):
		if self._N > 0:
			return self._get_sum(self._N - 1)
		else:
			raise IndexError("No probabilities stored!\n")
	
	def upper_bound(self, p1):
		if self._mp is not None:
			return bisect.bisect(self, p1, 0, self._N, key = self._get_sum)
		else:
			return bisect.bisect(self._cdf, p1)
	
	def set_probs(self, cb, N):
		"""
			Allocate a local array to store N probabilities, and fill it
			up by calling cb.
		"""
		self._cdf = list(cb(i) for i in range(N))
		for i in range(1, N):
			self._cdf[i] += self._cdf[i-1]
		self._N = N

	def set_cdf(self, mp, mp_off, N):
		"""
			Set stored CDFs from an external source
		"""
		self._cdf = []
		self._mp = mp
		self._mp_off = mp_off
		self._N = N
	
	def get_array(self):
		"""
			Get a reference to the underlying array storing the CDF.
			Note: the caller should not modify this. Also, this only works
			if the array was allocated in set_probs() and not if it is
			externally supplied in set_cdf().
		"""
		return self._cdf


class btprob_tree(btprob_base):
	"""
		Store probabilities in a binary tree, potentially allowing them
		to change dynamically. Storage can be allocated here or set
		externally to a memory area opened with mmap.
	"""
	
	def __init__(self, cb = None):
		btprob_base.__init__(self)
		self._get_prob_cb = cb # callback function that gets the probability for any of the indices
		self._sums = []
		self._mp = None
		self._N = 0
		self._mp_off = 0
		self._mps = struct.Struct('=d')
	
	def _get_sum(self, ix):
		if self._mp is not None:
			off1 = self._mp_off + 8 * ix
			tmp1 = self._mps.unpack_from(self._mp, off1)
			return tmp1[0]
		else:
			return self._sums[ix] # this will raise an exception if the index is out of range
	
	def _set_sum(self, ix, val):
		if self._mp is not None:
			off1 = self._mp_off + 8 * ix
			self._mps.pack_into(self._mp, off1, val) # this will throw an exception if mp is not writeable
		else:
			self._sums[ix] = val # this will raise an exception if the index is out of range
	
	
	def _get_prob(self, ix):
		return self._get_prob_cb(ix)
	
	def _get_sum_or_prob(self, ix):
		size1 = int(self._N / 2)
		return self._get_sum(ix) if ix < size1 else self._get_prob(ix)
	
	
	def total_sum(self):
		return self._get_sum(0)
	
	def upper_bound(self, p1):
		# choose a migration target based on the stored probability distribution for at ix1
		if p1 >= self.total_sum():
			raise BaseException('Invalid probability requested!\n')
		
		sum1 = 0.0
		i = 0
		last_i = 0
		psize = self._N
		
		while 2*i + 1 < psize: # while we have children
			tmp1 = self._get_sum_or_prob(2*i + 1) # sum in the left child
			tmp2 = self._get_prob(i) # value of i
			sum2 = sum1 + tmp1 + tmp2 # value of i in the sorted list of sums
			if p1 < sum2:
				# i might be a good candidate or we have to descend to the left
				if p1 >= sum1 + tmp1:
					break # i is the good candidate
				last_i = i
				i = 2*i + 1
			else:
				# in this case, we always have to go to the right
				if 2*i + 2 < psize:
					sum1 = sum2
					i = 2*i + 2
				else:
					# this probably should not happen, but is possible due to rounding errors?
					i = last_i
					break
		return i
	
	def _recalculate_sum(self, ix):
		left = 2*ix + 1
		right = 2*ix + 2
		tmp = self._get_prob(ix)
		if left < self._N:
			tmp += self._get_sum_or_prob(left)
		if right < self._N:
			tmp += self._get_sum_or_prob(right)
		self._set_sum(ix, tmp)
	
	def set_all_probs(self, N, cb = None):
		"""
			Set the callback giving the probabilities and recalculate all probabilities.
		"""
		if cb is not None:
			self._get_prob_cb = cb
		self._mp = None
		self._N = N
		size1 = int(N/2)
		self._sums = list(0 for i in range(size1))
		for i in reversed(range(0, size1)):
			self._recalculate_sum(i)
	
	def recalculate_sum(self, ix):
		"""
			Function to call (externally) after probability at ix has changed,
			this will (recursively) update the stored CDF.
		"""
		size1 = int(self._N / 2)
		if size1 == 0:
			return
		while True:
			if ix < size1:
				self._recalculate_sum(ix)
			if ix == 0:
				break
			ix = int((ix - 1) / 2)
			
	def set_external(self, mp, mp_off, N):
		self._sums = []
		self._mp = mp
		self._mp_off = mp_off
		self._N = N
	
	def get_array(self):
		"""
			Get a reference to the underlying array storing the partial sums.
			Note: the caller should not modify this. Also, this only works
			if the array was allocated in set_all_probs() and not if it is
			externally supplied in set_external().
		"""
		return self._sums

