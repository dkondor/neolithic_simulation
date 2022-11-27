#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
calmult.py -- simple method for calibrating C14 dates and calculate SPDs
	without binning, optimized for speed
	main functions separated here

Created on Tue Aug  2 17:34:12 2022

@author: Daniel Kondor -- kondor@csh.ac.at
"""

import numpy as np
import pandas as pd

# base directory for input data files and output
base_dir = '/home/dkondor/CSH/HoloSim/data'

# read the calibration curve -- TODO: use a flexible location for this!
calcurve = pd.read_csv(base_dir + '/population/intcal20.14c', comment='#')

# range of times (in calBP) that we will use -- TODO: make this a parameter!
crange = [3000, 11000]

calcurve = calcurve[calcurve.CALBP >= crange[0]]
calcurve = calcurve[calcurve.CALBP <= crange[1]]
calcurve.reset_index(inplace=True)

calBPrange = np.arange(*crange)
mu = np.zeros_like(calBPrange)
tau2 = np.zeros_like(calBPrange)

j = calcurve.shape[0] - 1
for i in range(calBPrange.shape[0]):
	while calcurve.CALBP[j-1] < calBPrange[i]:
		j -= 1
	x = (calBPrange[i] - calcurve.CALBP[j]) / (calcurve.CALBP[j-1] - calcurve.CALBP[j])
	mu[i] = calcurve.C14BP[j] + x * (calcurve.C14BP[j-1] - calcurve.C14BP[j])
	tau2[i] = np.power((calcurve.Error[j] + x * (calcurve.Error[j] - calcurve.Error[j-1])), 2)


def _check_range(trange):
	if trange[0] > trange[1]:
		tmp2 = trange[0]
		trange[0] = trange[1]
		trange[1] = tmp2
	if trange[0] < crange[0] or trange[1] > crange[1]:
		raise BaseException('Unsupported date range!\n')
	return trange


# calibrate a set of dates together
def calmult(x, sd, trange = None, normalize = False, eps = 1e-5):
	if trange is not None:
		trange = _check_range(trange)
		mu1 = mu[(calBPrange >= trange[0]) & (calBPrange < trange[1])]
		tau21 = tau2[(calBPrange >= trange[0]) & (calBPrange < trange[1])]
	else:
		mu1 = mu
		tau21 = tau2
	# note: numpy's outer operations work with lists and np.arrays, but not with
	# pandas dataframe columns, these we need to case explicitly to np.arrays
	tmp1 = np.subtract.outer(x.to_numpy() if type(x) is pd.core.series.Series else x, mu1)
	if ((type(sd) is np.ndarray) or (type(sd) is list) or type(sd) is pd.core.series.Series) and len(sd) > 1:
		if len(sd) != len(x):
			raise BaseException('Dates and errors have to be the same length!\n')
		if type(sd) is pd.core.series.Series:
			tmp2 = sd.to_numpy()
			tmp2 = np.add.outer(np.multiply(tmp2, tmp2), tau21)
		else:
			tmp2 = np.add.outer(np.multiply(sd, sd), tau21)
		tmp1 = -0.5 * tmp1 * tmp1 / tmp2
		ret = (1.0 / np.sqrt(2*np.pi*tmp2)) * np.exp(tmp1) 
	else:
		tmp1 = -0.5 * tmp1 * tmp1 / (tau21 + sd*sd)
		ret = (1.0 / np.sqrt(2*np.pi*(tau21 + sd*sd))) * np.exp(tmp1)
	ret[ret < eps] = 0
	if normalize:
		sum1 = np.sum(ret, 1)
		ix1 = sum1 > 0
		ret[ix1, :] = ret[ix1, :] / sum1[ix1, np.newaxis]
	return ret
	
def spd(x, sd, trange = None, normalize = False, bins = None, eps = 1e-5):
	if trange is not None:
		trange = _check_range(trange)
	tmp1 = calmult(x, sd, trange, normalize, eps)
	if bins is not None:
		if len(bins) != tmp1.shape[0]:
			raise BaseException('Unexpected bins size!\n')
		if type(bins) is not np.ndarray:
			bins = np.array(bins)
		ubins = np.unique(bins)
		tmp2 = np.zeros((len(ubins), tmp1.shape[1]))
		for i in range(len(ubins)):
			b = ubins[i]
			ix = np.where(bins == b)
			tmp2[i, :] = np.sum(tmp1[ix], 0) / len(ix[0])
		tmp1 = np.sum(tmp2, 0)
	else:
		tmp1 = np.sum(tmp1, 0)
	return pd.DataFrame({'calBP': calBPrange if trange is None else np.arange(trange[0], trange[1]),
					  'PrDens': tmp1})

def uncalibrate(x):
	if np.any(np.less(x, crange[0])) or np.any(np.greater_equal(x, crange[1])):
		raise BaseException('Given dates are out of the supported range!\n')
	ix = np.subtract(x, crange[0]) # index into the mu and tau arrays (we assume x are integers)
	return pd.DataFrame({
		'calBP': x,
		'uncalBP': mu[ix],
		'err2': tau2[ix]
		})
	
	
	

