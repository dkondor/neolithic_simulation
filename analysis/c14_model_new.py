#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 15:06:58 2022

c14_model_new.py -- process all C14 data using the new methodology:
	1. Aggregate in overlapping quasi-rectangular regions
	2. Calibrate and calculate SPDs (using a binning generated with rcarbon)
	3. Fit a logistic growth model to each of the SPDs
	4. Generate synthetic dates using the above models as the underlying
		probability distribution (process similar to rcarbon's calsample method)
	5. Uncalibrate the synthetic dates (i.e. convert their ages to radiocarbon
		"measurements")
	6. Calibrate and calculate SPDs for the synthetic datasets
	7. Calculate average SPD and 95% confidence intervals (similarly to
		modelTest())
	8. Calculate ACFs for each of the synthetic SPDs along with ACF for the
		original data after detrending with the mean synthetic SPD
	9. Calculate 95% confidence intervals for the ACFs

@author: Daniel Kondor -- kondor@csh.ac.at
"""

#%% 1. Dependencies

import os
# note: change to the directory that contains calmult.py
os.chdir('/home/dkondor/CSH/HoloSim/neolithic_simulation/analysis')

import numpy as np
import pandas as pd
import scipy.optimize as sco
import tqdm
import multiprocessing as mp
import calmult
import pickle as pkl
import gzip
import itertools


#%% 2. General parameters
tmin1 = 5000 # time range for main analysis, in calBP
tmax1 = 9000
tedge = 1000 # padding for edge effects
lmax = 3000 # maximum lag to consider for ACFs
mindates = 500 # only use regions with 500 dates (without binning)
nsamples = 100 # number of synthetic SPDs per region
norm = False # use unnormalized SPDs

# base directory for input data files and output
base_dir = '/home/dkondor/CSH/HoloSim/data'
datadir = base_dir + '/population'

#%% 3. Main functions used

# fit a logistic model
def fitlog(data):
	"""
	Fit a logistic model to the given data.

	Parameters
	----------
	data : Data to use. Should be a pandas DataFrame, with columns calBP and PrDens.

	Returns
	-------
	Dictionary with the parameters. Also adds a column ("pred") to data with
	the fitted values.
	Returns None if fitting the model failed.
	"""
	def f1(x, T1, x0, N1, N0):
		return N0 + N1 / (1.0 + np.exp((x-x0)/T1))
	
	tmin0 = data.calBP.min()
	T1 = 100 # timescale
	N = data.PrDens[data.calBP < tmin0 + 1000].mean()
	x0 = data.calBP[data.PrDens > N / 2].max()
	try:
		tmp1 = sco.curve_fit(f1, data.calBP, data.PrDens, p0 = [T1, x0, N, 0],
						  bounds = ([50, x0-2000, N/2, 0],
					 [5000, x0+2000, data.PrDens.max(), N/2]))
		
		Tf, xf, N1f, N0f = tmp1[0]
		f1v = np.vectorize(f1)
		data['pred'] = f1v(data.calBP, *tmp1[0])
		return {'N1': N1f, 'N0': N0f, 'T1': Tf, 'x0': xf}
	except:
		return None


def acf(x, base, lmax):
	"""
	Calculate the autocorrelation of a time series after subtracting a trend.
	
	Parameters
	----------
	x : Time series to consider, array-like or list.
	base : Trend to subtract before calculating the ACF.
	lmax : Maximum lag to consider

	Returns
	-------
	Autocorrelation function up to lmax.
	
	"""
	
	n = len(x)
	y = x - base
	v = y.var()
	y = y - y.mean()
	r = np.correlate(y, y, 'full')
	r = r[n-1:n-1+lmax]
	if r[0] > 0:
		return r / r[0]
	else:
		return np.zeros(lmax)
		


def get_acf_min(acf, filt = None):
	"""
	Get the location of the first minimum in the results of ACF. Note that
	the value of the ACF has to be negative and optionally less than the
	values given in filt.

	Parameters
	----------
	acf : Result of ACF to find the minimum in. Must be a numpy array.
	filt : Optionally a filter; is given, minimum locations are only
		considered in regions where acf < filt (note: not used in the results)

	Returns
	-------
	Location (0-based index) of the first minimum or None if no suitable
	minimum was found.
	"""
	
	ds = np.sign(np.diff(acf))
	if np.any(ds == 0):
		# get rid of zeros (i.e. cases where consecutive values are
		# exactly the same in the acf -- this is very unlikely to occur)
		for i in range(len(ds)-1, 0, -1):
			if ds[i-1] == 0:
				ds[i-1] = ds[i]
	min1 = np.nonzero(np.diff(ds) == 2)[0] + 1
	if len(min1) == 0:
		return None
	min1 = min1[acf[min1] < 0]
	if len(min1) == 0:
		return None
	if filt is not None:
		if type(filt) in [np.ndarray, pd.core.series.Series]:
			min1 = min1[acf[min1] < filt[min1]]
		else:
			min1 = min1[acf[min1] < filt]
		if len(min1) == 0:
			return None
	return min1[0]
	


def model_fit_all(pop2, regions1, tmin1, tmax1, tedge, mindates, norm, bins = "bin100"):
	"""
	Calculate SPDs and fit logistic model in a set of grid cells.

	Parameters
	----------
	pop2 : C14 data (TODO: format description!).
	regions1 : Grid cells to use. Should be a pandas DataFrame, where each
		row defines a region by its corners (columns "lonmin", "lonmax",
		"latmin", "latmax") and a numeric ID ("ID" column, should be unique).
	tmin1 : Lower time limit (in calBP). The SPD is cut to this time range and
		the fit is performed considering only the values in this range.
	tmax1 : Upper time limit (in calBP).
	tedge : Padding in time to deal with edge effect. Calibration is performed
		for the extended time range [tmin1 - tedge, tmax1 + tedge].
	mindates : Minimum number of dates required in a region to include in the
		analysis.
	norm : Whether to use normalized SPDs (i.e. normalize the probability
		distributions after calibration).
	bins : Name of the column in pop2 containing bins for the dates. If None,
		then binning is not used.

	Returns
	-------
	Pandas DataFrame containing the fitted parameters for each region (where
	there where enough dates and the fit was successful).

	"""
	res_all = list()
	for ID in regions1.ID:
		region2 = regions1[regions1.ID == ID].reset_index(drop=True)
		dates1 = pop2[(pop2.lon >= region2.lonmin[0]) & (pop2.lon <= region2.lonmax[0]) &
					  (pop2.lat >= region2.latmin[0]) & (pop2.lat <= region2.latmax[0])].copy()
		dates1 = dates1[(dates1.calBPmed >= tmin1 - tedge) & (dates1.calBPmed <= tmax1 + tedge)]
		if dates1.shape[0] < mindates:
			continue
		dates1.reset_index(drop = True, inplace = True)
		
		# original SPD
		spd1 = calmult.spd(dates1.c14age, dates1.c14std, trange = [tmin1 - tedge, tmax1 + tedge],
					 normalize=norm, bins = dates1[bins] if bins is not None else None)
		spd1 = spd1[(spd1.calBP >= tmin1) & (spd1.calBP < tmax1)]
		
		tmp1 = fitlog(spd1)
		if tmp1 is not None:
			tmp1['ID'] = ID
			res_all.append(tmp1)
	return pd.DataFrame(res_all)




def model_test_one(dates1, pred, tmin1, tmax1, tedge, lmax, nsamples, norm, rng0,
				   bins = "bin100", rollmean = None):
	"""
	Generate a set of synthetic dates and their SPDs for one region (and one
	set of dates), calculate their ACF and return summary statistics about
	the mean, median and 1%, 5%, 95% and 99% percentiles.

	Parameters
	----------
	dates1 : Filtered dates in the region of interest.
	pred : Fitted model to use as a basis for generating synthetic data.
		Should be a pandas DataFrame with columns "calBP" and "prob", the
		latter giving the probabilities for generating dates in each of the
		years (note: probabilities should sum to one).
	tmin1 : Lower time limit (in calBP). SPDs are cut to this time range and
		ACFs are calculated considering only the time series in this range.
	tmax1 : Upper time limit (in calBP).
	tedge : Padding in time to deal with edge effect. Synthetic dates are
		generated and calibration is performed for the extended time range
		[tmin1 - tedge, tmax1 + tedge].
	lmax : Maximum lag to consider in the ACFs.
	nsamples : Number of synthetic data sets to generate.
	norm : Whether to use normalized SPDs (i.e. normalize the probability
		distributions after calibration).
	rng0 : Random number generator (from numpy.random) to use.
	bins : Name of the column in dates1 containing bins for the dates. If None,
		then binning is not used.
	rollmean : if given (as an integer or as a list of integers), repeat all
		results with performing a moving average over it

	Raises
	------
	BaseException if internal sanity checks fail (should not happen).

	Returns
	-------
	A dictionary with two pandas DataFrames, "spdc" and "acfc", containing
	SPDs and ACFs for the real dataset and summary statistics for the
	synthetic datasets.

	"""
	
	if rollmean is None:
		rollmean = [0]
	elif type(rollmean) is int:
		if rollmean != 0:
			rollmean = [0, rollmean]
		else:
			rollmean = [0]
	elif type(rollmean) is list:
		rollmean = rollmean.copy()
		if 0 not in rollmean:
			rollmean.append(0)
	elif type(rollmean) in [np.ndarray, pd.core.series.Series]:
		rollmean = list(rollmean)
		if 0 not in rollmean:
			rollmean.append(0)
	
	# original SPD
	spd0 = calmult.spd(dates1.c14age, dates1.c14std, trange = [tmin1 - tedge, tmax1 + tedge],
				 normalize=norm, bins = dates1[bins] if bins is not None else None)
	spd0 = spd0[(spd0.calBP >= tmin1) & (spd0.calBP < tmax1)]
	PrScale = sum(pred.predlog[(pred.calBP >= tmin1) & (pred.calBP < tmax1)])
	
	ndates = len(dates1[bins].unique()) if bins is not None else dates1.shape[0]
	
	spd20 = np.zeros([nsamples, tmax1-tmin1], dtype=float)
	
	for j in range(nsamples):
		dates2 = rng0.choice(pred.calBP, ndates, True, pred.prob)
		err0 = rng0.choice(dates1.c14std, ndates, True)
		dates2 = calmult.uncalibrate(np.int32(dates2))
		spd3 = calmult.spd(np.array(dates2.uncalBP), err0, trange = [tmin1 - tedge, tmax1 + tedge], normalize=norm)
		spd3 = spd3[(spd3.calBP >= tmin1) & (spd3.calBP < tmax1)]
		# normalize the resulting SPD to have the same total weight as in
		# the model used for sampling
		# (note: the alternative would be to normalize with the observed SPD)
		spd3['PrDens'] = spd3['PrDens'] * PrScale / sum(spd3['PrDens'])
		# sanity check that the time ranges match
		if np.any(spd0.calBP != spd3.calBP):
			raise BaseException('SPD time ranges do not match!\n')
		spd20[j,] = spd3.PrDens
	
	ret = list()
	for r in rollmean:
		t0 = 0
		if r == 0:
			spd1 = spd0.copy()
			spd2 = spd20
		else:
			ctmp = np.ones(r) / r
			t0 = r // 2
			l1 = tmax1 - tmin1
			# note: tmp1 will be an array with length l1-r+1
			tmp1 = np.convolve(spd0.PrDens, ctmp, mode='valid')
			spd1 = spd0.iloc[t0:(t0+l1-r+1)].reset_index(drop=True)
			spd1.PrDens = tmp1
			
			spd2 = np.zeros([nsamples, l1-r+1], dtype=float)
			for j in range(nsamples):
				spd2[j,] = np.convolve(spd20[j,], ctmp, mode='valid')
		
		spd1['mavg'] = np.mean(spd2, 0)
		spd1['mmed'] = np.median(spd2, 0)
		spd1['m1']  = np.quantile(spd2, 0.005, 0)
		spd1['m5']  = np.quantile(spd2, 0.025, 0)
		spd1['m95'] = np.quantile(spd2, 0.975, 0)
		spd1['m99'] = np.quantile(spd2, 0.995, 0)
		spd1['sd'] = np.std(spd2, 0)
		spd1 = pd.merge(spd1, pred, on = 'calBP')
		
		tmp1 = {'rm': r}
		
		
		# calculate the ACF of the detrended time series
		acf1 = acf(spd1.PrDens, spd1.mavg, lmax)
		acf1 = pd.DataFrame({'lag': np.arange(lmax), 'c14': acf1})
		
		tmp1['acfc'] = acf1
	
		# calculate a p-value for the SPD
		# for this, we need to convert SPD values to z-scores
		zsim = (spd2 - spd1.mavg.to_numpy()) / spd1.sd.to_numpy()
		zlo = np.quantile(zsim, 0.025, 0)
		zhi = np.quantile(zsim, 0.975, 0)
		
		tmp21 = zsim - zlo
		tmp22 = zsim - zhi
		tmp21[tmp21 > 0] = 0
		tmp22[tmp22 < 0] = 0
		expected = -1 * np.sum(tmp21, 1) + np.sum(tmp22, 1)
		
		tmp2 = (spd1['PrDens'] - spd1['mavg']) / spd1['sd']
		tmp21 = tmp2 - zlo
		tmp22 = tmp2 - zhi
		observed = sum(abs(tmp21[tmp21 < 0])) + sum(tmp22[tmp22 > 0])
		tmp1['pval'] = (sum(expected > observed) + 1) / (nsamples + 1)
		tmp1['spdc'] = spd1
		
		ret.append(tmp1)
	
	return ret
	


def model_test_all_new(fit0, pop2, regions1, tmin1, tmax1, tedge, lmax, nsamples,
					   norm, rng0, bins = "bin100", rollmean = None):
	"""
	Perform the generation of synthetic dates based on previously fitted
	logistic models for a set of regions. Note: all parameters should be the
	same that were given to model_fit_all() when generating fit0.

	Parameters
	----------
	fit0 : Result of model_fit_all() for the same set of regions.
	pop2 : C14 data.
	regions1 : Regions used for aggregating (should contain the same regions
		that were used with model_fit_all())
	tmin1 : Lower time limit (in calBP). SPDs are cut to this time range and
		ACFs are calculated considering only the time series in this range.
	tmax1 : Upper time limit (in calBP).
	tedge : Padding in time to deal with edge effect. Synthetic dates are
		generated and calibration is performed for the extended time range
		[tmin1 - tedge, tmax1 + tedge].
	lmax : Maximum lag to consider in the ACFs.
	nsamples : Number of synthetic data sets to generate.
	norm : Whether to use normalized SPDs (i.e. normalize the probability
		distributions after calibration).
	rng0 : Random number generator (from numpy.random) to use.
	bins : Name of the column in dates1 containing bins for the dates. If None,
		then binning is not used.
	rollmean : if given (as an integer or as a list of integers), repeat all
		results with performing a moving average over it


	Returns
	-------
	List of list of dicts containing the region ID, rolling mean window and
	pandas DataFrames of SPDs and ACFs of the real dataset and summaries for
	the synthetic dates.
	
	"""
	res_all = list()
	for i in range(fit0.shape[0]):
		ID = fit0.ID[i]
		region2 = regions1[regions1.ID == ID].reset_index(drop=True)
		dates1 = pop2[(pop2.lon >= region2.lonmin[0]) & (pop2.lon <= region2.lonmax[0]) &
					  (pop2.lat >= region2.latmin[0]) & (pop2.lat <= region2.latmax[0])].copy()
		dates1 = dates1[(dates1.calBPmed >= tmin1 - tedge) & (dates1.calBPmed <= tmax1 + tedge)]
		dates1.reset_index(drop = True, inplace = True)
		
		pred = pd.DataFrame({'calBP': np.arange(tmin1-tedge, tmax1+tedge)})
		pred['predlog'] = fit0.N0[i] + fit0.N1[i] / (1 + np.exp((pred.calBP - fit0.x0[i]) / fit0.T1[i]))
		pred['prob'] = pred.predlog/ sum(pred.predlog)
		
		tmp1 = model_test_one(dates1, pred, tmin1, tmax1, tedge, lmax, nsamples,
						norm, rng0, bins, rollmean)
		for x in tmp1:
			x['ID'] = ID
		
		res_all.append(tmp1)
	return res_all
		

def do_one_run2(pars):
	"""
	Do one full analysis for one set of regions (selected by dx and dy) and
	save the results in a pickle file. Suitable to be called in a
	multiprocessing.Pool.

	Parameters
	----------
	pars : Tuple of parameters, containing the following members (in order):
		
		pop2 : C14 data.
		regions_all : Set of overlapping regions to use. The actual subset is
			selected based on the dx and dy parameters below.
		dx : Offset in the X (longitude) direction when selecting the grid to use
		dy : Offset in the Y (latitude) direction when selecting the grid to use
		tmin1 : Lower time limit (in calBP). SPDs are cut to this time range and
			ACFs are calculated considering only the time series in this range.
		tmax1 : Upper time limit (in calBP).
		tedge : Padding in time to deal with edge effect. Synthetic dates are
			generated and calibration is performed for the extended time range
			[tmin1 - tedge, tmax1 + tedge].
		mindates : Minimum number of dates needed in a region to consider for
			the analysis.
		lmax : Maximum lag to consider in the ACFs.
		nsamples : Number of synthetic data sets to generate.
		norm : Whether to use normalized SPDs (i.e. normalize the probability
			distributions after calibration).
		seed : Random seed to use (a new random number generator is initialized with this).
		bins : Name of the column in dates1 containing bins for the dates. If None,
			then binning is not used.
		out_base : Base filename for output (dx and dy parameters are appended).
		rollmean : if given (as an integer or as a list of integers), repeat all
			results with performing a moving average over it
	Returns
	-------
	None.

	"""
	pop2, regions_all, dx, dy, tmin1, tmax1, tedge, mindates, lmax, nsamples, norm, seed, bins, out_base, rollmean = pars
	
	rng1 = np.random.default_rng(seed)
	regions1 = regions_all[(regions_all.dx == dx) & (regions_all.dy == dy)]
	fit0 = model_fit_all(pop2, regions1, tmin1, tmax1, tedge, mindates, norm, bins)
	res = model_test_all_new(fit0, pop2, regions1, tmin1, tmax1, tedge, lmax,
						  nsamples, norm, rng1, bins, rollmean)
	fn1 = '{}_dx{}_dy{}.pkl'.format(out_base, dx, dy)
	with gzip.open(fn1, 'wb') as f:
		pkl.dump({'res': res, 'fit': fit0}, f)




#%% 4. Read the data
aggr_regions = pd.read_csv('{}/new_grids_all.csv'.format(base_dir))
regions_all = aggr_regions[aggr_regions.res == 500] # only use the grid linear size of 500 km
pop2 = pd.read_csv('{}/c14_data_bins1.csv'.format(datadir))
rng = np.random.default_rng()


##############################################################################
#%% 5. run the main analysis
outdir1 = 'new_analysis'
nsamples = 1000
pool = mp.Pool(5)


fnbase = '{}/{}/res'.format(datadir, outdir1)
all_pars = list(itertools.product(range(5), range(5)))

for _ in tqdm.tqdm(pool.imap_unordered(do_one_run2, ((pop2, regions_all, x, y,
				tmin1, tmax1, tedge, mindates, lmax, nsamples, norm,
				rng.integers(2**31), 'bin100', fnbase, [0])
			for x, y in all_pars)), total=len(all_pars)):
	pass

pool.close()
del pool
# runtime: ~9 hours


#%% save the results in CSV for use with R
acfm1 = list()
pval = list()
for dx in tqdm.tqdm(range(5)):
	for dy in tqdm.tqdm(range(5)):
		with gzip.open('{}_dx{}_dy{}.pkl'.format(fnbase, dx, dy), 'rb') as f:
			tmp1 = pkl.load(f)
			res1 = tmp1['res']
			for y in res1:
				for x in y:
					rm = x['rm']
					ID = x['ID']
					acfc = x['acfc']
					m1 = get_acf_min(acfc.c14)
					if m1 is not None:
						acfm1.append({'dx': dx, 'dy': dy, 'ID': ID, 'm': m1, 'rm': rm})
					pval.append({'dx': dx, 'dy': dy, 'ID': ID, 'rm': rm, 'pval': x['pval']})

# write out acf minimum locations
acfm12 = pd.DataFrame(acfm1)
acfm12.to_csv('{}_min.csv'.format(fnbase), index=False)

# write out p-values
pval = pd.DataFrame(pval)
pval.to_csv('{}_pval.csv'.format(fnbase), index=False)

# write out example results
dx = 0
dy = 0

tmp1 = None
with gzip.open('{}_dx{}_dy{}.pkl'.format(fnbase, dx, dy), 'rb') as f:
	tmp1 = pkl.load(f)
	model_res1 = tmp1['res']
	acf_all = pd.concat(pd.concat([model_res1[i][j]['acfc'], pd.DataFrame(
			{'ID': np.repeat(model_res1[i][j]['ID'], model_res1[i][j]['acfc'].shape[0]),
			 'rm': np.repeat(model_res1[i][j]['rm'], model_res1[i][j]['acfc'].shape[0])}
			) ], axis=1) for i in range(len(model_res1)) for j in range(len(model_res1[i])) )
	acf_all.to_csv('{}_acf_zero_all_dx{}_dy{}.csv'.format(fnbase, dx, dy), index=False)
	spd_all = pd.concat(pd.concat([model_res1[i][j]['spdc'], pd.DataFrame(
			{'ID': np.repeat(model_res1[i][j]['ID'], model_res1[i][j]['spdc'].shape[0]),
			 'rm': np.repeat(model_res1[i][j]['rm'], model_res1[i][j]['spdc'].shape[0])}
			) ], axis=1) for i in range(len(model_res1)) for j in range(len(model_res1[i])) )
	spd_all.to_csv('{}_spd_all_dx{}_dy{}.csv'.format(fnbase, dx, dy), index=False)


#%% 
# calculate coefficient of variation values after detrending
cv1 = list()
for dx in tqdm.tqdm(range(5)):
	for dy in tqdm.tqdm(range(5)):
		with gzip.open('{}_dx{}_dy{}.pkl'.format(fnbase, dx, dy), 'rb') as f:
			tmp1 = pkl.load(f)
			res1 = tmp1['res']
			for y in res1:
				for x in y:
					rm = x['rm']
					ID = x['ID']
					spdc = x['spdc']
					tmp2 = spdc.PrDens - spdc.mavg
					cv = tmp2.std() / spdc.mavg.mean()
					cv1.append({'dx': dx, 'dy': dy, 'ID': ID, 'rm': rm, 'cv': cv})

cv1 = pd.DataFrame(cv1)
cv1.to_csv('{}_cv.csv'.format(fnbase), index=False)
	

