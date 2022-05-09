#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gaez_convert.py -- load GAEZ yield estimate, do a pre-processing before
matching with the DGGRID cell data
"""

import numpy as np
import pandas as pd
import tqdm
import itertools
import geopandas as gpd
import shapely as shp
import zipfile


# base directory for data files -- need to adjust this to the correct paths!
base_dir = '/home/dkondor/CSH/HoloSim/data'
# data directory with data downloaded from GAEZ (zip files)
data_base_dir = base_dir + '/gaez_crop'


# base filename
base1 = 'gaez_yield'
# crop types
types = ['spring_wheat_low', 'winter_wheat_low']
"""
Note: this assumes that the "potential yield" datasets were saved with the
following file names:
gaez_yield_spring_wheat_low.zip
gaez_yield_winter_wheat_low.zip
(for spring and winter wheat respectively)
"""


class gaezdata:
    _ncols = 0
    _nrows = 0
    _xll = 0
    _yll = 0
    _csize = 0
    _values = None
    
    def _key_validate(self, k, v):
        """
        Validate one line of metadata read from the input file
        """
        if k == "ncols":
            self._ncols = int(v)
        elif k == "nrows":
            self._nrows = int(v)
        elif k == "xllcorner":
            self._xll = float(v)
        elif k == "yllcorner":
            self._yll = float(v)
        elif k == "cellsize":
            self._csize = float(v) # TODO: check for rounding??
        elif k == "NODATA_value":
            pass # should we handle this separately?
        else:
            raise BaseException("gaezdata._key_validate(): invalid metadata key!\n")
    
    def load(self, fn, zip_archive = None, dtype = np.int32):
        """
        Load data in a grid from the given file.
        If zip_archive is not none, it is the name of a compressed file, and fn is the name of the data file inside it (typically data.asc).
        The parameter dtype can specify the type of data in the file; it will be used to create the data matrices.
        """
        f1 = None
        if zip_archive is not None:
            zf = zipfile.ZipFile(zip_archive, "r")
            f1 = zf.open(fn, "r")
        else:
            f1 = open(fn, "r")
        
        # first 6 lines are metadata (TODO: check if this is always the case)
        i = 0
        while i < 6:
            l1 = f1.readline()
            if zip_archive is not None:
                # note: zip files are opened in "binary" mode
                l1 = l1.decode('utf-8')
            l2 = l1.split()
            if len(l2) == 0:
                continue
            if len(l2) != 2:
                raise BaseException("gaezdata.load(): invalid metadata!\n")
            self._key_validate(*l2)
            i += 1
        
        # check if metadata is valid
        if self._ncols <= 0 or self._nrows <= 0 or self._csize <= 0:
            raise BaseException("gaezdata.load(): invalid metadata!\n")
        
        self._values = np.zeros([self._nrows, self._ncols], dtype = dtype)
        # note: data starts from the top (north) -- we fill up the array in
        # the "opposite" direction
        for i in range(self._nrows):
            l1 = f1.readline()
            if zip_archive is not None:
                l1 = l1.decode('utf-8')
            self._values[self._nrows - i - 1, ] = [dtype(x) for x in l1.split()]
    
    def get_value_xy(self, x, y):
        if x >= self._ncols or y >= self._nrows or x < 0 or y < 0:
            raise BaseException("gaezdata.get_value_xy(): position out of range!\n")
        return self._values[y, x]
    
    def lonlat_to_xy(self, lon, lat):
        x = int((lon - self._xll) / self._csize)
        y = int((lat - self._yll) / self._csize)
        return (x, y)
    
    def xy_to_lonlat(self, x, y):
        lon = self._xll + x * self._csize
        lat = self._yll + y * self._csize
        return (lon, lat)
    
    def get_value_lonlat(self, lon, lat):
        x, y = self.lonlat_to_xy(lon, lat)
        return self.get_value_xy(x, y)
    
    def get_geom_xy(self, x, y):
        shp1 = shp.geometry.Polygon([
                self.xy_to_lonlat(x, y),
                self.xy_to_lonlat(x + 1, y),
                self.xy_to_lonlat(x + 1, y + 1),
                self.xy_to_lonlat(x, y + 1),
                self.xy_to_lonlat(x, y)
                ])
        return shp1
    
    def get_geom_lonlat(self, lon, lat):
        tmp = self.lonlat_to_xy(lon, lat)
        return self.get_geom_xy(*tmp)
    
    def get_geom_with_value_lonlat(self, lon, lat):
        return gpd.GeoDataFrame({'value': [self.get_value_lonlat(lon, lat)],
                                'geometry' : [self.get_geom_lonlat(lon, lat)]},
                                crs="EPSG:4326")
        
    def get_geom_with_value_xy(self, x, y):
        return gpd.GeoDataFrame({'value': [self.get_value_xy(x, y)],
                                'geometry' : [self.get_geom_xy(x, y)]},
                                crs="EPSG:4326")
    
    def get_geoms_box_lonlat(self, lon1, lat1, lon2, lat2):
        """
        Get a GeoDataFrame with all cell geometries and values in the given
        bounding box.
        """
        x1, y1 = self.lonlat_to_xy(min(lon1, lon2), min(lat1, lat2))
        x2, y2 = self.lonlat_to_xy(max(lon1, lon2), max(lat1, lat2))
        return gpd.GeoDataFrame(({
                'value': self.get_value_xy(x, y),
                'geometry': self.get_geom_xy(x, y)
                }
            for x, y in itertools.product(range(x1, x2 + 1), range(y1, y2 + 1))
            ), crs="EPSG:4326")
        
        

##############################################################################
# 1. Load the potential yield datasets
crop_data1 = gaezdata()
crop_data2 = gaezdata()

crop_data1.load('data.asc', '{}/{}_{}.zip'.format(data_base_dir, base1, types[0]),
                dtype=np.float)
crop_data2.load('data.asc', '{}/{}_{}.zip'.format(data_base_dir, base1, types[1]),
                dtype=np.float)

eu_x = -12 # lower left corner for EU
eu_y = 36.5 

off_x = int((eu_x - crop_data1._xll) / crop_data1._csize)
off_y = int((eu_y - crop_data1._yll) / crop_data1._csize)

xsize = 732 # size of the grid used in the current study
ysize = 258

# create a matrix of yield values only in the study area
crop_data_combined1 = np.zeros([xsize, ysize, 2], dtype = np.float)
for x in range(xsize):
    crop_data_combined1[x,:,0] = crop_data1._values[off_y:(off_y + ysize), x + off_x]
    crop_data_combined1[x,:,1] = crop_data2._values[off_y:(off_y + ysize), x + off_x]
crop_data_combined = np.max(crop_data_combined1, axis=2)

crop_data_combined_df = pd.DataFrame([
    {'x' : x, 'y' : y, 'v' : crop_data_combined[x, y]} for (x, y) in
    itertools.product(range(xsize), range(ysize)) if crop_data_combined[x, y] > 0
    ])
# 130745 valid rows

# add IDs to the grid cells that are used
crop_data_combined_df['id'] = crop_data_combined_df.y * xsize + crop_data_combined_df.x

# add coordinates -- note that these are the centers of the cells
crop_data_combined_df['lon'] = eu_x + (crop_data_combined_df.x + 0.5) / 12.0
crop_data_combined_df['lat'] = eu_y + (crop_data_combined_df.y + 0.5) / 12.0

# same in a CSV format that is easy to use later
crop_data_combined_df.to_csv('{}/eu_combined_estimate.csv'.format(data_base_dir),
                             index=False)


###############################################################################
# 2. Process slope and water ratio data from GAEZ
# note: filenames should be gaez_slope_30p.zip and gaez_water.zip

slope = gaezdata()
slope.load('data.asc', '{}/gaez_slope_30p.zip'.format(data_base_dir), dtype=np.float)

slope1 = np.zeros([xsize, ysize], dtype = np.float)
for x in range(xsize):
    slope1[x,:] = slope._values[off_y:(off_y + ysize), x + off_x]

slope_df = pd.DataFrame([
    {'id' : y * xsize + x, 'slope_pct' : slope1[x, y]} for (x, y) in
    itertools.product(range(xsize), range(ysize)) if slope1[x, y] >= 0
    ])


water = gaezdata()
water.load('data.asc', '{}/gaez_water.zip'.format(data_base_dir), dtype=np.float)

water1 = np.zeros([xsize, ysize], dtype = np.float)
for x in range(xsize):
    water1[x,:] = water._values[off_y:(off_y + ysize), x + off_x]

water_df = pd.DataFrame([
    {'id' : y * xsize + x, 'water_pct' : water1[x, y]} for (x, y) in
    itertools.product(range(xsize), range(ysize)) if slope1[x, y] >= 0
    ])


# calculate a combined scaling factor, save this together
ws_df = pd.merge(water_df, slope_df, on='id')
ws_df['factor'] = np.maximum(1.0 - (ws_df.water_pct + ws_df.slope_pct) / 100.0, 0.0)
ws_df2 = ws_df[['id', 'factor']]


# do a matching with coastline data
ne_10m = gpd.GeoDataFrame.from_file(data_base_dir + '/ne_10m_land.shp')

cells_per_deg = 12
ne_overlap_res = list()
for x in tqdm.tqdm(range(xsize)):
    for y in range(ysize):
        poly1 = shp.geometry.Polygon([
            (eu_x + x/cells_per_deg, eu_y + y/cells_per_deg),
            (eu_x + (x+1)/cells_per_deg, eu_y + y/cells_per_deg),
            (eu_x + (x+1)/cells_per_deg, eu_y + (y+1)/cells_per_deg),
            (eu_x + x/cells_per_deg, eu_y + (y+1)/cells_per_deg),
            (eu_x + x/cells_per_deg, eu_y + y/cells_per_deg)
                ])
        i1 = ne_10m.intersection(poly1)
        s1 = sum(i1.area)
        s2 = poly1.area
        ne_overlap_res.append({'id': y * xsize + x, 'land': s1/s2})

ne_df = pd.DataFrame(ne_overlap_res)

cmb_df1 = pd.merge(ne_df, ws_df2)
cmb_df1['factor_cmb'] = np.maximum(cmb_df1.land + cmb_df1.factor - 1.0, 0.0)

cmb_df2 = cmb_df1[['id', 'factor_cmb']]
cmb_df2.to_csv('{}/eu_land_share.csv'.format(data_base_dir), index=False)









