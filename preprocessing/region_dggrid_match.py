#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
region_dggrid_match.py -- classify the DGGRID hexagons into the new regions
used for spatial aggregation
"""

import pandas as pd
import geopandas as gpd
import shapely as shp

# base directory with data files -- need to be adjusted!
base_dir = '/home/dkondor/CSH/HoloSim/data'

new_regions = pd.read_csv('{}/new_grids_all.csv'.format(base_dir))

def get_poly(x):
    return [
        (x.lonmin, x.latmin),
        (x.lonmax, x.latmin),
        (x.lonmax, x.latmax),
        (x.lonmin, x.latmax),
        (x.lonmin, x.latmin),
    ]

new_regions2 = gpd.GeoSeries(list(shp.geometry.Polygon(
    get_poly(new_regions.iloc[i])
    ) for i in range(new_regions.shape[0])), crs = "EPSG:4326")

new_regions2.sindex


pts = pd.read_csv('{}/gaez_crop/eu_dggrid_coords_K_land4a.csv'.format(base_dir))
pts2 = gpd.GeoDataFrame({'dgid': pts.dgid, 'geometry': (shp.geometry.Point(
    pts.iloc[i].lon, pts.iloc[i].lat) for i in range(pts.shape[0]))},
    crs = "EPSG:4326")

pts2.sindex


# do the query
matches = pts2.sindex.query_bulk(new_regions2)
matches2 = pd.concat((new_regions.iloc[matches[0]].reset_index(drop = True),
                      pts.iloc[matches[1]].reset_index(drop = True)), axis=1)
matches2 = matches2[['ID', 'dx', 'dy', 'res', 'dgid']]
matches2.to_csv('{}/new_grid_aggr.csv'.format(base_dir), index = False)


# create a "flat" version of the above, where one numeric ID corresponds to
# each region on each aggregation type -- this will help with efficiently
# aggregating simulation results and limiting the output file sizes
new_regions_ix = new_regions[['ID', 'dx', 'dy', 'res']].reset_index()
new_regions_ix.to_csv('{}/new_grid_ids.csv'.format(base_dir), index = False)

matches2_ix = pd.merge(matches2, new_regions_ix, on = ['ID', 'dx', 'dy', 'res'])
matches2_ix = matches2_ix[['dgid', 'index']]
matches2_ix.to_csv('{}/new_grid_aggr_flat.csv'.format(base_dir), index = False)



###############################################################################
# next steps: filter out cells with missing climate data
climate_data = pd.read_table('{}/climate/crop_data_cmb.dat'.format(base_dir),
							 delim_whitespace=True, names = ['year','wid', 'factor'])
cell_years = climate_data.groupby('wid').count()
valid_cells = cell_years.index[cell_years.year == 7500]

cells_dg = pd.read_csv('{}/climate/cell_ids_dggrid.csv'.format(base_dir))
valid_dgid = cells_dg.dgid[cells_dg.cid.isin(valid_cells)]

# subset that is actually used in the simulations
valid_dgid = valid_dgid[valid_dgid.isin(pts.dgid)].unique()
# 31370 cells remain (out of 36400)
valid_dgid = pd.DataFrame({'dgid': valid_dgid})
valid_dgid.to_csv('{}/new_dgid_aggr_filter.dat'.format(base_dir), header = False, index = False)


###############################################################################
# create a CSV version of DGGRID polygons used in the simulation, to be used
# for creating videos of simulation results

shp1 = gpd.GeoDataFrame.from_file('{}/dggrid/land_isea3h12eu.shp'.format(base_dir))
shp1 = shp1[shp1.global_id.isin(pts.dgid)]
shp1.reset_index(inplace = True, drop = True)

shp2 = pd.concat((pd.DataFrame({'id': shp1.global_id[i],
								'x': shp1.geometry[i].boundary.coords[j][0],
								'y': shp1.geometry[i].boundary.coords[j][1]}
							   for j in range(len(shp1.geometry[i].boundary.coords)))
				  for i in range(shp1.shape[0])))
shp2.to_csv('{}/dggrid/dggrid_poly.csv'.format(base_dir), index = False, header = False)




