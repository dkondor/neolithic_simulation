#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dggrid_process.py -- process DGGRID hexagonal grid, estimate agricultural
productivitiy data

Main approach: we calculate intersections (in projected coordinate space)
with geographies of interest: land data, GAEZ data, etc.
"""


import geopandas as gpd
import shapely as shp
import shapely.strtree as strt
import pandas as pd
import tqdm
import itertools
import networkx as nx

# base directory for data files -- need to adjust this to the correct paths!
base_dir = '/home/dkondor/CSH/HoloSim/data'

###############################################################################
# 1. Identify which cells fall on land vs sea; make note of coastal cells
# separately as well

# load the shapefile with the hexagon grid
shp1 = gpd.GeoDataFrame.from_file('{}/dggrid/isea3h12eu.shp'.format(base_dir))
# note: convert IDs to integers (stored as strings in the shapefile to avoid overflows)
shp1['global_id'] = list(int(x) for x in shp1.global_id)

# load a shapefile with polygons corresponding to the world's landmass
# data from here: https://www.naturalearthdata.com/downloads/50m-physical-vectors/
land = gpd.GeoDataFrame.from_file('{}/gaez_crop/ne_50m_land.shp'.format(base_dir))


res1 = pd.DataFrame({'id': shp1['global_id'], 'land': False, 'inside': False})
for i in tqdm.tqdm(range(len(land))):
    tmp1 = shp1.intersects(land.iloc[i].geometry)
    res1.land = res1.land | tmp1
    tmp2 = tmp1 & shp1.within(land.iloc[i].geometry)
    res1.inside = res1.inside | tmp2
# 140525 land cells in total (runtime: ~30 min)

res1['coastal'] = res1.land & (~res1.inside)
# save this result
res1.to_csv('{}/dggrid/land_coast_isea3h12eu.csv'.format(base_dir), index = False)


# filter the shapefile to include only land cells
shp1 = pd.merge(shp1, res1, left_on='global_id', right_on='id')
shp1 = shp1[shp1['land']]
shp1 = shp1[['global_id', 'geometry']]
shp1.reset_index(inplace=True, drop=True)
shp1.to_file('{}/dggrid/land_isea3h12eu.shp'.format(base_dir))


###############################################################################
# 2. Match with the GAEZ rectangular grid, calculate weighted average of
# agricultural yield and related variables in each cell


# parameters for creating the grid
eu_x = -12 # lower left corner for EU
eu_y = 36.5 
xsize = 732
ysize = 258

cells_per_deg = 12

eu_cells = gpd.GeoDataFrame((
    {'id' : y * xsize + x,
     'geometry': shp.geometry.Polygon([
            (eu_x + x/cells_per_deg, eu_y + y/cells_per_deg),
            (eu_x + (x+1)/cells_per_deg, eu_y + y/cells_per_deg),
            (eu_x + (x+1)/cells_per_deg, eu_y + (y+1)/cells_per_deg),
            (eu_x + x/cells_per_deg, eu_y + (y+1)/cells_per_deg),
            (eu_x + x/cells_per_deg, eu_y + y/cells_per_deg)
                ])}
     for (x,y) in itertools.product(range(xsize), range(ysize))), crs = "EPSG:4326")


srtree = strt.STRtree(eu_cells.geometry)
geom_map = dict((id(eu_cells.geometry[i]), eu_cells.id[i]) for i in range(eu_cells.shape[0]))

res1 = list()
for i in tqdm.tqdm(range(shp1.shape[0])):
    tmp1 = srtree.query(shp1.geometry[i])
    area1 = shp1.geometry[i].area if len(tmp1) > 0 else 0
    for j in range(len(tmp1)):
        id1 = geom_map[id(tmp1[j])]
        tmp2 = tmp1[j].intersection(shp1.geometry[i])
        if tmp2.area > 0:
            res1.append({'dgid' : shp1.global_id[i], 'cid' : id1,
                         'match' : tmp2.area / area1})
        
gaez_match = pd.DataFrame(res1)

gaez_match.to_csv('{}/gaez_crop/dggrid_gaez_match.csv'.format(base_dir), index=False)


# further processing: for GAEZ data, do a match here (weighting values by the relative overlap)
# for weather data, do a separate matching with the larger cells used there

eu_est = pd.read_csv('{}/gaez_crop/eu_combined_estimate.csv'.format(base_dir))
eu_est = eu_est[['id','v']]
eu_est2 = pd.merge(eu_est, gaez_match, left_on = 'id', right_on = 'cid')
eu_est2['v2'] = eu_est2.v * eu_est2.match
eu_est2 = eu_est2[['dgid', 'v2']]

eu_est3 = eu_est2.groupby('dgid').sum()
# 79481 cells remain

# cut to the largest connected component
dgids = set(eu_est3.index)
g = nx.Graph()
nbr = open('{}/dggrid/isea3h12eun.nbr'.format(base_dir), 'r')
for l in nbr.readlines():
    l2 = l.split()
    id1 = int(l2[0])
    if id1 in dgids:
        for i in range(1,len(l2)):
            id2 = int(l2[i])
            if id2 in dgids:
                g.add_edge(id1, id2)
    

c1 = list(nx.connected_components(g))
len(c1) # 34
maxc = max(len(x) for x in c1) # 74797
ix1 = -1
for i in range(len(c1)):
    if len(c1[i]) == maxc:
        ix1 = i

eu_est4 = eu_est3[eu_est3.index.isin(c1[ix1])]
eu_est4.to_csv('{}/gaez_crop/eu_dggrid_estimate1.csv'.format(base_dir))

# match with center coordinates (needed for the simulation)
coords1 = pd.read_csv('{}/dggrid/isea3h12eup.txt'.format(base_dir),
                      names=['dgid', 'lon', 'lat'])
coords1 = coords1[~coords1.lon.isna()] # note: this is a workaround for DGGRID adding the string
coords1['dgid'] = pd.to_numeric(coords1['dgid']) # "END" to the last line of output in some cases
eu_est5 = pd.merge(eu_est4, coords1, left_index=True, right_on = 'dgid')
eu_est5 = eu_est5[['dgid', 'lon', 'lat', 'v2']]
eu_est5.to_csv('{}/gaez_crop/eu_dggrid_coords_K.csv'.format(base_dir), index=False)


# similar, but with other data: slope and loess soil
gaez_match2 = gaez_match
gaez_match2 = gaez_match2[gaez_match2.dgid.isin(c1[ix1])]


eu_slope = pd.read_csv('{}/gaez_crop/eu_land_share.csv'.format(base_dir))
eu_slope = pd.merge(eu_slope, gaez_match2, left_on = 'id', right_on = 'cid')
eu_slope['factor2'] = eu_slope.factor_cmb * eu_slope.match
eu_slope = eu_slope[['dgid', 'factor2']]
eu_slope = eu_slope.groupby('dgid').sum()
eu_slope.to_csv('{}/gaez_crop/eu_dggrid_land_share.csv'.format(base_dir))


###############################################################################
# 3. Match with the coarser grid (0.5x0.5 degree rectangles) used in the
# climate data

# identify the unique cell IDS
climate_cells = set()
with open('{}/climate/crop_data_cmb.dat'.format(base_dir), 'r') as cellsf:
    for l in cellsf.readlines():
        l2 = l.split()
        climate_cells.add(int(l2[1]))
# note: 2727 cells in total
climate_cells = pd.DataFrame({'id': list(climate_cells)})
climate_cells['y'] = list(int(x / xsize) for x in climate_cells.id)
climate_cells['x'] = climate_cells.id % xsize


cellsize = 6 # one cell in the climate data is 6x6 cells in the original dataset
climate_cells_shp = gpd.GeoDataFrame((
    {'id' : climate_cells.id[i],
     'geometry': shp.geometry.Polygon([
            (eu_x + climate_cells.x[i]/cells_per_deg,
             eu_y + climate_cells.y[i]/cells_per_deg),
            (eu_x + (climate_cells.x[i] + cellsize)/cells_per_deg,
             eu_y + climate_cells.y[i]/cells_per_deg),
            (eu_x + (climate_cells.x[i] + cellsize)/cells_per_deg,
             eu_y + (climate_cells.y[i] + cellsize)/cells_per_deg),
            (eu_x + climate_cells.x[i]/cells_per_deg,
             eu_y + (climate_cells.y[i] + cellsize)/cells_per_deg),
            (eu_x + climate_cells.x[i]/cells_per_deg,
             eu_y + climate_cells.y[i]/cells_per_deg)
                ])}
     for i in range(climate_cells.shape[0])), crs = "EPSG:4326")


srtree = strt.STRtree(climate_cells_shp.geometry)
geom_map = dict((id(climate_cells_shp.geometry[i]), climate_cells_shp.id[i])
                for i in range(climate_cells_shp.shape[0]))

res1 = list()
for i in tqdm.tqdm(range(shp1.shape[0])):
    tmp1 = srtree.query(shp1.geometry[i])
    area1 = shp1.geometry[i].area if len(tmp1) > 0 else 0
    for j in range(len(tmp1)):
        id1 = geom_map[id(tmp1[j])]
        tmp2 = tmp1[j].intersection(shp1.geometry[i])
        if tmp2.area > 0:
            res1.append({'dgid' : shp1.global_id[i], 'cid' : id1,
                         'match' : tmp2.area / area1})
        
climate_match = pd.DataFrame(res1)
# 87032 matching rows
sum_match = climate_match[['dgid', 'match']].groupby('dgid').sum()
# 56412 remaining (most cells match only one climate cell)
sum_match = sum_match.rename(columns={'match': 'smatch'})
climate_match = pd.merge(climate_match, sum_match, left_on = 'dgid', right_index=True)
climate_match['factor'] = climate_match.match / climate_match.smatch
climate_match = climate_match[['dgid', 'cid', 'factor']]
climate_match.to_csv('{}/climate/cell_ids_dggrid.csv'.format(base_dir), index=False)



###############################################################################
# 4. Cut the study area to the ones used in the main simulations and match
# with the data about coastal cells;
# note: this is done based on the center points of cells
pts_new = pd.read_csv('{}/gaez_crop/eu_dggrid_coords_K.csv'.format(base_dir))
pts_new2 = gpd.GeoDataFrame(({'id': int(pts_new.iloc[i].dgid),
                'keep': (pts_new.iloc[i].lon < 7.8 or pts_new.iloc[i].lat < 45.5),
                'geometry': shp.geometry.Point(pts_new.iloc[i][['lon','lat']])}
                            for i in range(pts_new.shape[0])), crs = "EPSG:4326")
filter_shp = gpd.GeoDataFrame.from_file('{}/neolithic_4000bc_cut_ne.shp'.format(base_dir))

tmpix = pts_new2.columns.get_loc('keep')
for i in tqdm.tqdm(range(pts_new2.shape[0])):
    if not pts_new2.iloc[i, tmpix]:
        pts_new2.iloc[i, tmpix] = pts_new2.iloc[i].geometry.within(filter_shp.iloc[0].geometry)

pts_new3 = pts_new2[['id', 'keep']]
pts_new3 = pd.merge(pts_new3, pts_new, left_on = 'id', right_on = 'dgid')
# keep only the cells in the selected area
pts_new3 = pts_new3[pts_new3.keep]
# additional filter -- leave out eastern Anatolia
pts_new3 = pts_new3[pts_new3.lon < 32]
# 36400 cells remain

# match with the status of coastal cells
coastal_status = pd.read_csv('{}/dggrid/land_coast_isea3h12eu.csv'.format(base_dir))
pts_new3 = pd.merge(pts_new3, coastal_status, on = 'id')
pts_new3.reset_index(inplace=True, drop=True)

# create a "type" field for land / coastal status (coastal: 1, land: 2)
# note: we do not use coastal status for cells on the Atlantic coast, only
# in the Mediterranean
pts_new3['type'] = 2
x1 = -8.9
y1 = 37
x2 = 5.3
y2 = 47.5
tmpix = pts_new3.columns.get_loc('type')
for i in range(pts_new3.shape[0]):
	if pts_new3.coastal[i]:
		if pts_new3.lat[i] <= 47.5:
			# cells are only considered south of latitude 47.5 and 
			# southeast from the line between (x1, y1) and (x2, y2)
			tmp1 = ((pts_new3.lon[i] - x1)*(y2 - y1) -
				(pts_new3.lat[i] - y1)*(x2 - x1))
			if tmp1 >= 0:
				pts_new3.iloc[i, tmpix] = 1
# note: 1810 coastal cells in total
# save the result that will be used as the input of the simulation
pts_new4 = pts_new3[['dgid', 'lon', 'lat', 'v2', 'type']]
pts_new4.to_csv('{}/gaez_crop/eu_dggrid_coords_K_land4a.csv'.format(base_dir),index=False)

