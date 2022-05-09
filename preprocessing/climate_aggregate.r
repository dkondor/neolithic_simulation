# climate_aggregate.r -- pre-process and aggregate the past and baseline
#  climate datasets

library(ncdf4)
library(reshape2)

# base directory for data files (downloaded datasets)
# this should include a 'climate' subdirectory with the following:
# 1. AgMERRA climate dataset (baseline climate between 1980 and 2010),
#   precipitation and temperature
#   (AgMERRA_{1980..2010}_{prate,tavg}.nc4 files)
# 2. Past climate dataset, precipitation and temperature, for the
#   time period between 2,500 and 10,000 BP
#   (bias_regrid_{tas,pr}_{2.5_5,5_7.5,7.5_10}kyr.nc)
# see the script "agri_data_download.fsh" for the source of these
# datasets
base_dir = '~/CSH/HoloSim/data/'
datadir = paste0(base_dir, 'climate')


# filter a dataframe (with lon and lat columns) to the study area
# (using a coarse bounding box)
df_filter_eu = function(df) {
  eu_coords = c(-12, 28, 40, 72)
  return(df[df$lon >= eu_coords[1] & df$lon <= eu_coords[3] &
              df$lat >= eu_coords[2] & df$lat <= eu_coords[4], ])
}



#######################################################################
# 1. pre-process the baseline climate dataset (AgMERRA, 1980-2010),
# calculate monthly averages in the study area
base_year = 1980
nyears = 31

# size of the grid used in the dataset
xsize = 1440
ysize = 720
nmonths = 12

avgT = array(0, c(xsize, ysize, nmonths))
avgW = array(0, c(xsize, ysize, nmonths))

filesT = list()
filesW = list()

# open the data files
for(i in 1:nyears) {
  year = base_year + i - 1
  filesW[[i]] = nc_open(paste0(datadir, '/AgMERRA_', year, '_prate.nc4'))
  filesT[[i]] = nc_open(paste0(datadir, '/AgMERRA_', year, '_tavg.nc4'))
}

month_len = c(31,28,31,30,31,30,31,31,30,31,30,31)

# 1.1. calculate monthly averages for the whole dataset
for(i in 1:nmonths) {
  start_ix = 1
  if(i > 1) start_ix = start_ix + sum(month_len[1:(i-1)])
  len1 = month_len[i]
  days_sum = 0
  for(j in 1:nyears) {
    start_ix2 = start_ix
    len2 = len1
    if(filesT[[j]]$dim$time$len == 366) {
      # if this is a leap year
      if(i > 2) start_ix2 = start_ix2 + 1
      if(i == 2) len2 = len2 + 1
    }
    days_sum = days_sum + len2 # count the actual days in sum
    
    # sum for temperature
    tmp1 = ncvar_get(filesT[[j]], 'tavg', start = c(1,1,start_ix2),
                     count = c(-1,-1,len2))
    avgT[,,i] = avgT[,,i] + rowSums(tmp1, na.rm = FALSE, dims = 2)
    
    # sum for precipitation
    tmp1 = ncvar_get(filesW[[j]], 'prate', start = c(1,1,start_ix2),
                     count = c(-1,-1,len2))
    avgW[,,i] = avgW[,,i] + rowSums(tmp1, na.rm = FALSE, dims = 2)
  }
  
  # divide by the number of days
  avgT[,,i] = avgT[,,i] / days_sum
  avgW[,,i] = avgW[,,i] / days_sum
}

for(i in 1:nyears) {
  nc_close(filesW[[i]])
  nc_close(filesT[[i]])
}
rm(filesW, filesT)


# 1.2. assign explicit coordinates and filter for the study area
avgT3 = melt(avgT, na.rm = TRUE)
avgW3 = melt(avgW, na.rm = TRUE)
rm(avgT, avgW)

names(avgT3) = c("row", "col", "month", "T")
names(avgW3) = c("row", "col", "month", "W")

avgcl = merge(avgT3, avgW3, by = c("row", "col", "month"))
rm(avgT3, avgW3)

# note: filter out cells that do not have results for all months
tmp1 = aggregate(avgcl["T"], by=avgcl[c("row", "col")], FUN=length)
tmp1 = tmp1[tmp1$T == 12, c("row", "col")] # 251632 remains out of 251644
avgcl = merge(avgcl, tmp1, by=c("row", "col"))

# add coordinates
# note: coordinates are floats, but in this dataset, they are
# 1/4 degree resolution, so they can be represented exactly
# (no rounding errors or similar issues)
avgcl$lon = (avgcl$row - 1) / 4.0 # do not adjust center position for
avgcl$lat = 90 - (avgcl$col - 1) / 4.0 # better aggregation later
# adjust longitudes: we use negative values for the western hemisphere
avgcl$lon[avgcl$lon >= 180] = avgcl$lon[avgcl$lon >= 180] - 360

# filter for the study area of western Eurasia
avgcl = df_filter_eu(avgcl)
# save separately
save(avgcl, file = paste0(datadir, '/AgMERRA_months_avg_eu_df.RData'))

rm(avgcl)


#####################################################################
# 2. pre-process the past climate dataset

# one file contains 2500 years
nyears = 2500

# data is in 2,500 year chunks
for(years_select in c("2.5_5", "5_7.5", "7.5_10")) {
  # open the data files
  tmp3 = nc_open(paste0(datadir, "/bias_regrid_pr_",
                        years_select, "kyr.nc"))
  tmp4 = nc_open(paste0(datadir, "/bias_regrid_tas_",
                        years_select, "kyr.nc"))
  
  climate1 = data.frame()
  
  for(i in 1:nyears) {
    i1 = (i-1)*12 + 1
    tas1 = ncvar_get(tmp4, "tas", start = c(1, 1, i1), count = c(-1, -1, 12))
    pr1 = ncvar_get(tmp3, "pr", start = c(1, 1, i1), count = c(-1, -1, 12))
    tas1 = rowMeans(tas1, na.rm = FALSE, dims = 2)
    pr1 = rowMeans(pr1, na.rm = FALSE, dims = 2)
    tas1 = melt(tas1, na.rm = TRUE)
    pr1 = melt(pr1, na.rm = TRUE)
    names(tas1) = c("row", "col", "T")
    names(pr1) = c("row", "col", "W")
    cl1 = merge(tas1, pr1, by=c("row", "col"))
    
    cl1$lon = (cl1$row - 1) * 0.5
    cl1$lat = (cl1$col - 1) * 0.5
    cl1$lon[cl1$lon >= 180] = cl1$lon[cl1$lon >= 180] - 360
    cl1 = df_filter_eu(cl1)
    
    cl1$year = i
    climate1 = rbind(climate1, cl1)
  }
  
  # save in binary format
  save(climate1,  file = paste0(datadir, "/climate_processed_",
                                years_select, "kyr.RData"))
  nc_close(tmp3)
  nc_close(tmp4)
  rm(climate1)
}






