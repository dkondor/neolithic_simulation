# agri_model_estimate.r -- script to estimate (relative) variations in
# Neolithic agriculture based on climate data using an agricultural
# yield emulator

library(reshape2)
library(ncdf4)


# base directory for data files -- need to adjust to the actual path!
# This should have a 'climate' subdirectory that contains
# the yield estimator data files:
#   LPJmL_spring_wheat_ggcmi_phase2_emulator_A0.nc4
#   LPJmL_winter_wheat_ggcmi_phase2_emulator_A0.nc4
# together with the outputs of the "climate_aggregate.r" script:
#   AgMERRA_months_avg_eu_df.RData
#   climate_processed_2.5_5kyr.RData
#   climate_processed_5_7.5kyr.RData
#   climate_processed_7.5_10kyr.RData
base_dir = '~/CSH/HoloSim/data/'
datadir = paste0(base_dir, 'climate')


# crop estimator types used
crop_types = c("LPJmL_spring_wheat", "LPJmL_winter_wheat")


######################################################################
# 1. Helper functions

# 1.1. CO2 data
# load the CO2 concentration data -- TODO: indicate the source of this file!
co2 = read.csv(paste0(base_dir, "/ghg-concentrations2.csv"))

# select whichever measurement is available
co2$co2 = co2$epica_vostok
co2$co2[is.na(co2$co2)] = co2$law_dome[is.na(co2$co2)]
co2$co2[is.na(co2$co2)] = co2$siple_station[is.na(co2$co2)]
# we only care about the recent past
co2 = co2[co2$year > -15000, c("year", "co2")]

# function to approximate the atmospheric CO2 concentration in a
# given year
get_co2 = function(year1) {
  if(year1 < min(co2$year) | year1 > max(co2$year)) stop("Year out of range for CO2 data!\n")
  i1 = which(co2$year == max(co2$year[co2$year < year1]))
  i2 = which(co2$year == min(co2$year[co2$year >= year1]))
  if(year1 == co2$year[i2]) return(co2$co2[i2])
  tmp1 = year1 - co2$year[i1]
  tmp2 = co2$year[i2] - year1
  tmp3 = tmp1 + tmp2
  return(co2$co2[i1] * (1-tmp1/tmp3) + co2$co2[i2] * (1-tmp2/tmp3))
}


# 1.2. yield emulator model
# calculate yields based on a vector of values, K is the matrix
# of coefficients loaded from the data files
calculate_yield_m = function(K, C, T, W, N) {
  return(
    K[, 1] + K[, 2] * C + K[, 3] * T + K[, 4] * W + K[, 5] * N + K[, 6] * C*C
    + K[, 7] * C*T + K[, 8] * C*W + K[, 9] * C*N + K[, 10] * T*T + K[, 11] * T*W
    + K[, 12] * T*N + K[, 13] * W*W + K[, 14] * W*N + K[, 15] * N*N
    + K[, 16] * C*C*C + K[, 17] * C*C*T + K[, 18] * C*C*W + K[, 19] * C*C*N
    + K[, 20] * C*T*T + K[, 21] * C*T*W + K[, 22] * C*T*N + K[, 23] * C*W*W
    + K[, 24] * C*W*N + K[, 25] * C*N*N + K[, 26] * T*T*T + K[, 27] * T*T*W
    + K[, 28] * T*N + K[, 29] * T*W + K[, 30] * T*W*N + K[, 31] * T*N*N
    + K[, 32] * W*W*W + K[, 33] * W*W*N + K[, 34] * W*N*N)
}


####################################################################
# 2. Preprocess the baseline climate data

# load the baseline climate data and aggregate it to the
# 0.5x0.5 degree grid used in the past climate data
# (originally, it is in a 0.25x0.25 degree grid)
load(paste0(datadir, '/AgMERRA_months_avg_eu_df.RData'))
# variable name: avgcl

avgcl$lon = floor(2 * avgcl$lon) / 2
avgcl$lat = floor(2 * avgcl$lat) / 2

# aggregate in space and calculate yearly averages
# note: this does averaging by space and time together
avgcl = aggregate(avgcl[c("T", "W")], by = avgcl[c("lon", "lat")], FUN=mean)

# limit to middle latitudes: up to 58 deg north (tip of Denmark)
# and down to 36.5 degrees (bottom of Turkey, Spain)
avgcl = avgcl[avgcl$lat < 58 & avgcl$lat > 36.5,]


#####################################################################
# 3. Main processing based on the past climate data

# note: we process each segment of the past climate data separately
# set the segments here
years_segments = c("2.5_5", "5_7.5", "7.5_10")
start_years_BP = c(5000, 7500, 10000)
# each segment has 2,500 years of data
nyears = 2500


for(i1 in 1:3) {
  years_select = years_segments[i1]
  year0 = 1950 - start_years_BP[i1]

  # load the past climate data -- variable name: climate1
  load(paste0(datadir, '/climate_processed_',years_select,'kyr.RData'))

  # limit to middle latitudes: up to 58 deg north (tip of Denmark)
  # and down to 36.5 degrees (bottom of Turkey, Spain)
  climate1 = climate1[climate1$lat < 58 & climate1$lat > 36.5,]

  # Do the estimation for both crop types
  for(crop_type in crop_types) {
    # load the model
    agri1 = nc_open(paste0(datadir, "/", crop_type, "_ggcmi_phase2_emulator_A0.nc4"))

    # extract the variable that corresponds to rain-fed agriculture
    K_rf = ncvar_get(agri1, "K_rf")

    # convert to a "sparse matrix" representation, only keeping cells with
    # nonzero coefficients

    # find out how many cells have nonzero result
    tmp1 = rowSums(K_rf == 0, dims = 2)
    # just deal with the nonzero entries
    tmp2 = which(tmp1 < 34, arr.ind = TRUE)

    # save the coordinate values and corresponding coordinate
    coords = data.frame(row = tmp2[,1], col = tmp2[,2])
    coords$lon = (coords$row - 1)*0.5 - 180
    coords$lat = 90 - (coords$col - 1)*0.5
    
    res1 = data.frame()
    for(i in 1:nyears) {
      year1 = year0 + nyears - i # note: time in climate data is reversed
      tmp1 = climate1[climate1$year == i, c("T", "W", "lon", "lat")]
      names(tmp1)[1:2] = c("T1", "W1")
      tmp1 = merge(tmp1, avgcl, by=c("lon", "lat"))
      tmp1$dT = tmp1$T1 - tmp1$T
      tmp1$rW = tmp1$W1 / (30 * tmp1$W)
      
      tmp1 = merge(tmp1, coords, by=c("lon", "lat"))
      # TODO: this should be done only once in the loop
      K_rf2 = matrix(0, nrow(tmp1), 34)
      for(j in 1:nrow(tmp1)) {
        K_rf2[j, ] = K_rf[tmp1$row[j], tmp1$col[j], ]
      }
      
      tmp1$Y = calculate_yield_m(K_rf2, get_co2(year1),
                                 tmp1$dT, tmp1$rW, 0)
      tmp1$year = year1
      res1 = rbind(res1, tmp1[c("lon", "lat", "year", "Y")])
      rm(tmp1)
    }
    
    save(res1, file=paste0(datadir, "/",crop_type,"_yield_estimate_",years_select,"kyr.RData"))
    rm(K_rf, coords)
    nc_close(agri1)
  }
  
  rm(climate1)
  gc()
}


#####################################################################
# 4. Calculate relative yields over the whole time period

# 4.1. estimate yields in the baseline period and use this to select
# the crop type with higher yield in each cell

baseline = list()

for(i in 1:2) {
  agri1 = nc_open(paste0(datadir, "/", crop_types[i], "_ggcmi_phase2_emulator_A0.nc4"))
  
  # extract the variable that corresponds to rain-fed agriculture
  K_rf = ncvar_get(agri1, "K_rf")
  tmp1 = rowSums(K_rf == 0, dims = 2)
  tmp2 = which(tmp1 < 34, arr.ind = TRUE)
  
  # save the coordinate values and corresponding coordinate
  coords = data.frame(row = tmp2[,1], col = tmp2[,2])
  coords$lon = (coords$row - 1)*0.5 - 180
  coords$lat = 90 - (coords$col - 1)*0.5
  
  Ktmp = matrix(0, nrow(coords), 34)
  for(j in 1:nrow(coords)) {
    Ktmp[j, ] = K_rf[coords$row[j], coords$col[j], ]
  }
  coords$Ybl = calculate_yield_m(Ktmp, 360, 0, 1, 200)
  
  baseline[[i]] = coords
  rm(coords, K_rf)
  nc_close(agri1)
}

# combine, select the better yield in each cell
baseline_cmb = merge(baseline[[1]][c("lon", "lat", "Ybl")],
                     baseline[[2]][c("lon", "lat", "Ybl")],
                     by = c("lon", "lat"), all = TRUE)
baseline_cmb$max_ix = 1
baseline_cmb$max_ix[is.na(baseline_cmb$Ybl.x)] = 2
baseline_cmb$max_ix[baseline_cmb$Ybl.y > baseline_cmb$Ybl.x] = 2
rm(baseline)


crop_data_cmb = data.frame()
for(years_select in c("7.5_10", "5_7.5", "2.5_5")) {
  crop_data = list()
  # load both crops
  for(i in 1:2) {
    load(paste0(datadir, "/",crop_types[i],"_yield_estimate_",years_select,"kyr.RData"))
    crop_data[[i]] = res1
    rm(res1)
  }
  # select the appropriate crop based on the baseline and
  # combine them to the common result
  for(i in 1:2) {
    crop_data[[i]] = merge(crop_data[[i]],
                           baseline_cmb[baseline_cmb$max_ix == i, c("lon", "lat")],
                           by = c("lon", "lat"))
    crop_data_cmb = rbind(crop_data_cmb, crop_data[[i]])
  }
  rm(crop_data)
  gc()
}

# remove negative values
crop_data_cmb$Y[crop_data_cmb$Y < 0] = 0
# calculate average for each cell
tmp1 = aggregate(crop_data_cmb["Y"],
                 by = crop_data_cmb[c("lon", "lat")], FUN=mean)
names(tmp1)[3] = "Yavg"
# calculate relative values by scaling with the average
crop_data_cmb = merge(crop_data_cmb, tmp1, by=c("lon", "lat"))
crop_data_cmb$Yr = crop_data_cmb$Y / crop_data_cmb$Yavg

# add IDs, based on the same scheme used for the GAEZ crop data
eu_x = -12 # lower left corner for the study area
eu_y = 36.5
xsize = 732 # number of cells in x direction (longitude)
ysize = 258 # number of cells in y direction (latitude)
scale = 12 # number of cells per degree (i.e. each cell has 1/12th degree linear size)

crop_data_cmb$x = (crop_data_cmb$lon - eu_x)*scale
crop_data_cmb$y = (crop_data_cmb$lat - eu_y)*scale
crop_data_cmb$cid = crop_data_cmb$y * xsize + crop_data_cmb$x
# save the result as a binary file (used only for evaluating ACFs later)
save(crop_data_cmb, file = paste0(datadir, "/crop_data_cmb.RData"))

# save only the relative yields and IDs, sorted by year,
# so that it can be used directly during the simulation
crop_data_cmb = crop_data_cmb[c("year", "cid", "Yr")]
crop_data_cmb = crop_data_cmb[order(crop_data_cmb$year,
                                    crop_data_cmb$cid),]
write.table(crop_data_cmb, paste0(datadir, "/crop_data_cmb.dat"),
            row.names = FALSE, col.names = FALSE)






