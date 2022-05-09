# c14_data_acf_spd_new.r -- aggregate C14 dates into regions of
# quasi-rectangular tilings of the study area, calculate SPD of
# dates in each of the tiles
library(rcarbon)
library(foreach)
library(doParallel)

# base directory of all data files
base_dir = '~/CSH/HoloSim/data'
datadir = paste0(base_dir, '/population') # directory with C14 data

# number of CPU cores used for parallel computations -- adjust this as
# necessary (by default it uses one less than all available cores)
ncores = max(detectCores() - 1, 1)

# minimum number of dates to require in a region to do an SPD
# (500 is a conservative choice)
thresh = 500

# time range to cut dates to (in BP)
# these can be specified as command line arguments or adjusted here
args = commandArgs(trailingOnly=TRUE)
tmin1 = 5000
tmax1 = 10000
if(length(args) >= 2) {
  tmax1 = as.integer(args[1])
  tmin1 = as.integer(args[2])
  if(is.na(tmin1) || is.na(tmax1)) {
    error('Invalid date range arguments!')
  }
  if(tmin1 >= tmax1) {
    error('Invalid date range arguments!')
  }
}



#########################################################
# 1. load the data

# C14 dates
popc14 = read.csv(paste0(datadir, "/2021-09-08-dataset-nodup-no-sd-filter.csv"))
popc14 = popc14[!is.na(popc14$X),]
rownames(popc14) = popc14$X
popc14 = popc14[,-1]
names(popc14)[11:12] = c("caldate5", "caldate95")
popc14$caldatemean1 = (popc14$caldate5 + popc14$caldate95)/2


# load the tilings used for aggregation
load(paste0(base_dir, "/new_grids_all.RData"))

all_comb = unique(grids_all[c("res", "dx", "dy")])



#####################################################################
# 2. run aggregation and SPD creation in parallel for all tiles

# note: makeForkCluster() does not work on Windows -- uncomment the
# version below if needed
cl = parallel::makeForkCluster(ncores)
# cl = parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

all_res = foreach(i = 1:nrow(all_comb)) %dopar% {
  r = all_comb$res[i]
  dx = all_comb$dx[i]
  dy = all_comb$dy[i]
  tmp1 = popc14
  tmp1 = tmp1[tmp1$caldatemean1 <= tmax1 & tmp1$caldatemean1 >= tmin1,]
  
  grid1 = grids_all[grids_all$res == r & grids_all$dx == dx &
                      grids_all$dy == dy,]
  # match all points in the grid -- we could use a more efficient
  # method here instead of looping, but there are not so many points
  # and grid regions as well
  tmp1$bin = NA
  for(i in 1:nrow(tmp1)) {
    x = tmp1$lon[i]
    y = tmp1$lat[i]
    tmp2 = grid1$ID[grid1$lonmin <= x & grid1$lonmax > x &
                       grid1$latmin <= y & grid1$latmax > y]
    # note: length(tmp2) > 1 is an error
    if(length(tmp2) > 0) tmp1$bin[i] = tmp2[1]
  }
  tmp1 = tmp1[!is.na(tmp1$bin),]
  tmp1$cnt = 1
  tmp2 = aggregate(tmp1["cnt"], by=tmp1["bin"], FUN=length)
  
  spd_cmb = data.frame()
  errors = NULL
  
  for(b in tmp2$bin[tmp2$cnt >= thresh]) {
    x = tryCatch( {
      tmp3 = tmp1[tmp1$bin == b, ]
      tmp3c = calibrate(tmp3$c14age, tmp3$c14std, normalised = FALSE)
      tmp3$sitef = factor(tmp3$site)
      cbins = binPrep(sites = as.integer(tmp3$sitef),
                      ages = tmp3$c14age, h = 100)
      spd1 = spd(tmp3c, bins = cbins, timeRange = c(tmax1,tmin1))
      tmp4 = spd1$grid
      tmp4$bin_ID = b
      spd_cmb = rbind(spd_cmb, tmp4)
      1
    }, error = function(x){return(0)})
    if(x == 0) errors = c(errors, b)
  }
  
  res12 = list()
  res12$r = r
  res12$dx = dx
  res12$dy = dy
  res12$calBPmax = tmax1
  res12$calBPmin = tmin1
  res12$spd_cmb = spd_cmb
  res12$errors = errors
  return(res12)
}

# save the result
save(all_res, file = paste0(datadir, "/c14_spd_new_", tmax1, "_", tmin1, ".RData"))

