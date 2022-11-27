# c14_data_preprocess_new.r -- new preprocessing for C14 data that will
# make it easier to work together in R and Python
# replace all strings with numeric IDs and perform a binning with rcarbon

library(rcarbon)
Sys.setlocale("LC_COLLATE", "C")

# base directory of all data files
base_dir = '~/CSH/HoloSim/data'
datadir = paste0(base_dir, '/population/') # directory with C14 data

popc14 = read.csv(paste0(datadir, "/2021-09-08-dataset-nodup-no-sd-filter.csv"))
popc14 = popc14[!is.na(popc14$X),]
names(popc14)[12:13] = c("calBP5", "calBP95")

# remove unnecessary columns
popc14 = popc14[c("X", "labnr", "site", "c14age", "c14std", "lon", "lat",
                  "calBP5", "calBP95")]
# create numeric IDs for all sites
sites = sort(unique(popc14$site))
sites = data.frame(site = sites, sitef = 1:length(sites))
popc14 = merge(popc14, sites, by="site")
popc14 = popc14[c("X", "labnr", "sitef", "c14age", "c14std", "lon", "lat",
                  "calBP5", "calBP95")]
popc14 = popc14[order(popc14$X),]

# save this data
write.csv(sites, paste0(datadir, "c14_sites.csv"), row.names = FALSE)
write.csv(popc14, paste0(datadir, "c14_data_ids.csv"), row.names=FALSE)


# calibrate all dates, estimate the median dates
popc14c = calibrate(popc14$c14age, popc14$c14std, ncores = 1)

popc14$calBPmed = medCal(popc14c)

# calculate possible general binning with thresholds of 50, 100 and 200
# years, based on the full (unfiltered) dataset
# (note: there might be slight differences compared to doing this after
# cutting the dataset to a shorter time interval, but these should not
# be very significant)
popc14$bin50 = binPrep(popc14$sitef, popc14$calBPmed, 50)
popc14$bin100 = binPrep(popc14$sitef, popc14$calBPmed, 100)
popc14$bin200 = binPrep(popc14$sitef, popc14$calBPmed, 200)

# save this dataset
write.csv(popc14, paste0(datadir, "c14_data_bins1.csv"), row.names=FALSE)



