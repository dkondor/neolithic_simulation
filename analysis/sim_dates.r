# sim_dates.r -- script to simulate a "sampling" process on simulation
# results, evaluate the resulting sampled set of dates

library(ggplot2)
library(foreach)
library(parallel)

base_dir = "~/CSH/HoloSim/data"


#######################################################################
# 1. load base data on spatial aggregation, define functions
# for calculating ACFs

regions = read.csv(paste0(base_dir, "/new_grid_ids.csv"))
aggr_ids = read.csv(paste0(base_dir, "/new_grid_aggr_flat.csv"))
# filter for the r = 500 km tiling -- only this is used in this script
regoins = regions[regions$res == 500,]
aggr_ids = aggr_ids[aggr_ids$index %in% regions$index,]

# indices of one example tiling
tmpix = regions$index[regions$res == 500 & regions$dx == 0 & regions$dy == 0]


# helper for the below functions
get_diff_sign = function(acf1) {
  tmps = sign(diff(acf1))
  for(i in length(tmps):2) {
    if(tmps[i-1] == 0) tmps[i-1] = tmps[i]
  }
  return(tmps)
}

# get the location of the first maximum (after changing sign twice)
get_acf_max = function(acf1) {
  tmp1p = which(sign(acf1) < 0)
  if(length(tmp1p) == 0) return(NULL)
  tmp1p = min(tmp1p)
  tmp2p = which(diff(get_diff_sign(acf1)) == -2) + 1
  tmp2p = tmp2p[tmp2p > tmp1p]
  tmp2p = tmp2p[acf1[tmp2p] > 0]
  return(head(tmp2p, 1))
}

# get the location of the first (negative) minimum
get_acf_min = function(acf1) {
  if(length(acf1) < 3) return(NULL)
  tmp2n = which(diff(get_diff_sign(acf1)) == 2) + 1
  tmp2n = tmp2n[acf1[tmp2n] < 0]
  return(head(tmp2n, 1))
}


# do an ACF for one regional time series
do_one_acf_base = function(tmp2, ix, two_step_fit = FALSE,
                           need_ts = FALSE, no_detrend = FALSE) {
  l1 = 4000
  
  res1pd = NULL
  res2pd = NULL
  res1nd = NULL
  res2nd = NULL
  
  tmp2 = tmp2[order(tmp2$year),]
  
  if(no_detrend) {
    tmp2$pred = mean(tmp2$spd)
    tmp2$diff = tmp2$spd - tmp2$pred
  } else {
    # detrending
    N = mean(tail(tmp2$spd, 1000))
    x0 = min(tmp2$year[tmp2$spd >= N/2])
    T1 = 100 # TODO: maybe adjust this?
    nls1 = tryCatch( {
      if(two_step_fit) {
        nls1 = nls(spd ~ N / (1 + exp(-(year - x0)/T1)) , tmp2,
                   list(N = N))
        N = nls1$m$getPars()[1]
        nls1 = nls(spd ~ N / (1 + exp(-(year - x0)/T1)) , tmp2,
                   list(x0 = x0, T1 = T1))
      } else {
        nls1 = nls(spd ~ N / (1 + exp(-(year - x0)/T1)) , tmp2,
                   list(N = N, x0 = x0, T1 = T1))
      }
      nls1
    }, error = function(x) NULL)
    if(is.null(nls1)) return(NULL)
    
    tmp2$pred = predict(nls1)
    tmp2$diff = tmp2$spd - tmp2$pred
  }
  
  # calculate a coefficient of variation
  cv = sd(tmp2$diff) / mean(tmp2$pred)
  
  # do the ACF
  acf1 = acf(tmp2$diff, lag.max = l1, plot = FALSE)
  
  if(sum(is.nan(acf1$acf)) < l1) {
    tmp2p = get_acf_max(acf1$acf)
    tmp2n = get_acf_min(acf1$acf)
    res1pd = c(res1pd, acf1$lag[tmp2p])
    res2pd = c(res2pd, acf1$acf[tmp2p])
    res1nd = c(res1nd, acf1$lag[tmp2n])
    res2nd = c(res2nd, acf1$acf[tmp2n])
  }
  
  sum_res1 = list()
  sum_res1$ix = ix
  if(length(res1pd) > 0) sum_res1$all_res_max_dt = 
    data.frame(lag = res1pd, val = res2pd, ix = ix)
  if(length(res1nd) > 0) sum_res1$all_res_min_dt = 
    data.frame(lag = res1nd, val = res2nd, ix = ix)
  sum_res1$acfdt = data.frame(lag=acf1$lag, acf=acf1$acf, ix = ix)
  if(need_ts) sum_res1$tsdt = tmp2 # save the actual detrended TS
  sum_res1$cv = cv
  return(sum_res1)
}

# do the full analysis for one aggregated simulation result,
# stored as a list
do_one_acf.list = function(res1, two_step_fit = FALSE,
                           need_ts = FALSE, no_detrend = FALSE) {
  sum_res1 = list()
  for(i in 1:length(res1)) {
    ix = res1[[i]]$ix
    j = length(sum_res1) + 1
    sum_res1[[j]] = do_one_acf_base(res1[[i]]$spd, ix, two_step_fit,
                                    need_ts, no_detrend)
  }
  return(sum_res1)
}

do_one_acf.data.frame = function(res1, two_step_fit = FALSE,
                                 need_ts = FALSE, no_detrend = FALSE) {
  sum_res1 = list()
  for(ix in unique(res1$ix)) {
    tmp1 = res1[res1$ix == ix,]
    j = length(sum_res1) + 1
    sum_res1[[j]] = do_one_acf_base(tmp1[c("year", "spd")],
                                    ix, two_step_fit, need_ts, no_detrend)
  }
  return(sum_res1)
}

do_one_acf = function(res1, two_step_fit = FALSE,
                      need_ts = FALSE, no_detrend = FALSE) {
  UseMethod("do_one_acf")
}

get_bin = function(x, bw) {
  return(bw * floor(x / bw) + bw / 2)
}

create_hist = function(tmp, var, aggr, binwidth) {
  tmp$bin = get_bin(tmp[,var], binwidth)
  tmp$count = 1
  tmp2 = aggregate(tmp["count"], by=tmp[c(aggr, "bin")], FUN=sum)
  scnt = aggregate(tmp2$count, by=tmp2[aggr], FUN=sum)
  names(scnt)[!names(scnt) %in% aggr] = "sum"
  tmp2 = merge(tmp2, scnt, by=aggr)
  tmp2$p = tmp2$count / tmp2$sum
  return(tmp2)
}


######################################################################
# 2. load results (for ACF and CV) for the C14 data

# base data -- acf
type1 = "log"
rr = 5
c14_peaks = read.csv(paste0(base_dir, "/population/acf_result/",
                            "acf_peaks_min_", type1, "_r", rr, ".csv"))
c14_peaks$bin = get_bin(c14_peaks$lag, 100)
c14_peaks$count = 1
c14_peaks = aggregate(c14_peaks["count"], by=c14_peaks["bin"], FUN=sum)
scnt1 = sum(c14_peaks$count)
c14_peaks$p = c14_peaks$count / scnt1


# base data -- CV
all_var_c14_log = read.csv(paste0(base_dir, "/population/acf_result/",
                                  "all_res_cv_", type1, "_r", rr, ".csv"))
all_var_c14_log$bin = get_bin(all_var_c14_log$cv, 0.04)
all_var_c14_log$count = 1
cv_c14 = aggregate(all_var_c14_log["count"], by=all_var_c14_log["bin"], FUN=sum)
scnt1 = sum(cv_c14$count)
cv_c14$p = cv_c14$count / scnt1




#####################################################################
# 3. use stdev from the real C14 data
popc14 = read.csv(paste0(base_dir, '/population/2021-09-08-dataset-nodup-no-sd-filter.csv'))
popc14 = popc14[!is.na(popc14$X),]
names(popc14)[12:13] = c("calBP5", "calBP95")
popc14$calBPmn = (popc14$calBP5 + popc14$calBP95) / 2

tmpsd1 = popc14$c14std[popc14$calBPmn >= 4950 & popc14$calBPmn <= 8950]
# 21904 dates are in the interval


######################################################################
# 4. main function to process the sampled dates
# for one simulation result
tmin1 = 0
tmax1 = 4000
thresh = 500

process_rand_one2 = function(fn1, aggr_ids_filt = NULL) {
  dates1 = read.table(fn1, header = FALSE)
  names(dates1) = c("dgid", "sim_year")
  years = tmin1:(tmax1-1)
  
  dates1$sd1 = sample(tmpsd1, nrow(dates1), replace = TRUE)
  dates1$date2 = dates1$sim_year + rnorm(nrow(dates1), 0, dates1$sd1)
  
  spd_all = data.frame()
  if(is.null(aggr_ids_filt)) aggr_ids_filt = unique(aggr_ids$index)
  for(ix in aggr_ids_filt) {
    tmpdates = dates1[dates1$dgid %in% aggr_ids$dgid[aggr_ids$index == ix],
                      c("date2", "sd1")]
    if(nrow(tmpdates) < thresh) next
    spd1 = vapply(years, FUN = function(x) {
      return(sum(dnorm(x, tmpdates$date2, tmpdates$sd1)))
    }, FUN.VALUE = 0)
    spd_all = rbind(spd_all, data.frame(ix = ix, year = years, spd = spd1))
  }
  return(spd_all)
}



# register helpers for parallel execution on 4 threads
# (comment out the following lines to use a single thread)
# cl = parallel::makeCluster(4) -- uncomment this on Windows
cl = parallel::makeForkCluster(4) # this likely only works on UNIX / Linux
doParallel::registerDoParallel(cl)



#####################################################################
# 5. process sampled results for simulation results with conflict
G = 80
s = 4
AA = c(1,2,5,10,20,50)
EE = c(1,5,10,20,100,200)
dir1 = paste0(base_dir, "/simulation/")
base1 = "sample_pop"

acf_peaksc = data.frame()
cv2c = data.frame()
tmp1 = foreach(A = AA) %:% foreach(E = EE) %dopar% {
  fn1 = paste0(dir1, base1, "_G", G, "_C10_E", E, "_A", A, "_s", s, ".dat")
  tmp1 = process_rand_one2(fn1)
  acf_all = do_one_acf(tmp1, TRUE)
  tmp2 = foreach(x = acf_all, .combine = rbind) %do% { return(x$all_res_min_dt) }
  tmp3 = foreach(x = acf_all, .combine = rbind) %do% {
    return(data.frame(ix = x$ix, cv = x$cv))
  }
  tmp2$A = A
  tmp2$E = E
  tmp3$A = A
  tmp3$E = E
  return(list(acf = tmp2, cv = tmp3))
}
for(i in 1:length(tmp1)) for(j in 1:length(tmp1[[i]])) {
  acf_peaksc = rbind(acf_peaksc, tmp1[[i]][[j]]$acf)
  cv2c = rbind(cv2c, tmp1[[i]][[j]]$cv)
}
save(acf_peaksc, cv2c, file=paste0(dir1, "sampled_resc.RData"))

aggr = c("A", "E")
acf_peaksc2 = create_hist(acf_peaksc, "lag", aggr, 100)
cv2c2 = create_hist(cv2c, "cv", aggr, 0.04)

acf_peaksc2$pa = 1 / acf_peaksc2$A
acf_peaksc2$pe = 1 / acf_peaksc2$E

w=7.6; h=4.5 # figure size -- 6x6 panels
p1 = ggplot(mapping = aes(x=bin, y=p))
p1 = p1 + geom_col(data = c14_peaks, color="grey", fill="grey")
p1 = p1 + geom_col(data = acf_peaksc2, color="red", fill="red", alpha=0.45)
p1 = p1 + theme_bw(6) + xlab("Lag [years]") + ylab("Relative frequency")
p1 = p1 + scale_x_continuous(limits=c(0,2200))
p1 = p1 + facet_grid(pa~pe, labeller = label_bquote(rows = p[A] == .(pa),
                                                    cols = p[E] == .(pe)))
fn1 = paste0(dir1, "acf_peaks_sample_pop")
ggsave(paste0(fn1, ".pdf"), p1, width=w, height=h)
ggsave(paste0(fn1, ".png"), p1, width=w, height=h, dpi=300)

# separate plot for the A = 5 and E = 5 case for the main figure
# in the paper
tmp1 = acf_peaksc2[acf_peaksc2$A == 5 & acf_peaksc2$E == 5,]
p1 = ggplot(mapping = aes(x=bin, y=p))
p1 = p1 + geom_col(data = c14_peaks, color="grey", fill="grey")
p1 = p1 + geom_col(data = tmp1, color="red", fill="red", alpha=0.45)
p1 = p1 + theme_bw(6) + xlab("Lag [years]") + ylab("Relative frequency")
p1 = p1 + scale_x_continuous(limits=c(0,2200))

ggsave(paste0(dir1, "acf_peaks_sample_pop_A5_E5.pdf"), p1,
       width=1.6, height=1.2)


cv2c2$pa = 1 / cv2c2$A
cv2c2$pe = 1 / cv2c2$E

p1 = ggplot(mapping = aes(x=bin, y=p))
p1 = p1 + geom_col(data = cv_c14, color="grey", fill="grey")
p1 = p1 + geom_col(data = cv2c2, color="red", fill="red", alpha=0.45)
p1 = p1 + theme_bw(6) + xlab("Coefficient of variation") +
  ylab("Relative frequency")
p1 = p1 + scale_x_continuous(limits=c(0,1.5))
p1 = p1 + facet_grid(pa~pe, labeller = label_bquote(rows = p[A] == .(pa),
                                                    cols = p[E] == .(pe)))
fn1 = paste0(dir1, "cv_sample_pop")
ggsave(paste0(fn1, ".pdf"), p1, width=w, height=h)
ggsave(paste0(fn1, ".png"), p1, width=w, height=h, dpi=300)

tmp1 = cv2c2[cv2c2$A == 5 & cv2c2$E == 5,]
p1 = ggplot(mapping = aes(x=bin, y=p))
p1 = p1 + geom_col(data = cv_c14, color="grey", fill="grey")
p1 = p1 + geom_col(data = tmp1, color="red", fill="red", alpha=0.45)
p1 = p1 + theme_bw(6) + xlab("Coefficient of variation") +
  ylab("Relative frequency")
p1 = p1 + scale_x_continuous(limits=c(0,1.5))

ggsave(paste0(dir1, "cv_sample_pop_A5_E5.pdf"), p1,
       width=1.6, height=1.2)


#####################################################################
# 6. process sampled results for simulation results with conflict
G = c(20,30,40)
s = 1:4
dir1 = paste0(base_dir, "/simulation/")
b = "sample_popn"
tmp1 = foreach(G1 = G) %:% foreach(s1 = s) %dopar% {
  fn1 = paste0(dir1, b, "_G", G1, "C10ws", s1, ".dat")
  tmp1 = process_rand_one2(fn1)
  acf_all = do_one_acf(tmp1, TRUE)
  tmp2 = foreach(x = acf_all, .combine = rbind) %do% { return(x$all_res_min_dt) }
  tmp3 = foreach(x = acf_all, .combine = rbind) %do% {
    return(data.frame(ix = x$ix, cv = x$cv))
  }
  tmp2$G = G1
  tmp2$s = s1
  tmp3$G = G1
  tmp3$s = s1
  return(list(acf = tmp2, cv = tmp3))
}
acf_peaksn = data.frame()
cv2n = data.frame()
for(i in 1:length(tmp1)) for(j in 1:length(tmp1[[i]])) {
  acf_peaksn = rbind(acf_peaksn, tmp1[[i]][[j]]$acf)
  cv2n = rbind(cv2n, tmp1[[i]][[j]]$cv)
}
save(acf_peaksn, cv2n, file=paste0(dir1, "sampled_resn.RData"))

aggr = c("G", "s")
acf_peaksn2 = create_hist(acf_peaksn, "lag", aggr, 100)
cv2n2 = create_hist(cv2n, "cv", aggr, 0.04)

w=4.7; h=2.37; # figure size -- 4x3 panels
p1 = ggplot(mapping = aes(x=bin, y=p))
p1 = p1 + geom_col(data = c14_peaks, color="grey", fill="grey")
p1 = p1 + geom_col(data = acf_peaksn2, color="red", fill="red", alpha=0.45)
p1 = p1 + theme_bw(6) + xlab("Lag [years]") + ylab("Relative frequency")
p1 = p1 + scale_x_continuous(limits=c(0,2200))
p1 = p1 + facet_grid(G~s, labeller = label_bquote(rows = G == .(G),
                                                  cols = s == .(s)))
fn1 = paste0(dir1, "acf_peaks_samplen_pop")
ggsave(paste0(fn1, ".pdf"), p1, width=w, height=h)
ggsave(paste0(fn1, ".png"), p1, width=w, height=h, dpi=300)

p1 = ggplot(mapping = aes(x=bin, y=p))
p1 = p1 + geom_col(data = cv_c14, color="grey", fill="grey")
p1 = p1 + geom_col(data = cv2n2, color="red", fill="red", alpha=0.45)
p1 = p1 + theme_bw(6) + xlab("Coefficient of variation") +
  ylab("Relative frequency")
p1 = p1 + scale_x_continuous(limits=c(0,1.5))
p1 = p1 + facet_grid(G~s, labeller = label_bquote(rows = G == .(G),
                                                  cols = s == .(s)))
fn1 = paste0(dir1, "cv_samplen_pop")
ggsave(paste0(fn1, ".pdf"), p1, width=w, height=h)
ggsave(paste0(fn1, ".png"), p1, width=w, height=h, dpi=300)

tmp1 = acf_peaksn2[acf_peaksn2$G == 40 & acf_peaksn2$s == 4,]
p1 = ggplot(mapping = aes(x=bin, y=p))
p1 = p1 + geom_col(data = c14_peaks, color="grey", fill="grey")
p1 = p1 + geom_col(data = tmp1, color="red", fill="red", alpha=0.45)
p1 = p1 + theme_bw(6) + xlab("Lag [years]") + ylab("Relative frequency")
p1 = p1 + scale_x_continuous(limits=c(0,2200))
ggsave(paste0(dir1, "acf_peaks_samplen_pop_G40_s4.pdf"), p1,
       width=1.6, height=1.2)

tmp1 = cv2n2[cv2n2$G == 40 & cv2n2$s == 4,]
p1 = ggplot(mapping = aes(x=bin, y=p))
p1 = p1 + geom_col(data = cv_c14, color="grey", fill="grey")
p1 = p1 + geom_col(data = tmp1, color="red", fill="red", alpha=0.45)
p1 = p1 + theme_bw(6) + xlab("Coefficient of variation") +
  ylab("Relative frequency")
p1 = p1 + scale_x_continuous(limits=c(0,1.5))
ggsave(paste0(dir1, "cv_samplen_pop_G40_s4.pdf"), p1,
       width=1.6, height=1.2)



#####################################################################
# 7. main figures for the paper (Fig. 1 -- Fig. 3)
# + SI Figs. S9-S10
# these show individual time series (SPDs and population numbers)
# in a set of tiles in one tiling scheme (dx = 0 and dy = 0 position)

# 7.1. preparation -- load the data sets, re-create the results needed
# in these cases

G = 80
s = 4
A = 5
E = 5
dir1 = paste0(base_dir, "/simulation/")
b1 = "sample_pop"

fn1 = paste0(dir1, b1, "_G", G, "_C10_E", E, "_A", A, "_s", s, ".dat")
tmp1 = process_rand_one2(fn1, tmpix)
acf_all = do_one_acf(tmp1, TRUE)
acf_all = foreach(x = acf_all, .combine = rbind) %do% { return(x$acfdt) }

spdc = tmp1
acfc = acf_all

# load the original results (i.e. without the sampling process)
fn1 = paste0(dir1, "/runmr_G", G, "_C10_E", E, "_A", A, "_s", s, "spaggr.out.gz")
res1 = read.table(gzfile(fn1))
names(res1) = c("year", "ix", "cnt", "spd")
res1 = res1[res1$ix %in% tmpix,]
spdc0 = res1
acfc0 = do_one_acf(res1, TRUE)
acfc0 = foreach(x = acfc0, .combine = rbind) %do% { return(x$acfdt) }


# simulation results without conflict
G = 40
s = 4
dir1 = paste0(base_dir, "/simulation/")
b1 = "sample_popn"

fn1 = paste0(dir1, b1, "_G", G, "C10ws", s, ".dat")
tmp1 = process_rand_one2(fn1, tmpix)
acf_all = do_one_acf(tmp1, TRUE)
acf_all = foreach(x = acf_all, .combine = rbind) %do% { return(x$acfdt) }

spdn = tmp1
acfn = acf_all

fn1 = paste0(dir1, "/runG", G, "C10ws", s, "spaggr.out.gz")
res1 = read.table(gzfile(fn1))
names(res1) = c("year", "ix", "cnt", "spd")
res1 = res1[res1$ix %in% tmpix,]
spdn0 = res1
acfn0 = do_one_acf(res1, TRUE)
acfn0 = foreach(x = acfn0, .combine = rbind) %do% { return(x$acfdt) }


# C14 data
load(paste0(base_dir, "/population/c14_spd_new_9000_5000.RData"))
res3 = Filter(function(x) (x$dx == 0 & x$dy == 0 & x$r == 500), all_res)
res3 = res3[[1]]$spd_cmb
names(res3) = c("year", "pop", "ID")
res3$year = 9000 - res3$year
tmpregions = regions[regions$res == 500 & regions$dx == 0 & regions$dy == 0,]
res3 = merge(res3, tmpregions, by.x = "ID", by.y="ID")

spdc14 = res3[c("index", "year", "pop")]
names(spdc14)[1] = "ix"
names(spdc14)[3] = "spd"

acf_all = do_one_acf(spdc14, TRUE)
acf_all = foreach(x = acf_all, .combine = rbind) %do% { return(x$acfdt) }
acfc14 = acf_all


ix1 = unique(acfc$ix) # 14
ix2 = unique(acfn$ix) # 19
ix3 = unique(acfc14$ix) # 17

sum(ix1 %in% ix2) # 14
sum(ix1 %in% ix3) # 10
sum(ix2 %in% ix3) # 13

tmpix2 = ix1[ix1 %in% ix3]

colors = c("#5aa732","#005c96","#e75216", "#009ace","#ffd500",
           "#8671a7", "#e78502","#db6aa2","#007939","#9a2846","#815234","#000000")


################################################################
# 7.2. Fig. 1 -- SPD time series and ACFs for the C14 data

tmp1 = spdc14[spdc14$ix %in% tmpix2,]
tmp1$ix = factor(tmp1$ix, levels = tmpix2)

# Fig 1a (SPD of C14 data)
p1 = ggplot(tmp1) + geom_line(aes(x=7050-year, y=spd, group=ix,
                                  color=ix), size=0.2)
p1 = p1 + scale_x_reverse() + theme_bw(6)
p1 = p1 + xlab("Date [BCE]") + ylab("SPD (arbitrary units)")
p1 = p1 + scale_color_manual(values=colors) + guides(color=FALSE)
fn1 = paste0(dir1, "spdc14")
ggsave(paste0(fn1, ".pdf"), width=1.6, height=1.2)


# Fig 1b (ACF of C14 data)
tmp1 = acfc14[acfc14$ix %in% tmpix2,]
tmp1$ix = factor(tmp1$ix, levels = tmpix2)

p1 = ggplot(tmp1) + geom_line(aes(x=lag, y=acf, group=ix,
                                  color=ix), size=0.2)
p1 = p1 + theme_bw(6) + xlab("Lag [years]") + ylab("ACF")
p1 = p1 + scale_color_manual(values=colors) + guides(color=FALSE)
fn1 = paste0(dir1, "acfc14")
ggsave(paste0(fn1, ".pdf"), width=1.6, height=1.2)



##################################################################
# 7.3. Fig. 2 -- sampled population and ACF for the simulation
# variant without conflict
tmp1 = spdn[spdn$ix %in% tmpix2,]
tmp1$ix = factor(tmp1$ix, levels = tmpix2)

# Fig. 2a
p1 = ggplot(tmp1) + geom_line(aes(x=7050-year, y=spd, group=ix,
                                  color=ix), size=0.2)
p1 = p1 + scale_x_reverse() + theme_bw(6)
p1 = p1 + xlab("Date [BCE]") + ylab("Sampled population")
p1 = p1 + scale_color_manual(values=colors) + guides(color=FALSE)
fn1 = paste0(dir1, "spdn")
ggsave(paste0(fn1, ".pdf"), width=1.6, height=1.2)


# Fig. 2b
tmp1 = acfn[acfn$ix %in% tmpix2,]
tmp1$ix = factor(tmp1$ix, levels = tmpix2)

p1 = ggplot(tmp1) + geom_line(aes(x=lag, y=acf, group=ix,
                                  color=ix), size=0.2)
p1 = p1 + theme_bw(6) + xlab("Lag [years]") + ylab("ACF")
p1 = p1 + scale_color_manual(values=colors) + guides(color=FALSE)
fn1 = paste0(dir1, "acfn")
ggsave(paste0(fn1, ".pdf"), width=1.6, height=1.2)



###################################################################
# 7.4. Fig 3 -- sampled population and ACF for the simulation
# variant including violent conflict

# Fig 3a
tmp1 = spdc[spdc$ix %in% tmpix2,]
tmp1$ix = factor(tmp1$ix, levels = tmpix2)

p1 = ggplot(tmp1) + geom_line(aes(x=7050-year, y=spd, group=ix,
                                  color=ix), size=0.2)
p1 = p1 + scale_x_reverse() + theme_bw(6)
p1 = p1 + xlab("Date [BCE]") + ylab("Sampled population")
p1 = p1 + scale_color_manual(values=colors) + guides(color=FALSE)
fn1 = paste0(dir1, "spdc")
ggsave(paste0(fn1, ".pdf"), width=1.6, height=1.2)


# Fig 3b
tmp1 = acfc[acfc$ix %in% tmpix2,]
tmp1$ix = factor(tmp1$ix, levels = tmpix2)

p1 = ggplot(tmp1) + geom_line(aes(x=lag, y=acf, group=ix,
                                  color=ix), size=0.2)
p1 = p1 + theme_bw(6) + xlab("Lag [years]") + ylab("ACF")
p1 = p1 + scale_color_manual(values=colors) + guides(color=FALSE)
fn1 = paste0(dir1, "acfc")
ggsave(paste0(fn1, ".pdf"), width=1.6, height=1.2)



#####################################################################
# 8. faceted plots (SI)
# For the simulation results, we display the resulting population
# numbers and the result of the sampling process together

ix1 = unique(spdc$ix)
ix12 = unique(spdc0$ix)
ix2 = unique(spdn$ix)
ix22 = unique(spdn0$ix)

ix13 = unique(acfc$ix)
ix14 = unique(acfc0$ix)
ix23 = unique(acfn$ix)
ix24 = unique(acfn0$ix)

create_pop_plot = function(tmp1, tmp2) {
  
  mn1 = aggregate(tmp1["spd"], by=tmp1["ix"], FUN=sum)
  mn2 = aggregate(tmp2["spd"], by=tmp2["ix"], FUN=sum)
  mn1 = merge(mn1, mn2, by="ix")
  mn1$r = mn1$spd.x / mn1$spd.y
  tmp2 = merge(tmp2, mn1[c("ix", "r")], by="ix")
  tmp2$spd = tmp2$spd * tmp2$r
  tmp1 = rbind(tmp1, tmp2[c("ix", "year", "spd", "type")])
  tmp1 = tmp1[tmp1$ix != 12,] # outlier region
  
  p1 = ggplot(tmp1) + geom_line(aes(x=7000-year, y=spd / 1000, group=type,
                                    color=type), size=0.2)
  p1 = p1 + facet_wrap(~ix, ncol = 4) + theme_bw(6) + scale_x_reverse()
  p1 = p1 + xlab("Date [BCE]") + ylab("Population [x1000]")
  p1 = p1 + theme(strip.background = element_blank(), 
                  strip.text.x = element_blank())
  p1 = p1 + scale_color_manual(values=c("black", "red"))
  p1 = p1 + theme(legend.position = c(0.84, 0.12))
  p1 = p1 + labs(color="")
  p1 = p1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  return(p1)
}

create_acf_plot = function(tmp1) {
  p1 = ggplot(tmp1) + geom_line(aes(x=lag, y=acf, group=type, color=type),
                                size=0.2)
  p1 = p1 + facet_wrap(~ix, ncol = 4) + theme_bw(6)
  p1 = p1 + xlab("Lag [years]") + ylab("ACF")
  p1 = p1 + theme(strip.background = element_blank(), 
                  strip.text.x = element_blank())
  p1 = p1 + scale_color_manual(values=c("black", "red"))
  p1 = p1 + theme(legend.position = c(0.84, 0.12))
  p1 = p1 + labs(color="")
  p1 = p1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  return(p1)
}

#####################################################################
# 8.1. Fig S10: simulation with conflict: 16 cases
tmp1 = spdc0[spdc0$ix %in% ix1, c("ix", "year", "spd")]
tmp1$type = "original"
tmp1 = tmp1[tmp1$year < 4000,]
tmp2 = spdc
tmp2$type = "sampled"

p1 = create_pop_plot(tmp1, tmp2)

fn1 = paste0(dir1, "sampled_pop_ex_G80_A5_E5_s4")
ggsave(paste0(fn1, ".pdf"), p1, width=3.2, height=2.6)
# adjusted legend position in PDF (with Inkscape): 240, 192

tmp1 = acfc0[acfc0$ix %in% ix13, c("ix", "lag", "acf")]
tmp1$type = "original"
tmp2 = acfc
tmp2$type = "sampled"
tmp1 = rbind(tmp1, tmp2)
tmp1 = tmp1[!(tmp1$ix %in% c(12,25)),]

p1 = create_acf_plot(tmp1)

fn1 = paste0(dir1, "sampled_acf_ex_G80_A5_E5_s4")
ggsave(paste0(fn1, ".pdf"), p1, width=3.2, height=2.6)


####################################################################
# 8.2. Fig S11: simulation without conflict
tmp1 = spdn0[spdn0$ix %in% ix2, c("ix", "year", "spd")]
tmp1$type = "original"
tmp1 = tmp1[tmp1$year < 4000,]
tmp2 = spdn
tmp2$type = "sampled"

p1 = create_pop_plot(tmp1, tmp2)

fn1 = paste0(dir1, "sampled_pop_exn_G40_s4")
ggsave(paste0(fn1, ".pdf"), p1, width=3.2, height=3.17)


tmp1 = acfn0[acfn0$ix %in% ix23, c("ix", "lag", "acf")]
tmp1$type = "original"
tmp2 = acfn
tmp2$type = "sampled"
tmp1 = rbind(tmp1, tmp2)
tmp1 = tmp1[tmp1$ix != 12,]

p1 = create_acf_plot(tmp1)

fn1 = paste0(dir1, "sampled_acfn_ex_G40_s4")
ggsave(paste0(fn1, ".pdf"), p1, width=3.2, height=3.17)



