# acf_new_aggr_plot.r -- process the results of ACF in the new 
#   aggregation, create plots

library(ggplot2)
library(foreach)
library(doParallel) # optionally, speed up some computations with parallel processing

base_dir = "~/CSH/HoloSim/data"


##################################################################
# 1. function definitions for calculating ACFs

# read base data for spatial aggregation used
aggr_key = read.csv(paste0(base_dir, '/new_grid_ids.csv'))
# note: we only process the 500 km linear tiling size in this script
# (omit the following line to change this)
aggr_key = aggr_key[aggr_key$res == 500,]
aggr_types = unique(aggr_key[c("res", "dx", "dy")])



# helper for the below functions
get_diff_sign = function(acf1) {
  tmps = sign(diff(acf1))
  for(i in length(tmps):2) {
    if(tmps[i-1] == 0) tmps[i-1] = tmps[i]
  }
  return(tmps)
}

# get the location of the first (negative) minimum
get_acf_min = function(acf1) {
  if(length(acf1) < 3) return(NULL)
  tmp2n = which(diff(get_diff_sign(acf1)) == 2) + 1
  tmp2n = tmp2n[acf1[tmp2n] < 0]
  return(head(tmp2n, 1))
}


# do the full analysis for one aggregated simulation result
# using a logistic detrending
do_one_acf = function(res1, tmin1, tmax1, two_step_fit = FALSE, T1 = 100) {
  sum_res1 = list()
  rr = unique(aggr_types$res)
  for(i in 1:length(rr)) { # note: length(rr) == 3
    sum_res1[[i]] = list()
    sum_res1[[i]]$r = rr[i]
    sum_res1[[i]]$all_res_min_dt = data.frame()
    sum_res1[[i]]$errors = data.frame()
  }
  
  for(i in 1:nrow(aggr_types)) {
    r = aggr_types$res[i]
    j = 1
    while(1) {
      if(j > length(sum_res1)) stop("Invalid r!")
      if(sum_res1[[j]]$r == r) break
      j = j + 1
    }
    
    dx = aggr_types$dx[i]
    dy = aggr_types$dy[i]
    regions = aggr_key[aggr_key$dx == dx & aggr_key$dy == dy & aggr_key$res == r,]
    
    tmp1 = merge(regions[c("index", "ID")], res1, by="index")
    
    res1nd = NULL
    res2nd = NULL
    
    all_years = tmin1:(tmax1-1)
    
    for(ID1 in unique(tmp1$ID)) {
      tmp2 = tmp1[tmp1$ID == ID1 & tmp1$year >= tmin1 & tmp1$year < tmax1,
                  c("year", "pop")]
      # fill in missing values
      tmpix = !(all_years %in% tmp2$year)
      if(sum(tmpix) > 0) {
        tmp2m = data.frame(year = all_years[!(all_years %in% tmp2$year)],
                           pop = 0)
        tmp2 = rbind(tmp2, tmp2m)
      }
      tmp2 = tmp2[order(tmp2$year),]
      
      # detrending
      N = mean(tail(tmp2$pop, 1000))
      x0 = min(tmp2$year[tmp2$pop >= N/2])
      nls1 = tryCatch( {
        if(two_step_fit) {
          nls1 = nls(pop ~ N / (1 + exp(-(year - x0)/T1)) , tmp2,
                     list(N = N))
          N = nls1$m$getPars()[1]
          nls1 = nls(pop ~ N / (1 + exp(-(year - x0)/T1)) , tmp2,
                     list(x0 = x0, T1 = T1))
        } else {
          nls1 = nls(pop ~ N / (1 + exp(-(year - x0)/T1)) , tmp2,
                     list(N = N, x0 = x0, T1 = T1))
        }
        nls1
      }, error = function(x) {
        tmperr = data.frame(i = i, res = r, dx = dx, dy = dy, ID = ID1)
        sum_res1[[j]]$errors = rbind(sum_res1[[j]]$errors, tmperr)
        NULL
      }
      )
      if(is.null(nls1)) next
      
      tmp2$pred = predict(nls1)
      tmp2$diff = tmp2$pop - tmp2$pred
      
      # do the ACF
      acf1 = acf(tmp2$diff, lag.max = 4999, plot = FALSE)
      
      if(sum(is.nan(acf1$acf)) < (tmax1 - tmin1)) {
        tmp2n = get_acf_min(acf1$acf)
        res1nd = c(res1nd, acf1$lag[tmp2n])
        res2nd = c(res2nd, acf1$acf[tmp2n])
      }
      
    }
    
    if(length(res1nd) > 0) {
      tmp1 = data.frame(lag = res1nd, val = res2nd, dx = dx, dy = dy)
      sum_res1[[j]]$all_res_min_dt = rbind(sum_res1[[j]]$all_res_min_dt, tmp1)
    }
  }
  
  return(sum_res1)
}

cnts_thresh = 10 # mininum number of cells to include in the results
# read one simulation result, already aggregated to the tiles used here
read_one_file = function(fn1) {
  res1 = read.table(gzfile(fn1))
  names(res1) = c("year", "index", "cnt", "pop")
  
  # filter out regions with too few cells
  tmp1 = aggregate(res1["cnt"], by=res1["index"], FUN=max)
  tmp1 = tmp1[tmp1$cnt >= cnts_thresh,]
  res1 = res1[res1$index %in% tmp1$index, c("year", "index", "pop")]
  return(res1)
}

# simple helper to create histograms
binwidth = 100
get_bin = function(x, bw) {
  return(bw * floor(x / bw) + bw / 2)
}

# aggregate in a data.frame by a set of variables
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


##################################################################
# 2. read the base data (C14 data results, r = 500 km tiling size)
rr1 = 5
type1 = "log"
c14_peaks = read.csv(paste0(base_dir, "/population/acf_result/",
                              "acf_peaks_min_", type1, "_r", rr1, ".csv"))
c14_peaks$bin = get_bin(c14_peaks$lag, binwidth)
c14_peaks$count = 1
c14_peaks = aggregate(c14_peaks["count"], by=c14_peaks["bin"], FUN=sum)
  
# create a normalized version as well
scnt1 = sum(c14_peaks$count)
c14_peaks$p = c14_peaks$count / scnt1
  

######################################################################
# 3. ACFs for two example cases
tmin1 = 0 # note: we filter for the first 4000 year of the simulation
tmax1 = 4000 # which covers 7000 to 3000 BCE

# 3.1. simulation results with conflict
G1 = 80
s1 = 4
dir1 = "simulation"
base_fn = "runmr"
E1 = 20
A1 = 5

fn01 = paste0(base_dir, "/", dir1, "/", base_fn, "_G", G1, "_C10_E", E1,
             "_A", A1, "_s", s1)

# 3.2. simulation results without conflict
G1 = 40
base_fn = "run"
fn02 = paste0(base_dir, "/", dir1, "/", base_fn, "G", G1, "C10ws", s1)

# process both and create plots
for(fn1 in c(fn01, fn02)) {
  res1 = read_one_file(paste0(fn1, "spaggr.out.gz"))
  acf_res = do_one_acf(res1, tmin1, tmax1, TRUE)
  
  # plot the results together with the C14 data baseline
  tmp1 = acf_res[[1]]$all_res_min_dt
  tmp1$bin = get_bin(tmp1$lag, binwidth)
  tmp1$count = 1
  tmp1 = aggregate(tmp1["count"], by=tmp1["bin"], FUN=sum)
  scnt1 = sum(tmp1$count)
  tmp1$p = tmp1$count / scnt1
  
  p1 = ggplot(mapping = aes(x=bin, y=p))
  p1 = p1 + geom_col(data=c14_peaks, color="grey", fill="grey")
  p1 = p1 + geom_col(data=tmp1, color="red", fill="red", alpha=0.45)
  p1 = p1 + theme_bw(6) + xlab("Lag [years]") + ylab("Relative frequency")
  p1 = p1 + scale_x_continuous(limits=c(0,2200))
  ggsave(paste0(fn1, "_acf.pdf"), p1, width=3.2, height=2)
  ggsave(paste0(fn1, "_acf.png"), p1, width=3.2, height=2, dpi=300)
}


#######################################################################
# 4. faceted plots for a set of cases together
# 4.1. simulation without conflict, 3x4 cases (G = 20,30,40; s = 1-4)
# 4.2. simulation with conflict, G = 80, s = 4, 6x6 cases
#  (variation in p_A and p_E)

# register helpers for parallel execution on 4 threads
# (comment out the following lines to use a single thread)
# cl = parallel::makeCluster(4) -- uncomment this on Windows
cl = parallel::makeForkCluster(4) # this likely only works on UNIX / Linux
doParallel::registerDoParallel(cl)



# 4.1. simulation without conflict
G1 = c(20,30,40)
s1 = 1:4
base_fn = "run"

all_res_n = foreach(G = G1, .combine = rbind) %:%
  foreach(s = s1, .combine = rbind) %dopar% {
    fn1 = paste0(base_dir, "/", dir1, "/", base_fn, "G", G, "C10ws", s)
    res1 = read_one_file(paste0(fn1, "spaggr.out.gz"))
    acf_res = do_one_acf(res1, tmin1, tmax1, TRUE)
    # note: we only calculate results for the r = 500 km tiling resolution
    # (filtering is done in do_one_acf()); this could be modified if needed
    tmp1 = acf_res[[1]]$all_res_min_dt
    tmp1$G = G
    tmp1$s = s
    return(tmp1)
}

all_res_nh = create_hist(all_res_n, "lag", c("G", "s"), 100)

p1 = ggplot(mapping = aes(x=bin, y=p))
p1 = p1 + geom_col(data=c14_peaks, color="grey", fill="grey")
p1 = p1 + geom_col(data=all_res_nh, aes(x=bin, y=p), color="red",
                   fill="red", alpha=0.45)
p1 = p1 + theme_bw(6) + xlab("Lag [years]") + ylab("Relative frequency")
p1 = p1 + scale_x_continuous(limits=c(0,2200))
p1 = p1 + facet_grid(G~s, labeller = labeller(s = function(x){
  return(paste0("s =  ", x))}, G = function(x) {
    return(paste0("G = ", x, " km"))}))

fn1 = paste0(base_dir, "/", dir1, "/acf_cmb_n")
ggsave(paste0(fn1, ".pdf"), p1, width=4.7, height=2.37)
ggsave(paste0(fn1, ".png"), p1, width=4.7, height=2.37, dpi=300)



# 4.2. simulation with conflict
G1 = 80
s1 = 4
dir1 = "simulation"
base_fn = "runmr"
E1 = c(1,5,10,20,100,200) # 1 / p_E parameter
A1 = c(1,2,5,10,20,50) # 1 / p_a parameter

all_res_c = foreach(E = E1, .combine = rbind) %:%
  foreach(A = A1, .combine = rbind) %dopar% {
    fn1 = paste0(base_dir, "/", dir1, "/", base_fn, "_G", G1, "_C10_E",
                 E, "_A", A, "_s", s1)
    res1 = read_one_file(paste0(fn1, "spaggr.out.gz"))
    acf_res = do_one_acf(res1, tmin1, tmax1, TRUE)
    # note: we only calculate results for the r = 500 km tiling resolution
    # (filtering is done in do_one_acf()); this could be modified if needed
    tmp1 = acf_res[[1]]$all_res_min_dt
    tmp1$A = A
    tmp1$E = E
    return(tmp1)
}

all_res_ch = create_hist(all_res_c, "lag", c("A", "E"), 100)
all_res_ch$pa = 1 / all_res_ch$A
all_res_ch$pe = 1 / all_res_ch$E

p1 = ggplot(mapping = aes(x=bin, y=p))
p1 = p1 + geom_col(data=c14_peaks, color="grey", fill="grey")
p1 = p1 + geom_col(data=all_res_ch, aes(x=bin, y=p), color="red",
                   fill="red", alpha=0.45)
p1 = p1 + theme_bw(6) + xlab("Lag [years]") + ylab("Relative frequency")
p1 = p1 + scale_x_continuous(limits=c(0,2200))
p1 = p1 + facet_grid(pa~pe, labeller = label_bquote(rows = p[A] == .(pa),
                                                    cols = p[E] == .(pe)))

fn1 = paste0(base_dir, "/", dir1, "/acf_cmb_c")
ggsave(paste0(fn1, ".pdf"), p1, width=7.6, height=4.5)
ggsave(paste0(fn1, ".png"), p1, width=7.6, height=4.5, dpi=300)




