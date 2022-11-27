# acf_new_aggr_rep.r -- parallel analysis for processing new aggregated
#  simulation results; version for repeated simulation realizations

library(foreach)
library(doParallel)


setwd('/mnt/4TBSSD/kondor')

aggr_key = read.csv('new_grid_ids.csv')
aggr_regions = read.csv('new_grids_all.csv')

aggr_types = unique(aggr_key[c("res", "dx", "dy")])

ncores = 50 # maximum number of CPU cores to use

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


# do the full analysis for one aggregated simulation result
do_one_acf = function(res1, tmin1, tmax1, two_step_fit = FALSE) {
  sum_res1 = list()
  rr = unique(aggr_types$res)
  for(i in 1:length(rr)) { # note: length(rr) == 3
    sum_res1[[i]] = list()
    sum_res1[[i]]$r = rr[i]
    sum_res1[[i]]$all_res_min = data.frame()
    sum_res1[[i]]$all_res_max = data.frame()
    sum_res1[[i]]$all_res_min_dt = data.frame()
    sum_res1[[i]]$all_res_max_dt = data.frame()
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
    
    res1p = NULL
    res2p = NULL
    res1n = NULL
    res2n = NULL
    
    res1pd = NULL
    res2pd = NULL
    res1nd = NULL
    res2nd = NULL
    
    cv1 = data.frame()
    
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
      
      acf1 = acf(tmp2$pop, lag.max = 4999, plot = FALSE)
      
      if(sum(is.nan(acf1$acf)) < (tmax1 - tmin1)) {
        tmp2p = get_acf_max(acf1$acf)
        tmp2n = get_acf_min(acf1$acf)
        res1p = c(res1p, acf1$lag[tmp2p])
        res2p = c(res2p, acf1$acf[tmp2p])
        res1n = c(res1n, acf1$lag[tmp2n])
        res2n = c(res2n, acf1$acf[tmp2n])
      }
      
      # detrending
      N = mean(tail(tmp2$pop, 1000))
      x0 = min(tmp2$year[tmp2$pop >= N/2])
      T1 = 100 # TODO: maybe adjust this?
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
      cv0 = sd(tmp2$diff) / mean(tmp2$pred)
      cv1 = rbind(cv1, data.frame(ID = ID1, cv = cv0, dx = dx, dy = dy))
      
      # do the ACF
      acf1 = acf(tmp2$diff, lag.max = 4999, plot = FALSE)
      
      if(sum(is.nan(acf1$acf)) < (tmax1 - tmin1)) {
        tmp2p = get_acf_max(acf1$acf)
        tmp2n = get_acf_min(acf1$acf)
        res1pd = c(res1pd, acf1$lag[tmp2p])
        res2pd = c(res2pd, acf1$acf[tmp2p])
        res1nd = c(res1nd, acf1$lag[tmp2n])
        res2nd = c(res2nd, acf1$acf[tmp2n])
      }
      
    }
    
    if(length(res1p) > 0) {
      tmp1 = data.frame(lag = res1p, val = res2p, dx = dx, dy = dy)
      sum_res1[[j]]$all_res_max = rbind(sum_res1[[j]]$all_res_max, tmp1)
    }
    if(length(res1n) > 0) {
      tmp1 = data.frame(lag = res1n, val = res2n, dx = dx, dy = dy)
      sum_res1[[j]]$all_res_min = rbind(sum_res1[[j]]$all_res_min, tmp1)
    }
    
    if(length(res1pd) > 0) {
      tmp1 = data.frame(lag = res1pd, val = res2pd, dx = dx, dy = dy)
      sum_res1[[j]]$all_res_max_dt = rbind(sum_res1[[j]]$all_res_max_dt, tmp1)
    }
    if(length(res1nd) > 0) {
      tmp1 = data.frame(lag = res1nd, val = res2nd, dx = dx, dy = dy)
      sum_res1[[j]]$all_res_min_dt = rbind(sum_res1[[j]]$all_res_min_dt, tmp1)
    }
    
    sum_res1[[j]]$cv = rbind(sum_res1[[j]]$cv, cv1)
  }
  
  return(sum_res1)
}



# command line arguments: base filename, number of realizations
args = commandArgs(trailingOnly=TRUE)
fnbase = args[1]
N = as.integer(args[2])


cl <- parallel::makeForkCluster(max(ncores,N))
doParallel::registerDoParallel(cl)

tmin1 = 0
tmax1 = 4000 # match for start date of 7000 BCE
y1 = 7000 # only one start date
cnts_thresh = 10 # mininum number of cells to include in the results

all_res = foreach(i = 1:N) %dopar% {
  fn1 = paste0(fnbase, "_", i, ".out.gz")
  res1 = read.table(gzfile(fn1))
  names(res1) = c("year", "index", "cnt", "pop")
  
  # filter out regions with too few cells
  tmp1 = aggregate(res1["cnt"], by=res1["index"], FUN=max)
  tmp1 = tmp1[tmp1$cnt >= cnts_thresh,]
  res1 = res1[res1$index %in% tmp1$index, c("year", "index", "pop")]
  
  sum_res1 = tryCatch( {
    do_one_acf(res1, tmin1, tmax1, TRUE)
  }, error = function(x) { return(NULL) })
    
  all_res1 = list()
  all_res1$res = sum_res1
  all_res1$i = i
  return(all_res1)
}
save(all_res, file = paste0(fnbase, "_all_res.RData"))
rm(all_res)



