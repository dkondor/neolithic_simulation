# c14_data_acf3.r -- analysis of SPDs of C14 data:
# detrending, calculation of ACFs, identification minimum locations,
# and summary plots

library(ggplot2)
library(foreach)

# base directory -- need to be set to the location of the input data!
base_dir = "~/CSH/HoloSim/data"
datadir = paste0(base_dir, "/population")
# output directory for the identified peaks and plots
outdir = paste0(datadir, "/acf_result")

# read the tiling used for spatial aggregation
aggr_regions = read.csv(paste0(base_dir, '/new_grids_all.csv'))


###################################################################
# 1. Function definitions

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

# create a histogram of all results
do_one_acf = function(all_res, two_step_fit = FALSE, T1 = 300) {
  # all resolutions in the input
  rr = unique(unlist(lapply(all_res, FUN=function(x) x$r)))
  sum_res1 = list()
  for(i in 1:length(rr)) { # note: length(rr) == 3
    sum_res1[[i]] = list()
    sum_res1[[i]]$r = rr[i]
    sum_res1[[i]]$errors = data.frame()
  }
  
  for(i in 1:length(rr)) {
    r = sum_res1[[i]]$r # current resolution
    
    # results without detrending
    all_res_min = data.frame()
    all_res_max = data.frame()
    
    # results with exponential detrending
    all_res_exp_min = data.frame()
    all_res_exp_max = data.frame()
    
    # results with logistic detrending
    all_res_log_min = data.frame()
    all_res_log_max = data.frame()
    
    # results for the coefficient of variation
    all_res_cv_exp = data.frame()
    all_res_cv_log = data.frame()
    
    for(j in 1:length(all_res)) if(all_res[[j]]$r == r) {
      tmpres1 = all_res[[j]]$spd_cmb
      
      res1p = data.frame()
      res1n = data.frame()
      
      res1p_exp = data.frame()
      res1n_exp = data.frame()
      
      res1p_log = data.frame()
      res1n_log = data.frame()
      
      # coefficient of variation (after detrending)
      cv_exp = data.frame()
      cv_log = data.frame()
      
      for(b in unique(tmpres1$bin_ID)) {
        ID1 = b # for compatibility
        tmpgrid1 = tmpres1[tmpres1$bin_ID == b,]
        tmpgrid1 = tmpgrid1[order(tmpgrid1$calBP, decreasing = TRUE),]
        
        # ACF without detrending
        acf1 = acf(tmpgrid1$PrDens, lag.max = 4999, plot = FALSE)
        tmp2p = get_acf_max(acf1$acf)
        tmp2n = get_acf_min(acf1$acf)
        if(length(tmp2p) == 1) res1p = rbind(res1p, data.frame(ID = b,
            lag = acf1$lag[tmp2p], val = acf1$acf[tmp2p]))
        if(length(tmp2n) == 1) res1n = rbind(res1n, data.frame(ID = b,
            lag = acf1$lag[tmp2n], val = acf1$acf[tmp2n]))
        
        # detrending -- exponential
        nls1 = nls(PrDens ~ exp(a*calBP + b), tmpgrid1, list(a = -0.0001, b = 0))
        tmpgrid1$fit = predict(nls1)
        tmpgrid1$d1 = tmpgrid1$PrDens - tmpgrid1$fit
        
        # repeat acf
        acf1 = acf(tmpgrid1$d1, lag.max = 4999, plot = FALSE)
        tmp2p = get_acf_max(acf1$acf)
        tmp2n = get_acf_min(acf1$acf)
        
        if(length(tmp2p) == 1) res1p_exp = rbind(res1p_exp, data.frame(ID = b,
                      lag = acf1$lag[tmp2p], val = acf1$acf[tmp2p]))
        if(length(tmp2n) == 1) res1n_exp = rbind(res1n_exp, data.frame(ID = b,
                      lag = acf1$lag[tmp2n], val = acf1$acf[tmp2n]))
        
        # calculate a coefficient of variation
        cv = sd(tmpgrid1$d1) / mean(tmpgrid1$fit)
        cv_exp = rbind(cv_exp, data.frame(ID = b, cv = cv))
        
        # try a logistic fit
        N = mean(tail(tmpgrid1$PrDens, 1000))
        x0 = max(tmpgrid1$calBP[tmpgrid1$PrDens >= N/2])
        nls1 = tryCatch( {
          if(two_step_fit) {
            nls1 = nls(PrDens ~ N / (1 + exp((calBP - x0)/T1)) , tmpgrid1,
                       list(N = N))
            N = nls1$m$getPars()[1]
            nls1 = nls(PrDens ~ N / (1 + exp((calBP - x0)/T1)) , tmpgrid1,
                       list(x0 = x0, T1 = T1))
          } else {
            nls1 = nls(PrDens ~ N / (1 + exp((calBP - x0)/T1)) , tmpgrid1,
                       list(N = N, x0 = x0, T1 = T1))
          }
          nls1
        }, error = function(x) NULL)
        if(is.null(nls1)) {
          tmperr = data.frame(i = i, j = j, res = r, ID = ID1,
                              dx = all_res[[j]]$dx, dy = all_res[[j]]$dy)
          sum_res1[[i]]$errors = rbind(sum_res1[[i]]$errors, tmperr)
          next
        }
        tmpgrid1$pred = predict(nls1)
        tmpgrid1$diff = tmpgrid1$PrDens - tmpgrid1$pred
        
        # do the ACF
        acf1 = acf(tmpgrid1$diff, lag.max = 4999, plot = FALSE)
        
        if(sum(is.nan(acf1$acf)) < (4999)) {
          tmp2p = get_acf_max(acf1$acf)
          tmp2n = get_acf_min(acf1$acf)
          
          if(length(tmp2p) == 1) res1p_log = rbind(res1p_log, data.frame(ID = b,
                        lag = acf1$lag[tmp2p], val = acf1$acf[tmp2p]))
          if(length(tmp2n) == 1) res1n_log = rbind(res1n_log, data.frame(ID = b,
                        lag = acf1$lag[tmp2n], val = acf1$acf[tmp2n]))
        }
        
        # calculate a coefficient of variation
        cv = sd(tmpgrid1$diff) / mean(tmpgrid1$pred)
        cv_log = rbind(cv_log, data.frame(ID = b, cv = cv))
      }
      
      if(nrow(res1p) > 0) {
        res1p$dx = all_res[[j]]$dx
        res1p$dy = all_res[[j]]$dy
        all_res_max = rbind(all_res_max, res1p)
      }
      if(nrow(res1n) > 0) {
        res1n$dx = all_res[[j]]$dx
        res1n$dy = all_res[[j]]$dy
        all_res_min = rbind(all_res_min, res1n)
      }
      
      if(nrow(res1p_exp) > 0) {
        res1p_exp$dx = all_res[[j]]$dx
        res1p_exp$dy = all_res[[j]]$dy
        all_res_exp_max = rbind(all_res_exp_max, res1p_exp)
      }
      if(nrow(res1n_exp) > 0) {
        res1n_exp$dx = all_res[[j]]$dx
        res1n_exp$dy = all_res[[j]]$dy
        all_res_exp_min = rbind(all_res_exp_min, res1n_exp)
      }
      
      if(nrow(res1p_log) > 0) {
        res1p_log$dx = all_res[[j]]$dx
        res1p_log$dy = all_res[[j]]$dy
        all_res_log_max = rbind(all_res_log_max, res1p_log)
      }
      if(nrow(res1n_log) > 0) {
        res1n_log$dx = all_res[[j]]$dx
        res1n_log$dy = all_res[[j]]$dy
        all_res_log_min = rbind(all_res_log_min, res1n_log)
      }
      if(nrow(cv_exp) > 0) {
        cv_exp$dx = all_res[[j]]$dx
        cv_exp$dy = all_res[[j]]$dy
        all_res_cv_exp = rbind(all_res_cv_exp, cv_exp)
      }
      if(nrow(cv_log) > 0) {
        cv_log$dx = all_res[[j]]$dx
        cv_log$dy = all_res[[j]]$dy
        all_res_cv_log = rbind(all_res_cv_log, cv_log)
      }
    }
    
    sum_res1[[i]]$all_res_min = all_res_min
    sum_res1[[i]]$all_res_max = all_res_max
    sum_res1[[i]]$all_res_exp_min = all_res_exp_min
    sum_res1[[i]]$all_res_exp_max = all_res_exp_max
    sum_res1[[i]]$all_res_log_min = all_res_log_min
    sum_res1[[i]]$all_res_log_max = all_res_log_max
    sum_res1[[i]]$all_res_cv_exp = all_res_cv_exp
    sum_res1[[i]]$all_res_cv_log = all_res_cv_log
  }
  
  return(sum_res1)
}


##################################################################
# 2. Load SPD results, calculate ACFs (main case)

tmin = 5000
tmax = 9000

tr1 = paste0(tmax, "_", tmin)
load(paste0(datadir, "/c14_spd_new_", tr1, ".RData"))
# (variable name: all_res)

# calculate ACFs, identify first minimum locations
sum_res1 = do_one_acf(all_res, TRUE, T1 = 100)

# variable names and captions
names = c("all_res_min", "all_res_max" ,"all_res_exp_min",
          "all_res_exp_max", "all_res_log_min", "all_res_log_max") 
desc = c("First minimum, no detrending", "First maximum, no detrending",
         "First minimum, exponential detrending",
         "First maximum, exponential detrending",
         "First minimum, logistic detrending",
         "First maximum, logistic detrending")
vars = c("First minimum", "First maximum", "First minimum",
         "First maximum", "First minimum", "First maximum")

# save the results and create plots of histograms
dir.create(outdir, recursive = TRUE)

for(i in 1:length(sum_res1)) {
  r = sum_res1[[i]]$r
  for(j in 1:length(names)) {
    tmp1 = sum_res1[[i]][[names[j]]]
    p1 = ggplot(tmp1) + geom_histogram(aes(x=lag),
                                       binwidth = 100, boundary = 0)
    p1 = p1 + theme_bw(6) + scale_x_continuous(limits=c(0,2200))
    p1 = p1 + xlab("Lag [years]") + ylab("Frequency")
    p1 = p1 + ggtitle(paste0(vars[j], ", resolution: ", r, " km"))
    fn1 = paste0(outdir, "/", names[j], "_r", r / 100)
    ggsave(paste0(fn1, ".pdf"), width=2.3, height=1.6)
    ggsave(paste0(fn1, ".png"), width=2.3, height=1.6, dpi=300)
  }
  # save the distribution of minimums with expenontial and logistic detrending
  r = r / 100
  tmp1 = sum_res1[[i]]$all_res_exp_min
  write.csv(tmp1, paste0(outdir, "/acf_peaks_min_exp_r", r, ".csv"), row.names=FALSE)
  tmp1 = sum_res1[[i]]$all_res_log_min
  write.csv(tmp1, paste0(outdir, "/acf_peaks_min_log_r", r, ".csv"), row.names=FALSE)
}


######################################################################
# 3. create a faceted plot based on different tiling locations
# (only for the logistic detrending, minimum locations)

# size of figures set manually
w1 = c(90, 122, 171)
h1 = c(75, 100, 137)
rr1 = c(500, 700, 1000)

for(i in 1:length(rr1)) {
  r = rr1[i]
  tmp1 = NA
  for(j in 1:length(sum_res1)) if(sum_res1[[j]]$r == r) {
    tmp1 = sum_res1[[j]]$all_res_log_min
  }
  
  dx1 = unique(tmp1$dx)
  xlab1 = paste0("dx = ", 100*dx1, " km")
  names(xlab1) = dx1
  dy1 = unique(tmp1$dy)
  ylab1 = paste0("dy = ", 100*dy1, " km")
  names(ylab1) = dy1
  
  p1 = ggplot(tmp1) + geom_histogram(aes(x=lag),
                                     binwidth = 100, boundary = 0)
  p1 = p1 + theme_bw(6) + scale_x_continuous(limits=c(0,NA))
  p1 = p1 + xlab("Lag [years]") + ylab("Frequency")
  p1 = p1 + ggtitle(paste0("Tiling linear size: ", r, " km"))
  p1 = p1 + facet_grid(factor(dy, levels = sort(unique(dy), decreasing = TRUE))~dx,
            labeller = labeller(.rows = ylab1, .cols = xlab1))
  p1 = p1 + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  fn1 = paste0(outdir, "/all_res_grid_log_min_r", r / 100, "_grid")
  ggsave(paste0(fn1, ".pdf"), width=w1[i], height=h1[i], units = "mm")
  ggsave(paste0(fn1, ".png"), width=w1[i], height=h1[i], units = "mm", dpi=300)
}


###################################################################
# 4. Create plot of the distribution of coefficients of variation
names = c("all_res_cv_exp", "all_res_cv_log") 
desc = c("Exponential detrending", "Logistic detrending")
for(i in 1:length(sum_res1)) {
  r = sum_res1[[i]]$r
  for(j in 1:length(names)) {
    tmp1 = sum_res1[[i]][[names[j]]]
    p1 = ggplot(tmp1) + geom_histogram(aes(x=cv),
                                       binwidth = 0.04, boundary = 0)
    p1 = p1 + theme_bw(6) + scale_x_continuous(limits=c(0,NA))
    p1 = p1 + xlab("Coefficient of variation") + ylab("Frequency")
    p1 = p1 + ggtitle(paste0(desc[j], ", resolution: ", r, " km"))
    fn1 = paste0(outdir, "/", names[j], "_r", r / 100)
    ggsave(paste0(fn1, ".pdf"), width=2.3, height=1.6)
    ggsave(paste0(fn1, ".png"), width=2.3, height=1.6, dpi=300)
  }
}

for(i in 1:length(sum_res1)) {
  r = sum_res1[[i]]$r / 100
  for(j in 1:length(names)) {
    tmp1 = sum_res1[[i]][[names[j]]]
    write.csv(tmp1, paste0(outdir, "/", names[j], "_r", r, ".csv"), row.names=FALSE)
  }
}


