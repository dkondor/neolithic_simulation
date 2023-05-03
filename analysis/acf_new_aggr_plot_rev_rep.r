# acf_new_aggr_plot_rep.r -- process the results of ACF for the repeated
# simulation realizations, create plots of histograms


stop("Do not run the whole script!")

library(ggplot2)
library(foreach)
library(doParallel)

base_dir = "/home/dkondor/CSH/HoloSim/data"

binwidth = 100 # binwidth for ACF peak distributions
cv_bw = 0.04 # binwidth for CV distributions
get_bin = function(x, bw) {
  return(bw * floor(x / bw) + bw / 2)
}
create_hist = function(peaks1, bw, var) {
  tmp1 = peaks1
  tmp1["bin"] = get_bin(tmp1[var], bw)
  tmp1$count = 1
  tmp1 = aggregate(tmp1["count"], by=tmp1["bin"], FUN=sum)
  scnt = sum(tmp1$count)
  tmp1$p = tmp1$count / scnt
  return(tmp1)
}


# load the "baseline" data (results from processing C14 dates)
# ACF minimum distribution
c14_peaks = read.csv(paste0(base_dir, '/population/new_analysis/res_min.csv'))
# confirm that we're using the results without a moving average
c14_peaks = c14_peaks[c14_peaks$rm == 0,]
c14_peaks = create_hist(c14_peaks, binwidth, "m")

# CV distribution
c14_cv = read.csv(paste0(base_dir, '/population/new_analysis/res_cv.csv'))
c14_cv = c14_cv[c14_cv$rm == 0,]
c14_cv = create_hist(c14_cv, cv_bw, "cv")


# read the results of running multiple realizations of the simulation
# (distribution of ACF minima or CV)
read_histograms = function(fn, r = 500, varname1 = "all_res_min_dt",
                           varname2 = "lag", bw = binwidth) {
  load(fn) # loads variable all_res
  N = length(all_res)
  
  histograms = data.frame()
  for(i in 1:N) {
    tmp1 = all_res[[i]]$res
    peaks = data.frame()
    for(j in 1:length(tmp1)) {
      if(tmp1[[j]]$r %in% r) {
        tmp2 = tmp1[[j]][[varname1]]
        if(length(r) > 1) tmp2$r = tmp1[[j]]$r
        peaks = rbind(peaks, tmp2)
      }
    }
    if(is.null(peaks)) stop('Cannot find result!\n')
    if(nrow(peaks) == 0) next
    peaks$bin = get_bin(peaks[[varname2]], bw)
    peaks$count = 1
    if(length(r) > 1) {
      peaks = aggregate(peaks["count"], by=peaks[c("bin", "r")], FUN=sum)
      sum1 = aggregate(peaks["count"], by=peaks["r"], FUN=sum)
      names(sum1)[names(sum1) == "count"] = "scnt"
      peaks = merge(peaks, sum1, by="r")
      peaks$p = peaks$count / peaks$scnt
    } else {
      peaks = aggregate(peaks["count"], by=peaks["bin"], FUN=sum)
      sum1 = sum(peaks$count)
      peaks$p = peaks$count / sum1
    }
    peaks$i = i
    histograms = rbind(histograms, peaks)
  }
  rm(all_res)
  return(histograms)
}

# aggregate histograms from multiple simulation realizations
aggr_hist = function(histograms) {
  all_i = unique(histograms$i)
  all_bins = unique(histograms$bin)
  tmp1 = data.frame(i = rep(all_i, each=length(all_bins)),
                    bin = rep(all_bins, length(all_i)))
  tmp1 = merge(tmp1, histograms, by=c("i", "bin"), all.x = TRUE)
  tmp1[is.na(tmp1)] = 0
  tmp1 = aggregate(tmp1["p"], by=tmp1["bin"], FUN=function(x) {
    return(c(mean = mean(x), sd = sd(x)))
  })
  tmp1 = do.call(data.frame, tmp1)
  return(tmp1)
}



######################################################################
# main results: stationary aggressors variant
GG = c(50,60,70,80)
s1 = 0:4
base_dir2 = paste0(base_dir, "/simulation/runs_rep")
base_fn = "res"
AA = c(1,2,5,10,20,50)
EE = c(1,5,10,20,100,200)
w=6.9; h=3.75 # figure size

cl = parallel::makeForkCluster(4)
doParallel::registerDoParallel(cl)

foreach(G = GG) %:% foreach(s = s1) %dopar% {
  all_aggr = data.frame()
  for(A in AA) for(E in EE) {
    fn1 = paste0(base_dir2, "/", base_fn, "_G", G, "_C10_E", E, "_A", A,
                 "_s", s, "_all_res.RData")
    tmp1 = read_histograms(fn1)
    aggr1 = aggr_hist(tmp1)
    aggr1$pe = 1 / E
    aggr1$pa = 1 / A
    all_aggr = rbind(all_aggr, aggr1)
  }
  
  p1 = ggplot(all_aggr) + geom_col(data=c14_peaks, aes(x=bin, y=p), color="grey", fill="grey")
  p1 = p1 + geom_col(aes(x=bin, y=p.mean), color="red", fill="red", alpha = 0.45)
  
  p1 = p1 + theme_bw(6) + xlab("Lag [years]") + ylab("Relative frequency")
  p1 = p1 + scale_x_continuous(limits=c(0,2200))
  p1 = p1 + facet_grid(pa~pe, labeller = label_bquote(rows = p[A] == .(pa),
                                                      cols = p[E] == .(pe)))
  
  fnbase1 = paste0(base_dir2, "/acfrep_G", G, "_s", s)
  
  ggsave(paste0(fnbase1, ".pdf"), p1, width=w, height=h)
  ggsave(paste0(fnbase1, ".png"), p1, width=w, height=h, dpi=300)
}

# CV values -- only for G = 80
base_dir2 = paste0(base_dir,"/simulation/runs_rep")
G = 80

foreach(s = s1) %dopar% {
  all_aggr = data.frame()
  for(A in AA) for(E in EE) {
    fn1 = paste0(base_dir2, "/", base_fn, "_G", G, "_C10_E", E, "_A", A,
                 "_s", s, "_all_res.RData")
    tmp1 = read_histograms(fn1, 500, "cv", "cv", cv_bw)
    aggr1 = aggr_hist(tmp1)
    aggr1$pe = 1 / E
    aggr1$pa = 1 / A
    all_aggr = rbind(all_aggr, aggr1)
  }
  
  p1 = ggplot(all_aggr) + geom_col(data=c14_cv, aes(x=bin, y=p), color="grey", fill="grey")
  p1 = p1 + geom_col(aes(x=bin, y=p.mean), color="red", fill="red", alpha = 0.45)
  
  p1 = p1 + theme_bw(6) + xlab("CV") + ylab("Relative frequency")
  p1 = p1 + scale_x_continuous(limits=c(0,1.5))
  p1 = p1 + facet_grid(pa~pe, scales = "free_y",
                       labeller = label_bquote(rows = p[A] == .(pa),
                                               cols = p[E] == .(pe)))
  
  fnbase1 = paste0(base_dir2, "/cvrep_G", G, "_s", s)
  
  ggsave(paste0(fnbase1, ".pdf"), p1, width=w, height=h)
  ggsave(paste0(fnbase1, ".png"), p1, width=w, height=h, dpi=300)
}


# alternate presentation: only E = 1 and A = 20 (main pars),
# but all G,s pairs in one figure
GG = c(50,60,70,80)
s1 = 0:4
base_dir2 = paste0(base_dir, "/simulation/runs_rep")
base_fn = "res"

A = 20
E = 1
all_aggr = data.frame()
all_cv = data.frame()
for(G in GG) for(s in s1) {
  fn1 = paste0(base_dir2, "/", base_fn, "_G", G, "_C10_E", E, "_A", A,
                 "_s", s, "_all_res.RData")
  tmp1 = read_histograms(fn1)
  aggr1 = aggr_hist(tmp1)
  aggr1$s = s
  aggr1$G = G
  all_aggr = rbind(all_aggr, aggr1)
  cv1 = read_histograms(fn1, 500, "cv", "cv", cv_bw)
  cv1 = aggr_hist(cv1)
  cv1$s = s
  cv1$G = G
  all_cv = rbind(all_cv, cv1)
}

w=6.9; h=3.75 # figure size

p1 = ggplot(all_aggr) + geom_col(data=c14_peaks, aes(x=bin, y=p), color="grey", fill="grey")
p1 = p1 + geom_col(aes(x=bin, y=p.mean), color="red", fill="red", alpha = 0.45)

p1 = p1 + theme_bw(6) + xlab("Lag [years]") + ylab("Relative frequency")
p1 = p1 + scale_x_continuous(limits=c(0,2200))
p1 = p1 + facet_grid(G~s, labeller = labeller(s = function(x){
  return(paste0("s = ", x))}, G = function(x){
    return(paste0("G = ", x, " km"))}))

fnbase1 = paste0(base_dir2, "/acfrep_A", A, "_E", E)

ggsave(paste0(fnbase1, ".pdf"), p1, width=w, height=h)
ggsave(paste0(fnbase1, ".png"), p1, width=w, height=h, dpi=300)


p1 = ggplot(all_cv) + geom_col(data=c14_cv, aes(x=bin, y=p), color="grey", fill="grey")
p1 = p1 + geom_col(aes(x=bin, y=p.mean), color="red", fill="red", alpha = 0.45)

p1 = p1 + theme_bw(6) + xlab("CV") + ylab("Relative frequency")
p1 = p1 + scale_x_continuous(limits=c(0,1.5))
p1 = p1 + facet_grid(G~s, labeller = labeller(s = function(x){
  return(paste0("s = ", x))}, G = function(x){
    return(paste0("G = ", x, " km"))}))

fnbase1 = paste0(base_dir2, "/cvrep_A", A, "_E", E)

ggsave(paste0(fnbase1, ".pdf"), p1, width=w, height=h)
ggsave(paste0(fnbase1, ".png"), p1, width=w, height=h, dpi=300)



######################################################################
# model without conflict
w = 4.7; h = 2.4
base_dir2 = paste0(base_dir, "/simulation/runs_rep")
base_fn = "resn"
GG = c(20,30,40)
s1 = 1:4

all_aggr = data.frame()
for(G in GG) for(s in s1) {
  fn1 = paste0(base_dir2, "/", base_fn, "_G", G, "_s", s, "_all_res.RData")
  tmp1 = read_histograms(fn1)
  aggr1 = aggr_hist(tmp1)
  aggr1$G = G
  aggr1$s = s
  all_aggr = rbind(all_aggr, aggr1)
}
  
p1 = ggplot(all_aggr) + geom_col(data=c14_peaks, aes(x=bin, y=p), color="grey", fill="grey")
p1 = p1 + geom_col(aes(x=bin, y=p.mean), color="red", fill="red", alpha = 0.45)

p1 = p1 + theme_bw(6) + xlab("Lag [years]") + ylab("Relative frequency")
p1 = p1 + scale_x_continuous(limits=c(0,2200))
p1 = p1 + facet_grid(G~s, labeller = label_bquote(rows = G == .(G),
                                                    cols = s == .(s)))
p1 = p1 + facet_grid(G~s, labeller = labeller(s = function(x){
  return(paste0("s = ", x))}, G = function(x){
    return(paste0("G = ", x, " km"))}))

fnbase1 = paste0(base_dir2, "/acfrepn_all")

ggsave(paste0(fnbase1, ".pdf"), p1, width=w, height=h)
ggsave(paste0(fnbase1, ".png"), p1, width=w, height=h, dpi=300)


#######################################################################
# same for cv
w = 5; h = 2.4
base_dir2 = paste0(base_dir, "/simulation/runs_rep")
base_fn = "resn"
GG = c(20,30,40)
s1 = 1:4

all_aggr = data.frame()
for(G in GG) for(s in s1) {
  fn1 = paste0(base_dir2, "/", base_fn, "_G", G, "_s", s, "_all_res.RData")
  tmp1 = read_histograms(fn1, 500, "cv", "cv", cv_bw)
  aggr1 = aggr_hist(tmp1)
  aggr1$G = G
  aggr1$s = s
  all_aggr = rbind(all_aggr, aggr1)
}

p1 = ggplot(all_aggr) + geom_col(data=c14_cv, aes(x=bin, y=p), color="grey", fill="grey")
p1 = p1 + geom_col(aes(x=bin, y=p.mean), color="red", fill="red", alpha = 0.45)

p1 = p1 + theme_bw(6) + xlab("CV") + ylab("Relative frequency")
p1 = p1 + scale_x_continuous(limits=c(0,1.5))
p1 = p1 + facet_grid(G~s, labeller = labeller(s = function(x){
  return(paste0("s = ", x))}, G = function(x){
    return(paste0("G = ", x, " km"))}))

fnbase1 = paste0(base_dir2, "/cvrepn_all")

ggsave(paste0(fnbase1, ".pdf"), p1, width=w, height=h)
ggsave(paste0(fnbase1, ".png"), p1, width=w, height=h, dpi=300)


p1 = p1 + facet_wrap(G~s, scales="free_y", labeller = labeller(s = function(x){
  return(paste0("s = ", x))}, G = function(x){
    return(paste0("G = ", x, " km"))}))

ggsave(paste0(fnbase1, "_2.pdf"), p1, width=5, height=3.2)




###################################################################
# individual panels for the main paper figures
# Fig. 4 -- simulation without conflict
G = 40
s = 2
fn1 = paste0(base_dir2, "/", base_fn, "_G", G, "_s", s, "_all_res.RData")
cv1 = read_histograms(fn1, 500, "cv", "cv", cv_bw)
cv1 = aggr_hist(cv1)
acf1 = read_histograms(fn1)
acf1 = aggr_hist(acf1)

p1 = ggplot(cv1) + geom_col(data=c14_cv, aes(x=bin, y=p), color="grey", fill="grey")
p1 = p1 + geom_col(aes(x=bin, y=p.mean), color="red", fill="red", alpha = 0.45)

p1 = p1 + theme_bw(6) + xlab("Coefficient of variation") + ylab("Relative frequency")
p1 = p1 + scale_x_continuous(limits=c(0,1.5))

fn1 = paste0(base_dir2, "/fig2d.pdf")
ggsave(fn1, p1, width=1.6, height=1.2)


p1 = ggplot(acf1) + geom_col(data=c14_peaks, aes(x=bin, y=p), color="grey", fill="grey")
p1 = p1 + geom_col(aes(x=bin, y=p.mean), color="red", fill="red", alpha = 0.45)

p1 = p1 + theme_bw(6) + xlab("Lag [years]") + ylab("Relative frequency")
p1 = p1 + scale_x_continuous(limits=c(0,2200))

fn1 = paste0(base_dir2, "/fig2c.pdf")
ggsave(fn1, p1, width=1.6, height=1.2)


# Fig. 5 main variant with stationary aggressors
G = 80
s = 2
E = 1
A = 20

base_dir2 = "simulation/runs_rep"
base_fn = "res"
fn1 = paste0(base_dir2, "/", base_fn, "_G", G, "_C10_E", E, "_A", A,
             "_s", s, "_all_res.RData")
acf1 = read_histograms(fn1)
acf1 = aggr_hist(acf1)
cv1 = read_histograms(fn1, 500, "cv", "cv", cv_bw)
cv1 = aggr_hist(cv1)

p1 = ggplot(acf1) + geom_col(data=c14_peaks, aes(x=bin, y=p), color="grey", fill="grey")
p1 = p1 + geom_col(aes(x=bin, y=p.mean), color="red", fill="red", alpha = 0.45)

p1 = p1 + theme_bw(6) + xlab("Lag [years]") + ylab("Relative frequency")
p1 = p1 + scale_x_continuous(limits=c(0,2200))

fn1 = paste0(base_dir2, "/fig3c.pdf")
ggsave(fn1, p1, width=1.6, height=1.2)

p1 = ggplot(cv1) + geom_col(data=c14_cv, aes(x=bin, y=p), color="grey", fill="grey")
p1 = p1 + geom_col(aes(x=bin, y=p.mean), color="red", fill="red", alpha = 0.45)

p1 = p1 + theme_bw(6) + xlab("Coefficient of variation") + ylab("Relative frequency")
p1 = p1 + scale_x_continuous(limits=c(0,1.5))

fn1 = paste0(base_dir2, "/fig3d.pdf")
ggsave(fn1, p1, width=1.6, height=1.2)


##################################################################
# model with static conflict (non-dynamic variable)
GG = c(50,60,70,80)
s1 = 0:4
base_dir2 = paste0(base_dir, "/simulation/runs_rep")
base_fn = "resnr"
EE = c(1,5,10,20,100,200)
w=6.9; h=3.75 # figure size


cl = parallel::makeForkCluster(4)
doParallel::registerDoParallel(cl)

foreach(G = GG) %dopar% {
  all_aggr = data.frame()
  for(E in EE) for(s in s1) {
    fn1 = paste0(base_dir2, "/", base_fn, "_G", G, "_C10_E", E,
                 "_s", s, "_all_res.RData")
    tmp1 = read_histograms(fn1)
    aggr1 = aggr_hist(tmp1)
    aggr1$pe = 1 / E
    aggr1$s = s
    all_aggr = rbind(all_aggr, aggr1)
  }
  
  p1 = ggplot(all_aggr) + geom_col(data=c14_peaks, aes(x=bin, y=p), color="grey", fill="grey")
  p1 = p1 + geom_col(aes(x=bin, y=p.mean), color="red", fill="red", alpha = 0.45)
  
  p1 = p1 + theme_bw(6) + xlab("Lag [years]") + ylab("Relative frequency")
  p1 = p1 + scale_x_continuous(limits=c(0,2200))
  p1 = p1 + facet_grid(s~pe, labeller = label_bquote(rows = s == .(s),
                                                      cols = p[E] == .(pe)))
  
  fnbase1 = paste0(base_dir2, "/acfrepnr_G", G)
  
  ggsave(paste0(fnbase1, ".pdf"), p1, width=w, height=h)
  ggsave(paste0(fnbase1, ".png"), p1, width=w, height=h, dpi=300)
}


# CV plots (only G = 80)
G = 80
s1 = 0:4
base_dir2 = paste0(base_dir, "/simulation/runs_rep")
base_fn = "resnr"
EE = c(1,5,10,20,100,200)
w=6.9; h=3.75 # figure size


all_aggr = data.frame()
all_cv = data.frame()
for(E in EE) for(s in s1) {
  fn1 = paste0(base_dir2, "/", base_fn, "_G", G, "_C10_E", E,
               "_s", s, "_all_res.RData")
  tmp1 = read_histograms(fn1)
  aggr1 = aggr_hist(tmp1)
  aggr1$pe = 1 / E
  aggr1$s = s
  all_aggr = rbind(all_aggr, aggr1)
  tmp1 = read_histograms(fn1, 500, "cv", "cv", cv_bw)
  cv1 = aggr_hist(tmp1)
  cv1$pe = 1 / E
  cv1$s = s
  all_cv = rbind(all_cv, cv1)
}

p1 = ggplot(all_aggr) + geom_col(data=c14_peaks, aes(x=bin, y=p), color="grey", fill="grey")
p1 = p1 + geom_col(aes(x=bin, y=p.mean), color="red", fill="red", alpha = 0.45)

p1 = p1 + theme_bw(6) + xlab("Lag [years]") + ylab("Relative frequency")
p1 = p1 + scale_x_continuous(limits=c(0,2200))
p1 = p1 + facet_grid(s~pe, labeller = label_bquote(rows = s == .(s),
                                                   cols = p[E] == .(pe)))

fnbase1 = paste0(base_dir2, "/acfrepnr_G", G)

ggsave(paste0(fnbase1, ".pdf"), p1, width=w, height=h)
ggsave(paste0(fnbase1, ".png"), p1, width=w, height=h, dpi=300)


p1 = ggplot(all_cv) + geom_col(data=c14_cv, aes(x=bin, y=p), color="grey", fill="grey")
p1 = p1 + geom_col(aes(x=bin, y=p.mean), color="red", fill="red", alpha = 0.45)

p1 = p1 + theme_bw(6) + xlab("CV") + ylab("Relative frequency")
p1 = p1 + scale_x_continuous(limits=c(0,1.5))
p1 = p1 + facet_grid(s~pe, scales = "free_y", labeller = label_bquote(rows = s == .(s),
                                                   cols = p[E] == .(pe)))

fnbase1 = paste0(base_dir2, "/cvrepnr_G", G)

ggsave(paste0(fnbase1, ".pdf"), p1, width=w, height=h)
ggsave(paste0(fnbase1, ".png"), p1, width=w, height=h, dpi=300)


######################################################################
# 3. additional model variants (only for a limited set of parameters)
G = 80
s = 2
A = 10
EE = c(1,5)

w = 3.4; h = 2.2

var_P = c("0", "P")
names_P = c("E", "F")
var_S = c("0", "S")
names_S = c("C", "D")
var_L = c("0", "L")
names_L = c("A", "B")

all_aggr = data.frame()
all_cv = data.frame()
base_dir2 = paste0(base_dir, "/simulation/runs_rep_new")
base_fn = "res"
for(E in EE) for(i in 1:2) for(j in 1:2) for(k in 1:2) {
  P = var_P[i]
  S = var_S[j]
  L = var_L[k]
  name1 = paste0(names_L[k], names_S[j], names_P[i])
  fn1 = paste0(base_dir2, "/", base_fn, P, S, L,
               "_G", G, "_C10_E", E, "_A", A,
               "_s", s, "_all_res.RData")
  tmp1 = read_histograms(fn1)
  aggr1 = aggr_hist(tmp1)
  aggr1$E = E
  aggr1$pe = 1 / E
  aggr1$name = name1
  aggr1$P = P
  aggr1$S = S
  aggr1$L = L
  all_aggr = rbind(all_aggr, aggr1)
  
  tmp1 = read_histograms(fn1, 500, "cv", "cv", cv_bw)
  aggr1 = aggr_hist(tmp1)
  aggr1$E = E
  aggr1$pe = 1 / E
  aggr1$name = name1
  aggr1$P = P
  aggr1$S = S
  aggr1$L = L
  all_cv = rbind(all_cv, aggr1)
}


# do 2x2 figures (by P and S parameters -> 4 figs. in total)
for(E in EE) for(L in c("0", "L")) {
  tmp1 = all_aggr[all_aggr$E == E & all_aggr$L == L,]
  tmp2 = all_cv[all_cv$E == E & all_cv$L == L,]
  
  p1 = ggplot(tmp1) + geom_col(data=c14_peaks, aes(x=bin, y=p), color="grey", fill="grey")
  p1 = p1 + geom_col(aes(x=bin, y=p.mean), color="red", fill="red", alpha = 0.45)
  p1 = p1 + theme_bw(6) + xlab("Lag [years]") + ylab("Relative frequency")
  p1 = p1 + scale_x_continuous(limits=c(0,2200))
  p1 = p1 + facet_grid(P~S)
  
  fnbase1 = paste0(base_dir2, "/acfrep", L, "_E", E)
  ggsave(paste0(fnbase1, ".pdf"), p1, width=w, height=h)
  ggsave(paste0(fnbase1, ".png"), p1, width=w, height=h, dpi=300)
  
  p1 = ggplot(tmp2) + geom_col(data=c14_cv, aes(x=bin, y=p), color="grey", fill="grey")
  p1 = p1 + geom_col(aes(x=bin, y=p.mean), color="red", fill="red", alpha = 0.45)
  p1 = p1 + theme_bw(6) + xlab("CV") + ylab("Relative frequency")
  p1 = p1 + scale_x_continuous(limits=c(0,1.5))
  p1 = p1 + facet_grid(P~S)
  
  fnbase1 = paste0(base_dir2, "/cvrep", L, "_E", E)
  ggsave(paste0(fnbase1, ".pdf"), p1, width=w, height=h)
  ggsave(paste0(fnbase1, ".png"), p1, width=w, height=h, dpi=300)
}


# do figures together, based on model variant "names"
w = 4.5; h = 2.4
for(E in EE) {
  tmp1 = all_aggr[all_aggr$E == E,]
  tmp2 = all_cv[all_cv$E == E,]
  tmp1$lab = paste0(tmp1$name, " model variant")
  tmp2$lab = paste0(tmp2$name, " model variant")
  
  p1 = ggplot(tmp1) + geom_col(data=c14_peaks, aes(x=bin, y=p), color="grey", fill="grey")
  p1 = p1 + geom_col(aes(x=bin, y=p.mean), color="red", fill="red", alpha = 0.45)
  p1 = p1 + theme_bw(6) + xlab("Lag [years]") + ylab("Relative frequency")
  p1 = p1 + scale_x_continuous(limits=c(0,2200))
  p1 = p1 + facet_wrap(~lab, nrow = 2, ncol = 4)
  
  fnbase1 = paste0(base_dir2, "/acfrep_cmb_E", E)
  ggsave(paste0(fnbase1, ".pdf"), p1, width=w, height=h)
  ggsave(paste0(fnbase1, ".png"), p1, width=w, height=h, dpi=300)
  
  
  p1 = ggplot(tmp2) + geom_col(data=c14_cv, aes(x=bin, y=p), color="grey", fill="grey")
  p1 = p1 + geom_col(aes(x=bin, y=p.mean), color="red", fill="red", alpha = 0.45)
  p1 = p1 + theme_bw(6) + xlab("CV") + ylab("Relative frequency")
  p1 = p1 + scale_x_continuous(limits=c(0,1.5))
  p1 = p1 + facet_wrap(~lab, nrow = 2, ncol = 4)
  
  fnbase1 = paste0(base_dir2, "/cvrep_cmb_E", E)
  ggsave(paste0(fnbase1, ".pdf"), p1, width=w, height=h)
  ggsave(paste0(fnbase1, ".png"), p1, width=w, height=h, dpi=300)
}




