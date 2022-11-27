# acf_new_aggr_ex.r -- example plots of regional time series
# for the main paper figures; this uses the results from the script
# neolithic_cpp/simulation_runs.fsh

library(ggplot2)

# base directory, change this as appropriate!
base_dir = "~/CSH/HoloSim/data"
base_dir2 = paste0(base_dir, "/simulation/") # directory with simulation results

colors = c("#5aa732","#005c96","#e75216", "#009ace","#ffd500",
           "#8671a7", "#e78502","#db6aa2","#007939","#9a2846","#815234","#000000")

# time range for ACFs (relative to simulation start date)
tmin1 = 0
tmax1 = 4000

# fit logistic growth and calculate the ACF of one time series
do_one_acf = function(tmp2) {
  tmp2 = tmp2[tmp2$year >= tmin1 & tmp2$year < tmax1,]
  all_years = tmin1:(tmax1-1)
  
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
  T1 = 100 # TODO: maybe adjust this?
  nls1 = tryCatch( {
    nls1 = nls(pop ~ N / (1 + exp(-(year - x0)/T1)) , tmp2,
               list(N = N))
    N = nls1$m$getPars()[1]
    nls1 = nls(pop ~ N / (1 + exp(-(year - x0)/T1)) , tmp2,
               list(x0 = x0, T1 = T1), algorithm = "port")
    nls1
  }, error = function(x) {
    NULL
  }
  )
  if(is.null(nls1)) return(NULL)
  
  tmp2$pred = predict(nls1)
  tmp2$diff = tmp2$pop - tmp2$pred
  
  # do the ACF
  acf1 = acf(tmp2$diff, lag.max = 4999, plot = FALSE)
  
  tmp3 = data.frame(lag = acf1$lag, acf = acf1$acf)
  return(tmp3)
}



aggr_key = read.csv(paste0(base_dir, '/new_grid_ids.csv'))
aggr_regions = read.csv(paste0(base_dir, '/new_grids_all.csv'))

aggr_types = unique(aggr_key[c("res", "dx", "dy")])
aggr_types = aggr_types[aggr_types$res == 500,]


y1 = 7000 # simulation start date
cnts_thresh = 10 # mininum number of cells to include in the results

# main parameters
s = 2
G = 80
E = 1
A = 20

fn1 = paste0(base_dir2, "runs_G", G, "_C10_E", E, "_A",
             A, "_s", s, "spaggr.out.gz")
res1 = read.table(gzfile(fn1))

names(res1) = c("year", "index", "cnt", "pop")

# filter out regions with too few cells
tmp1 = aggregate(res1["cnt"], by=res1["index"], FUN=max)
tmp1 = tmp1[tmp1$cnt >= cnts_thresh,]
res1 = res1[res1$index %in% tmp1$index, c("year", "index", "pop")]

# use only one grid position
r = 500
dx = 0
dy = 0
regions = aggr_key[aggr_key$dx == dx & aggr_key$dy == dy & aggr_key$res == r,]

tmp1 = merge(regions[c("index", "ID")], res1, by="index")

all_years = tmin1:(tmax1-1)

err_ID = NULL

# calculate ACFs
all_acf = data.frame()
for(ID1 in unique(tmp1$ID)) {
  tmp2 = tmp1[tmp1$ID == ID1 & tmp1$year >= tmin1 & tmp1$year < tmax1,
              c("year", "pop")]
  tmp3 = do_one_acf(tmp2)
  if(is.null(tmp3)) err_ID = c(err_ID, ID1)
  else {
    tmp3$ID = ID1
    all_acf = rbind(all_acf, tmp3)
  }
}


# select the regions with the top 10 population as examples
# (main text figures)
tmp2 = aggregate(tmp1["pop"], by=tmp1["ID"], FUN=mean)
IDs1 = tmp2$ID[order(tmp2$pop, decreasing = TRUE)]
IDs1 = IDs1[1:10]

tmp12 = tmp1[tmp1$ID %in% IDs1,]

p1 = ggplot(tmp12) + geom_line(size = 0.2,
  aes(x=-1*(year-y1), y=pop/1000, group=ID, color=factor(ID)))
p1 = p1 + scale_color_manual(values = colors) + theme_bw(6)
p1 = p1 + xlab("Date [BCE]") + ylab("Population [x1000]")
p1 = p1 + scale_x_reverse(limits=c(7000, 3000))
p1 = p1 + guides(color="none")

fn1 = paste0(base_dir, "/fig3a.pdf")
ggsave(fn1, p1, width=1.6, height=1.2)


tmp22 = all_acf[all_acf$ID %in% IDs1,]
p1 = ggplot(tmp22) + geom_line(size = 0.2,
  aes(x=lag, y=acf, group=ID, color=factor(ID)))
p1 = p1 + scale_color_manual(values = colors) + theme_bw(6)
p1 = p1 + xlab("Lag [years]") + ylab("ACF")
p1 = p1 + guides(color="none")

fn1 = paste0(base_dir, "/fig3b.pdf")
ggsave(fn1, p1, width=1.6, height=1.2)



###################################################################
# Fig 2. -- simulation without conflict
G = 40
fn1 = paste0(base_dir2, "runG", G, "C10ws", s, "spaggr.out.gz")
res1 = read.table(gzfile(fn1))
names(res1) = c("year", "index", "cnt", "pop")
# filter for the regions in one grid position
res1 = merge(regions[c("index", "ID")], res1, by="index")

# select the regions with the top 10 population after 3000 years as examples
# tmp2 = aggregate(res1["pop"], by=res1["index"], FUN=mean)
tmp2 = res1[res1$year == 3000, c("index", "pop")]
IDs1 = tmp2$index[order(tmp2$pop, decreasing = TRUE)]
IDs1 = IDs1[1:10]

res2 = res1[res1$index %in% IDs1,]
tmp1 = res2[1:4]
names(tmp1)[2] = "ID"

# acfs
all_years = tmin1:(tmax1-1)
err_ID = NULL
all_acf = data.frame()
for(ID1 in unique(tmp1$ID)) {
  tmp2 = tmp1[tmp1$ID == ID1 & tmp1$year >= tmin1 & tmp1$year < tmax1,
              c("year", "pop")]
  tmp3 = do_one_acf(tmp2)
  if(is.null(tmp3)) err_ID = c(err_ID, ID1)
  else {
    tmp3$ID = ID1
    all_acf = rbind(all_acf, tmp3)
  }
}


p1 = ggplot(tmp1) + geom_line(size = 0.2, aes(x=-1*(year-y1),
                          y=pop/1000, group=ID, color=factor(ID)))
p1 = p1 + scale_color_manual(values = colors) + theme_bw(6)
p1 = p1 + xlab("Date [BCE]") + ylab("Population [x1000]")
p1 = p1 + scale_x_reverse(limits=c(7000, 3000))
p1 = p1 + guides(color="none")
p1 = p1 + scale_y_continuous(limits=c(0,400))

fn1 = paste0(base_dir, "/fig2a.pdf")
ggsave(fn1, p1, width=1.6, height=1.2)


tmp22 = all_acf[all_acf$ID %in% IDs1,]
p1 = ggplot(tmp22) + geom_line(size = 0.2, aes(x=lag, y=acf, group=ID,
                                               color=factor(ID)))
p1 = p1 + scale_color_manual(values = colors) + theme_bw(6)
p1 = p1 + xlab("Lag [years]") + ylab("ACF")
p1 = p1 + guides(color="none")

fn1 = paste0(base_dir, "/fig2b.pdf")
ggsave(fn1, p1, width=1.6, height=1.2)



####################################################################
# Figures for the SI -- multiple s values together
G = 40
IDs1 = NULL
all_res1 = data.frame()
for(s in 1:4) {
  fn1 = paste0(base_dir2, "runG", G, "C10ws", s, "spaggr.out.gz")
  res1 = read.table(gzfile(fn1))
  names(res1) = c("year", "index", "cnt", "pop")
  # filter for the regions in one grid position
  res1 = merge(regions[c("index", "ID")], res1, by="index")

  tmp2 = res1[res1$year >= 2000, c("index", "pop")]
  tmp2 = aggregate(tmp2["pop"], by=tmp2["index"], FUN=mean)
  if(s == 1) IDs1 = tmp2$index[tmp2$pop > 10000]
  
  res1 = res1[res1$index %in% IDs1,]
  res1 = res1[res1$year <= 4000,]
  res1$s = s
  all_res1 = rbind(all_res1, res1)
}


colors = c("#5aa732","#005c96","#e78502","#db6aa2")

all_res12 = all_res1[all_res1$year%%10 == 0,]
p1 = ggplot(all_res12) + geom_line(aes(x=-1*(year-y1), y=pop/1000,
                          group=interaction(index,s), color=factor(s)))
p1 = p1 + theme_bw(6)
p1 = p1 + xlab("Date [BCE]") + ylab("Population [x1000]")
p1 = p1 + scale_x_reverse(limits=c(7000, 3000))
p1 = p1 + facet_wrap(~index, ncol=4, nrow=6, scales="free_y")
p1 = p1 + scale_color_manual(values=colors)
ggsave(paste0(base_dir, "/simn_regions.pdf"), width=210, height=210, units="mm")

# create ACFs of these
all_acf = data.frame()
all_err = data.frame()
for(s in 1:4) for(ID in unique(all_res1$index)) {
  tmp1 = all_res1[all_res1$s == s & all_res1$index == ID, c("year", "pop")]
  tmp2 = do_one_acf(tmp1)
  if(is.null(tmp2)) all_err = rbind(all_err, data.frame(index=ID, s=s))
  else {
    tmp2$s = s
    tmp2$index = ID
    all_acf = rbind(all_acf, tmp2)
  }
}

all_acf2 = all_acf[all_acf$lag%%10 == 0,]
p1 = ggplot(all_acf2) + geom_line(aes(x=lag, y=acf, color=factor(s),
                                     group=interaction(index,s)))
p1 = p1 + theme_bw(6) + xlab("Lag [years]") + ylab("ACF")
p1 = p1 + facet_wrap(~index, ncol=4, nrow=6, scales="free_y")
p1 = p1 + scale_color_manual(values=colors)
ggsave(paste0(base_dir, "/simn_regions_acf.pdf"), width=210, height=210, units="mm")


# same for the simulation with conflict (stationary aggressors)
G = 80
E1 = c(1,5)
A1 = c(10,20)

r = 500
dx = 0
dy = 0
regions = aggr_key[aggr_key$dx == dx & aggr_key$dy == dy & aggr_key$res == r,]

all_res2 = data.frame()
for(s in 0:4) for(A in A1) for(E in E1) {
  fn1 = paste0(base_dir2, "runs_s", s,  "/runs_G", G, "_C10_E", E, "_A",
             A, "_s", s, "spaggr.out.gz")
  res1 = read.table(gzfile(fn1))
  names(res1) = c("year", "index", "cnt", "pop")
  res1 = res1[res1$year <= tmax1,]
  res1 = res1[res1$index %in% regions$index,]
  res1$s = s
  res1$A = A
  res1$E = E
  all_res2 = rbind(all_res2, res1)
}
write.csv(all_res2, paste0(base_dir, "/all_res_ex.csv"), row.names = FALSE)
# all_res2 = read.csv(paste0(base_dir, "/all_res_ex.csv"))

tmp2 = all_res2[all_res2$year >= 2000 & all_res2$s == 1,]
tmp2 = aggregate(tmp2["pop"], by=tmp2[c("index", "s", "A", "E")], FUN=mean)
tmp2$ind = tmp2$pop >= 700

IDs1 = tmp2$index[tmp2$s == 1 & tmp2$A == 20 & tmp2$E == 5 & tmp2$pop >= 700]

all_res2 = all_res2[all_res2$index %in% IDs1,]

colors = c("#5aa732","#005c96","#e78502","#db6aa2","#9a2846")

# only plot every 10th year to decrease figure sizes
# this makes the plots slightly less "hairy", but there is no
# big (noticeable) difference
all_res21 = all_res2[all_res2$year %% 10 == 0,]

for(A in A1) for(E in E1) {
  tmp1 = all_res21[all_res21$A == A & all_res21$E == E,]
  p1 = ggplot(tmp1)
  p1 = p1 + geom_line(aes(x=-1*(year-y1), y=pop/1000, color=factor(s),
                          group=interaction(index,s),))
  p1 = p1 + theme_bw(6)
  p1 = p1 + xlab("Date [BCE]") + ylab("Population [x1000]")
  p1 = p1 + scale_x_reverse(limits=c(7000, 3000))
  p1 = p1 + facet_wrap(~index, ncol=4, nrow=6, scales="free_y")
  p1 = p1 + scale_color_manual(values=colors)
  ggsave(paste0(base_dir, "/simc_regions_A", A, "_E", E, ".pdf"), width=210, height=210, units="mm")
}


# create ACFs of these
all_acf = data.frame()
all_err = data.frame()
for(s in 0:4) for(ID in unique(all_res2$index)) for(A in A1) for(E in E1) {
  tmp1 = all_res2[all_res2$s == s & all_res2$index == ID &
                    all_res2$A == A & all_res2$E == E, c("year", "pop")]
  tmp2 = do_one_acf(tmp1)
  if(is.null(tmp2)) all_err = rbind(all_err, data.frame(index=ID, s=s))
  else {
    tmp2$s = s
    tmp2$index = ID
    tmp2$A = A
    tmp2$E = E
    all_acf = rbind(all_acf, tmp2)
  }
}

all_acf2 = all_acf[all_acf$lag%%10 == 0,]
for(A in A1) for(E in E1) {
  tmp1 = all_acf2[all_acf2$A == A & all_acf2$E == E,]
  p1 = ggplot(tmp1) + geom_line(aes(x=lag, y=acf, color=factor(s),
                                       group=interaction(index,s)))
  p1 = p1 + theme_bw(6) + xlab("Lag [years]") + ylab("ACF")
  p1 = p1 + facet_wrap(~index, ncol=4, nrow=6, scales="free_y")
  p1 = p1 + scale_color_manual(values=colors)
  ggsave(paste0(base_dir, "/simc_regions_A", A, "_E", E, "_acf.pdf"), width=210, height=210, units="mm")
}

