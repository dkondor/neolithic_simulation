# c14_model_plots1.r -- create plots of the histogram of ACF minima

library(ggplot2)

# base directory of all data files
base_dir = '~/CSH/HoloSim/data'
datadir = paste0(base_dir, '/population') # directory with C14 data
outdir1 = paste0(datadir, '/new_analysis')


#######################################################################
# set of figures of regional time series and ACFs
# create plots of the results
plot_one_spd = function(spdc, pval = NULL) {
  p1 = ggplot(spdc) + scale_x_reverse()
  p1 = p1 + geom_line(aes(x=calBP, y=PrDens, color="C14 data"))
  p1 = p1 + geom_line(aes(x=calBP, y=predlog, color="model"))
  p1 = p1 + geom_line(aes(x=calBP, y=mavg, color="synthetic"))
  p1 = p1 + geom_ribbon(aes(x=calBP, ymin=m5, ymax=m95), fill="red", alpha=0.45)
  p1 = p1 + scale_color_manual(values=c("black", "blue", "red"))
  p1 = p1 + theme_bw(5) + xlab("Date [calBP]") + ylab("Density")
  if(!is.null(pval)) p1 = p1 + ggtitle(paste0("p-value: ", format(pval, digits=4)))
  p1 = p1 + guides(color="none")
  return(p1)
}

plot_one_acf = function(acfc) {
  p2 = ggplot(acfc)
  p2 = p2 + geom_line(aes(x=lag, y=c14, color="C14 data"))
  p2 = p2 + theme_bw(5) + xlab("Lag [years]") + ylab("ACF")
  p2 = p2 + guides(color="none")
  return(p2)
}

plot_one_res = function(spdc, acfc, ID, outdir, base_fn, pval = NULL) {
  p1 = plot_one_spd(spdc, pval)
  p2 = plot_one_acf(acfc)
  fn1 = paste0(outdir, "/", base_fn, "_spd_", ID)
  ggsave(paste0(fn1, ".pdf"), p1, width=1.6, height=1.2)
  ggsave(paste0(fn1, ".png"), p1, width=1.6, height=1.2, dpi=300)
  fn1 = paste0(outdir, "/", base_fn, "_acf_", ID)
  ggsave(paste0(fn1, ".pdf"), p2, width=1.6, height=1.2)
  ggsave(paste0(fn1, ".png"), p2, width=1.6, height=1.2, dpi=300)
}

plot_all_res_list = function(all_res, outdir, base_fn) {
  for(i in 1:length(all_res)) {
    plot_one_res(all_res[[i]]$spdc, all_res[[i]]$acfc, all_res[[i]]$ID,
                 outdir, base_fn)
  }
}

plot_all_res_dfs = function(all_spd, all_acf, outdir, base_fn, pval = NULL) {
  IDs = unique(all_spd$ID)
  for(ID in IDs) {
    spd = all_spd[all_spd$ID == ID,]
    acf = all_acf[all_acf$ID == ID,]
    pval1 = NULL
    if(!is.null(pval)) pval1 = pval$pval[pval$ID == ID]
    if(nrow(acf) == 0) error("Missing ACF!")
    plot_one_res(spd, acf, ID, outdir, base_fn, pval1)
  }
}




###################################################################
# new main results: normalized SPDs
for(outdir in outdir1) {
  for(base1 in c('res')) {
    acfm12 = read.csv(paste0(outdir, '/', base1, '_min.csv'))
    
    m1 = 2200 # maximum of x-axis on the histogram plots
    tmp1 = m1 %% 500
    tmp2 = m1 %/% 500
    if(tmp1 > 0) tmp2 = tmp2 + 1
    
    l1 = c(0,m1) # limits for the plots
    br1 = (0:tmp2)*500 # breaks for the plots
    
    for(rm1 in unique(acfm12$rm)) {
      acfm1 = acfm12[acfm12$rm == rm1,]
      acfm1$dx2 = 100*acfm1$dx
      acfm1$dy2 = 100*acfm1$dy
      
      p1 = ggplot(acfm1) + geom_histogram(aes(x=m), binwidth=100, boundary=0)
      p1 = p1 + theme_bw(6) + xlab("ACF first minimum location [years]") +
        ylab("Frequency")
      p1 = p1 + scale_x_continuous(limits=l1, breaks = br1)
      
      fn0 = paste0(outdir, '/acf_dist')
      fn0 = paste0(fn0, '_rm', rm1)
      
      ggsave(paste0(fn0, '.pdf'), p1, width=1.6, height=1.2)
      ggsave(paste0(fn0, '.png'), p1, width=1.6, height=1.2, dpi=300)
      
      # disaggregated results
      p1 = p1 + facet_grid(dx2~dy2, labeller = labeller(dx2 = function(x){
        return(paste0("dx = ", x, " km"))}, dy2 = function(x){
          return(paste0("dy = ", x, " km"))}))
      p1 = p1 + theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
      ggsave(paste0(outdir, '/c14_acf_grid_rm', rm1, '.pdf'), p1, width=90, height=80, units = "mm")
    }
  }
  
  # example plots
  base_fn = "res"
  dx = 0
  dy = 0
  all_spd = read.csv(paste0(outdir, '/', base_fn, '_spd_all_dx',
                            dx, '_dy', dy, '.csv'))
  pval = read.csv(paste0(outdir, '/', base_fn, '_pval.csv'))
  pval = pval[pval$dx == dx & pval$dy == dy,]
  outdir2 = paste0(outdir, '/region_figs_rm')
  dir.create(outdir2)
  
  for(base2 in c('zero')) {
    all_acf = read.csv(paste0(outdir, '/', base_fn, '_acf_', base2,
                              '_all_dx', dx, '_dy', dy, '.csv'))
    for(rm in c(0,100)) {
      acf1 = all_acf[all_acf$rm == rm,]
      spd1 = all_spd[all_spd$rm == rm,]
      pval1 = pval[pval$rm == rm,]
      plot_all_res_dfs(spd1, acf1, outdir2, paste0('res_', base2, '_rm',rm), pval1)
    }
  }
}


###############################################################
# results for coefficient of variation
for(outdir in outdir1) {
  tmp1 = read.csv(paste0(outdir, '/res_cv.csv'))
  for(rm in unique(tmp1$rm)) {
    tmp2 = tmp1[tmp1$rm == rm,]
    p1 = ggplot(tmp2)
    p1 = p1 + geom_histogram(aes(x=cv), binwidth = 0.04, boundary = 0)
    p1 = p1 + theme_bw(6) + xlab("Coefficient of variation") +
      ylab("Relative frequency") + scale_x_continuous(limits=c(0,1.2))
    fn1 = paste0(outdir, '/res_cv_rm', rm)
    ggsave(paste0(fn1, '.pdf'), p1, width=1.6, height=1.2)
    ggsave(paste0(fn1, '.png'), p1, width=1.6, height=1.2, dpi=300)
    
    # disaggregated results
    tmp2$dx2 = 100 * tmp2$dx
    tmp2$dy2 = 100 * tmp2$dy
    p1 = ggplot(tmp2)
    p1 = p1 + geom_histogram(aes(x=cv), binwidth = 0.04, boundary = 0)
    p1 = p1 + theme_bw(6) + xlab("Coefficient of variation") +
      ylab("Relative frequency")
    p1 = p1 + scale_x_continuous(limits=c(0,1.12), breaks = (0:4)*0.25)
    p1 = p1 + facet_grid(dx2~dy2, labeller = labeller(dx2 = function(x){
      return(paste0("dx = ", x, " km"))}, dy2 = function(x){
        return(paste0("dy = ", x, " km"))}))
    p1 = p1 + theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
    ggsave(paste0(outdir, '/c14_cv_grid_rm', rm, '.pdf'), p1, width=90, height=80, units = "mm")
  }
}



#####################################################################
# p-values
for(outdir in outdir1) {
  pval1 = read.csv(paste0(outdir, '/res_pval.csv'))
  for(rm in unique(pval1$rm)) {
    pval2 = pval1[pval1$rm == rm,]
  
    p1 = ggplot(pval2) + geom_histogram(aes(x=pval), binwidth = 0.025,
                    boundary = 0, fill="red", color="red", alpha=0.4)
    p1 = p1 + theme_bw(6) + scale_x_continuous(limits=c(0,1))
    p1 = p1 + xlab("p-value") + ylab("Frequency")
    fn1 = paste0(outdir, '/pval_rm', rm)
    ggsave(paste0(fn1, '.pdf'), p1, width=3.2,height=2.2)
    ggsave(paste0(fn1, '.png'), p1, width=3.2,height=2.2,dpi=300)
  }
}


####################################################################
# plot example SPDs and ACFs in only 10 regions together for the
# figures in the main text

colors = c("#5aa732","#005c96","#e75216", "#009ace","#ffd500",
           "#8671a7", "#e78502","#db6aa2","#007939","#9a2846","#815234","#000000")

base_fn = "res"
outdir = outdir1
dx = 0
dy = 0
all_spd = read.csv(paste0(outdir, '/', base_fn, '_spd_all_dx',
                          dx, '_dy', dy, '.csv'))
all_spd = all_spd[all_spd$rm == 0,]

base2 = 'zero'
all_acf = read.csv(paste0(outdir, '/', base_fn, '_acf_', base2,
                            '_all_dx', dx, '_dy', dy, '.csv'))
all_acf = all_acf[all_acf$rm == 0,]


# select 10 regions randomly
IDs1 = sample(unique(all_spd$ID), 10)
# note: use the IDs below to reproduce the result of the selection
# used in the paper
# IDs1 = c(2, 32, 9, 28, 39, 8, 33, 25, 10, 19)

tmp1 = all_spd[all_spd$ID %in% IDs1,]
p1 = ggplot(tmp1) + geom_line(size = 0.2, aes(x=calBP - 1950, y=PrDens,
                                        group=ID, color=factor(ID)))
p1 = p1 + scale_color_manual(values = colors) + theme_bw(6)
p1 = p1 + xlab("Date [BCE]") + ylab("SPD (arbitrary units)")
p1 = p1 + scale_x_reverse(limits=c(7000, 3000))
p1 = p1 + guides(color="none")
ggsave(paste0(outdir, "/fig2a.pdf"), p1, width=1.6, height=1.2)

tmp2 = all_acf[all_acf$ID %in% IDs1,]
p2 = ggplot(tmp2) + geom_line(size = 0.2, aes(x=lag, y=c14,
                                      group=ID, color=factor(ID)))
p2 = p2 + scale_color_manual(values = colors) + theme_bw(6)
p2 = p2 + xlab("Lag [years]") + ylab("ACF")
p2 = p2 + guides(color="none")
ggsave(paste0(outdir, "/fig2b.pdf"), p2, width=1.6, height=1.2)

