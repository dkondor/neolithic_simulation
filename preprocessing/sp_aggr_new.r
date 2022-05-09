# sp_aggr_new.r -- create overlapping tilings of approximately
# square regions over the study area

do_plot = TRUE # note: set this to FALSE to avoid creating example plots
if(do_plot) {
  library(ggplot2)
  library(ggmap)
  library(rworldmap)
}

base_dir = '~/CSH/HoloSim/data/'

# sptial extents of the study area
xmin1 = -10.5
xmax1 = 33.5
ymin1 = 33.5
ymax1 = 62


###############################################################
# 1. function definitions

# create a set of quasi-rectangular regions covering the area
# given by the coordinates, with res being the side in kms
create_grid = function(xmin, ymin, xmax, ymax, res) {
  R = 6371 # radius of the Earth in kms
  D = 2*pi*R # approximate circumference
  
  res1 = data.frame()
  ystep = 360.0 * (res / D) # step in latitude
  id1 = 0 # cell ID
  
  y = ymin
  while(y < ymax) {
    y2 = y + ystep
    ym = y + ystep / 2.0
    xstep = ystep / cos(pi*ym/180.0)
    x1 = seq(xmin, xmax, by=xstep)
    n1 = length(x1)
    tmp1 = data.frame(lonmin = x1, lonmax = x1 + xstep,
                      latmin = y, latmax = y2, ID = id1 + 1:n1)
    res1 = rbind(res1, tmp1)
    y = y2
    id1 = id1 + n1
  }
  return(res1)
}

# create a series of shifted tilings
create_grids = function(xmin, ymin, xmax, ymax, res, delta) {
  R = 6371 # radius of the Earth in kms
  D = 2*pi*R # approximate circumference
  ystep = 360.0 * (delta / D) # step in latitude
  xstep = ystep / cos(pi*ymin/180)
  nsteps = ceiling(res / delta) - 1
  
  grids = data.frame()
  for(i in 0:nsteps) for(j in 0:nsteps) {
    xmin1 = xmin - i*xstep
    ymin1 = ymin - j*ystep
    tmpgrid = create_grid(xmin1, ymin1, xmax, ymax, res)
    tmpgrid$dx = i
    tmpgrid$dy = j
    tmpgrid$xmin = xmin1
    tmpgrid$ymin = ymin1
    grids = rbind(grids, tmpgrid)
  }
  
  return(grids)
}


##################################################################
# 2. create and save a set of tilings

grids1000 = create_grids(xmin1, ymin1, xmax1, ymax1, 1000, 100)
grids700 = create_grids(xmin1, ymin1, xmax1, ymax1, 700, 100)
grids500 = create_grids(xmin1, ymin1, xmax1, ymax1, 500, 100)

grids_all = data.frame()
tmp1 = grids500
tmp1$res = 500
grids_all = rbind(grids_all, tmp1)
tmp1 = grids700
tmp1$res = 700
grids_all = rbind(grids_all, tmp1)
tmp1 = grids1000
tmp1$res = 1000
grids_all = rbind(grids_all, tmp1)

save(grids_all, file = paste0(base_dir, "new_grids_all.RData"))
write.csv(grids_all, paste0(base_dir, "new_grids_all.csv"), row.names=FALSE)


#######################################################################
# 3. plot an example of the tiling

# helper function
fortify_grid = function(grid1) {
  res1 = data.frame()
  for(i in 1:nrow(grid1)) {
    x = grid1[i,]
    tmp1 = data.frame(ID=x$ID, group=x$ID, order = 1:5,
                      x = c(x$lonmin, x$lonmax, x$lonmax, x$lonmin, x$lonmin),
                      y = c(x$latmin, x$latmin, x$latmax, x$latmax, x$latmin))
    res1 = rbind(res1, tmp1)
  }
  return(res1)
}


# get the base map
base = getMap(resolution="low") #extract basemap
base2 = fortify(base)
base2 = base2[base2$long >= -50 & base2$long < 80 & base2$lat > 25,]
base2 = base2[order(base2$group, base2$order),]


tmp1 = grids500[grids500$dx == 0 & grids500$dy == 0,]

regions1 = fortify_grid(tmp1)
p1 = ggplot(base2) + geom_polygon(aes(x=long, y=lat, group=group),
                                  color="grey", fill="black", size=0.1)
p1 = p1 + geom_polygon(data = regions1, aes(x=x, y=y, group=group),
                       color="red", fill="red", alpha=0.4)
p1 = p1 + theme_nothing() + coord_map(xlim=c(-10.8, 44), ylim=c(34,64),
                                      clip="on")
fn1 = paste0(base_dir, "grid_map1")
ggsave(paste0(fn1, ".png"), p1, width=3.2, height=3.2, dpi=300)
ggsave(paste0(fn1, ".pdf"), p1, width=2.4, height=2.2)


p1 = p1 + coord_map(projection = "albers", lat0 = 34, lat1 = 64,
                    xlim=c(-10, 33), ylim=c(34,64), clip="on")
fn1 = paste0(base_dir, "grid_map2")
ggsave(paste0(fn1, ".png"), p1, width=3.2, height=3.2, dpi=300)
ggsave(paste0(fn1, ".pdf"), p1, width=2.4, height=2.2)




