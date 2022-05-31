library(raster)
library(gstat)


# input and output dirs:
# if on laptop
if (strsplit(getwd(), '/')[[1]][2] == 'home'){
  on.laptop=T
  data.dir = '/media/deth/SLAB/seasonality/other/rf_vars'
  analysis.dir = '/media/deth/SLAB/seasonality/other/rf_vars'
  # if on Savio
} else {
  on.laptop=F
  data.dir = '/global/scratch/users/drewhart/seasonality/'
  # data.dir = '/global/home/groups/fc_landgen/' # directory for Lauren 
  # analysis.dir = '/global/home/users/ldimaggio/ondemand/' # directory for Lauren 
  analysis.dir = '/global/scratch/users/drewhart/seasonality/'
}


# rename all bands 
names = c('phn.asy', 'tmp.min.asy', 'tmp.max.asy',
          'tmp.min.mea', 'ppt.asy',
          'ppt.sea', 'def.asy', 'cld.asy', 'vrm.med',
          'riv.dis', 'eco.dis')

# load rasters of prepped variables
vars = brick(paste0(data.dir, "/asynch_model_all_vars.tif"))
names(vars) = names

# convert to SpatialPixelsDataFrame, and drop NAs
df = as(vars, "SpatialPixelsDataFrame")
df = df[complete.cases(df@data),]

# only take every 1500 rows
# (should give a sample of about 1800, and thus (n*(n-1))/2=~1.6M variogram points to be calculated)
df = df[seq(1, nrow(df), 1500),]

# fit variogram
v = variogram(phn.asy~1, df)
(f = fit.variogram(v, vgm("Sph"))) # NOTE: change 'Sph' to fit diff functional form
cat('SILL: ', f$psill[2])
cat('RANGE: ', f$range[2])
cat('NUGGET: ', f$psill[1])

# plot it
jpeg('asynch_variogram.jpg', width = 900, height = 700, units='px', quality=90)
plot(v,f)
dev.off()

