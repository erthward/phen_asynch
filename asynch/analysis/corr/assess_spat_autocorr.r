library(raster)
library(gstat)


# use all vars, or just asynch?
use_all_vars = F


# input and output dirs:
# if on laptop
if (strsplit(getwd(), '/')[[1]][2] == 'home'){
  on.laptop=T
  data.dir = '/home/deth/Desktop/CAL/research/projects/seasonality/results/maps'
  analysis.dir = '/media/deth/SLAB/seasonality/other/rf_vars'
  # if on Savio
} else {
  on.laptop=F
  data.dir = '/global/scratch/users/drewhart/seasonality/'
  # data.dir = '/global/home/groups/fc_landgen/' # directory for Lauren 
  # analysis.dir = '/global/home/users/ldimaggio/ondemand/' # directory for Lauren 
  analysis.dir = '/global/scratch/users/drewhart/seasonality/'
}


if (use_all_vars){
  # load rasters of prepped variables
  rast = brick(paste0(data.dir, "/asynch_model_all_vars.tif"))
  names = c('phn.asy', 'tmp.min.asy', 'tmp.max.asy',
            'tmp.min.mea', 'ppt.asy',
            'ppt.sea', 'def.asy', 'cld.asy', 'vrm.med',
            'riv.dis', 'eco.dis')
  names(rast) = names
  } else {
   rast = brick(paste0(data.dir, '/NIRv_global_asynch.tif'))
   names(rast) = c('phn.asy.EMPTY', 'phn.asy.sd.EMPTY', 'phn.asy', 'phn.asy.sd', 'n')
}

# get cells that are not missing
notna <- which(!is.na(values(rast[['phn.asy']])))

# save all the ranges
ranges = c()

for (i in seq(1000)){
  cat('ITERATION NUMBER ', i, '...\n\n')
  # draw a random non-NA cell in the raster
  samp = sample(notna, 1)

  # get its coordinates
  coords = as.data.frame(xyFromCell(rast[['phn.asy']], samp))
  coordinates(coords) = ~ x + y
  proj4string(coords) <- "+init=epsg:4326"

  # create a circle centered on that cell, then crop and mask to it
  coords_utm = spTransform(coords, CRS("+init=epsg:8857"))
  buff_utm = rgeos::gBuffer(coords_utm, width = 5e5, quadsegs = 250L)
  buff = spTransform(buff_utm, CRS(proj4string(coords)))
  subrast = mask(crop(rast, buff), buff)

  # convert to SpatialPixelsDataFrame, and drop NAs
  df = as(subrast, "SpatialPixelsDataFrame")
  df = df[,3]
  df = df[complete.cases(df@data), ]

  # calculate range of semivariogram for that area
  v = variogram(phn.asy~1, spTransform(df, CRS("+init=epsg:8857")))
  (f = fit.variogram(v, vgm("Sph"))) # NOTE: change 'Sph' to fit diff functional form
  ranges = c(ranges, f$range[2])

  if (i == 1){
    # plot it
    jpeg('asynch_variogram.jpg', width = 900, height = 700, units='px', quality=90)
    plot(v,f)
    dev.off()
  }

}
cat('---------------------------------------\n\n')
cat('MEDIAN RANGE: ', median(ranges)/1000, ' km\n')
cat('MEAN RANGE: ', mean(ranges)/1000, ' km\n')
cat('STD RANGE: ', sd(ranges)/1000, ' km\n')

jpeg('hist_asynch_ranges.jpg', width=900, heigh=900, units='px', quality=90)
hist(ranges, breaks=100)
dev.off()

write.csv(as.data.frame(ranges), 'estimated_ranges.csv', row.names=F)
