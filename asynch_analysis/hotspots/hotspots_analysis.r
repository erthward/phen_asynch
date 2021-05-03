library(sf)
library(raster)
library(tmap)
library(RColorBrewer)

tmap_mode('view')

# load data
spots = st_read('../../../data/hotspots/hotspots_2016_1.shp')
asynch = brick('../data/global_asynch_result.tif')
asynch_rast = asynch[[3]][,,drop=FALSE]

# drop invalid hotspot geometries
valids = c()
for (i in seq(nrow(spots))){
    valids = c(valids, st_is_valid(spots[i, c('geometry')])) 
}
valids[is.na(valids)] = FALSE
spots = spots[valids, ]

# reproject hotspots
spots = st_transform(spots, crs(asynch_rast))

# overlay
#tm_shape(asynch_rast) + 
#    tm_raster() +
#tm_shape(spots) + 
#    tm_polygons()
plot(asynch_rast, col=colorRampPalette(c("white", "yellow", "orange",
                                         "darkorange", "darkorange2",
                                         "darkorange3","red", "darkred"))(20))
plot(spots$geometry, add=T)


# compare asynch values inside and outside hotspots
