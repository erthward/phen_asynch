library(raster)
library(ggplot2)
library(RColorBrewer)
library(nlme)
library(ape)
library(sp)
library(sf)
library(tmap)
library(gtools)
library(randomForest)
library(spatialEco)
library(SpatialML)
library(wesanderson)
library(ggthemes)


data_dir = "/home/deth/Desktop/UCB/research/projects/seasonality/seasonal_asynchrony/asynch_analysis/data/"


# load the spatial data, transforming everything to Equal Earth
asynch = brick(paste0(data_dir, 'global_asynch_result.tif'))
names(asynch) = c('asynch', 'R2', 'asynch_euc', 'R2_euc', 'n')

targ_crs = st_crs(proj4string(crs(asynch)))

countries = st_transform(st_read(paste0(data_dir, 'countries/countries.shp')), targ_crs)
countries = st_crop(countries, asynch)

vrm = raster(paste0(data_dir, 'vrm_10kmneigh_reproj_export.tif'))
names(vrm) = c('vrm')

dist2riv = raster(paste0(data_dir, 'dist_to_rivers_SIF.tif'))
names(dist2riv) = c('dist2riv')

alt = raster(paste0(data_dir, 'alt_SRTM_reproj_export.tif'))
names(alt) = c('alt')

lonlat = brick(paste0(data_dir, 'lonlat_export.tif'))
names(lonlat) = c('lon', 'lat')

#bioclim = stack(paste0(data_dir, 'bioclim/', list.files(paste0(data_dir, 'bioclim/'))))
#names(bioclim) = sub('wc2.1_2.5m_', '', names(bioclim))

for (var in c(vrm, alt, lonlat, dist2riv)){ #, bioclim)){
  if (!(st_crs(var) == st_crs(asynch))){
    var_name = deparse(substitute(var))
    assign(var_name, projectRaster(var, projectExtent(var, crs(asynch))))
  }
}

# resample the bioclim vars
#bioclim = resample(bioclim, asynch)

# write out, then read in, the resampled bioclim
#for (lyr in names(bioclim)){
#    writeRaster(bioclim[[lyr]],
#                paste0(data_dir, 'bioclim_resamp/', lyr, '.tif'),
#                format='GTiff')
#}

bioclim = stack(paste0(data_dir, 'bioclim_resamp/', list.files(paste0(data_dir, 'bioclim_resamp/'))))
bioclim = bioclim[[mixedsort(names(bioclim))]]

# crop and mask all rasters together
distriv = crop(dist2riv, asynch)
asynch = crop(asynch, vrm)
asynch = crop(asynch, dist2riv)
vrm = crop(vrm, asynch)
alt = crop(alt, asynch)
lonlat = crop(lonlat, asynch)
distriv = crop(dist2riv, asynch)
bioclim = crop(bioclim, asynch)

# resample dist2riv
dist2riv = resample(dist2riv, asynch)

# mask out the oceans and other spots
vrm = mask(vrm, asynch[[3]])
alt = mask(alt, asynch[[3]])
lonlat = mask(lonlat, asynch[[3]])
dist2riv = mask(dist2riv, asynch[[3]])
bioclim = crop(bioclim, asynch[[3]])

# get the lat*vrm intxn
latxvrm = (90 + lonlat[['lat']]) * vrm
names(latxvrm) = c('latxvrm')

# gather the covars into a stack
vars = stack(asynch, lonlat, vrm, alt, latxvrm, dist2riv, bioclim)





# grab just the Euclidean asynchrony layer
#euc_asynch = asynch[[3]]

#corr.hex <- raster.modified.ttest(euc_asynch, vars, sub.sample = TRUE, type = "hex")   
#corr.hex <- raster.modified.ttest(as(euc_asynch, 'SpatialPixelsDataFrame'),
#                                  as(vrm, 'SpatialPixelsDataFrame'),
#                                  sub.sample = TRUE, type = "hex")   

#head(corr.hex@data)
#bubble(corr.hex, "corr")






# get sample points
pts = st_as_sf(spsample(as_Spatial(countries), n=10000, type='hexagonal'))
df <- na.omit(data.frame(extract(vars, pts)))

# get training and test data
train_frac = 0.5
indices = sample(1:nrow(df), size = round(train_frac * nrow(df)))
trn = df[indices,]
tst = df[-indices,]

# prepare inverse-distance weighting matrix for Moran's I tests
geoMat<-pointDistance(SpatialPoints(tst[,c('lon', 'lat')]),longlat=TRUE) # pointDistance calculates pairwise geographic distances
geoMat[upper.tri(geoMat)] <- geoMat[lower.tri(geoMat)] # store these distances in a symmetrical matrix
spdata.IDW <- 1/geoMat # find the inverse distances
diag(spdata.IDW) <- 0 # fill in the diagonal elements



########################
# random forest modeling

# set formula for RF
formrf = asynch ~ lon + lat + alt + vrm

# build the model
rf = randomForest(formrf, data=trn,
                  ntree=500, mtry=2, importance=TRUE)

# take a look at the result
varImpPlot(rf)

# make predictions
preds1 = predict(rf, tst[,-c(1,2,3)])
tst$err1 = tst[,1] - preds1

# bubble plot
ggplot(tst) +
    geom_point(aes(x=lon, y=lat, col=err1), size=1) +
    #scale_color_gradientn(colours=wes_palette("Zissou1", 50, type='continuous'))
    scale_color_gradient2(low='#cc003d', mid='#dbdbdb', high='#009e64') +
    coord_map() + 
    theme_bw()

# and a not-so-manual bubbleplot
bubble(as_Spatial(st_as_sf(tst, coords = c('lon', 'lat'), crs=st_crs(countries))),
       'err1')

# make predictions for full dataset (to map as raster)
full_preds1 = predict(rf1, df[,-c(1,2,3)])
df$preds1 = full_preds1
df$err1 = df[,1] - full_preds1


#############################
# run same model but with the explcit interaction term
formrf2 = asynch ~ lon + lat + alt + vrm + latxvrm
rf2 = randomForest(formrf2, data=trn,
                   ntree=500, mtry=2, importance=TRUE)

varImpPlot(rf2)

preds2 = predict(rf2, tst[,-c(1,2,3)])
tst$err2 = tst[,1] - preds2
ggplot(tst) +
    geom_point(aes(x=lon, y=lat, col=err2), size=1) +
    #scale_color_gradientn(colours=wes_palette("Zissou1", 50, type='continuous'))
    scale_color_gradient2(low='#cc003d', mid='#dbdbdb', high='#009e64')
# make predicitons for full dataset (to map as raster)
full_preds2 = predict(rf2, df[,-c(1,2,3)])
df$preds2 = full_preds2
df$err2 = df[,1] - full_preds2


# Perform a global test of Moran's I on the point observations
Moran.I(tst$err2, spdata.IDW)


######################################
# TODO: CONTINUE SUSSING OUT RF RESULTS' SENSITIVITY
#       TO mtry, ntree, etc...
# run same model but without the lon and lat terms
formrf3 = asynch ~ alt + vrm + latxvrm
rf3 = randomForest(formrf3, data=trn,
                   ntree=5000, importance=TRUE)

varImpPlot(rf3)

preds3 = predict(rf3, tst[,-c(1,2,3)])
tst$err3 = tst[,1] - preds3
ggplot(tst) +
    geom_point(x=tst$lon, y=tst$lat, size=1, colour='black') +
    geom_point(aes(x=lon, y=lat, col=err3), size=1) +
    #scale_color_gradientn(colours=wes_palette("Zissou1", 50, type='continuous'))
    scale_color_gradient2(low='#cc003d', mid='#dbdbdb', high='#009e64')

# make predicitons for full dataset (to map as raster)
full_preds3 = predict(rf3, df[,-c(1,2,3)])
df$preds3 = full_preds3
df$err3 = df[,1] - full_preds3
ggplot(df) +
    geom_point(x=df$lon, y=df$lat, size=2, colour='black') +
    geom_point(aes(x=lon, y=lat, col=err3), size=2) +
    #scale_color_gradientn(colours=wes_palette("Zissou1", 50, type='continuous'))
    scale_color_gradient2(low='#cc003d', mid='#dbdbdb', high='#009e64')


# Perform a global test of Moran's I on the point observations
Moran.I(tst$err3, spdata.IDW)



# geographical random forest
#grf3 = grf(formrf3, trn, 400, 'adaptive', trn[,c('lon', 'lat')])



######################################
# include bioclim vars, then run rfs again
formrf4 = formula(paste('asynch ~ alt + vrm + latxvrm + ', paste(names(bioclim), collapse=' + ')))
rf4 = randomForest(formrf4, data=trn,
                   ntree=5000, importance=TRUE)

varImpPlot(rf4)

preds4 = predict(rf4, tst[,-c(1,2,3)])
tst$err4 = tst[,1] - preds4
ggplot(tst) +
    geom_point(x=tst$lon, y=tst$lat, size=1, colour='black') +
    geom_point(aes(x=lon, y=lat, col=err4), size=1) +
    #scale_color_gradientn(colours=wes_palette("Zissou1", 50, type='continuous'))
    scale_color_gradient2(low='#cc003d', mid='#dbdbdb', high='#009e64')

# make predicitons for full dataset (to map as raster)
full_preds4 = predict(rf4, df[,-c(1,2,3)])
df$preds4 = full_preds4
df$err4 = df[,1] - full_preds4
ggplot(df) +
    geom_point(x=df$lon, y=df$lat, size=2, colour='black') +
    geom_point(aes(x=lon, y=lat, col=err4), size=2) +
    #scale_color_gradientn(colours=wes_palette("Zissou1", 50, type='continuous'))
    scale_color_gradient2(low='#cc003d', mid='#dbdbdb', high='#009e64')
# and a not-so-manual bubbleplot
bubble(as_Spatial(st_as_sf(tst, coords = c('lon', 'lat'), crs=st_crs(countries))),
       'err4')



# Perform a global test of Moran's I on the point observations
Moran.I(tst$err4, spdata.IDW)


# geographical random forest
#grf4 = grf(formrf4, trn, 400, 'adaptive', trn[,c('lon', 'lat')])

######################################
# include all vars, then run rfs again
# (except drop intxn; seems stupid to include one in RF anyhow)
formrf5 = formula(paste('asynch ~ lon + lat + alt + vrm + dist2riv + ',
                        paste(names(bioclim), collapse=' + ')))
rf5 = randomForest(formrf5, data=trn,
                   ntree=5000, importance=TRUE)

varImpPlot(rf5)

preds5 = predict(rf5, tst[,-c(1,2,3)])
tst$err5 = tst[,1] - preds5
ggplot(tst) +
    geom_point(x=tst$lon, y=tst$lat, size=1, colour='black') +
    geom_point(aes(x=lon, y=lat, col=err5), size=1) +
    #scale_color_gradientn(colours=wes_palette("Zissou1", 50, type='continuous'))
    scale_color_gradient2(low='#cc003d', mid='#dbdbdb', high='#009e64')

# make predicitons for full dataset (to map as raster)
full_preds5 = predict(rf5, df[,-c(1,2,3)])
df$preds5 = full_preds5
df$err5 = df[,1] - full_preds5
ggplot(df) +
    geom_point(x=df$lon, y=df$lat, size=2, colour='black') +
    geom_point(aes(x=lon, y=lat, col=err5), size=2) +
    #scale_color_gradientn(colours=wes_palette("Zissou1", 50, type='continuous'))
    scale_color_gradient2(low='#cc003d', mid='#dbdbdb', high='#009e64')

# and a not-so-manual bubbleplot
bubble(as_Spatial(st_as_sf(tst, coords = c('lon', 'lat'), crs=st_crs(countries))),
       'err5')


# Perform a global test of Moran's I on the point observations
Moran.I(tst$err5, spdata.IDW)


# geographical random forest
grf5 = grf(formrf5, trn, 100, 'adaptive', trn[,c('lon', 'lat')], ntree=50)

med_rsqs = c()
for (i in seq(2136)){
    med_rsqs = c(med_rsqs, median(grf5[['Forests']][[i]][['rsq']]))
}
trn$grf_med_rsq = med_rsqs
ggplot(trn) +
    geom_point(x=trn$lon, y=trn$lat, size=2, colour='black') +
    geom_point(aes(x=lon, y=lat, col=grf_med_rsq), size=2) +
    #scale_color_gradientn(colours=wes_palette("Zissou1", 50, type='continuous'))
    scale_color_gradient2(low='#cc003d', mid='#dbdbdb', high='#009e64')





######################################
# hand-select bioclim vars, then run rfs again
# bio_7 instead of bio_4
# bio_18 and bio_19
# bio_13 instead of bio_16 (and instaed of bio_12)
# bio_14 instead of bio_17
# bio_5 and bio_6 instead of 11, 1, 3, etc
# bio_8 and bio_9
# bio_3

formrf6 = asynch ~ lon + lat + alt + vrm + dist2riv + bio_7 + bio_18 + bio_13 + bio_14 + bio_5 + bio_6 + bio_8 + bio_9 + bio_3
rf6 = randomForest(formrf6, data=trn,
                   ntree=5000, importance=TRUE)

varImpPlot(rf6)

preds6 = predict(rf6, tst[,-c(1,2,3)])
tst$err6 = tst[,1] - preds6
ggplot(tst) +
    geom_point(x=tst$lon, y=tst$lat, size=1, colour='black') +
    geom_point(aes(x=lon, y=lat, col=err6), size=1) +
    #scale_color_gradientn(colours=wes_palette("Zissou1", 50, type='continuous'))
    scale_color_gradient2(low='#cc003d', mid='#dbdbdb', high='#009e64')

# make predicitons for full dataset (to map as raster)
full_preds6 = predict(rf6, df[,-c(1,2,3)])
df$preds6 = full_preds6
df$err6 = df[,1] - full_preds6
ggplot(df) +
    geom_point(x=df$lon, y=df$lat, size=2, colour='black') +
    geom_point(aes(x=lon, y=lat, col=err6), size=2) +
    #scale_color_gradientn(colours=wes_palette("Zissou1", 50, type='continuous'))
    scale_color_gradient2(low='#cc003d', mid='#dbdbdb', high='#009e64')

# and a not-so-manual bubbleplot
bubble(as_Spatial(st_as_sf(tst, coords = c('lon', 'lat'), crs=st_crs(countries))),
       'err6')


# Perform a global test of Moran's I on the point observations
Moran.I(tst$err6, spdata.IDW)


# geographical random forest
grf6 = grf(formrf6, trn, 100, 'adaptive', trn[,c('lon', 'lat')], ntree=50)

med_rsqs = c()
for (i in seq(2138)){
    med_rsqs = c(med_rsqs, median(grf6[['Forests']][[i]][['rsq']]))
}
trn$grf_med_rsq6 = med_rsqs
ggplot(trn) +
    geom_point(x=trn$lon, y=trn$lat, size=2, colour='black') +
    geom_point(aes(x=lon, y=lat, col=grf_med_rsq6), size=2) +
    #scale_color_gradientn(colours=wes_palette("Zissou1", 50, type='continuous'))
    scale_color_gradient2(low='#cc003d', mid='#dbdbdb', high='#009e64')






