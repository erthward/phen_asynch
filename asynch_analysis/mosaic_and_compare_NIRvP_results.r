library(raster)

# load the coeffs and asynch results for both the NA and SA regions
na_c = brick('./NA_agg/NA_agg_coeffs.tif')
na_a = brick('./NA_agg/NA_agg_asynch.tif')
sa_c = brick('./SA_agg/SA_agg_coeffs.tif')
sa_a = brick('./SA_agg/SA_agg_asynch.tif')

# mosaic them
nirvp_c = mosaic(na_c, sa_c, fun=mean)
nirvp_a = mosaic(na_a, sa_a, fun=mean)

# write them out
writeRaster(nirvp_c, 'NA_SA_agg_coeffs.tif', format='GTiff')
writeRaster(nirvp_a, 'NA_SA_agg_asynch.tif', format='GTiff')

# load the SIF asynch map
sif_a = brick('../seasonal_asynchrony/asynch_analysis/data/global_asynch_result.tif')

# crop SIF asynch
sif_a = crop(sif_a, nirvp_a)

# reproject and resample the SIF asynch
sif_a = projectRaster(sif_a, projectExtent(sif_a, nirvp_a@crs))
sif_a = resample(sif_a, nirvp_a, 'bilinear')

# mutually mask missing values
sif_a[is.na(nirvp_a)] = NA
nirvp_a[is.na(sif_a)] = NA

# scatter and assess their asynch, R^2, and n values
par(mfrow=c(1,3))
plot(sif_a[[3]][,], nirvp_a[[3]][,],
     main='asynch', xlab='SIF', ylab='NIRvP',
     cex.lab=3, cex.axis=3, cex.main=3, cex.sub=3)
plot(sif_a[[4]][,], nirvp_a[[4]][,],
     main='asynch R^2s', xlab='SIF', ylab='NIRvP',
     cex.lab=3, cex.axis=3, cex.main=3, cex.sub=3)
plot(sif_a[[5]][,], nirvp_a[[5]][,],
     main='asynch sample sizes', xlab='SIF', ylab='NIRvP',
     cex.lab=3, cex.axis=3, cex.main=3, cex.sub=3)

mod = lm(nirvp_a[[3]] ~ sif_a[[3]])
summary(mod)


