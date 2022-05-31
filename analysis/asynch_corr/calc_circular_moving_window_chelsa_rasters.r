library(raster)

min_t = raster('/global/scratch/users/drewhart/seasonality/CHELSA_bio6_1981-2010_V.2.1.tif')
ppt_s = raster('/global/scratch/users/drewhart/seasonality/CHELSA_bio15_1981-2010_V.2.1.tif')

#set the focal weights (NOTE: radius of 10 cells gives ~50km at equator)
fw_min_t <- focalWeight(min_t, 10, type='circle') 
fw_ppt_s <- focalWeight(ppt_s, 10, type='circle') 

# apply moving window, then save result
min_t_nmean <-focal(min_t, w=fw_min_t, fun=mean, na.rm=TRUE) 
writeRaster(min_t_nmean, '/global/scratch/users/drewhart/seasonality/CHELSA_bio6_1981-2010_V.2.1_10CELLRAD_NEIGHMEAN.tif', 'GTiff')

ppt_s_nstd <-focal(ppt_s, w=fw_min_t, fun=std, na.rm=TRUE) 
writeRaster(ppt_s_nstd, '/global/scratch/users/drewhart/seasonality/CHELSA_bio15-2010_V.2.1_10CELLRAD_NEIGHSTD.tif', 'GTiff')
