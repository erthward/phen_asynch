library(raster)

min_t = raster('/global/scratch/users/drewhart/seasonality/CHELSA_bio6_1981-2010_V.2.1_5km.tif')
ppt_s = raster('/global/scratch/users/drewhart/seasonality/CHELSA_bio15_1981-2010_V.2.1_5km.tif')

#set the focal weights (NOTE: radius of 10 cells gives ~50km at equator)
fw_min_t <- focalWeight(min_t, 10, type='circle') 
fw_min_t[fw_min_t > 0] <- 1
fw_ppt_s <- focalWeight(ppt_s, 10, type='circle') 
fw_ppt_s[fw_ppt_s > 0] <- 1

# apply moving window, then save result
min_t_nmean <-focal(min_t, w=fw_min_t, fun=mean, na.rm=TRUE) 
writeRaster(min_t_nmean, '/global/scratch/users/drewhart/seasonality/CHELSA_bio6_1981-2010_V.2.1_5km_10CELLRAD_NEIGHMEAN.tif', 'GTiff', overwrite=T)

ppt_s_nsd <-focal(ppt_s, w=fw_min_t, fun=sd, na.rm=TRUE) 
writeRaster(ppt_s_nsd, '/global/scratch/users/drewhart/seasonality/CHELSA_bio15-2010_V.2.1_5km_10CELLRAD_NEIGHSD.tif', 'GTiff', overwrite=T)
