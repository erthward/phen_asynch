library(raster)
library(sf)
library(tmap)
tmap_mode('plot')

# load asynch map
asynch = raster('/home/deth/Desktop/tmp_GPP_seasonality_validation/NIRvP_euc_asynch.tif')
# floor at 0
asynch[asynch<0] = 0

# load fluxnet sites
fluxnet = read.csv('./FLUXNET_sites.csv')
# drop rows missing coords
# and drop rows missing FLUXNET2015 data
fluxnet = fluxnet[!is.na(fluxnet$LOCATION_LONG),] 
fluxnet = fluxnet[!is.na(fluxnet$FLUXNET2015),] 

# coerce fluxnet to sf object and check CRS's line up
flux_sf = st_as_sf(fluxnet, coords = c("LOCATION_LONG", "LOCATION_LAT"), crs=4326)
st_crs(asynch) == st_crs(flux_sf)

# overlay
pdf('asynch_fluxnet_overaly.pdf', width=20, height=5)
overlay = tm_shape(asynch) + 
  tm_raster(style = "quantile", n = 12, title = "seasonal asynchrony",
            palette = colorRampPalette( c("darkolivegreen4","yellow", "brown"))(12),
            legend.hist = TRUE)+
  tm_legend(outside = TRUE, hist.width = 2) + 
tm_shape(flux_sf) + 
    tm_dots(size=0.0025)
overlay
dev.off()
