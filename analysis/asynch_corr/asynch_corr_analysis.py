import numpy as np
import rasterio as rio
from rasterio.crs import CRS
import rioxarray as rxr
import geopandas as gpd
import matplotlib.pyplot as plt

# load the spatial data, reproject to 8857 (Equal Earth)
data_dir = "/home/deth/Desktop/UCB/research/projects/seasonality/seasonal_asynchrony/asynch_analysis/data/"
asynch = rio.open(data_dir + '/global_asynch_result.tif')
vrm = rio.open(data_dir + 'vrm_10kmneigh_reproj_export.tif')
alt = rio.open(data_dir + 'alt_SRTM_reproj_export.tif')
countries = gpd.read_file(data_dir + 'countries/countries.shp')

# load the spatial data, transforming everything to Equal Earth


latlon = brick('./lonlat_export.tif')
alt = raster('./alt_SRTM_reproj_export.tif')

# mask out the oceans and other spots
vrm = mask(vrm, asynch[['asynch']])
latlon = mask(latlon, asynch[['asynch']])
alt = mask(alt, asynch[['asynch']])

# gather into a stack
vars = stack(asynch, latlon, alt, vrm)

# need to reproj?

# get sample points
pts = st_as_sf(spsample(as_Spatial(countries), n=1000, type='hexagonal'))

# coerce to a data.frame
df <- data.frame(na.omit(values(vars)))
colnames(df) = c('asynch', 'n', 'R2', 'lon', 'lat', 'alt', 'vrm')
head(df)

# non-spatial global regression
mod = lm(asynch ~ lon + lat + alt + vrm + alt*lat + vrm*lat, data=df)
summary(mod)

# NOTE: definitely need to account for spatial autocorrelation
# map residuals
map_resids = function(mod, df){
    resids = residuals(mod)
    breaks <- c(min(resids), min(resids)/3, 0, max(resids)/3, max(resids))  
    BlueRed <- colorRampPalette(c("darkblue", "white", "darkred"))
    map.resids <- SpatialPointsDataFrame(data=data.frame(resids), coords=cbind(df$lon,df$lat))
    map.resids <- st_as_sf(data.frame(resids=resids, lon=df$lon, lat=df$lat), coords=c('lon', 'lat'), crs=4326)
    print(map.resids)
    tm_shape(map.resids) + tm_dots(col='resids', size=0.5)
    #spplot(map.resids, cuts=breaks, col.regions=BlueRed(5), cex=1) 
}

map_resids(mod, df)
    
# NOTE: HOW TO RUN THIS FOR FULL DATASET?
# set row-step for subsetting df (for now)
step = 100
subdf = df[seq(1,nrow(df),step),]

# define the formula
formula = asynch ~ lon + lat + alt + vrm + alt*lat + vrm*lat

# build a null model
gls.null = gls(formula, data=subdf, method='ML')

# then build with spatial covariance structures
gls.gau = gls(formula, correlation=corGaus(form = ~ lat + lon, nugget=T),
              method='ML', data=subdf)

gls.sph = gls(formula, correlation=corSpher(form = ~ lat + lon, nugget=T),
              method='ML', data=subdf)

summary(gls.null)
map_resids(gls.null, subdf)

summary(gls.gau)
map_resids(gls.gau, subdf)

summary(gls.sph)
map_resids(gls.sph, df)



########################
# random forest modeling

# create explicit interaction var
df$vrm.x.lat = df$vrm * df$lat

# get training and test data
train_frac = 0.01
indices = sample(1:nrow(df), size = round(train_frac * nrow(df)))
trn = df[indices,]
tst = df[-indices,]

# set formula for RF
formrf = asynch ~ lon + lat + alt + vrm

# build the model
rf = randomForest(formrf, data=trn,
                  ntree=500, mtry=2, importance=TRUE)

# take a look at the result
varImpPlot(rf)

# make predictions
preds = predict(rf, tst[,-c(1,2,3)])
tst$err = tst[,1] - preds

ggplot(tst) +
    geom_point(aes(x=lon, y=lat, col=err), size=1) +
    #scale_color_gradientn(colours=wes_palette("Zissou1", 50, type='continuous'))
    scale_color_gradient2(low='#cc003d', mid='#dbdbdb', high='#009e64') +
    coord_map() + 
    theme_bw()

# make predicitons for full dataset (to map as raster)
full_preds = predict(rf, df[,-c(1,2,3)])
df$preds = full_preds
df$err = df[,1] - full_preds
dfrast <- rasterFromXYZ(df[, c('lon', 'lat', 'err')])
re_df = as.data.frame(dfrast, xy=T)
colnames(re_df) = c('lon', 'lat', 'error')

out = ggplot() +
        geom_raster(data=re_df, aes(x=lon, y=lat, fill=error)) +
        #scale_fill_gradient2(low='#940500', mid='#f5f5f5', high='#00138f') +
        scale_fill_cmocean(name='curl', direction=-1, end = 0.68) + 
        coord_quickmap() +
        theme_bw() +
        theme(legend.key.size = unit(4, 'cm'),
              legend.key.width = unit(3, 'cm'),
              legend.text = element_text(size=13),
              legend.title = element_text(size=16))
out

ggsave(out, file='rf_pred_error_map.png', width=35, height=45, units='cm', dpi=1000)


# run same model but with the explcit interaction term
formrf2 = asynch ~ lon + lat + alt + vrm + vrm.x.lat
rf2 = randomForest(formrf2, data=trn,
                   ntree=500, mtry=2, importance=TRUE)
varImpPlot(rf2)
preds2 = predict(rf2, tst[,-c(1,2,3)])
tst$err2 = tst[,1] - preds2
ggplot(tst) +
    geom_point(aes(x=lon, y=lat, col=err2), size=1) +
    #scale_color_gradientn(colours=wes_palette("Zissou1", 50, type='continuous'))
    scale_color_gradient2(low='red', mid='white', high='blue')


