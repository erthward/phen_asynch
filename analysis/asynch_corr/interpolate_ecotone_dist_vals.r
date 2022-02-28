library(raster)
library(akima)

dist_to_ecotone = raster('/media/deth/SLAB/seasonality/other/rf_vars/dist_to_ecotone_raw.tif')

df = as.data.frame(dist_to_ecotone, xy=T)

yes =  df[!is.na(df$dist_to_ecotone_raw),]
no = df[is.na(df$dist_to_ecotone_raw),]

x = yes$x
y = yes$y
z = yes$dist_to_ecotone_raw
xo = df$x
yo = df$y

all_vals = interp(x, y, z, xo, yo, extrap=T)
