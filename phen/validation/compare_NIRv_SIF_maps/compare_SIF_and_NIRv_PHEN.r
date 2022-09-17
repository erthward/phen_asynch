library(sf)
library(raster)
librbary(reshape)
library(RColorBrewer)
library(latex2exp)
library(ggplot2)

# simple diverging palette
pal = colorRampPalette(c('red', 'white', 'blue'))

# data directory
data.dir = "/home/deth/Desktop/CAL/research/projects/seasonality/results/maps/"

# load country boundaries
world = st_read(paste0(data.dir, 'NewWorldFile_2020.shp'))

# load both rasters
sif = stack(paste0(data.dir, 'SIF_global_coeffs.tif'))
nirv = stack(paste0(data.dir, 'NIRv_global_coeffs.tif'))

# coerce to vectors of single values
sif_vals = melt(as.data.frame(sif)[, 2:5])$value
nirv_vals = melt(as.data.frame(nirv)[, 2:5])$value

# print result of regression of SIF on NIRV
# (minus coefficient, because two variables should be have same values
#  to begin with)
mod = lm(nirv_vals ~ sif_vals)
print(summary(mod))
