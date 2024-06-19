library(sp)                   # spatial data
library(raster)               # raster data
library(terra)                # newer raster data
library(sf)                   # newer spatial data
library(spdep)                # spatial autocorrelation
library(rsample)              # function for stratified random sampling
library(randomForest)         # global RFs
library(RRF)                  # fast RF var selection (w/ conservative results in Bag et al. 2022)
library(SpatialML)            # GWRFs
library(GWmodel)              # GW models
library(vip)                  # var importance plots
library(pdp)                  # partial depend. plots (& ICE curves)
library(DALEX)                # feature importance
library(ggplot2)              # plotting
library(ggthemes)             # plot themes
library(grid)                 # textGrob for plot grids
library(cowplot)              # easy plot gridding
library(tmap)                 # mapping
library(maps)                 # countries map data
library(RColorBrewer)         # colors
library(cmocean)              # cmocean palettes
library(dplyr)                # reshaping dfs
library(caret)                # Recursive Feature selection
library(rfUtilities)          # Jeff Evans R package for model selection




##########################################################################
# SETUP 


#######################
# SET BEHAVIORAL PARAMS
#######################

args = commandArgs(trailingOnly=T)

# variable for which to prep data
asynch.var = args[1]
#asynch.var = 'NIRv'
cat(paste0('\nVAR: ', asynch.var, '\n'))

# set neighborhood radius (in km) to use for asynch analysis
neigh.rad = args[2]
cat(paste0('\nNEIGH RAD: ', neigh.rad, '\n'))

# input and output dirs:
  # if on laptop
if (strsplit(getwd(), '/')[[1]][2] == 'home'){
       on.laptop=T
       data.dir = '/media/deth/SLAB/diss/3-phn/final_maps_and_results/'
  # if on Savio
} else {
       on.laptop=F
       data.dir = '/global/scratch/users/drewhart/seasonality/rf_data/'
}

# seed
seed.num = 12345
set.seed(seed.num)

# number of strata and proportion of each stratum
# to use for stratified random sampling
n.strata = 5
strat.samp.prop = 0.25

# training data fraction
# NOTE: use 60% for training, 40% for testing, which should be plenty,
#       given how large a dataset we have
train.frac = 0.6


############
# HELPER FNS
############

align.rast.x.to.y = function(rast.x, rast.y, mask.it=F){
   # resample
   rast.x.resamp = raster::resample(rast.x, rast.y, "bilinear")
   # crop
   ext = extent(rast.y)
   rast.x.out = raster::crop(rast.x.resamp, ext)
   # mask?
   if (mask.it){
      rast.x.out = raster::mask(rast.x.out, rast.y)
   }
   return(rast.x.out)
}

read.file = function(var.str, asynch.file=T, align.to=NA, mask.it=F){
   # candidate files
   files = list.files(data.dir)
   # find right file
   if (asynch.file){
      file.name = files[grep(paste0(var.str, '_asynch_', as.character(neigh.rad), 'km'), files)]
      stopifnot(length(file.name) == 1)
      # NOTE: TAKING 3RD LYR, THE EUC DIST-BASED ASYNCHRONY VALUE
      data = brick(paste0(data.dir, '/', file.name))[[1]]
   } else {
      file.name = files[grep(var.str, files)]
      stopifnot(length(file.name) == 1)
      data = raster(paste0(data.dir, '/', file.name))
   }
   # align to the align.to target raster, if align.to is a raster
   if (class(align.to)[1] == 'RasterLayer'){
      data = align.rast.x.to.y(data, align.to, mask.it=mask.it) 
   }
   return(data)
}



##########################################################################
# DATA PREP 

  
####################
# LOAD RESPONSE VARS
####################

# load phenological asynch data
phn.asy = read.file(paste0(asynch.var, '_STRICT'), T)


#################
# LOAD PREDICTORS
#################

# asynchrony in min and max temperatures
tmp.min.asy = read.file('tmmn', asynch.file=T,
                    align.to=phn.asy, mask.it=F)
tmp.max.asy = read.file('tmmx', asynch.file=T,
                    align.to=phn.asy, mask.it=F)

# asynchrony in precipitation
ppt.asy = read.file('pr', asynch.file=T,
                    align.to=phn.asy, mask.it=F)

# asynchrony in climatic water deficit
def.asy = read.file('def', asynch.file=T,
                    align.to=phn.asy, mask.it=F)

# asynchrony in cloud cover
cld.asy = read.file('cloud', asynch.file=T,
                    align.to=phn.asy, mask.it=F)

# vector ruggedness metric, ~100km agg med and sd
# NOTE: calculated as fixed pixels, not moving windows
# NOTE: masking isn't necessary because incomplete rows of the stacked
#       variables are dropped later on, but makes it easier to inspect
#       the individual layers manually now because of the crazy blocky artefacts
#       their dataset has surrounding the continents
vrm.med = read.file('vrm_100KMmd', asynch.file=F,
                    align.to=phn.asy, mask.it=T)

# 100 km neighborood entropy in MODIS IGBP vegetation type
# (reclassed to forest, shrubland, savanna, grassland, permanent wetland, or invalid)
veg.ent = read.file('MODIS_IGBP_veg_entropy', asynch.file=F,
                    align.to=phn.asy, mask.it=F)
# mask to other datasets
veg.ent = raster::mask(veg.ent, phn.asy)



###########
# PREP DATA
###########

# gather into a stack
vars = stack(phn.asy,
             tmp.min.asy,
             tmp.max.asy,
             ppt.asy,
             def.asy,
             cld.asy,
             vrm.med,
             veg.ent)

par(mfrow=c(3,3))
for (lyr in c(phn.asy,
             tmp.min.asy,
             tmp.max.asy,
             ppt.asy,
             def.asy,
             cld.asy,
             vrm.med,
             veg.ent)){
   plot(lyr)#@extent@xmax)
}

par(mfrow=c(1,1))

# rename all bands 
names = c('phn.asy',
          'tmp.min.asy',
          'tmp.max.asy',
          'ppt.asy',
          'def.asy',
          'cld.asy',
          'vrm.med',
          'veg.ent')

names(vars) = names

# write stack to file
raster::writeRaster(vars, paste0(data.dir, "/asynch_model_all_vars_",
                                 asynch.var, '_',
                                 as.character(neigh.rad), "km.tif"), overwrite=T)

# coerce to a data.frame
df = as.data.frame(vars, xy=T)

# drop NAs
df = df[complete.cases(df),]

# scale variables
df[,3:ncol(df)] = scale(df[,3:ncol(df)])

# take a stratified random sample, stratifying by phn.asy values
# (it has a somewhat right-skewed distribution, and at any rate,
#  important that model is trained on data including high asynch values)
df.strat = df %>%
   mutate(strata = base::cut(phn.asy, breaks=n.strata, labels=seq(n.strata))) %>%
   group_by(strata) %>%
   slice_sample(prop=strat.samp.prop) %>%
   ungroup() %>%
   as.data.frame()
df.strat = df.strat[, !colnames(df.strat) %in% c('strata')]

# write data.frame to file (for reload in 'analyze' mode),
write.csv(df, paste0(data.dir, "/asynch_model_all_vars_prepped_",
                     asynch.var, '_',
                     as.character(neigh.rad), "km.csv"), row.names=F)
write.csv(df.strat, paste0(data.dir, "/asynch_model_all_vars_prepped_strat_",
                           asynch.var, '_',
                           as.character(neigh.rad), "km.csv"), row.names=F)

