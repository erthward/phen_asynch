library(sp)                   # spatial data
library(raster)               # raster data
library(terra)                # newer raster data
library(sf)                   # newer spatial data
library(spdep)                # spatial autocorrelation
library(rsample)              # function for stratified random sampling
library(randomForest)         # global RFs
library(RRF)                  # fast RF var selection (w/ conservative results in Bag et al. 2022)
#library(ranger)               # faster RFs
#library(h2o)                  # distributed RFs (on cloud)
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
#library(wesanderson)          # wes palettes
library(dplyr)                # reshaping dfs
library(caret)                # Recursive Feature selection
library(rfUtilities)          # Jeff Evans R package for model selection




##########################################################################
# SETUP 


#######################
# SET BEHAVIORAL PARAMS
#######################

# input and output dirs:
  # if on laptop
if (strsplit(getwd(), '/')[[1]][2] == 'home'){
       on.laptop=T
       data.dir = '/media/deth/SLAB/seasonality/other/rf_vars'
       analysis.dir = '/media/deth/SLAB/seasonality/other/rf_vars'
  # if on Savio
} else {
       on.laptop=F
       data.dir = '/global/scratch/users/drewhart/seasonality/'
       # data.dir = '/global/home/groups/fc_landgen/' # directory for Lauren 
       # analysis.dir = '/global/home/users/ldimaggio/ondemand/' # directory for Lauren 
       analysis.dir = '/global/scratch/users/drewhart/seasonality/'
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
      file.name = files[grep(paste0(var.str, '_global_asynch'), files)]
      stopifnot(length(file.name) == 1)
      # NOTE: TAKING 3RD LYR, THE EUC DIST-BASED ASYNCHRONY VALUE
      data = brick(paste0(data.dir, '/', file.name))[[3]]
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

# SIF-based asynch
#phn.asy = read.file('SIF', asynch.file=T)

# NIRv-based asynch
phn.asy = read.file('NIRv', T)


#################
# LOAD PREDICTORS
#################

# asynchrony in min temperature
tmp.min.asy = read.file('tmmn', asynch.file=T,
                    align.to=phn.asy, mask.it=F)

# 50km neighborhood mean daily min temp, coldest month
tmp.min.nmn = read.file('CHELSA_bio6_1981-2010_V.2.1_5km_10CELLRAD_NEIGHMEAN',
                        asynch.file=F, align.to=phn.asy, mask.it=F)

# asynchrony in precipitation
ppt.asy = read.file('pr', asynch.file=T,
                    align.to=phn.asy, mask.it=F)

# 50km neighborhood standard deviation in precipitation seasonality
ppt.sea.nsd = read.file('CHELSA_bio15_1981-2010_V.2.1_5km_10CELLRAD_NEIGHSD',
                    asynch.file=F, align.to=phn.asy, mask.it=F)

# asynchrony in climatic water deficit
def.asy = read.file('def', asynch.file=T,
                    align.to=phn.asy, mask.it=F)

# asynchrony in cloud cover
cld.asy = read.file('cloud', asynch.file=T,
                    align.to=phn.asy, mask.it=F)

# vector ruggedness metric, ~50km agg med and sd
# NOTE: CALCULATED AS FIXED PIXELS, NOT MOVING WINDOWS!
# NOTE: masking isn't necessary because incomplete rows of the stacked
#       variables are dropped later on, but makes it easier to inspect
#       the individual layers now because of the crazy blocky artifacts
#       surrounding continents in their dataset
vrm.med = read.file('vrm_50KMmd', asynch.file=F,
                    align.to=phn.asy, mask.it=T)

# distance from rivers
# NOTE: INCLUDE? IF SO, NEIGH MEAN AND/OR SD?
riv.dis = read.file('dist_to_rivers', asynch.file=F,
                    align.to=phn.asy, mask.it=F)

# distance from ecotones
# NOTE: has gaps in it because GEE didn't calculate distance beyond a max,
#       and too computationally intensive to quickly extrapolate into them,
#       so for now just backfilling with the max
#       TODO: come back to this? even though it's likely not that influential...
eco.dis = read.file('dist_to_ecotone', asynch.file=F,
                    align.to=phn.asy, mask.it=F)
# NOTE: backfill NAs with max val
eco.dis[is.na(eco.dis)] = cellStats(eco.dis, max)
# then mask to other datasets
eco.dis = raster::mask(eco.dis, phn.asy)



###########
# PREP DATA
###########

# gather into a stack
vars = stack(phn.asy, tmp.min.asy, tmp.min.nmn,
             ppt.asy, ppt.sea.nsd, def.asy, cld.asy, vrm.med,
             riv.dis, eco.dis)

# aggregate to coarser res, if working on laptop
if (on.laptop){
   vars <- aggregate(vars, fact=12)
}

# rename all bands 
names = c('phn.asy', 'tmp.min.asy',
          'tmp.min.nmn', 'ppt.asy',
          'ppt.sea.nsd', 'def.asy', 'cld.asy', 'vrm.med',
          'riv.dis', 'eco.dis')

names(vars) = names

# write stack to file
raster::writeRaster(vars, paste0(data.dir, "/asynch_model_all_vars.tif"), overwrite=T)
#vars = rast(vars)
#terra::writeRaster(vars, paste0(data.dir, "/asynch_model_all_vars.tif"), overwrite=T)

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
write.csv(df, paste0(data.dir, "/asynch_model_all_vars_prepped.csv"), row.names=F)
write.csv(df.strat, paste0(data.dir, "/asynch_model_all_vars_prepped_strat.csv"), row.names=F)

