# CHRIS PACIOREK'S WORKAROUND TO FIX THE rgdal/sp/sf ISSUE I RAN INTO:
#Sys.setenv(GDAL_DATA = "/global/home/groups/consultsw/sl-7.x86_64/modules/gdal/2.2.3/share/gdal")

library(rgdal)                # GDAL bindings
library(sp)                   # spatial data
library(raster)               # raster data
library(rsample)              # function for stratified random sampling
library(dplyr)                # defines the %>% pipe function, among other things


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
       data.dir = '/global/scratch/users/drewhart/seasonality/asynch_drivers_analysis/'
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
cat('\nreading LSP asynchrony file...\n')
phn.asy = read.file(paste0(asynch.var, '_STRICT'), T)


#################
# LOAD PREDICTORS
#################

cat('\nreading climate asynchrony files...\n')
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

cat('\nreading VRM file...\n')
# vector ruggedness metric, ~100km agg med and sd
# NOTE: calculated as fixed pixels, not moving windows
# NOTE: masking isn't necessary because incomplete rows of the stacked
#       variables are dropped later on, but makes it easier to inspect
#       the individual layers manually now because of the crazy blocky artefacts
#       their dataset has surrounding the continents
vrm.med = read.file('vrm_100KMmd', asynch.file=F,
                    align.to=phn.asy, mask.it=T)

cat('\nreading veg entropy file...\n')
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
cat('\nwriting raster stack to disk...\n')
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
cat('\nwriting CSVs of extracted raster values to disk...\n')
write.csv(df, paste0(data.dir, "/asynch_model_all_vars_prepped_",
                     asynch.var, '_',
                     as.character(neigh.rad), "km.csv"), row.names=F)
write.csv(df.strat, paste0(data.dir, "/asynch_model_all_vars_prepped_strat_",
                           asynch.var, '_',
                           as.character(neigh.rad), "km.csv"), row.names=F)
cat('\ndata prep complete.\n\n')

