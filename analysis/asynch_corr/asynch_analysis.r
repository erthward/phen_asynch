library(raster)               # raster data
library(terra)                # newer raster data
library(sp)                   # spatial data
library(sf)                   # newer spatial data
library(spdep)                # spatial autocorrelation
library(randomForest)         # global RFs
library(h2o)                  # distributed RFs (on cloud)
library(SpatialML)            # GWRFs
library(GWmodel)              # GW models
library(vip)                  # var importance plots
library(pdp)                  # partial depend. plots (& ICE curves)
library(DALEX)                # feature importance
library(ggplot2)              # plotting
library(ggthemes)             # plot themes
library(cowplot)              # easy plot gridding
library(tmap)                 # mapping
library(RColorBrewer)         # colors
library(cmocean)              # cmocean palettes
library(wesanderson)          # wes palettes

# TODO:
# 1. when/why is the extent being truncated so that Australia and such get lost??


#######################
# SET BEHAVIORAL PARAMS
#######################

# input and output dirs:
  # if on laptop
if (strsplit(getwd(), '/')[[1]][2] == 'home'){
       on.laptop=T
       data.dir = '/home/deth/Desktop/CAL/research/projects/seasonality/results/maps'
    analysis.dir = '/home/deth/Desktop/CAL/research/projects/seasonality/results/maps'
  # if on Savio
} else {
       on.laptop=F
       data.dir = '/global/scratch/users/drewhart/seasonality/'
       analysis.dir = '/global/scratch/users/drewhart/seasonality/'
}

# save plots?
save.plots=F

# RF params
ntree=1000
mtry=NULL


############
# HELPER FNS
############

align.rast.x.to.y = function(rast.x, rast.y, mask.it=F){
   # resample
   rast.x.resamp = raster::resample(rast.x, rast.y, "bilinear")
   # crop
   ext = extent(rast.y)
   rast.x.out = crop(rast.x.resamp, ext)
   # mask?
   if (mask.it){
      rast.x.out = mask(rast.x.out, rast.y)
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
   if (class(align.to)[1] == 'RasterLayer'){
      data = align.rast.x.to.y(data, align.to, mask.it=mask.it) 
   }
   return(data)
}



####################
# LOAD RESPONSE VARS
####################

# NIRv-based asynch
# NOTE: ADD ME
#phn.asy = read.file('NIRv', T)

# SIF-based asynch
phn.asy = read.file('SIF', asynch.file=T)


#################
# LOAD PREDICTORS
#################

# asynchrony in mean temperature
tmp.mea.asy = read.file('tmmean', asynch.file=T,
                    align.to=phn.asy, mask.it=F)

# asynchrony in min temperature
tmp.min.asy = read.file('tmmn', asynch.file=T,
                    align.to=phn.asy, mask.it=F)

# asynchrony in max temperature
tmp.max.asy = read.file('tmmx', asynch.file=T,
                    align.to=phn.asy, mask.it=F)

# mean daily min temp, coldest month
# NOTE: CALC AS NEIGH MEAN AND/OR SD?
tmp.min.mea = read.file('CHELSA_bio6', asynch.file=F,
                    align.to=phn.asy, mask.it=F)

# mean daily max temp, warmest month
# NOTE: CALC AS NEIGH MEAN AND/OR SD?
tmp.max.mea = read.file('CHELSA_bio5', asynch.file=F,
                    align.to=phn.asy, mask.it=F)

# temperature seasonality
# NOTE: CALC AS NEIGH MEAN AND/OR SD?
tmp.sea = read.file('CHELSA_bio4', asynch.file=F,
                    align.to=phn.asy, mask.it=F)

# number of growing degree days (above 0deg C)
# NOTE: CALC AS NEIGH MEAN AND/OR SD?
ngd = read.file('CHELSA_ngd0', asynch.file=F,
                    align.to=phn.asy, mask.it=F)

# asynchrony in precipitation
ppt.asy = read.file('pr', asynch.file=T,
                    align.to=phn.asy, mask.it=F)

# precipitation seasonality
# NOTE: CALC AS NEIGH MEAN AND/OR SD?
ppt.sea = read.file('CHELSA_bio15', asynch.file=F,
                    align.to=phn.asy, mask.it=F)

# asynchrony in climatic water deficit
def.asy = read.file('def', asynch.file=T,
                    align.to=phn.asy, mask.it=F)

# asynchrony in cloud cover
cld.asy = read.file('cloud', asynch.file=T,
                    align.to=phn.asy, mask.it=F)

# asynchrony in wind speed
wds.asy = read.file('vs', asynch.file=T,
                    align.to=phn.asy, mask.it=F)

# vector ruggedness metric
# NOTE: CALC AS NEIGH MEAN AND/OR SD?
vrm = read.file('vrm', asynch.file=F,
                    align.to=phn.asy, mask.it=F)

# neighborhood Shannon entropy of land cover
# NOTE: CALCULATE THIS!
#lc.ent = read.file('MODIS_IGBP', asynch.file=F,
                    #align.to=phn.asy, mask.it=F)

# distance from rivers
# NOTE: INCLUDE? IF SO, NEIGH MEAN AND/OR SD?
riv.dis = read.file('dist_to_rivers', asynch.file=F,
                    align.to=phn.asy, mask.it=F)



###########
# PREP DATA
###########

# coproject and align extent

# gather into a stack
vars = stack(phn.asy, tmp.mea.asy, tmp.min.asy, tmp.max.asy, tmp.min.mea, tmp.max.mea,
             tmp.sea, ngd, ppt.asy, ppt.sea, def.asy, cld.asy, wds.asy, vrm, #lc.ent
             riv.dis)

# aggregate to coarser res, if working on laptop
if (on.laptop){
   vars <- aggregate(vars, fact=12)
}
names(vars) = c('phn.asy', 'tmp.mea.asy', 'tmp.min.asy', 'tmp.max.asy',
                'tmp.min.mea', 'tmp.max.mea', 'tmp.sea', 'ngd', 'ppt.asy',
                'ppt.sea', 'def.asy', 'cld.asy', 'wds.asy', 'vrm', #'lc.ent',
                'riv.dis')

# write stack to file (for quicker reload later),
# converting first with terra to preserve band names
#vars = rast(vars)
#terra::writeRaster(vars, "asynch_model_all_vars.tif")
#vars = brick("asynch_model_all_vars.tif")

# coerce to a data.frame
df = as.data.frame(vars, xy=T)

# drop NAs
# TODO: FIGURE OUT WHY WINDSPEED LAYER IS EMPTY, THEN DELETE NEXT LINE
df = df[,!names(df)%in%c('wds.asy')]
df = df[complete.cases(df),]

# scale variables
df[,2:ncol(df)] = scale(df[,2:ncol(df)])

# write data.frame to file (for quicker reload lateR)
#write.csv(df, "asynch_model_all_vars_prepped.csv", row.names=F)
#df = read.csv("asynch_model_all_vars_prepped.csv")

#####################
# RUN GLOBAL RF MODEL
#####################


# get training and test data
train_frac = 0.25
indices = sample(1:nrow(df), size = round(train_frac * nrow(df)))
trn = df[indices,]
tst = df[-indices,]

# set formula for RF
formrf = parse(text=paste0("phn.asy ~ ", paste(colnames(df)[4:ncol(df)], collapse=" + ")))

# build the model
# NOTE: leaving lat and lon out of model
rf = randomForest(phn.asy ~ ., data=trn[,3:ncol(trn)],
#rf = randomForest(phn.asy ~ ., data=trn,
                  ntree=ntree, mtry=mtry, importance=TRUE)

# take a look at the result
rf
varImpPlot(rf)

# make predictions
preds = predict(rf, tst[,3:ncol(tst)])
#preds = predict(rf, tst)
tst$err = tst[,'phn.asy'] - preds

ggplot(tst) +
    geom_point(aes(x=x, y=y, col=err), size=1) +
    scale_color_gradient2(low='#cc003d', mid='#dbdbdb', high='#009e64') +
    #coord_map() + 
    theme_bw()

# make predicitons for full dataset (to map as raster)
full_preds = predict(rf, df[,3:ncol(df)])
df$preds = full_preds
df$err = df[,'phn.asy'] - full_preds
dfrast <- rasterFromXYZ(df[, c('x', 'y', 'phn.asy', 'preds', 'err')])
re_df = as.data.frame(dfrast, xy=T)
colnames(re_df) = c('x', 'y', 'phn.asy', 'preds', 'err')
cbar.limits = c(min(re_df[,c('phn.asy', 'preds')], na.rm=T),
           max(re_df[,c('phn.asy', 'preds')], na.rm=T))
obs_map = ggplot() + 
        geom_raster(data=re_df, aes(x=x, y=y, fill=phn.asy)) +
        #scale_fill_gradient2(low='#940500', mid='#f5f5f5', high='#00138f') +
        scale_fill_cmocean(name='dense', direction=-1, limits=cbar.limits) +
        #coord_quickmap() +
        theme_bw() +
        theme(legend.key.size = unit(2, 'cm'),
              legend.key.width = unit(1, 'cm'),
              legend.text = element_text(size=13),
              legend.title = element_text(size=16))
preds_map = ggplot() +
        geom_raster(data=re_df, aes(x=x, y=y, fill=preds)) +
        #scale_fill_gradient2(low='#940500', mid='#f5f5f5', high='#00138f') +
        scale_fill_cmocean(name='dense', direction=-1, limits=cbar.limits) +
        #coord_quickmap() +
        theme_bw() +
        theme(legend.key.size = unit(2, 'cm'),
              legend.key.width = unit(1, 'cm'),
              legend.text = element_text(size=13),
              legend.title = element_text(size=16))
err_map = ggplot() +
        geom_raster(data=re_df, aes(x=x, y=y, fill=err)) +
        #scale_fill_gradient2(low='#940500', mid='#f5f5f5', high='#00138f') +
        scale_fill_cmocean(name='curl', direction=-1, end = 0.88) + 
        #coord_quickmap() +
        theme_bw() +
        theme(legend.key.size = unit(2, 'cm'),
              legend.key.width = unit(1, 'cm'),
              legend.text = element_text(size=13),
              legend.title = element_text(size=16))
obs_vs_pred = ggplot() +
        geom_point(data=re_df, aes(x=phn.asy, y=preds))
plots = cowplot::plot_grid(obs_map, preds_map, err_map, obs_vs_pred, nrow=2, ncol=2)

plots




if (save.plots){
   ggsave(plots, file='global_rf_results.png', width=35, height=45, units='cm', dpi=1000)
}



##########################################################
##########################################################
#############DELETE ME####################################
##########################################################
##########################################################

# NOTE: some of this code pulled and tweaked from:
# https://zia207.github.io/geospatial-r-github.io/geographically-wighted-random-forest.html

# grid search for hyperparameters
h2o.init(nthreads=-1, max_mem_size="48g", enable_assertions=F)

# define h2o dfs
test.mf = df[seq(1, nrow(df), 3), 3:ncol(df)] 
valid.mf = df[seq(2, nrow(df), 3), 3:ncol(df)]
train.mf = df[seq(3, nrow(df), 3), 3:ncol(df)]
test.hex<-  as.h2o(test.mf)
valid.hex<-  as.h2o(valid.mf)
train.hex<-  as.h2o(train.mf)

# define response and predictors
response <- "phn.asy"
predictors <- setdiff(names(train.hex), response)

# hyperparameters
drf_hyper_params <-list(
          ntrees  = seq(10, 5000, by = 10),
          max_depth=c(10,20,30,40,50),
          sample_rate=c(0.7, 0.8, 0.9, 1.0)
          )

#  search criteria
drf_search_criteria <- list(
          strategy = "RandomDiscrete", 
          max_models = 200,
          max_runtime_secs = 900,
          stopping_tolerance = 0.001,
          stopping_rounds = 2,
          seed = 1345767
          )

# Grid Search
drf_grid <- h2o.grid(
          algorithm="randomForest",
          grid_id = "drf_grid_IDy",
          x= predictors,
          y = response,
          training_frame = train.hex,
          validation_frame = valid.hex,
          stopping_metric = "RMSE",
          nfolds=10,
          keep_cross_validation_predictions = TRUE,
          hyper_params = drf_hyper_params,
          search_criteria = drf_search_criteria,
          seed = 42
          )

# get RF grid params
drf_get_grid <- h2o.getGrid("drf_grid_IDx",sort_by="RMSE",decreasing=FALSE)
drf_get_grid@summary_table[1,]

# get the best DRF model
best_drf <- h2o.getModel(drf_get_grid@model_ids[[1]]) 
#capture.output(print(summary(best_drf)),file =  "DRF_summary_N_RY.txt")
best_drf

# get the DRF CV result
cv.drf<-best_drf@model$cross_validation_metrics_summary%>%.[,c(1,2)]
cv.drf





#################################################################
#################################################################
#############DELETE ME###########################################
#################################################################
#################################################################


######################
# RUN LOCAL GWRF MODEL
######################

coords = trn[,1:2]
grf.model <- grf(formula=phn.asy ~ tmp.mea.asy + tmp.min.asy + tmp.min.mea + tmp.max.mea +
                               tmp.sea + ngd + ppt.asy + ppt.sea + def.asy + cld.asy + 
                               vrm + riv.dis,
                 dframe=trn[, 3:ncol(trn)], 
                 bw=100,
                 kernel="fixed",
                 ntree=ntree, 
                 mtry=mtry,
                 forests = FALSE,
                 coords=coords)




###############################
# ASSESS PREDICTIVE PERFORMANCE
###############################

# use K-fold CV or bootstrapping? matter much?


