library(sp)                   # spatial data
library(raster)               # raster data
#library(terra)                # newer raster data
library(sf)                   # newer spatial data
library(spdep)                # spatial autocorrelation
library(rsample)              # function for stratified random sampling
library(RRF)                  # fast RF var selection (w/ conservative results in Bag et al. 2022)
library(ranger)               # faster RFs
library(randomForest)         # regular RFs
#library(h2o)                  # distributed RFs (on cloud)
library(Boruta)               # feature selection w/ boruta algo (i.e., shadow features)
library(fastshap)             # SHAP values 
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



# TODO:
# 1. if runtime is slow then use ranger instead (but then have to figure out how to not make range::importance clobber randomForest::importance that grf depends on!)
# 2. check that no factor/non-numeric variables being fed into RF



# WORKFLOW:

# X. draw subsample (COMPARE RUNNNING WITH AND WITHOUT, LATER)
  # X. check that it's small enough that independent bootstraps can be drawn from it (how??) (OR IS IT JUST MEASURE VARIANCE ACROSS BOOTSTRAPS?)
# 1. RF model (or GB?)
  #1.5 MAYBE LATER CONSIDER GRADIENT-BOOSTING AS WELL/INSTEAD?
# 2. boruta to select vars
# 3. rerun parsimonious model (use Lauren's and HOML approach to choosing reasonable starting param values and to tuning)
# 4. plot and assess var importance
# 5. map SHAP values (and how to get globally integrated ones??)
# 6. run geo-weighted RF and map top vars



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
       analysis.dir = '/global/scratch/users/drewhart/seasonality/'
       # data.dir = '/global/home/groups/fc_landgen/' # directory for Lauren 
       # analysis.dir = '/global/home/users/ldimaggio/ondemand/' # directory for Lauren 
}

# save plots?
save.plots = T
save.rasters = T

# verbose?
verbose = T

# seed
seed.num = 12345
set.seed(seed.num)

# subset the data before building the RF?
subset = T
subset.by.range = F
subset.fracs = c(0.005, 0.05)

# training data fraction
# NOTE: use 60% for training, 40% for testing, which should be plenty,
#       given how large a dataset we have
train.frac = 0.6

# include geo coords in RF?
# (if so, will be incorporated as 3rd-order polynom of projected coords)
include.coords = T


############
# HELPER FNS
############

# Function to buffer points in XY space:

# TAKEN FROM https://davidrroberts.wordpress.com/2015/09/25/spatial-buffering-of-points-in-r-while-retaining-maximum-sample-size/

# Returns the original data table with buffered points removed.
# Runs numerous iterations, as the random point selection can result in more/fewer output points.
# 1) Randomly select a single point
# 2) Remove points within 50km of that point
# 3) Randomly select of the remaining points
# 4) ...
# foo - a data.frame to select from with columns x, y
# buffer - the minimum distance between output points
# reps - the number of repetitions for the points selection
buffer.f <- function(foo, buffer, reps){
     # Make list of suitable vectors
     suitable <- list()
  for(k in 1:reps){
         # Make the output vector
         outvec <- as.numeric(c())
      # Make the vector of dropped (buffered out) points
      dropvec <- c()
          for(i in 1:nrow(foo)){
                   # Stop running when all points exhausted
                   if(length(dropvec)<nrow(foo)){
                              # Set the rows to sample from
                              if(i>1){
                                           rowsleft <- (1:nrow(foo))[-c(dropvec)]
              } else {
                           rowsleft <- 1:nrow(foo)
                      }
              # Randomly select point
              outpoint <- as.numeric(sample(as.character(rowsleft),1))
                      outvec[i] <- outpoint
                      # Remove points within buffer
                      outcoord <- foo[outpoint,c("x","y")]
                              dropvec <- c(dropvec, which(sqrt((foo$x-outcoord$x)^2 + (foo$y-outcoord$y)^2)<buffer))
                              # Remove unnecessary duplicates in the buffered points
                              dropvec <- dropvec[!duplicated(dropvec)]
                                    } 
          } 
          # Populate the suitable points list
          suitable[[k]] <- outvec
            }
    # Go through the iterations and pick a list with the most data
    best <- unlist(suitable[which.max(lapply(suitable,length))])
    foo[best,]
}


# function to plot variable importance for a ranger RF object
plot.ranger.import = function(ranger_rf){
  import = as.data.frame(ranger::importance(ranger_rf))
  import$var = rownames(import)
  colnames(import) = c('import', 'var')
  imp_plot = ggplot(import, aes(x=reorder(var,import), y=import, fill=import))+ 
    geom_bar(stat="identity", position="dodge")+ coord_flip()+
    ylab("Variable Importance")+
    xlab("")+
    ggtitle("Information Value Summary")+
    guides(fill="none")+
    scale_fill_gradient(low="red", high="blue")
  return(imp_plot)
}


# function to make partial dependence plot
# (adapted from: https://zia207.github.io/geospatial-r-github.io/geographically-wighted-random-forest.html)
make.pdp = function(rf, trn, response, varn1, varn2, seed.num=NA){
   # set seed, if requested
   if (!is.na(seed.num)){
      set.seed(seed.num)
   }
   # custom prediction functions
   pred <- function(object, newdata)  {
        results <- as.vector(predict(object, newdata))
     return(results$predictions)
       }
   pdp_pred <- function(object, newdata)  { 
      results <- mean(as.vector(predict(object, newdata))) 
      return(results$predictions) 
   }              
   # get variable imporatnce (by % inc MSE when removed)
   per.var.imp<-vip( 
                    rf, 
                    train = trn,
                    method = "permute", 
                    target = response,
                    metric = "RMSE", 
                    nsim = 5, 
                    sample_frac = 0.5, 
                    pred_wrapper = pred
   )$data
   print(per.var.imp)
   # grab the two variables of interest (numbered by decreased importance by % inc MSE)
   var.names=per.var.imp[1] 
   var_1=var.names[varn1,]$Variable 
   var_2=var.names[varn2,]$Variable 
   # get partial dependence of both vars
   pd_var_1<- pdp::partial(
                          rf, 
                          train = trn,
                          pred.var = var_1, 
                          pred.fun = pdp_pred, 
                          parallel = T,
                          grid.resolution = 10 
   )
   pd_var_2<- pdp::partial(
                          rf, 
                          train = trn,
                          pred.var = var_2, 
                          pred.fun = pdp_pred, 
                          parallel = T,
                          grid.resolution = 10 
   )
   # create their autoplots
   pd_1<-autoplot(pd_var_1, 
                  rug = TRUE, 
                  train=trn) +
        theme(text = element_text(size=15)) 
   pd_2<-autoplot(pd_var_2, 
                  rug = TRUE, 
                  train=trn) +
        theme(text = element_text(size=15)) 
   # bivaraite partial dependence and its autplot
   pd_var_1_v_2<- pdp::partial( 
                                 rf, 
                                 train = trn, 
                                 pred.var = c(var_1, var_2), 
                                 pred.fun = pdp_pred, 
                                 parallel = F,
                                 grid.resolution = 10 
   ) 
  pd_1_v_2<-autoplot(pd_var_1_v_2, contour = TRUE) + 
      theme(text = element_text(size=15)) 
  # arrange plots in grid
  table_grob = grid.arrange(pd_1,pd_2,pd_1_v_2,  
               ncol= 3, 
               heights = c(40,8), 
               top = textGrob(paste0("Partial Dependence Plot: ",
                                     var_1,
                                     " vs. ",
                                     var_2), gp=gpar(fontsize=25))) 
  return(table_grob)
}





######################
# LOAD AND SUBSET DATA
######################

# load countries polygons (for use as a simple basemap)
world = map_data('world')

# get band names
names = c('phn.asy', 'tmp.min.asy', 'tmp.max.asy',
          'tmp.min.mea', 'ppt.asy',
          'ppt.sea', 'def.asy', 'cld.asy', 'vrm.med',
          'riv.dis', 'eco.dis')

# load rasters of prepped variables
vars = brick(paste0(data.dir, "/asynch_model_all_vars.tif"))
names(vars) = names

# load data frame of prepped variables
df_full_unproj = read.csv(paste0(data.dir, "/asynch_model_all_vars_prepped.csv"))

# transform to projected coordinate system, then add cols for polynomial of coords
#df_full_unproj_sf = st_as_sf(df_full_unproj, coords=c('x', 'y'), crs=4326)
#df_full_proj_sf = st_transform(df_full_unproj_sf, 8857)
#TODO: GET TRANSFORMED X AND Y VALS AS COLUMNS AGAIN
#df_full = df_full_proj_sf
#df_full['x2'] = df_full['x']^2
#df_full['y2'] = df_full['y']^2
#df_full['x3'] = df_full['x']^3
#df_full['y3'] = df_full['y']^3

# TODO: FIX BROKEN SOFTWARE THAT CAN'T FIND GDAL_DATA SO CAN'T ASSIGN CRS!!!
#       IN MEANTIME, JUST OMITTING COORDS...
df_full = df_full_unproj






#########################################
# TUNE HYPERPARAMETERS OF GLOBAL RF MODEL
#########################################

# NOTE: ADAPTED FROM https://bradleyboehmke.github.io/HOML/random-forest.html#hyperparameters

n.features = (ncol(df_full)-3)

# for each of two subset fractions:
for (subset.frac in subset.fracs){
  cat('=============\n')
  cat('SUBSET FRAC: ', subset.frac, '\n\n')
  # subset the data
  if (subset.by.range){
      # TODO: DEBUG, IF USING
      # get cell coordinates in global equal-area projected CRS
      xy = df_full[, c('x', 'y')]
      coordinates(xy) = c('x', 'y')
      xy_proj = spTransform(xy, proj4str='+init=epsg:8857')
      # select a subset with all points >= 150km from other points
      # (based on my previous estimation of the range of the spatial autocorrelation in the data)
      df = buffer.f(as.data.frame(xy_proj@coords), buffer=150000, reps=1)
  }
  else {
      # take a stratified random sample, stratified by the response var
      split_subset  <- initial_split(df_full, prop = subset.frac, strata = "phn.asy")
      df = training(split_subset)
  }
  # get training and test data, using a stratified random sample
  split_strat  <- initial_split(df, prop = train.frac, 
                                strata = "phn.asy")
  trn  <- training(split_strat)
  tst  <- testing(split_strat)
    #build a default ranger model
  default.fit <- ranger(
    formula = phn.asy ~ ., 
    data = trn[, 3:ncol(trn)],
    mtry = floor(n.features / 3),
    respect.unordered.factors = "order",
    seed = 123
  )
  # and get its OOB and test-set RMSEs
  default.rmse <- sqrt(default.fit$prediction.error)
  default.preds = predict(default.fit, tst[,3:ncol(tst)])$predictions
  default.errs = default.preds - tst[,'phn.asy'] 
  default.tst.rmse = sqrt(mean(default.errs^2))
  hyper_grid <- expand.grid(
    mtry = floor(n.features * c(.05, .25, .4)),
    ntree = c(150, 200, 250, 300),
    min.node.size = c(1, 3, 5, 10), 
    replace = c(TRUE, FALSE),        
    sample.fraction = c(.3, .55, .8),                       
    rmse = NA,
    tst.rmse = NA
  )
  for(i in seq_len(nrow(hyper_grid))) {
      # fit model for ith hyperparameter combination
    fit <- ranger(
      formula         = phn.asy ~ ., 
      data            = trn[, 3:ncol(trn)], 
      num.trees       = hyper_grid$ntree[i],
      mtry            = hyper_grid$mtry[i],
      min.node.size   = hyper_grid$min.node.size[i],
      replace         = hyper_grid$replace[i],
      sample.fraction = hyper_grid$sample.fraction[i],
      verbose         = FALSE,
      seed            = 123,
      respect.unordered.factors = 'order',
    )
    hyper_grid$rmse[i] <- sqrt(fit$prediction.error)
    preds = predict(fit, tst[,3:ncol(tst)])$predictions
    errs = preds - tst[,'phn.asy'] 
    tst.rmse = sqrt(mean(errs^2))
    hyper_grid$tst.rmse[i] = tst.rmse
    if (i%%60 == 0){
      cat(i/nrow(hyper_grid)*100, '% complete...\n')
    }
  }
  hyper_grid_complete = hyper_grid %>%
    arrange(rmse) %>%
    mutate(perc_gain = (default.rmse - rmse) / default.rmse * 100,
           #perc_tst_gain = (default.tst.rmse - tst.rmse) / default.tst.rmse * 100,
           )
  cat(head(hyper_grid_complete, 50))
  hyper_grids[[subset.frac]] = hyper_grid_complete
}


# choose global RF params based on output above
ntree = 100
replace = T
rf.samp.frac = 0.75
mtry = round(n.features/3)

# and subset data based on output above
subset.frac = 0.05
if (subset.by.range){
  # TODO: DEBUG, IF USING
  # get cell coordinates in global equal-area projected CRS
  xy = df_full[, c('x', 'y')]
  coordinates(xy) = c('x', 'y')
  xy_proj = spTransform(xy, proj4str='+init=epsg:8857')
  # select a subset with all points >= 150km from other points
  # (based on my previous estimation of the range of the spatial autocorrelation in the data)
  df = buffer.f(as.data.frame(xy_proj@coords), buffer=150000, reps=1)
} else {
  # take a stratified random sample, stratified by the response var
  split_subset  <- initial_split(df_full, prop = subset.frac, strata = "phn.asy")
  df = training(split_subset)
}
# get training and test data, using a stratified random sample
split_strat  <- initial_split(df, prop = train.frac, 
                              strata = "phn.asy")
trn  <- training(split_strat)
tst  <- testing(split_strat)
print(paste0('TRAINING DATA: ', nrow(trn), ' ROWS'))
print(paste0('TEST DATA: ', nrow(tst), ' ROWS'))

# scatter the training and test datasets
scat_trn = ggplot(tst) +
  geom_polygon(data=world, aes(x=long, y=lat, group=group), color="black", fill="white" ) +
  geom_point(data=trn, aes(x=x, y=y), col='blue', alpha=0.05, size=0.1)
scat_tst = ggplot(tst) +
  geom_polygon(data=world, aes(x=long, y=lat, group=group), color="black", fill="white" ) +
  geom_point(data=tst, aes(x=x, y=y), col='red', alpha=0.05, size=0.1)
grid.arrange(scat_trn, scat_tst)




##########################
# BORUTA FEATURE SELECTION
##########################

bor_res = Boruta(phn.asy ~ .,
                 data=trn[,3:ncol(trn)],
                 doTrace=verbose * 3,
                 num.trees=ntree,
                 mtry=mtry,
                 replace=replace,
                 sample.fraction=rf.samp.frac
)
print(attStats(bor_res))

if (save.plots){
  jpeg(paste0(data.dir, '/boruta_boxplot.jpg'),
       width=900,
       height=400,
       units='px',
       quality=90,
  )
  par(mfrow=c(1,2))
  plot(bor_res)
  plotImpHistory(bor_res)
  dev.off()
}
par(mfrow=c(1,1))

# remove any vars rejected by boruta
# from trn, tst, df, and df_full

# ...


###########################################
# BUILD TUNED, PARSIMONIOUS GLOBAL RF MODEL
###########################################

rf_global = ranger(phn.asy ~ .,
                   data=trn[,3:ncol(trn)],
                   num.trees=ntree,
                   mtry=mtry,
                   importance='permutation',
                   verbose=verbose,
                   replace=replace,
                   sample.fraction=rf.samp.frac,
)

# take a look at the results
print(rf_global) 

# var importance plots, with permutation based metric...
p_imp_permut = plot.ranger.import(rf_permut)
pfun <- function(object, newdata) {
  predict(object, data = newdata)$predictions
}
# ... and with SHAP values
shap <- fastshap::explain(rf_permut, X = df[, 4:ncol(df)], pred_wrapper = pfun, nsim = 10)
shap_imp <- data.frame(
  Variable = names(shap),
  Importance = apply(shap, MARGIN = 2, FUN = function(x) sum(abs(x)))
)
p_imp_shap = ggplot(shap_imp, aes(reorder(Variable, Importance), Importance)) +
  geom_col() +
  coord_flip() +
  xlab("") +
  ylab("mean(|Shapley value|)")
varimp_grob = grid.arrange(shap_imp, permut_imp, ncol=2)

if (save.plots){
   ggsave(varimp_grop, file='var_import_plots_permut_and_SHAP.png',
          width=45, height=35, units='cm', dpi=600)
}

# partial dependence plots
# ppt.asy vs ppt.sea
#pdp12 = make.pdp(rf_permut, trn, 'phn.asy', 1, 2, seed.num=seed.num)
#plot(pdp12)

#
#pdp13 = make.pdp(rf_permut, trn, 'phn.asy', 1, 3, seed.num=seed.num)
#plot(pdp13)

#pdp23 = make.pdp(rf?permut, trn, 'phn.asy', 2, 3, seed.num=seed.num)
#plot(pdp23)

#if (save.plots){
#   ggsave(pdp12, file='pdp12.png',
#          width=22, height=12, units='cm', dpi=500)
#   ggsave(pdp13, file='pdp13.png',
#          width=22, height=12, units='cm', dpi=500)
#   ggsave(pdp23, file='pdp23.png',
#          width=22, height=12, units='cm', dpi=500)
#}


# assess model externally using withheld test data
preds = predict(rf_permut, tst[,3:ncol(tst)])$predictions
tst$err = preds - tst[,'phn.asy'] 
preds_plot = ggplot(tst) +
    geom_polygon(data=world, aes(x=long, y=lat, group=group), color="black", fill="white" ) +
    geom_point(aes(x=x, y=y, col=err, alpha=abs(err)/max(abs(err))), size=1) +
    scale_color_gradient2(low='#cc003d', mid='#dbdbdb', high='#009e64') +
    #coord_map() + 
    theme_bw()
preds_plot

if (save.plots){
   ggsave(preds_plot, file='preds_plot.png',
          width=30, height=22, units='cm', dpi=500)
  }


# make predicitions for full dataset (to map as raster)
full_preds = predict(rf_permut, df_full[,3:ncol(df_full)])$predictions
df.res = df_full %>% mutate(preds = full_preds, err = full_preds - df_full[,'phn.asy'])
dfrast <- rasterFromXYZ(df.res[, c('x', 'y', 'phn.asy', 'preds', 'err')])
re_df = as.data.frame(dfrast, xy=T)
colnames(re_df) = c('x', 'y', 'phn.asy', 'preds', 'err')
asy.cbar.limits = quantile(re_df$phn.asy, c(0.01, 0.99), na.rm=T, type=8)
err.cbar.limits = rep(max(abs(quantile(re_df$err, c(0.01, 0.99), na.rm=T, type=8))),2) * c(-1,1)
obs_map = ggplot() + 
        geom_polygon(data=world, aes(x=long, y=lat, group=group), color="black", fill="white" ) +
        geom_raster(data=re_df, aes(x=x, y=y, fill=phn.asy)) +
        #scale_fill_gradient2(low='#940500', mid='#f5f5f5', high='#00138f') +
        scale_fill_cmocean(name='dense', direction=-1, limits=asy.cbar.limits) +
        #coord_quickmap() +
        theme_bw() +
        theme(legend.key.size = unit(2, 'cm'),
              legend.key.width = unit(1, 'cm'),
              legend.text = element_text(size=13),
              legend.title = element_text(size=16))
preds_map = ggplot() +
        geom_polygon(data=world, aes(x=long, y=lat, group=group), color="black", fill="white" ) +
        geom_raster(data=re_df, aes(x=x, y=y, fill=preds)) +
        #scale_fill_gradient2(low='#940500', mid='#f5f5f5', high='#00138f') +
        scale_fill_cmocean(name='dense', direction=-1, limits=asy.cbar.limits) +
        #coord_quickmap() +
        theme_bw() +
        theme(legend.key.size = unit(2, 'cm'),
              legend.key.width = unit(1, 'cm'),
              legend.text = element_text(size=13),
              legend.title = element_text(size=16))
err_map = ggplot() +
        geom_polygon(data=world, aes(x=long, y=lat, group=group), color="black", fill="white" ) +
        geom_raster(data=re_df, aes(x=x, y=y, fill=err)) +
        #scale_fill_gradient2(low='#940500', mid='#f5f5f5', high='#00138f') +
        scale_fill_cmocean(name='curl', direction=-1, limits=err.cbar.limits) +
        #coord_quickmap() +
        theme_bw() +
        theme(legend.key.size = unit(2, 'cm'),
              legend.key.width = unit(1, 'cm'),
              legend.text = element_text(size=13),
              legend.title = element_text(size=16))
obs_vs_pred = ggplot() +
        geom_point(data=re_df, aes(x=phn.asy, y=preds), alpha=0.05) +
        geom_abline(intercept=0, slope=1)
global_main_plots = cowplot::plot_grid(obs_map, preds_map, err_map, obs_vs_pred, nrow=2, ncol=2)
global_main_plots

if (save.plots){
   ggsave(global_grrf_main_plots, file='global_grrf_main_plot.png',
          width=75, height=55, units='cm', dpi=500)
  }
if (save.rasters){
  writeRaster(dfrast,
              paste0(data.dir, '\global_rf_map_results.tif'),
              'GTiff',
              overwrite = T
  )
}


# plot SHAP values
p1 <- autoplot(shap)
p2 <- autoplot(shap, type = "dependence", feature = "ppt.asy", X = df[,3:ncol(df)], alpha = 0.5,
               color_by = "ppt.sea", smooth = TRUE, smooth_color = "black") +
  scale_color_viridis_c()
gridExtra::grid.arrange(p1, p2, nrow = 1)
p2 <- autoplot(shap, type = "dependence", feature = "tmp.min.asy", X = df[,3:ncol(df)], alpha = 0.5,
               color_by = "tmp.min.mea", smooth = TRUE, smooth_color = "black") +
  scale_color_viridis_c()
gridExtra::grid.arrange(p1, p2, nrow = 1)

# map SHAP values
shap_full = fastshap::explain(rf_permut, X = df_full[, 4:ncol(df_full)], pred_wrapper = pfun, nsim = 10)
df_shap_full = cbind(df_full[,c('x', 'y')], shap_full)
df_shap_full_rast <- rasterFromXYZ(df_shap_full)
names(df_shap_full_rast) = colnames(df_shap_full)[3:ncol(df_shap_full)]

if (save.rasters){
  for (lyr in names(df_shap_rast_full)){
    writeRaster(df_shap_full_rast,
                paste0(data.dir, '\rf_SHAP_vals_', lyr, '.tif'),
                'GTiff',
                overwrite = T
    )
  }
}

######################
# RUN LOCAL GWRF MODEL
######################

# code adapted from https://zia207.github.io/geospatial-r-github.io/geographically-wighted-random-forest.html

# local RF params
bw.local = 150
ntree.local = ntree
replace.local = replace
rf.samp.frac.local = rf.samp.frac
mtry.local = mtry

coords = trn[,c('x', 'y')]

grf.model <- SpatialML::grf(formula=phn.asy ~ tmp.min.asy + tmp.max.asy + tmp.min.mea +
                               ppt.asy + ppt.sea + def.asy + cld.asy + 
                               vrm.med + riv.dis + eco.dis,
                 dframe=trn[, 3:ncol(trn)], 
                 bw=bw.local,
                 kernel="adaptive",
                 ntree=ntree.local,
                 mtry=mtry.local,
                 forests = TRUE,
                 coords=coords,
                 importance = TRUE
                 )

# global model summary
glob.imp = as.data.frame(randomForest::importance(grf.model$Global.Model, type=1)) # %IncMSE
colnames(glob.imp) = c('pct.inc.mse')
glob.imp = glob.imp %>% arrange(desc(pct.inc.mse))

print(paste0('GLOBAL MODEL MSE: ', mean(grf.model$Global.Model$mse)))

print(paste0('GLOBAL MODEL Rsq: ', mean(grf.model$Global.Model$rsq)))


# local model summary and variable importance
print(grf.model$LocalModelSummary)

#imp.var.to.use = "Local.Pc.IncMSE" 
imp.var.to.use = "Local.IncNodePurity"
var.imp.loc = grf.model[[imp.var.to.use]]

# plot local var importance, with plots in order of decreasing global var imp
plots = lapply(seq(nrow(glob.imp)), function(n){
  var = rownames(glob.imp)[n]
  glob.imp.val = glob.imp[var,]
  cbar.limits = rep(max(abs(quantile(var.imp.loc[,var], c(0.01, 0.99), na.rm=T, type=8))),2) * c(-1,1)
  p = ggplot() +
     geom_polygon(data=world, aes(x=long, y=lat, group=group),
                                  color="black", fill="white" ) +
     geom_point(aes(x=trn$x, y=trn$y, col=var.imp.loc[,var]), alpha=0.5, size=0.75) +
     scale_color_cmocean(name='curl', direction=-1, limits=cbar.limits) +
     ggtitle(paste0(var, ": Global ", imp.var.to.use, ": ", glob.imp.val)) +
     labs(col=var)
  return(p)
})
local_rf_main_plots = grid.arrange(plots[[1]], plots[[2]], plots[[3]],
             plots[[4]], plots[[5]], plots[[6]], 
             plots[[7]], plots[[8]], plots[[9]],
             plots[[10]],
             ncol=4)

if (save.plots){
   ggsave(local_rf_main_plots, file='local_rf_main_plot.png',
          width=75, height=55, units='cm', dpi=1000)
  }



# map goodness of fit 
gof.loc = grf.model$LGofFit

# plot local var importance, with plots in order of decreasing global var imp
plots = lapply(seq(ncol(gof.loc)), function(n){
  var = colnames(gof.loc)[n]
  cbar.limits = rep(max(abs(quantile(gof.loc[,var], c(0.01, 0.99), na.rm=T, type=8))),2) * c(-1,1)
  p = ggplot() +
     geom_polygon(data=world, aes(x=long, y=lat, group=group),
                                  color="black", fill="white" ) +
     geom_point(aes(x=grf.model$Locations[,1], y=grf.model$Locations[,2],
                    col=gof.loc[,var]), alpha=0.5, size=0.75) +
     scale_color_cmocean(name='curl', direction=-1, limits=cbar.limits) +
     ggtitle(paste0("GoF metric: ", var)) +
     labs(col=var)
  return(p)
})
local_rf_gof_plots = grid.arrange(plots[[1]], plots[[2]], plots[[3]],
             plots[[4]], plots[[5]], plots[[6]], plots[[7]],
             ncol=4)

if (save.plots){
   ggsave(local_rf_main_plots, file='local_rf_gof_plot.png',
          width=75, height=55, units='cm', dpi=1000)
  }


####################################
# MODEL/VAR SELECTION AND EVALUATION
####################################

# what to do here exactly?

# and should I use K-fold CV or bootstrapping? does this choice matter much?


# gather response and predictors for training and test data 
  # response - training
  y_trn = trn$phn.asy

  # features - training
  x_trn = trn[,4:ncol(trn)]

  # response - test
  ytest = tst$phn.asy

  # features - test
  xtest = tst[,4:(ncol(tst)-1)]


#
# FINE TUNE THE GLOBAL MODEL
# 

# FIRST FIND THE OPTIMAL ntree VALUE
  # SET mtry TO THE DEFAULT VALUE floor(p/3)

  # get the number of predictors
  nvars = ncol(x_trn)

  # set the default mtry value to floor(p/3)
  mtry.default = floor(nvars/3)
  
  # set max ntree value to try 
  max.ntree = 75
  
  # look at plot of error vs different ntree values to determine which ntree values to test
  test.ntree.grrf <- RRF(phn.asy ~ .,
                         data=trn[,3:ncol(trn)],
                         xtest = xtest,
                         ytest = ytest,
                         ntree = max.ntree,
                         mtry=mtry.default,
                         coefReg=coefReg,
                         flagReg=1)
    # Generate an error plot for the different ntree values
    ntree.err.plot.grrf = plot(test.ntree.grrf)
      # based on the plot, ntree values ranging from 50 to 75 appear to have the lowest error 

  # create a vector of different ntree values to test
  grrf.ntree.vals <- c(seq(ntree, max.ntree, 5))

  # Test Different ntree values
    # create empty df to store MSE and R-squared values 
    grrf.ntree.stats <- as.data.frame(matrix(nrow = length(grrf.ntree.vals), ncol = 2,
                                  dimnames = list(grrf.ntree.vals, # ntree values
                                                  c('MSE', "RSQ") 
                                                  )
                                  )
                                )

    # test with different ntree values
    for (i in 1:length(grrf.ntree.vals)) {
      grrf.test.ntree.model <- RRF(phn.asy ~ .,
                        data=trn[,3:ncol(trn)],
                        xtest = xtest,
                        ytest = ytest,
                        ntree = grrf.ntree.vals[i],
                        mtry=mtry.default,
                        coefReg=coefReg,
                        flagReg=1)
        # gather mse & r-squared values from models with different ntree values into a df
        grrf.ntree.stats[i,1] <- grrf.test.ntree.model$mse[length(grrf.test.ntree.model$mse)]
        grrf.ntree.stats[i,2] <- grrf.test.ntree.model$rsq[length(grrf.test.ntree.model$rsq)]
      }
    # select the model whose ntree value gives the largest r-squared and lowest MSE
      # select the max r-squared since MSE typically decreases as ntree increases
      opt.ntree.grrf = as.numeric(rownames(grrf.ntree.stats[which.max(grrf.ntree.stats[ ,2]), ]))
        # it looks like the MSE & RSQ stabilizes at ntree = 75 
        # the accuracy improves by 0.00124 from ntree = 70 to ntree = 75



# FINDING THE OPTIMAL mtry VALUE 

  # number of times to run a test that searches for the optimal mtry value
  ntests.mtry = 50

  # run test to find optimal mtry value
    # store the mtry values in an empty df
    grrf.mtry.vals <- as.data.frame(matrix(nrow = ntests.mtry, ncol=1,
                                           dimnames = list(
                                                      c(paste("test", as.character(1:ntests.mtry),sep="_")), # test number
                                                      'mtry.value'
                                                      )
                                           )
                                    )

    # tune the RRF model for the optimal mtry value
    for (i in 1:ntests.mtry) {
      # find optimal mtry value
      grrf.tune.mtry <- tuneRRF(x = x_trn,
                                y = y_trn,
                                mtryStart = mtry.default, #floor(nvars/2)
                                ntreeTry = opt.ntree.grrf,
                                stepFactor = 1,
                                improve = 0.0001, # amount of improvement in OOB error required to continue searching for better mtry values
                                trace = F,
                                plot = F)
      # pull the mtry values from the 50 simulations with the smallest OOB error and store in df 
      grrf.mtry.vals[i, ] = grrf.tune.mtry[grrf.tune.mtry[, 2] == min(grrf.tune.mtry[, 2]), 1]
    }

    # find the frequency of mtry values to determine which value was selected the most
    grrf.mtry.val.freq <- as.data.frame(table(grrf.mtry.vals), stringsAsFactors = F)

    # find the mtry value with the largest frequency (AKA the optimal mtry value)
    opt.mtry.grrf = as.numeric(grrf.mtry.val.freq[which.max(grrf.mtry.val.freq[ ,2]), 1])


# Run GRRF with optimal mtry and ntree values
    opt.grrf <- RRF(phn.asy ~ .,
                    data = trn[,3:ncol(trn)],
                    xtest = xtest,
                    ytest = ytest,
                    ntree = opt.ntree.grrf,
                    mtry = opt.mtry.grrf,
                    coefReg=coefReg,
                    flagReg=1, 
                    importance = T)
    # Mean of squared residuals: 0.6893117
    # % Var explained: 30.76  






# VARIABLE SELECTION FOR GLOBAL MODEL
    
    # look at variable importance plots for the parameter optimized GRRF model defined above 
    RRF::varImpPlot(opt.grrf)
    importance(opt.grrf, type = 1)

  ## K-FOLD CV

    # adapted from https://towardsdatascience.com/effective-feature-selection-recursive-feature-elimination-using-r-148ff998e4f7#:~:text=Recursive%20Feature%20Elimination%C2%B2%2C%20or%20shortly,the%20optimal%20combination%20of%20features. 

    # Define the control using a random forest selection function
    control <- rfeControl(functions = rfFuncs, # random forest
                          method = "repeatedcv", # repeated cv
                          repeats = 10, # number of repeats 
                          # 10 is sufficient according to https://www.listendata.com/2014/11/random-forest-with-r.html
                          number = 10) # number of folds
    

    # Run RFE (recursive feature selection) 
      # ISSUE: mtry defaults to floor(p/3) and I cannot find a way to change this.  I'm not sure if this impacts the eature selection results
    rfe.results <- rfe(x = x_trn, 
                       y = y_trn, 
                       sizes = c(8:10), # number of features to try for each fold
                       rfeControl = control)

    # list the chosen features
    predictors(rfe.results)
    # plot the results
    plot(rfe.re