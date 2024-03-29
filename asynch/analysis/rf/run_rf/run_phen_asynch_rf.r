# CHRIS PACIOREK'S WORKAROUND TO FIX THE rgdal/sp/sf ISSUE I RAN INTO:
Sys.setenv(GDAL_DATA = "/global/home/groups/consultsw/sl-7.x86_64/modules/gdal/2.2.3/share/gdal")

library(sp)                   # spatial data
library(raster)               # raster data
#library(terra)                # newer raster data
library(sf)                   # newer spatial data
library(spdep)                # spatial autocorrelation
library(rsample)              # function for stratified random sampling
library(RRF)                  # fast RF var selection (w/ conservative results in Bag et al. 2022)
library(ranger)               # faster RFs
library(randomForest)         # regular RFs
library(spatialRF)            # toolkit for running RFs on spatial data
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


##########################################################################
# SETUP


#######################
# SET BEHAVIORAL PARAMS
#######################

args = commandArgs(trailingOnly=T)

# phen-asynch var to use
asynch.var = args[1]
cat('\nVAR: ', asynch.var, '\n')

# asynchrony neighborhood radius to use (in km)
neigh.rad = args[2]
cat('\nNEIGH RAD: ', neigh.rad, '\n')

# include coordinates in RF?
coords.as.covars = args[3]
cat('\nCOORDS AS COVARS? ', coords.as.covars, '\n')


# input and output dirs:
# if on laptop
if (strsplit(getwd(), '/')[[1]][2] == 'home'){
  on.laptop=T
  data.dir = '/media/deth/SLAB/diss/3-phn/other/rf_vars'
  # if on Savio
} else {
  on.laptop=F
  data.dir = '/global/scratch/users/drewhart/seasonality/rf_data/'
}

# verbose?
verbose = T

# seed
seed.num = 12345
set.seed(seed.num)

# subset the data before building the RF?
subset = T
subset.fracs = c(0.005, 0.05)

# training data fraction
# NOTE: use 60% for training, 40% for testing, which should be plenty,
#       given how large a dataset we have
train.frac = 0.6


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





####################
# LOAD AND PREP DATA
####################

# load countries polygons (for use as a simple basemap)
world = map_data('world')

# get band names
names = c('phn.asy',
          'tmp.min.asy',
          'tmp.max.asy',
          'ppt.asy',
          'def.asy',
          'cld.asy',
          'vrm.med',
          'veg.ent')

# load rasters of prepped variables
vars = brick(paste0(data.dir, "/asynch_model_all_vars_",
                    asynch.var, '_',
                    as.character(neigh.rad), "km.tif"))
names(vars) = names

# load data frame of prepped variables
df_full_unproj = read.csv(paste0(data.dir, "/asynch_model_all_vars_prepped_",
                                 asynch.var, '_',
                                 as.character(neigh.rad), "km.csv"))

# transform to projected coordinate system, then add cols for metric-CRS coords
# NOTE: NOT DIRECTLY CONVERTING TO SF BECAUSE EPSG SUPPORT BROKEN FOR SOME REASON...
# Warning message:
#In CPL_crs_from_input(x) :
#  GDAL Error 4: Unable to open EPSG support file gcs.csv.  Try setting the GDAL_DATA environment variable to point to the directory containing EPSG csv files.
df_full_unproj_sf = st_as_sf(df_full_unproj, coords=c('x', 'y'), crs=4326)
df_full = st_transform(df_full_unproj_sf, 3857)

# get metric coordinates for inclusion as RF covariates
# to 'subsume' spatial process (even though overall results suggest
# their inclusion/exclusion makes little difference)
coords = st_coordinates(df_full)
df_full['x'] = coords[,'X']
df_full['y'] = coords[,'Y']

# drop the df geometry columns
df_full = st_drop_geometry(df_full)


#########################################
# TUNE HYPERPARAMETERS OF GLOBAL RF MODEL
#########################################

# NOTE: ADAPTED FROM https://bradleyboehmke.github.io/HOML/random-forest.html#hyperparameters

n.features = (ncol(df_full)-1)

# for each of two subset fractions:
# NOTE: only runs when interactive
if (F){
  for (subset.frac in subset.fracs){
    cat('=============\n')
    cat('SUBSET FRAC: ', subset.frac, '\n\n')

    # take a stratified random sample, stratified by the response var
    split_subset  <- initial_split(df_full, prop = subset.frac, strata = "phn.asy")
    df = training(split_subset)

    # get training and test data, using a stratified random sample
    split_strat  <- initial_split(df, prop = train.frac,
                                  strata = "phn.asy")
    trn  <- training(split_strat)
    tst  <- testing(split_strat)
   
    #build a default ranger model
    default.fit <- ranger(
      formula = phn.asy ~ .,
      data = trn,
      mtry = floor(n.features / 3),
      respect.unordered.factors = "order",
      seed = 123
    )
    # and get its OOB and test-set RMSEs
    default.rmse <- sqrt(default.fit$prediction.error)
    default.preds = predict(default.fit, tst)$predictions
    default.errs = default.preds - tst[,'phn.asy']
    default.tst.rmse = sqrt(mean(default.errs^2))
    hyper_grid <- expand.grid(
      # NOTE: DECIDED TO REPLACE THIS WITH HARD-CODED NUMBERS
      #       BECAUSE CHANGE IN COVAR NUM CAUSED THIS TO BE 0, 2, 4 INSTEAD!
      #mtry = floor(n.features * c(.1, .33, .5)),
      mtry = c(1, 3, 5)
      ntree = c(150, 200, 250, 300),
      min.node.size = c(1, 3, 5, 10),
      replace = c(TRUE, FALSE),
      sample.fraction = c(.3, .55, .8),
      r2 = NA,
      rmse = NA,
      tst.rmse = NA
    )
    for(i in seq_len(nrow(hyper_grid))) {
      # fit model for with hyperparameter combination
      fit <- ranger(
        formula         = phn.asy ~ .,
        data            = trn,
        num.trees       = hyper_grid$ntree[i],
        mtry            = hyper_grid$mtry[i],
        min.node.size   = hyper_grid$min.node.size[i],
        replace         = hyper_grid$replace[i],
        sample.fraction = hyper_grid$sample.fraction[i],
        verbose         = FALSE,
        seed            = 123,
        respect.unordered.factors = 'order',
      )
      hyper_grid$r2[i] <- fit$r.squared
      hyper_grid$rmse[i] <- sqrt(fit$prediction.error)
      preds = predict(fit, tst)$predictions
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
             perc_tst_gain = (default.tst.rmse - tst.rmse) / default.tst.rmse * 100,
             perc_rsq_inc = (default.fit$r.squared - r2) / default.fit$r.squared * 100,
      )
    print(head(hyper_grid_complete, 50))
    # save tuning results
    write.csv(as.data.frame(hyper_grid),
              paste0(data.dir,
                     'tuning_results_subset_frac_',
                     subset.frac,
                     '_',
                     asynch.var,
                     '_',
                     as.character(neigh.rad),
                     'km.csv'),
              row.names=F)
  }
}

# set RF hyperparams based on output above
ntree = 300
replace = F
rf.sample.fraction = 0.8
# NOTE: CHANGED BELOW FROM 4 TO 3 WHEN FIXING mtry ISSUE NOTED ABOVE
mtry = 3
min.node.size = 1

# and choose data subset based on output above
subset.frac = 0.05

# take a stratified random sample, stratified by the response var
split_subset  <- initial_split(df_full, prop = subset.frac, strata = "phn.asy")
df = training(split_subset)

# get training and test data, using a stratified random sample
split_strat  <- initial_split(df, prop = train.frac,
                              strata = "phn.asy")
trn  <- training(split_strat)
tst  <- testing(split_strat)
print(paste0('TRAINING DATA: ', nrow(trn), ' ROWS'))
print(paste0('TEST DATA: ', nrow(tst), ' ROWS'))

# scatter the training and test datasets
# NOTE: only for interactive
if (F){
  scat_trn = ggplot(tst) +
    geom_polygon(data=world, aes(x=long, y=lat, group=group), color="black", fill="white" ) +
    geom_point(data=trn, aes(x=x, y=y), col='blue', alpha=0.05, size=0.1)
  scat_tst = ggplot(tst) +
    geom_polygon(data=world, aes(x=long, y=lat, group=group), color="black", fill="white" ) +
    geom_point(data=tst, aes(x=x, y=y), col='red', alpha=0.05, size=0.1)
  grid.arrange(scat_trn, scat_tst)
}



##########################
# BORUTA FEATURE SELECTION
##########################

# NOTE: only for interactive
if (F){
  bor_res = Boruta(phn.asy ~ .,
                   data=trn,
                   doTrace=verbose * 3,
                   num.trees=ntree,
                   mtry=mtry,
                   replace=replace,
                   sample.fraction=rf.sample.fraction,
                   min.node.size=min.node.size
  )
  print(attStats(bor_res))
  
  
  jpeg(paste0(data.dir, '/boruta_boxplot_', asynch.var, '_', as.character(neigh.rad), 'km.jpg'),
       width=900,
       height=400,
       units='px',
       quality=90,
  )
  par(mfrow=c(1,2))
  plot(bor_res)
  plotImpHistory(bor_res)
  dev.off()
  par(mfrow=c(1,1))
  
  # remove any vars rejected by boruta
  # from trn, tst, df, and df_full
  
  # ... NONE DROPPED!
}



###########################
# DROP COLLINEAR COVARIATES
###########################

# filter covariates using correlation and variance inflation factor (VIF) thresholds
predictor.var.names = spatialRF::auto_cor(
         x = trn,
         cor.threshold = 0.75,
)  %>% spatialRF::auto_vif(
         vif.threshold=5,
)

# ...NONE COLLINEAR!


###########################################
# BUILD TUNED, PARSIMONIOUS GLOBAL RF MODEL
###########################################
if (coords.as.covars == 'y'){
  rf_final = ranger(phn.asy ~ .,
                    data=trn,
                    num.trees=ntree,
                    mtry=mtry,
                    importance='permutation',
                    verbose=verbose,
                    replace=replace,
                    sample.fraction=rf.sample.fraction,
                    min.node.size=min.node.size,
  )
} else {
  rf_final = ranger(phn.asy ~ tmp.min.asy + tmp.max.asy + ppt.asy +
                              def.asy + cld.asy +
                              vrm.med + veg.ent,
                    data=trn,
                    num.trees=ntree,
                    mtry=mtry,
                    importance='permutation',
                    verbose=verbose,
                    replace=replace,
                    sample.fraction=rf.sample.fraction,
                    min.node.size=min.node.size,
  )
}

# take a look at the results
print(rf_final)

# var importance plots, with permutation based metric...
p_imp_permut = plot.ranger.import(rf_final)
pfun <- function(object, newdata) {
  predict(object, data = newdata)$predictions
}
# and save permutation-based importance values
permut_imp <- as.data.frame(ranger::importance(rf_final))
write.csv(permut_imp, paste0(data.dir, 'rf_permut_importance_',
                             coords.as.covars, 'COORDS_',
                             asynch.var, '_',
                             as.character(neigh.rad), 'km.csv'), row.names=T)
# ... and with SHAP values
shap <- fastshap::explain(rf_final, X = df[, 2:ncol(df)], pred_wrapper = pfun, nsim = 10)
shap_imp <- data.frame(
  Variable = names(shap),
  Importance = apply(shap, MARGIN = 2, FUN = function(x) sum(abs(x)))
)
write.csv(shap_imp, paste0(data.dir, 'rf_SHAP_importance_',
                             coords.as.covars, 'COORDS_',
                            asynch.var, '_',
                            as.character(neigh.rad), 'km.csv'), row.names=F)
p_imp_shap = ggplot(shap_imp, aes(reorder(Variable, Importance), Importance)) +
  geom_col() +
  coord_flip() +
  xlab("") +
  ylab("mean(|Shapley value|)")
varimp_grob = grid.arrange(p_imp_shap, p_imp_permut, ncol=2)


ggsave(varimp_grob, file=paste0(data.dir, 'var_import_plots_permut_and_SHAP_',
                                coords.as.covars, 'COORDS_',
                                asynch.var, '_',
                                as.character(neigh.rad), 'km.png'),
       width=45, height=35, units='cm', dpi=600)


# assess model externally using withheld test data
preds = predict(rf_final, tst)$predictions
tst$err = preds - tst[,'phn.asy']
preds_plot = ggplot(tst) +
  geom_polygon(data=world, aes(x=long, y=lat, group=group), color="black", fill="white" ) +
  geom_point(aes(x=df_full_unproj[rownames(tst), 'x'], y=df_full_unproj[rownames(tst), 'y'],
                 col=err, alpha=abs(err)/max(abs(err))), size=1) +
  scale_color_gradient2(low='#cc003d', mid='#dbdbdb', high='#009e64') +
  #coord_map() +
  theme_bw()
preds_plot

ggsave(preds_plot, file=paste0(data.dir, 'preds_plot_', coords.as.covars, 'COORDS_', asynch.var, '_', as.character(neigh.rad), 'km.png'),
       width=30, height=22, units='cm', dpi=500)


# make predictions for full dataset (to map as raster)
# NOTE: only if interactive
#if (F){
  full_preds = predict(rf_final, df_full)$predictions
  df.res = df_full %>% mutate(preds = full_preds, err = full_preds - df_full[,'phn.asy'])
  write.csv(df.res, paste0(data.dir, 'rf_full_preds_',
                            coords.as.covars, 'COORDS_',
                            asynch.var, '_',
                            as.character(neigh.rad), 'km.csv'), row.names=F)
  
# map SHAP values
cat('\n\n\nNOW CALCULATING FULL SHAPLEY VALUES...\n\n\n')
shap_full = fastshap::explain(rf_final, X = df_full[, 2:ncol(df_full)], pred_wrapper = pfun, nsim = 10)
write.csv(shap_full, paste0(data.dir, 'rf_SHAP_vals_',
                            coords.as.covars, 'COORDS_',
                            asynch.var, '_',
                            as.character(neigh.rad), 'km.csv'), row.names=F)
shap_full_w_coords = cbind(df_full_unproj[, c('x', 'y')], shap_full)
write.csv(shap_full_w_coords, paste0(data.dir, 'rf_SHAP_vals_w_coords_',
                                     coords.as.covars, 'COORDS_',
                                     asynch.var, '_',
                                     as.character(neigh.rad), 'km.csv'), row.names=F)

