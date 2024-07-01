# CHRIS PACIOREK'S WORKAROUND TO FIX THE rgdal/sp/sf ISSUE I RAN INTO:
Sys.setenv(GDAL_DATA = "/global/home/groups/consultsw/sl-7.x86_64/modules/gdal/2.2.3/share/gdal")

library(rgdal)                # GDAL bindings
library(sp)                   # spatial data
library(raster)               # raster data
library(sf)                   # newer spatial data
library(maps)                 # countries map data
library(rsample)              # function for stratified random sampling
library(ranger)               # faster RFs
library(spatialRF)            # toolkit for running RFs on spatial data
library(Boruta)               # feature selection w/ boruta algo (i.e., shadow features)
library(fastshap)             # SHAP values
library(vip)                  # var importance plots
library(pdp)                  # partial depend. plots (& ICE curves)
library(ggplot2)              # plotting
library(ggthemes)             # plot themes
library(grid)                 # textGrob for plot grids
library(RColorBrewer)         # colors
library(dplyr)                # reshaping dfs


##########################################################################
# SETUP


#######################
# SET BEHAVIORAL PARAMS
#######################

args = commandArgs(trailingOnly=T)

# phen-asynch var to use
asynch.var = args[1]
# COMMENT PRIOR LINE AND UNCOMMENT NEXT LINE FOR INTERACTIVE HYPERPARAMETER TUNING!
#asynch.var = 'NIRv'
cat('\nVAR: ', asynch.var, '\n')

# asynchrony neighborhood radius to use (in km)
neigh.rad = args[2]
# COMMENT PRIOR LINE AND UNCOMMENT NEXT LINE FOR INTERACTIVE HYPERPARAMETER TUNING!
#neigh.rad = '100'
cat('\nNEIGH RAD: ', neigh.rad, '\n')

# include coordinates in RF?
coords.as.covars = args[3]
# COMMENT PRIOR LINE AND UNCOMMENT NEXT LINE FOR INTERACTIVE HYPERPARAMETER TUNING!
#coords.as.covars = 'y'
cat('\nCOORDS AS COVARS? ', coords.as.covars, '\n')


# input and output dirs:
# if on laptop
if (strsplit(getwd(), '/')[[1]][2] == 'home'){
  on.laptop=T
  data.dir = '/media/deth/SLAB/diss/3-phn/other/rf_vars/'
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

cat('\nReading prepped variables as a raster brick...\n')
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
vars = brick(paste0(data.dir, "asynch_model_all_vars_",
                    asynch.var, '_',
                    as.character(neigh.rad), "km.tif"))
names(vars) = names

# load data frame of prepped variables
cat('\nReading CSVs of prepped, extracted raster values...\n')
df_full_unproj = read.csv(paste0(data.dir, "asynch_model_all_vars_prepped_",
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
      mtry = c(1, 3, 5),
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
mtry = 5
min.node.size = 3

cat('\nSplitting input data into training and test subsets...\n')
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
  
  
  jpeg(paste0(data.dir, 'boruta_boxplot_', asynch.var, '_', as.character(neigh.rad), 'km.jpg'),
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
cat('\nFitting random forest...\n')
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
cat('\nCalculating and saving variable importance values...\n')
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
cat('\nUsing test data to assess model...\n')
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
cat('\nCalculating and saving model predictions for full global raster...\n')
full_preds = predict(rf_final, df_full)$predictions
df.res = df_full %>% mutate(preds = full_preds, err = full_preds - df_full[,'phn.asy'])
write.csv(df.res, paste0(data.dir, 'rf_full_preds_',
                          coords.as.covars, 'COORDS_',
                          asynch.var, '_',
                          as.character(neigh.rad), 'km.csv'), row.names=F)

# map SHAP values
cat('\nCalculating and saving SHAP values for full global dataset...')
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
cat('\nModeling complete.\n')
