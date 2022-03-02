library(raster)               # raster data
library(terra)                # newer raster data
library(sp)                   # spatial data
library(sf)                   # newer spatial data
library(spdep)                # spatial autocorrelation
library(randomForest)         # global RFs
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

# TODO:
# 2. if runtime is slow then use ranger instead (but then have to figure out how to not make range::importance clobber randomForest::importance that grf depends on!)
# 4. get running top to bottom
# 5. set up to run on Savio, then run it


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
}

# mode (either prep data, or analyze data)
#mode = 'prep'
#mode = 'analyze'
mode = 'both'

# save plots?
save.plots = T

# seed
seed.num = 12345
set.seed(seed.num)

# number of strata and proportion of each stratum
# to use for stratified random sampling
n.strata = 5
strat.samp.prop = 0.1
# NOTE: full raster should have ~15925248 pixels, so
#       this would use ~1/15th of the full dataset
strat.samp.n = 1000000/n.strata 

# training data fraction
train.frac = 0.3 # use 70% for training, 30% for testing


# global RF params
ntree = 50
mtry = 5

# local RF params
bw.local = 150
ntree.local = ntree
mtry.local = mtry


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
   # align to the align.to target raster, if align.to is a raster
   if (class(align.to)[1] == 'RasterLayer'){
      data = align.rast.x.to.y(data, align.to, mask.it=mask.it) 
   }
   return(data)
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
     return(results)
       }
   pdp_pred <- function(object, newdata)  { 
      results <- mean(as.vector(predict(object, newdata))) 
      return(results) 
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
                          parallel = F,
                          grid.resolution = 10 
   )
   pd_var_2<- pdp::partial(
                          rf, 
                          train = trn,
                          pred.var = var_2, 
                          pred.fun = pdp_pred, 
                          parallel = F,
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

if (mode %in% c('prep', 'both')){
  
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
  
  # vector ruggedness metric, ~50km agg med and sd
  # NOTE: CALCULATED AS FIXED PIXELS, NOT MOVING WINDOWS!
  # NOTE: masking isn't necessary because incomplete rows of the stacked
  #       variables are dropped later on, but makes it easier to inspect
  #       the individual layers now because of the crazy blocky artifacts
  #       surrounding continents in their dataset
  vrm.med = read.file('vrm_50KMmd', asynch.file=F,
                      align.to=phn.asy, mask.it=T)
  vrm.std = read.file('vrm_50KMsd', asynch.file=F,
                      align.to=phn.asy, mask.it=T)
  
  
  # neighborhood Shannon entropy of land cover
  # NOTE: CALCULATE THIS!
  #lc.ent = read.file('MODIS_IGBP', asynch.file=F,
                      #align.to=phn.asy, mask.it=F)
  # NOTE: for now, using 25km dissimilarity of EVI between neigh pixels
  # from earthenv.org instead
  # NOTE: masking isn't necessary because incomplete rows of the stacked
  #       variables are dropped later on, but makes it easier to inspect
  #       the individual layer now because of the ocean values
  hab.div = read.file('hab_shannon', asynch.file=F,
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
  eco.dis = mask(eco.dis, phn.asy)
  
  
  
  ###########
  # PREP DATA
  ###########
  
  # gather into a stack
  vars = stack(phn.asy, tmp.mea.asy, tmp.min.asy, tmp.max.asy, tmp.min.mea, tmp.max.mea,
               tmp.sea, ngd, ppt.asy, ppt.sea, def.asy, cld.asy, wds.asy, vrm.med, vrm.std,
               hab.div, riv.dis, eco.dis)

  # aggregate to coarser res, if working on laptop
  if (on.laptop){
     vars <- aggregate(vars, fact=12)
  }
 
  # rename all bands 
  names = c('phn.asy', 'tmp.mea.asy', 'tmp.min.asy', 'tmp.max.asy',
            'tmp.min.mea', 'tmp.max.mea', 'tmp.sea', 'ngd', 'ppt.asy',
            'ppt.sea', 'def.asy', 'cld.asy', 'wds.asy', 'vrm.med', 'vrm.std',
            'hab.div', 'riv.dis', 'eco.dis')

  names(vars) = names

  # write stack to file (for reload in 'analyze' mode),
  #raster::writeRaster(vars, paste0(data.dir, "/asynch_model_all_vars.tif"), overwrite=T)
  vars = rast(vars)
  terra::writeRaster(vars, paste0(data.dir, "/asynch_model_all_vars.tif"), overwrite=T)
  
  # coerce to a data.frame
  df = as.data.frame(vars, xy=T)
  
  # TODO: FIGURE OUT WHY WINDSPEED LAYER IS EMPTY, THEN DELETE NEXT LINE
  df = df[,!names(df)%in%c('wds.asy')]

  # drop NAs
  df = df[complete.cases(df),]
  
  # take a stratified random sample, stratifying by phn.asy values
  # (it has a somewhat right-skewed distribution, and at any rate,
  #  important that model is trained on data including high asynch values)
  df.strat = df %>%
     mutate(strata = base::cut(phn.asy, breaks=n.strata, labels=seq(n.strata))) %>%
     group_by(strata) %>%
     slice_sample(n=strat.samp.n) %>%
     ungroup() %>%
     as.data.frame()
  df.strat = df.strat[, !colnames(df.strat) %in% c('strata')]
  
  
  # scale variables
  df[,3:ncol(df)] = scale(df[,3:ncol(df)])
  df.strat[,3:ncol(df.strat)] = scale(df.strat[,3:ncol(df.strat)])
  
  # write data.frame to file (for reload in 'analyze' mode),
  write.csv(df, paste0(data.dir, "/asynch_model_all_vars_prepped.csv"), row.names=F)
  write.csv(df.strat, paste0(data.dir, "/asynch_model_all_vars_prepped_strat.csv"), row.names=F)
 

} 

if (mode %in% c('analyze', 'both')){

  # load rasters of prepped variables
  vars = brick(paste0(data.dir, "/asynch_model_all_vars.tif"))
  names(vars) = names

  # load data frames of prepped variables
  df = read.csv(paste0(data.dir, "/asynch_model_all_vars_prepped.csv"))
  df.strat = read.csv(paste0(data.dir, "/asynch_model_all_vars_prepped_strat.csv"))
 
  # load countries polygons (for use as a simple basemap)
  world = map_data('world')
  
  
  #####################
  # RUN GLOBAL RF MODEL
  #####################
  
  
  # get training and test data
  trn.indices = sample(1:nrow(df.strat), size = round(train.frac * nrow(df.strat)))
  trn = df.strat[trn.indices,]
  tst = df.strat[-trn.indices,]

  print(paste0('TRAINING DATA: ', nrow(trn), ' ROWS'))
  print(paste0('TEST DATA: ', nrow(tst), ' ROWS'))
  
  # build the model
  # NOTE: leaving lat and lon out of model
  rf = randomForest(phn.asy ~ ., data=trn[,3:ncol(trn)],
  #rf = randomForest(phn.asy ~ ., data=trn,
                    ntree=ntree, importance=TRUE)
  
  # take a look at the result
  print(rf)
  
  # quick assessment plots
  par(mfrow=c(1,2))
  p1 = vip(rf)
  p2 = ggplot() +
     geom_point(aes(x=trn$phn.asy, y=predict(rf)), alpha=0.2) +
     geom_abline(intercept = 0, slope = 1)
  quick_assess_grob = grid.arrange(p1, p2, ncol=2)
  plot(quick_assess_grob)

  if (save.plots){
     ggsave(quick_assess_grob, file='quick_assess_plots.png',
            width=45, height=35, units='cm', dpi=1000)
  }
  

  
  # partial dependence plots
  pdp12 = make.pdp(rf, trn, 'phn.asy', 1, 2, seed.num=seed.num)
  plot(pdp12)
  
  pdp13 = make.pdp(rf, trn, 'phn.asy', 1, 3, seed.num=seed.num)
  plot(pdp13)
  
  pdp23 = make.pdp(rf, trn, 'phn.asy', 2, 3, seed.num=seed.num)
  plot(pdp23)

  if (save.plots){
     ggsave(pdp12, file='pdp12.png',
            width=45, height=25, units='cm', dpi=1000)
     ggsave(pdp13, file='pdp13.png',
            width=45, height=25, units='cm', dpi=1000)
     ggsave(pdp23, file='pdp23.png',
            width=45, height=25, units='cm', dpi=1000)
  }
  
  
  # make predictions
  preds = predict(rf, tst[,3:ncol(tst)])
  #preds = predict(rf, tst)
  tst$err = preds - tst[,'phn.asy'] 
  preds_plot = ggplot(tst) +
      geom_polygon(data=world, aes(x=long, y=lat, group=group), color="black", fill="white" ) +
      geom_point(aes(x=x, y=y, col=err), size=1) +
      scale_color_gradient2(low='#cc003d', mid='#dbdbdb', high='#009e64') +
      #coord_map() + 
      theme_bw()
  preds_plot

  if (save.plots){
     ggsave(preds_plot, file='preds_plot.png',
            width=75, height=45, units='cm', dpi=1000)
    }

  
  # make predicitons for full dataset (to map as raster)
  full_preds = predict(rf, df[,3:ncol(df)])
  df.res = df %>% mutate(preds = full_preds, err = full_preds - df[,'phn.asy'])
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
  global_rf_main_plots = cowplot::plot_grid(obs_map, preds_map, err_map, obs_vs_pred, nrow=2, ncol=2)
  global_rf_main_plots

  if (save.plots){
     ggsave(global_rf_main_plots, file='global_rf_main_plot.png',
            width=75, height=55, units='cm', dpi=1000)
    }
  
  
  ######################
  # RUN LOCAL GWRF MODEL
  ######################
  
  # code adapted from https://zia207.github.io/geospatial-r-github.io/geographically-wighted-random-forest.html
  coords = trn[,1:2]
  grf.model <- grf(formula=phn.asy ~ tmp.mea.asy + tmp.min.asy + tmp.max.asy + tmp.min.mea + tmp.max.mea +
                                 tmp.sea + ngd + ppt.asy + ppt.sea + def.asy + cld.asy + 
                                 vrm.med + vrm.std + riv.dis + hab.div + eco.dis,
                   dframe=trn[, 3:ncol(trn)], 
                   bw=bw.local,
                   kernel="adaptive",
                   ntree=ntree.local,
                   mtry=mtry.local,
                   forests = TRUE,
                   coords=coords)
  
  # global model summary
  glob.imp = as.data.frame(importance(grf.model$Global.Model, type=1)) # %IncMSE
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
               plots[[10]], plots[[11]], plots[[12]],
               plots[[13]], plots[[14]], plots[[15]], plots[[16]],
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
               ncol=3)
 
  if (save.plots){
     ggsave(local_rf_main_plots, file='local_rf_gof_plot.png',
            width=75, height=55, units='cm', dpi=1000)
    }
 
  
  ####################################
  # MODEL/VAR SELECTION AND EVALUATION
  ####################################
  
  # what to do here exactly?
  
  # and should I use K-fold CV or bootstrapping? does this choice matter much?
  
}
