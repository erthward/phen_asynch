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



# TODO:
# 2. if runtime is slow then use ranger instead (but then have to figure out how to not make range::importance clobber randomForest::importance that grf depends on!)
# 4. get running top to bottom
# 5. set up to run on Savio, then run it



# workflow:
# 1. draw subsample
# 2. check that it's small enough that independent bootstraps can be drawn from it (how??)
# 3. RF model (or GB?)
# 4. boruta to select vars
# 5. rerun parsimonious model
# 6. Lauren's approach to tune params?
# 7. plot and assess var importance
# 8. map SHAP values (and how to get globally integrated ones??)
# 9. run geo-weighted RF and map top vars





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

# save plots?
save.plots = F

# seed
seed.num = 12345
set.seed(seed.num)

# training data fraction
# NOTE: use 60% for training, 40% for testing, which should be plenty,
#       given how large a dataset we have
train.frac = 0.6


# global RF params
ntree = 75
mtry = 5

# local RF params
bw.local = 150
ntree.local = ntree
mtry.local = mtry


############
# HELPER FNS
############

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



##########################################################################
# ANALYSIS 

# rename all bands 
names = c('phn.asy', 'tmp.min.asy', 'tmp.max.asy',
          'tmp.min.mea', 'ppt.asy',
          'ppt.sea', 'def.asy', 'cld.asy', 'vrm.med',
          'riv.dis', 'eco.dis')

# load rasters of prepped variables
vars = brick(paste0(data.dir, "/asynch_model_all_vars.tif"))
names(vars) = names

# load data frames of prepped variables
df = read.csv(paste0(data.dir, "/asynch_model_all_vars_prepped.csv"))
#df.strat = read.csv(paste0(data.dir, "/asynch_model_all_vars_prepped_strat.csv"))

# load countries polygons (for use as a simple basemap)
world = map_data('world')


#####################
# RUN GLOBAL RF MODEL
#####################


# get training and test data, using a stratified random sample
# get training and test data, using a stratified random sample

split_strat  <- initial_split(df, prop = train.frac, 
                              strata = "phn.asy")

trn  <- training(split_strat)
tst  <- testing(split_strat)
trn = df.strat[trn.indices,]
tst = df.strat[-trn.indices,]

print(paste0('TRAINING DATA: ', nrow(trn), ' ROWS'))
print(paste0('TEST DATA: ', nrow(tst), ' ROWS'))

# build the model
# NOTE: leaving lat and lon out of model

# standard RF
rf = RRF(phn.asy ~ .,
         data=trn[,3:ncol(trn)],
         ntree=ntree,
         mtry=mtry,
         importance=T,
         flagReg=0)

# regularized RF
rrf = RRF(phn.asy ~ .,
         data=trn[,3:ncol(trn)],
         ntree=ntree,
         mtry=mtry,
         importance=T,
         flagReg=1,
         )

# guided regularized RF
impRF <- rf$importance[,"IncNodePurity"]
imp <- impRF/(max(impRF))#normalize the importance score
gamma <- 0.5
coefReg <- (1-gamma)+gamma*imp #weighted average
grrf <- RRF(phn.asy ~ .,
            data=trn[,3:ncol(trn)],
            ntree=ntree,
            mtry=mtry,
            importance=T,
            coefReg=coefReg,
            flagReg=1)


# take a look at the results
# standard RF
print(rf) 
  # Mean of squared residuals: 0.7034528
    # low MSE
  # % Var explained: 29.76
# regularized RF
print(rrf) 
  # Mean of squared residuals: 0.7073498
    # low MSE
  # % Var explained: 29.38
# guided regularized RF
print(grrf)
  # Mean of squared residuals: 0.6945594
    # slightly lower MSE
  # % Var explained: 30.65

# quick assessment plots
par(mfrow=c(1,2))
p1 = varImpPlot(grrf)
# TODO: WHY STRANGE VERTICAL SUDDENLY SHOWED UP IN FOLLOWING PLOT?
p2 = ggplot() +
   geom_point(aes(x=trn$phn.asy, y=predict(grrf)), alpha=0.2) +
   geom_abline(intercept = 0, slope = 1)
# error below
quick_assess_grob = grid.arrange(p1, p2, ncol=2)
plot(quick_assess_grob)

if (save.plots){
   ggsave(quick_assess_grob, file='quick_assess_plots.png',
          width=45, height=35, units='cm', dpi=1000)
}

# partial dependence plots
# ppt.asy vs ppt.sea
pdp12 = make.pdp(grrf, trn, 'phn.asy', 1, 2, seed.num=seed.num)
plot(pdp12)

#
pdp13 = make.pdp(grrf, trn, 'phn.asy', 1, 3, seed.num=seed.num)
plot(pdp13)

pdp23 = make.pdp(grrf, trn, 'phn.asy', 2, 3, seed.num=seed.num)
plot(pdp23)

if (save.plots){
   ggsave(pdp12, file='pdp12.png',
          width=22, height=12, units='cm', dpi=500)
   ggsave(pdp13, file='pdp13.png',
          width=22, height=12, units='cm', dpi=500)
   ggsave(pdp23, file='pdp23.png',
          width=22, height=12, units='cm', dpi=500)
}


# make predictions
preds = predict(grrf, tst[,3:ncol(tst)])
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
          width=30, height=22, units='cm', dpi=500)
  }


# make predicitions for full dataset (to map as raster)
full_preds = predict(grrf, df[,3:ncol(df)])
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
global_grrf_main_plots = cowplot::plot_grid(obs_map, preds_map, err_map, obs_vs_pred, nrow=2, ncol=2)
global_grrf_main_plots

if (save.plots){
   ggsave(global_grrf_main_plots, file='global_grrf_main_plot.png',
          width=75, height=55, units='cm', dpi=500)
  }


######################
# RUN LOCAL GWRF MODEL
######################

# code adapted from https://zia207.github.io/geospatial-r-github.io/geographically-wighted-random-forest.html
coords = trn[,1:2]
grf.model <- SpatialML::grf(formula=phn.asy ~ tmp.min.asy + tmp.max.asy + tmp.min.mea +
                               ppt.asy + ppt.sea + def.asy + cld.asy + 
                               vrm.med + riv.dis + eco.dis,
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
    plot(rfe.results, type=c("g", "o"))
    
    # variable importance plot
    randomForest::varImpPlot(rfe.results$fit)



  ## Murphy et al., (2010) RF MODEL SELECTION (JEFF EVANS R PACKAGE)
    rf_model_sel <- rf.modelSel(xdata = x_trn, 
                                ydata = y_trn, 
                                final.model = T, 
                                seed = seed.num, 
                                ntree = opt.ntree.grrf, 
                                mtry = opt.mtry.grrf)
    rf_model_sel

    # variable importance plot
      randomForest::varImpPlot(rf_model_sel$rf.final)
