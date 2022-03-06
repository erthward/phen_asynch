library(ggplot2)
library(lme4)

# TODO:
  # change Python script to sample randomly globally,
  # then assign regions in a data-driven way

  # account for spatial autocorrelation somehow

  # plot fitted results

df = read.csv('./clim_dist_results.csv')

mod = lmer(seas_dist ~ geo_dist + (clim_dist|reg), data=df)
summary(mod)
coef(mod)
# KEY INSIGHT:
# seasonal distance strongly predicted by climatic distance
# in temperate high-asynch regions,
# but no better predicted by climatic distance in 
# tropical high-asynchrony than in any low-asynchrony regions
