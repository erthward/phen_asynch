library(ggplot2)
library(lme4)
library(ncf)

# TODO:

  # should this actually just be a multi-Mantel test of some sort, since they're distance vals so non-independent?

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

# HOWEVER: need to do this with a set of partial Mantels/MMRRs to be valid
