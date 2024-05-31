library(plotbiomes)
library(ggplot2)
library(viridis)
library(betareg)

df = read.csv('./FLUXNET_validation_results.csv')
df$map_cm = df$map/10
df$one_minus_r2 = 1 - df$r2
df$cell_dist = as.factor(df$cell_dist)

# plot Euclidean distance between RS-fitted seasonality and
# FLUXNET GPP seasonality on top of Whittaker biomes

#png('flux_val_results_biome.png',
#    width=1000, height=600)
#
#ggplot() +
#   geom_polygon(data = Whittaker_biomes, aes(x = temp_c,
#                                             y = precp_cm,
#                                             fill = biome,
#                                             alpha=-1.01),
#                colour = "gray97", # colour of polygon border
#                size   = -1.5) +    # thickness of polygon border
#   geom_point(data=df, aes(x=mat, y=map_cm, col=dist, size=dist, alpha=0.7)) + 
#   scale_color_viridis(option = "plasma")
#
#dev.off()

# build regression model of goodness of fit
beta_mod = betareg(r2 ~ lon + lat + mat + map + cell_dist + gpp_ts_len, data=df)
print(summary(beta_mod))

beta_mod_reduced = lm(r2 ~ mat + cell_dist + gpp_ts_len, data=df)
print(summary(beta_mod_reduced))

