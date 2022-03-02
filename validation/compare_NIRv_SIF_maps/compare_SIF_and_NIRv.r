library(raster)
library(RColorBrewer)
library(latex2exp)
library(ggplot2)

# simple diverging palette
pal = colorRampPalette(c('red', 'white', 'blue'))

# load both
sif = stack('SIF_global_asynch.tif')[[3]]
nirv = stack('NIRv_global_asynch.tif')[[3]]

# center and scale each one ((x-mu)/sigma)
sif_scale = scale(sif)
nirv_scale = scale(nirv)

# get difference and plot
diff_scale = nirv_scale - sif_scale

world_reproj = st_transform(world, st_crs(diff_scale))
png("scaled_NIRv-scaled_SIF_map.png", units="in", height=4.5, width=9, res=1000)
plot(world_reproj, col='black')
plot(diff_scale,
     col=pal(20),
     main=TeX("$asynch_{NIR_{V}}\\ -\\ $asynch_{SIF}$\\ (both\\ centered\\ and\\ scaled"),
     add=T)
dev.off()

# histogram
hist(diff_scale)

# scatter samples against one another and fit SLR
both = stack(nirv_scale, sif_scale)
samp = sampleRandom(both, size=500000, na.rm=T, sp=T)
names(samp) = c('scaled_NIRv', 'scaled_SIF')

mod = lm(scaled_NIRv ~ scaled_SIF, data=samp)
summary(mod)
min_lim = floor(min(samp@data))
max_lim  = ceiling(max(samp@data))
slr_plot = ggplot(samp@data) +
              geom_abline(intercept=0, slope=1) + 
              geom_point(aes(x=scaled_NIRv, y=scaled_SIF), alpha=0.01) +
              geom_smooth(method='lm', se=F, data=samp@data, aes(x=scaled_NIRv, y=scaled_SIF)) +
              annotate("text", x=10, y=-7.5, size=8,
                       label= TeX(paste0("$R^{2}=", round(summary(mod)$r.squared, 3), "$"))) +
              xlab(TeX("$NIR_{V}, scaled")) +
              ylab(TeX("$SIF, scaled"))
slr_plot = slr_plot + theme(axis.title.x = element_text(size=20, face="bold"),
                            axis.title.y = element_text(size=20, face="bold"),
                            axis.text.x = element_text(size=14),
                            axis.text.y = element_text(size=14)) +
                      coord_cartesian(xlim = c(min_lim, max_lim), ylim=c(min_lim, max_lim)) 
slr_plot

ggsave(slr_plot, file='SLR_scaled_SIF_asynch_vs_scaled_NIRv_asynch.png',
       width=45, height=45, units='cm', dpi=1000)


# save output
writeRaster(diff_scale, 'scaled_NIRv-scaled_SIF.tif', 'GTiff', overwrite=T)
