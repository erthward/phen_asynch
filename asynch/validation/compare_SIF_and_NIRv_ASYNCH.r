library(sf)
library(raster)
library(RColorBrewer)
library(latex2exp)
library(ggplot2)

# simple diverging palette
pal = colorRampPalette(c('red', 'white', 'blue'))

# data directory
data.dir = "/media/deth/SLAB/diss/3-phn/GEE_outputs/final/"
countries.data.dir = "/home/deth/Desktop/CAL/research/projects/seasonality/seasonal_asynchrony/data/bounds/"

# load country boundaries
world = st_read(paste0(countries.data.dir, 'NewWorldFile_2020.shp'))

# loop over neighborhood radii (in km)
for (neigh.rad in c('50', '100', '150')){

   cat ('\n\nRUNNING COMPARISON FOR ', neigh.rad, ' KM-RADIUS NEIGHBORHOOD...\n\n')

   # load both rasters
   sif = stack(paste0(data.dir, 'SIF_STRICT_asynch_', neigh.rad, 'km.tif'))[[1]]
   nirv = stack(paste0(data.dir, 'NIRv_STRICT_asynch_', neigh.rad, 'km.tif'))[[1]]
   
   # center and scale each one ((x-mu)/sigma)
   sif_scale = raster::scale(sif)
   nirv_scale = raster::scale(nirv)
   
   # get difference and plot
   diff_scale = nirv_scale - sif_scale
   max_abs_val=max(abs(c(cellStats(diff_scale, min), cellStats(diff_scale, max))))
   diff_scale[diff_scale>max_abs_val] = max_abs_val
   diff_scale[diff_scale<-max_abs_val] = -max_abs_val
   diff_scale_df = as.data.frame(diff_scale, xy=T)

   ggplot() +
      #geom_sf(data=world, col='black') + 
      geom_raster(data=diff_scale_df, aes(x, y, fill=layer)) +
          scale_fill_gradientn(
              colours = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(3),
              values = c(1.0, (0 - min(diff_scale_df$layer)) / (max(diff_scale_df$layer) - min(diff_scale_df$layer)), 0))


   png(paste0("scaled_NIRv-scaled_SIF_map_", neigh.rad, "km.png"), units="in", height=4.5, width=9, res=600)

   nbreaks=12
   breaks=seq(-max_abs_val, max_abs_val, length.out=nbreaks)
   world_reproj = st_transform(world, st_crs(diff_scale))
   plot(world_reproj$geometry, col='black')
   plot(diff_scale,
        #col=pal(nbreaks),
        col='RdBu;,
        breaks=breaks,
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
                 geom_point(aes(y=scaled_NIRv, x=scaled_SIF), alpha=0.01) +
                 geom_smooth(method='lm', se=F, data=samp@data, aes(y=scaled_NIRv, x=scaled_SIF)) +
                 annotate("text", x=7.5, y=-7.5, size=8, col='red',
                          label= TeX(paste0("$R^{2}=", round(summary(mod)$r.squared, 3), "$"))) +
                 ylab(TeX("$NIR_{V}, scaled")) +
                 xlab(TeX("$SIF, scaled"))
   slr_plot = slr_plot + theme(axis.title.x = element_text(size=20, face="bold"),
                               axis.title.y = element_text(size=20, face="bold"),
                               axis.text.x = element_text(size=14),
                               axis.text.y = element_text(size=14)) +
                         coord_cartesian(xlim = c(-10, 10), ylim=c(-10, 10)) 
   slr_plot
   
   ggsave(slr_plot,
          file=paste0('SLR_scaled_SIF_asynch_vs_scaled_NIRv_asynch_', neigh.rad, 'km.png'),
          width=45, height=45, units='cm', dpi=600)
   
   
   # save output
   writeRaster(diff_scale, paste0('scaled_NIRv-scaled_SIF_', neigh.rad, 'km.tif'),
               'GTiff', overwrite=T)
}
