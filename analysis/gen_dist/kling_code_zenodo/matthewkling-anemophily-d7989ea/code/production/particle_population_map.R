
ppm <- function(){
   
   library(raster)
   library(tidyverse)
   library(sf)
   library(lwgeom) 
   library(rworldmap)
   library(rgeos)
   
   
   # generate populations csv, on PC
   if(F){
      loadd(annotated)
      pops <- annotated %>%
         select(expected_wind, component, from, lon_from, lat_from) %>%
         filter(is.finite(lon_from)) %>%
         distinct()
      
      saveRDS(pops, "data/wind/particles/populations.rds")
      saveRDS(pops, "data_github/wind/particles/populations.rds")
   }
   
   
   ###############
   
   f <- readRDS("~/documents/wind/particle_figures/data_uniform_length_tmean_light.rds")
   
   pops <- readRDS("data_github/wind/particles/populations.rds") %>% 
      mutate(expected_wind = ifelse(expected_wind == .75, 1, expected_wind)) %>% 
      mutate(expected_wind = factor(expected_wind, labels = c("none", "partial", "full"))) %>%
      sample_n(nrow(.))
   
   world <- getMap(resolution = "low") %>% gBuffer(width = 0) %>% gUnaryUnion() %>% st_as_sf(world)
   md <- world %>% as_Spatial() %>% fortify()
   
   
   fs <- f %>% 
      select(id, time, x, y) %>%
      group_by(id) %>%
      mutate(brk = x - lag(x),
             brk = abs(brk) > 300) %>%
      filter(is.finite(brk)) %>%
      mutate(nbrk = sum(brk)) %>%
      filter(nbrk < 2) %>%
      mutate(brktime = ifelse(any(brk), time[brk], 0),
             group = ifelse(time >= brktime, "a", "b"),
             group = paste(id, group)) %>%
      select(id, x, y, time, group)
   
   fs <- fs %>% filter(time %% 2 == 0)
   
   ##### arrows ######
   
   # arrows
   arw <- expand_grid(x = seq(-165, 165, 30),
                      y = seq(-75, 75, 15))
   ffa <- f %>%
      mutate(xr = round(x),
             yr = round(y))
   
   arw <- ffa %>%
      filter(paste(xr, yr) %in% paste(arw$x, arw$y)) %>%
      group_by(id) %>%
      mutate(n = n()) %>%
      #filter(n > 1) %>%
      filter(time != 50) %>%
      group_by(xr, yr) %>%
      arrange(id) %>%
      slice(1) %>%
      right_join(ffa) %>%
      group_by(id) %>%
      filter(any(is.finite(n))) %>%
      filter(is.finite(n) |
                time == time[is.finite(n)]+1) %>%
      mutate(n = n()) %>%
      filter(n == 2) %>%
      mutate(x = c(x[1], x[1] + (x[2] - x[1])/1000),
             y = c(y[1], y[1] + (y[2] - y[1])/1000))
   
   
   
   ########
   
   
   p <- ggplot() + 
      geom_polygon(data = md, aes(long, lat, group = group), color = NA, fill="black") +
      geom_path(data = fs, aes(x, y, group=group), color="white", size = .025) + # wind
      geom_polygon(data = md, aes(long, lat, group = group), color = "black", fill=NA, size = .25) +
      geom_path(data = arw, aes(x, y, group = id), 
                arrow = arrow(type = "closed", length = unit(4, "mm"), angle = 10), 
                color = "black", size = .25) +
      geom_point(data = pops, aes(lon_from, lat_from, color = expected_wind), size=2) + # pops
      scale_color_manual(values = c("red", "purple", "dodgerblue")) +
      guides(colour = guide_legend(override.aes = list(size=6))) +
      coord_map(ylim = c(-52, 68), xlim = c(-163, 163)) +
      theme_void(base_size = 30) +
      theme(legend.position="top",
            panel.background = element_rect(fill = "gray40")) +
      labs(color = "genome wind dispersal level   ")
   
   ggsave("figures/manuscript/figure_2a.png", width=20, height=11, units="in")
}
