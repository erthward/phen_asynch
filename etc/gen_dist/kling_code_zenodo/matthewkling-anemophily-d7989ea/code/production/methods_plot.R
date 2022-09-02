
unwrap <- function(x, width=20){
        x <- x %>% crop(extent(180-width, 180, -90, 90)) %>% shift(-360) %>% merge(x)
        x <- x %>% crop(extent(-180, -180+width, -90, 90)) %>% shift(360) %>% merge(x)
        x
}

methods_figure <- function(annotated, traits, wind_graphs, h1_rand, h2_rand, h3_rand, div_rand){
        
        comp <- "17d"
        i <- 14
        
        d <- annotated %>% mutate(mig = dmig) %>%
                filter(component == comp) %>%
                filter(is.finite(clim), is.finite(distance), is.finite(wind), 
                       is.finite(mig), is.finite(fst))
        
        
        ### maps ###
        
        pops <- d %>%
                select(lon_from, lat_from) %>%
                mutate(lon_from_original = lon_from) %>%
                distinct() %>%
                mutate(d = sqrt((lon_from - mean(lon_from))^2 + (lon_from - mean(lon_from))^2)) %>%
                arrange(d)
        coordinates(pops) <- c("lon_from", "lat_from")
        crs(pops) <- crs('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
        pops <- pops %>% spTransform(crs(readRDS(wind_graphs[1])))
        
        pps <- as.data.frame(pops)
        pps$id <- 1:nrow(pps)
        focal <- filter(pps, id == i)
        zzz <- focal$lon_from_original
        
        d <- d %>% mutate(focal = lon_from == zzz | lon_to == zzz )
        
        ext <- extent(pops) * 1.25
        ext@xmin <- -35
        
        land <- raster("f:/cfsr/land.tif") %>%  
                rotate() %>% 
                unwrap(180) %>%
                reclassify(c(NA, NA, 0)) %>%
                crop(ext)
        
        
        if(F){ # only do this once
                for(p in 0:3){
                        rose <- paste0("data/wind/roses/p", p, "/3_4_5.tif")
                        g <- make_wind_graph(rose,
                                             paste0("data/wind/roses_sandbox/out_p", p, "_sum.rds"),
                                             direction = "downwind")
                        g <- make_wind_graph(rose,
                                             paste0("data/wind/roses_sandbox/in_p", p, "_sum.rds"),
                                             direction = "upwind")
                        g <- make_wind_graph(rose, 
                                             paste0("data/wind/roses_sandbox/out_p", p, "_prod.rds"), 
                                             direction = "downwind", logtrans = TRUE)
                        g <- make_wind_graph(rose, 
                                             paste0("data/wind/roses_sandbox/in_p", p, "_prod.rds"), 
                                             direction = "upwind", logtrans = TRUE)
                }
        }
        
        p = 1
        lcf = "_sum" 
        
        dp <- filter(d, wind_p == paste0("wind_p", p))
        
        graph_out <- readRDS(paste0("data/wind/roses_sandbox/out_p", p, lcf, ".rds"))
        graph_in <- readRDS(paste0("data/wind/roses_sandbox/in_p", p, lcf, ".rds"))
        
        wind_out <- (accCost(graph_out, pops[i,]) / 3600) %>% crop(ext)
        wind_in <- (accCost(graph_in, pops[i,]) / 3600) %>% crop(ext)
        distance <- distanceFromPoints(wind_out, coordinates(pops)[i, 1:2]) / 1000
        
        r <- stack(wind_out, wind_in, distance, land) %>%
                rasterToPoints() %>%
                as.data.frame() %>%
                rename(wind_out = layer.1,
                       wind_in = layer.2,
                       distance = layer.3,
                       land = layer.4) %>%
                mutate(wind_ratio = wind_out / wind_in,
                       wind_mean = (wind_out + wind_in) / 2)
        r$land[r$x < -10] <- 0
        r$land[r$x < 0 & r$y > 60] <- 0
        
        rr <- r %>% gather(stat, value, -x, -y) %>%
                group_by(stat) %>% mutate(value = scale(value))
        
        
        
        ptcol <- "red"
        pal <- c("#40004a", "darkblue", "dodgerblue", "lightblue")
        pal <- c("darkblue", "dodgerblue", "lightblue") %>% rev()
        bgcol <- pal[2]
        
        map_theme <- theme_bw(base_size = 16) +
                theme(#panel.border = element_rect(color = "black", fill = NA),
                        #legend.title = element_text(color = "white"),
                        #legend.text = element_text(color = "white"),
                        axis.title = element_blank(),
                        axis.text = element_blank(),
                        axis.ticks = element_blank(),
                        strip.background = element_rect(fill = "black"),
                        strip.text = element_text(color = "white"),
                        legend.justification = c(0,1),
                        legend.position = c(0, 1),
                        legend.background = element_rect(fill = NA, color = NA),
                        legend.box.background = element_rect(fill = NA, color = NA),
                        legend.key = element_rect(#fill = "cadetblue", 
                                fill = bgcol, #"#41a374",
                                color = NA),
                        legend.title = element_text(size = 15),
                        legend.text = element_text(size = 15))
        ip = i
        
        bwf <- diff(range(c(r$wind_out[r$land == 1], r$wind_in[r$land == 1]))) / 30
        mf <- ggplot() +
                #facet_grid(. ~"flow", switch = "y") +
                #geom_raster(data = r %>% filter(land == 1), aes(x, y), fill = "cadetblue") +
                geom_raster(data = r %>% filter(land == 1), aes(x, y), fill = bgcol) +
                geom_contour(data = r, aes(x, y, z = wind_in, color = "inbound\ntravel time (h)"), binwidth =  bwf) +
                geom_contour(data = r, aes(x, y, z = wind_out, color = "outbound\ntravel time (h)"), binwidth = bwf) +
                #geom_raster(data = r %>% filter(land == 0), aes(x, y), fill = "white") +
                geom_raster(data = r %>% filter(land == 0), aes(x, y), fill = "white") +
                geom_point(data = pps, aes(lon_from, lat_from), color = "black", size = 2) +
                geom_point(data = filter(pps, id == i), aes(lon_from, lat_from), color = ptcol, size = 2) +
                scale_color_manual(values = c("white", "black")) +
                coord_cartesian(xlim = range(r$x), ylim = range(r$y), expand = c(0, 0)) +
                #coord_map(xlim = range(r$x), ylim = range(r$y)) +
                labs(color = NULL) +
                map_theme + theme(legend.spacing.y = unit(0, "cm"))
        
        bwc <- diff(range(r$wind_mean[r$land == 1])) / 30
        mc <- ggplot() +
                #facet_grid(. ~ "connectivity", switch = "y") +
                geom_raster(data = r %>% filter(land == 1) %>% mutate(z = distance/wind_mean,
                                                                      z = ifelse(z < .8, .8, z)), 
                            aes(x, y, fill = z)) +
                geom_contour(data = r, aes(x, y, z = distance, color = "distance"), binwidth = 200) +
                geom_contour(data = r, aes(x, y, z = wind_mean, color = "wind connectivity"), binwidth = bwc) +
                geom_raster(data = r %>% filter(land == 0), aes(x, y), fill = "white") +
                geom_point(data = pps, aes(lon_from, lat_from), color = "black", size = 2) +
                geom_point(data = filter(pps, id == i), aes(lon_from, lat_from), color = ptcol, size = 2) +
                scale_color_manual(values = c("black", "white")) +
                # scale_fill_gradientn(colors = c("black", "cyan"), breaks = c(1, 2)) +
                # scale_fill_viridis_c(trans = "log10") +
                scale_fill_gradientn(trans = "log10", breaks = 1:2, colors = pal) +
                guides(color = guide_legend(order = 1), 
                       fill = guide_colorbar(order = 0, barheight = 4.5)) +
                coord_cartesian(xlim = range(r$x), ylim = range(r$y), expand = c(0, 0)) +
                labs(fill = "normalized\nconnectivity\n(km/h)\n",
                     color = NULL) +
                map_theme + theme(legend.spacing.y = unit(0, "cm"))
        
        ms <- ggplot() +
                #facet_grid(. ~ "directionality", switch = "y") +
                geom_raster(data = r %>% filter(land == 1), aes(x, y, fill = wind_ratio)) +
                geom_contour(data = r, aes(x, y, z = log(wind_ratio)), color = "white", binwidth = .15) +
                geom_raster(data = r %>% filter(land == 0), aes(x, y), fill = "white") +
                geom_point(data = pps, aes(lon_from, lat_from), color = "black", size = 2) +
                geom_point(data = filter(pps, id == i), aes(lon_from, lat_from), color = ptcol, size = 2) +
                guides(fill = guide_colorbar(barheight = 9)) +
                # scale_fill_gradientn(colors = c("cyan", "black"), trans = "log10") +
                # scale_fill_viridis_c(trans = "log10") +
                scale_fill_gradientn(trans = "log10", colors = pal) +
                coord_cartesian(xlim = range(r$x), ylim = range(r$y), expand = c(0, 0)) +
                labs(fill = "outbound travel time\n/ inbound travel time\n") +
                map_theme + theme(legend.spacing.y = unit(0, "cm"))
        
        
        ### network cartoons ###
        
        library(ggraph)
        library(igraph)
        
        
        bi <- tibble(from = c(1, 2, 2, 3, 1, 3),
                     to = c(2, 1, 3, 2, 3, 1),
                     pair = c(12, 12, 23, 23, 13, 13),
                     flow = c(5, 3, 3, 2.6, 2, 1))
        uni <- bi %>%
                group_by(pair) %>%
                summarize(from = from[1], to = to[1],
                          diff = sum(flow),
                          dir = flow[1] / flow[2]) %>%
                select(from, to, pair, diff, dir)
        
        gb <- graph_from_data_frame(bi)
        gu <- graph_from_data_frame(uni)
        
        y <- sqrt(1^2 - .5^2)
        np <- data.frame(x = c(0, .5, 1),
                         y = c(0, y, 0))
        
        arw <- arrow(angle = 10, type = "closed", length = unit(8, 'mm'))
        pad <- .2
        xl <- xlim(-pad, 1 + pad)
        yl <- ylim(-pad*.5, y + pad*.5)
        
        center_label <- function(txt, y = .22, fontface = "bold", ...)  ggplot2::annotate(geom = "text", x = .5, y = y, label = txt, 
                                                                                          fontface = fontface, lineheight = .9, ...)
        
        
        edge_weights <- scale_edge_width_continuous(range = c(.25, 2))
        
        
        layout <- create_layout(gb, layout = 'kk')
        layout[,1:2] <- np
        nf <- ggraph(layout) + 
                geom_edge_arc(aes(edge_width = flow),
                              arrow = arw, end_cap = circle(10, 'mm'), curvature = .15) + 
                geom_node_point(size = 10, color = "black", shape = 21, fill = "white") + 
                coord_fixed() + xl + yl + 
                theme_void() + theme(legend.position = "none") +
                center_label("flow", size = 7) +
                edge_weights
        
        layout <- create_layout(gu, layout = 'kk')
        layout[,1:2] <- np
        nc <- ggraph(layout) + 
                geom_edge_arc(aes(edge_width = diff),
                              curvature = c(.15, -.15, .15)) + 
                geom_node_point(size = 10, color = "black", shape = 21, fill = "white") + 
                coord_fixed() + xl + yl + 
                theme_void() + theme(legend.position = "none") +
                center_label("isolation", size = 7) +
                edge_weights
        
        ns <- ggraph(layout) + 
                geom_edge_arc(aes(edge_width = dir),
                              arrow = arw, end_cap = circle(10, 'mm'), 
                              curvature = c(.15, -.15, .15)) + 
                geom_node_point(size = 10, color = "black", shape = 21, fill = "white") + 
                coord_fixed() + xl + yl + 
                theme_void() + theme(legend.position = "none") +
                center_label("asymmetry", size = 7) +
                edge_weights
        
        nd <- ggraph(layout) + 
                geom_edge_arc(aes(edge_width = dir),
                              arrow = arw, end_cap = circle(10, 'mm'), 
                              curvature = c(.15, -.15, .15)) + 
                geom_node_point(size = c(5, 7, 10), color = "black", shape = 21, fill = "black") + 
                coord_fixed() + xl + yl + 
                theme_void() + theme(legend.position = "none") +
                center_label("diversity", size = 7) +
                edge_weights
        
        
        
        ##### wind-genetic scatterplots ######
        
        rand <- list(h1_rand, h2_rand, h3_rand, div_rand) %>%
                map(function(x) map_df(x, "s13")) 
        for(i in 1:4) rand[[i]]$h <- c("flow", "symmetry", "connectivity", "diversity")[i] 
        rand <- bind_rows(rand)
        rand <- rand %>%
                filter(component == "17d", var == "wind", mig_var %in% c("sim", "mig", "div_ar"),
                       wind_p == "wind_p1", method == "pearson", dist_filter == "alldist")
        rp <- rand %>%
                group_by(h) %>%
                summarize(p = 1 - ecdf(spc[i != 0])(spc[i == 0]),
                          r = spc[i == 0]) %>%
                mutate(txt = paste0("r = ", round(r, 2), "\np = ", round(p, 3)))
        
        ####
        
        x <- annotated %>% mutate(mig = dmig) %>%
                filter(component == "17d",
                       wind_p == "wind_p1") %>%
                mutate(fst = ifelse(fst < 0, 0, fst),
                       mig = ifelse(mig < 0, 0, mig)) %>%
                group_by(component, edge) %>%
                arrange(from) %>%
                mutate(n = n()) %>%
                filter(n == 2) %>%
                mutate(mean_mig = mean(mig),
                       max_mig = max(mig),
                       mean_wind = mean(wind),
                       symmetry = c(log10(mig[1] / mig[2]),
                                    log10(mig[2] / mig[1])),
                       div = log10(div_ar_to / div_ar_from),
                       diff = c((mig[1] - mig[2]),
                                (mig[2] - mig[1])),
                       symmetry_wind = c(log10(wind[1] / wind[2]),
                                         log10(wind[2] / wind[1])),
                       diff_wind = c((wind[1] - wind[2]),
                                     (wind[2] - wind[1]))) %>%
                ungroup() %>%
                mutate(fst = 1 / fst,
                       symmetry = 10 ^ symmetry,
                       div = 10 ^ div,
                       symmetry_wind = 10 ^ symmetry_wind,
                       focal = from == i | to == i) %>%
                select(from, to, edge, focal,
                       flow_gen = mig, conn_gen = fst, symm_gen = symmetry, div_gen = div,
                       flow_wind = wind, conn_wind = mean_wind, symm_wind = symmetry_wind)
        
        gl <- geom_path(aes(group = edge), color = "gray80")
        gp <- function(...) geom_point(size = 1, ...)
        gpf <- geom_point(data = filter(x, focal), color = "red", size = 1)
        gs <- geom_smooth(method = lm, se = F, color = "black")
        gsf <- geom_smooth(data = filter(x, focal), method = lm, se = F, color = "red")
        lx <- scale_x_log10()
        ly <- scale_y_log10()
        
        
        
        #######
        
        th <- theme_bw(base_size = 16) +
                theme(#axis.title = element_blank(),
                        #axis.text = element_blank(),
                        #axis.ticks = element_blank(),
                        panel.grid = element_blank(),
                        strip.background = element_rect(fill = "black", color = "black"),
                        strip.text = element_text(color = "white"))
        
        annotate <- ggplot2::annotate
        z <- plot_spacer()
        gcol <- "orange"
        wcol <- "cadetblue"
        
        gcf <- x %>% ggplot(aes(conn_gen, flow_gen)) + gl + gp(color = gcol) + gpf + gs + th +
                scale_x_log10() +
                labs(x = "genetic similarity (1/Fst)", y = "gene flow")
        gcs <- x %>% ggplot(aes(conn_gen, symm_gen)) + gl + gp(color = gcol) + gpf + gs + th +
                scale_x_log10() +
                scale_y_log10() +
                labs(x = "genetic similarity (1/Fst)", y = "gene flow asymmetry ratio")
        gcd <- x %>% ggplot(aes(conn_gen, div_gen)) + gl + gp(color = gcol) + gpf + gs + th +
                scale_x_log10() +
                scale_y_log10(breaks = seq(0, 2, .1)) +
                labs(x = "genetic similarity (1/Fst)", y = "genetic diversity ratio")
        gsf <- x %>% ggplot(aes(symm_gen, flow_gen)) + gl + gp(color = gcol) + gpf + gs + th +
                scale_x_log10() +
                labs(x = "gene flow asymmetry ratio", y = "gene flow")
        gfd <- x %>% ggplot(aes(flow_gen, div_gen)) + gl + gp(color = gcol) + gpf + gs + th +
                scale_y_log10(breaks = seq(0, 2, .1)) +
                labs(y = "genetic diversity ratio", x = "gene flow")
        gsd <- x %>% ggplot(aes(symm_gen, div_gen)) + gl + gp(color = gcol) + gpf + gs + th +
                scale_x_log10() +
                labs(x = "gene flow asymmetry ratio", y = "genetic diversity ratio")
        
        
        wcf <- x %>% ggplot(aes(conn_wind, flow_wind)) + gl + gp(color = wcol) + gpf + gs + th +
                scale_x_log10() +
                scale_y_log10() +
                labs(x = "wind connectivity (1/h)", y = "wind flow (1/h)")
        wcs <- x %>% ggplot(aes(conn_wind, symm_wind)) + gl + gp(color = wcol) + gpf + gs + th +
                scale_x_log10() +
                scale_y_log10() +
                labs(x = "wind connectivity (1/h)", y = "wind flow asymmetry ratio")
        wsf <- x %>% ggplot(aes(symm_wind, flow_wind)) + gl + gp(color = wcol) + gpf + gs + th +
                scale_x_log10() +
                scale_y_log10() +
                labs(x = "wind flow asymmetry ratio", y = "wind flow (1/h)")
        
        plt <- wcf + wsf + wcs + 
                gcf + gsf + gcs +  
                gcd + gfd + gsd
        
        ggsave("figures/manuscript/figure_s1.pdf", plt, width = 12, height = 12, units = "in")
        
        
        
        #########
        
        th <- theme_bw(base_size = 16) +
                theme(panel.grid = element_blank(),
                      strip.background = element_rect(fill = "black", color = "black"),
                      strip.text = element_text(color = "white"))
        
        fwg <- x %>% ggplot(aes(flow_wind, flow_gen)) + gl + gp(color = "gray30") + gpf + gs + th +
                geom_text(data = filter(rp, h == "flow"), aes(x = max(x$flow_wind), y = min(x$flow_gen), label = txt), 
                          hjust = 1, vjust = 0, lineheight = .85) +
                scale_x_log10(breaks = c(.001, .01)) +
                scale_y_continuous(breaks = c(0, .5, 1)) +
                labs(x = "wind flow (1/h)", y = "gene flow")
        
        cwg <- x %>% ggplot(aes(conn_wind, conn_gen)) + gl + gp(color = "gray30") + gpf + gs + th +
                geom_text(data = filter(rp, h == "connectivity"), aes(x = max(x$conn_wind), y = min(x$conn_gen), label = txt), 
                          hjust = 1, vjust = 0, lineheight = .85) +
                scale_x_log10(breaks = c(.001, .01)) +
                scale_y_log10() +
                labs(x = "wind connectivity (1/h)", y = "genetic similarity")
        
        swg <- x %>% ggplot(aes(symm_wind, symm_gen)) + gl + gp(color = "gray30") + gpf + gs + th +
                geom_text(data = filter(rp, h == "symmetry"), aes(x = max(x$symm_wind), y = min(x$symm_gen), label = txt), 
                          hjust = 1, vjust = 0, lineheight = .85) +
                scale_x_log10(breaks = c(.33, 1, 3)) +
                scale_y_log10(breaks = c(.33, 1, 3)) +
                labs(x = "wind asymmetry ratio", y = "gene flow asymmetry ratio")
        
        dwg <- x %>% ggplot(aes(symm_wind, div_gen)) + gl + gp(color = "gray30") + gpf + gs + th +
                geom_text(data = filter(rp, h == "diversity"), aes(x = max(x$symm_wind), y = min(x$div_gen, na.rm = T), label = txt), 
                          hjust = 1, vjust = 0, lineheight = .85) +
                scale_x_log10(breaks = c(.33, 1, 3)) +
                scale_y_log10(breaks = c(.85, 1, 1.2)) +
                labs(x = "wind asymmetry ratio", y = "genetic diversity ratio")
        
        
        
        ### final composite plot ###
        
        lo <- 
                "ABBC
        DEEF
        GHHI
        JKKL"
        plt <- nf + mf + fwg +
                nc + mc + cwg +
                ns + ms + swg +
                nd + plot_spacer() + dwg +
                plot_layout(design = lo, heights = c(1,1,1,1)) & 
                theme(strip.text = element_text(size = 16))
        
        x <- ggs("figures/manuscript/figure_1.pdf", plt, width=14, height=14,
                 add = grid.text(LETTERS[1:4],
                                 x=.03,
                                 y=c(.96, .71, .46, .21),
                                 gp=gpar(fontsize=24, fontface="bold", col="black")))
        
        x <- ggs("figures/manuscript/figure_1.png", plt, 
                 width=14, height=14, units = "in", dpi = 1500,
                 add = grid.text(LETTERS[1:4],
                                 x=.03,
                                 y=c(.96, .71, .46, .21),
                                 gp=gpar(fontsize=24, fontface="bold", col="black")))
}

