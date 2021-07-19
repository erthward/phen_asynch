
# code for appendix 3: 
# analyses of spatial structure in wind-genetic correlations

autocorelation <- function(){
        
        library(tidyverse)
        library(geosphere)
        library(vegan)
        
        select <- dplyr::select
        
        
        #loadd(h4_rand, h1_rand, h2_rand, h3_rand, annotated, nreps, pgs)
        
        
        # s13
        rand <- list(h1_rand, h2_rand, h3_rand, h4_rand) %>%
                map(function(x) map_df(x, "s13")) 
        for(i in 1:4) rand[[i]]$h <- c("flow", "asymmetry", "isolation", "diversity")[i] 
        rand <- bind_rows(rand)
        rand <- filter(rand, wind_p == "wind_p1")
        
        specs <- rand %>% select(h, mig_var) %>% distinct()
        specs <- expand_grid(flow = specs$mig_var[specs$h == "flow"],
                             asymmetry = specs$mig_var[specs$h == "asymmetry"],
                             isolation = specs$mig_var[specs$h == "isolation"],
                             diversity = specs$mig_var[specs$h == "diversity"],
                             method = unique(rand$method),
                             dist_filter = unique(rand$dist_filter))
        
        spec <- specs %>%
                filter(flow == "mig",
                       asymmetry == "mig",
                       isolation == "sim",
                       diversity == "div_ar",
                       method == "pearson",
                       dist_filter == "alldist") %>%
                gather(h, mig_var, flow:diversity)
        
        sum_fun <- median
        
        s <- left_join(spec, rand) %>% 
                filter(var %in% c("wind", "swind")) %>%
                mutate(h = factor(h, levels = c("flow", "isolation", "asymmetry", "diversity")),
                       syndrome = expected_wind) %>%
                filter(!is.na(component)) %>%
                group_by(h, syndrome, component) %>%
                mutate(n = n()) %>% filter(n == nreps+1)
        
        spr <- s %>% 
                ungroup() %>%
                group_by(h, component, syndrome) %>%
                summarize(p = 1 - ecdf(spc[i != 0])(spc[i == 0]),
                          spc = spc[i == 0]) %>%
                ungroup() %>%
                gather(stat, value, spc, p) %>%
                mutate(x = case_when(syndrome == "0" ~ 0,
                                     syndrome == "0.5" ~ .5,
                                     syndrome == "1" ~ 1,
                                     syndrome == "1~0" ~ 1.5))
        
        centroids <- annotated %>%
                group_by(component) %>%
                summarize(lon = mean(lon_from, na.rm = T),
                          lat = mean(lat_from, na.rm = T)) %>%
                arrange(component)
        
        fx <- function(r, g){
                r <- filter(r, is.finite(value))
                g <- filter(g, component %in% r$component)
                
                # mantel test
                rd <- dist(r$value)
                gd <- distm(g[,2:3]) %>% as.dist()
                m <- vegan::mantel(rd, gd)
                m <- tibble(stat = r$group[1], test = "mantel", 
                            r = m$statistic, r2 = r^2, p = m$signif)
                
                # moran's i test
                gd <- 1 / as.matrix(gd)
                diag(gd) <- 0
                gd[is.infinite(gd)] <- 0
                i <- ape::Moran.I(r$value, gd)
                i <- tibble(stat = r$group[1], test = "moran",
                            r = i$observed, r2 = r^2, p = i$p.value)
                
                bind_rows(m, i)
        }
        
        corrs <- spr %>%
                arrange(syndrome, h, stat, component) %>%
                mutate(group = paste(syndrome, h, stat)) %>%
                split(.$group)
        
        m <- map_df(corrs, fx, g = centroids) %>%
                separate(stat, c("syndrome", "h", "stat"), sep = " ")
        
        m <- m %>% mutate(r = signif(r, 3),
                          r2 = signif(r2, 3),
                          p = signif(p, 3))
        
        write_csv(m, "data_github/autocorrelation_tests.csv")
        
        
        #####################
        
        # after accounting for space using partial Mantel, 
        # does wind-genetic correlation still increase as expected?
        
        fx <- function(r, g){
                r <- filter(r, is.finite(value))
                g <- filter(g, component %in% r$component)
                
                # first-order and partial mantel tests
                rd <- dist(r$value)
                sd <- dist(r$syndrome)
                gd <- distm(g[,2:3]) %>% as.dist()
                m <- vegan::mantel(sd, rd)
                m <- tibble(stat = r$group[1], test = "mantel", 
                            r = m$statistic, r2 = r^2, p = m$signif)
                mp <- vegan::mantel.partial(sd, rd, gd)
                mp <- tibble(stat = r$group[1], test = "partial mantel", 
                             r = mp$statistic, r2 = r^2, p = mp$signif)
                
                bind_rows(m, mp)
        }
        
        corrs <- spr %>%
                arrange(h, stat, component) %>%
                mutate(group = paste(h, stat)) %>%
                split(.$group)
        
        m <- map_df(corrs, fx, g = centroids) %>%
                separate(stat, c("h", "stat"), sep = " ")
        
        m %>% select(h, stat, test, r) %>%
                mutate(stat = recode(stat, "spc" = "r"),
                       test = recode(test, "partial mantel" = "partial"),
                       r = signif(r, 3)) %>%
                spread(test, r) %>%
                mutate(difference = partial - mantel) %>%
                write_csv("data_github/autocorrelation_controls.csv")
}

