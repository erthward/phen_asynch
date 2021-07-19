


intergenome_data <- function(d, covariates){
        
        # identify relevant datasets
        igd <- d %>% 
                filter(wind_pollination + wind_dispersal > 0) %>%
                select(species, dataset_id, genome) %>%
                distinct() %>%
                count(species, dataset_id) %>%
                filter(n == 2) %>%
                left_join(d) 
        
        if(any(!grepl("Quercus", igd$species))) stop("this function assumes Quercus biology")
        
        # restructure
        igd <- igd %>%
                select(species, dataset_id, genome, 
                       edge, from, to, lat_from, lon_from, lat_to, lon_to,
                       y, covariates) %>%
                arrange(dataset_id, edge, from, to) %>%
                group_by(dataset_id, edge, from, to) %>% 
                mutate(n = n()) %>%
                filter(n == 2)
        
        # only the wind from the nuclear genome is relevant, for oaks 
        if("wind" %in% covariates) igd <- igd %>% mutate(wind = wind[genome == "nu"])
        if("swind" %in% covariates) igd <- igd %>% mutate(swind = swind[genome == "nu"]) 
        
        igd %>% 
                spread(genome, y) %>%
                filter(is.finite(nu), is.finite(cp))
}

remix <- function(x, i = NULL){
        y <- x %>% 
                mutate(from = as.character(from)) %>%
                gather(to, value, -from)
        
        pops <- unique(c(y$to, y$from))
        
        e <- expand_grid(to = pops, from = pops) %>%
                left_join(y) %>%
                spread(to, value)
        
        if(is.null(i)) i <- sample(length(pops))
        
        from <- select(e, from)
        e <- select(e, -from)
        z <- e[i, i]
        colnames(z) <- colnames(e)
        
        bind_cols(from, z)
}

mantel <- function(x, nreps, method, covariates){
        require(ppcor)
        message(x$component[1])
        pw <- x %>%
                select(to, from, y) %>% 
                spread(to, y)
        
        y <- tibble(component = x$component[1],
                    i = 0:nreps)
        for(v in covariates) y[,v] <- NA
        
        for(i in 0:nreps){
                set.seed(i)
                
                if(i == 0){pwi <- pw} else{
                        pwi <- try(remix(pw))
                        if(class(pwi) == "try-error") return(y)
                }
                
                f <- pwi %>%
                        gather(to, y, -from) %>% 
                        mutate(to = as.integer(to),
                               from = as.integer(from)) %>%
                        left_join(select(x, -y), ., by = c("to", "from"))
                
                s <- f %>% select(c("y", covariates)) %>% na.omit()
                s <- try(pcor(s, method = method))
                if(class(s) == "try-error") next()
                s <- s$estimate[2:(length(covariates)+1), 1]
                
                y[i + 1, 3:(2+length(s))] <- s
        }
        
        y
}

mantel_intergenome <- function(x, nreps, method, covariates){
        require(ppcor)
        
        pw_nu <- x %>%
                select(to, from, nu) %>%
                spread(to, nu)
        pw_cp <- x %>%
                select(to, from, cp) %>%
                spread(to, cp)
        
        y <- tibble(dataset_id = x$dataset_id[1],
                    i = 0:nreps)
        for(v in covariates) y[,v] <- NA
        
        npops <- length(unique(c(x$to, x$from)))
        
        for(i in 0:nreps){
                set.seed(i)
                
                if(i == 0){
                        nui <- pw_nu
                        cpi <- pw_cp
                }else{
                        # keep nu and cp together during resampling
                        ii = sample(npops)
                        nui <- remix(pw_nu, ii)
                        cpi <- remix(pw_cp, ii)
                }
                
                
                nuf <- nui %>%
                        gather(to, nu, -from) %>% 
                        mutate(to = as.integer(to), from = as.integer(from))
                cpf <- cpi %>%
                        gather(to, cp, -from) %>% 
                        mutate(to = as.integer(to), from = as.integer(from))
                
                f <- left_join(nuf, cpf, by = c("to", "from")) %>%
                        left_join(select(x, -nu, -cp), ., by = c("to", "from"))
                
                s <- f %>% select(c("nu", covariates)) %>% na.omit()
                s <- try(pcor(s, method = method))
                if(class(s) == "try-error") next()
                s <- s$estimate[2:(length(covariates)+1), 1]
                
                y[i + 1, 3:(2+length(s))] <- s
        }
        
        y
}

perms <- function(mig_var = "mig", method = "spearman", dist_filter = "alldist", 
                  p = NULL, nreps = 1000, covariates, d){
        
        # test 1: sign test
        # test 2: functional group test
        # test 3: intergenome test
        
        
        ### admin and randomizations ###
        
        tag <- paste0("_", mig_var, "_", method, "_", dist_filter, "_", p)
        message(tag)
        
        future::plan(multiprocess)
        
        d$y <- d[[mig_var]]
        
        if(!is.null(p)) d <- filter(d, wind_p == p)
        
        if(dist_filter == "shortdist"){
                d <- d %>% group_by(component) %>% 
                        mutate(y = ifelse(distance > median(distance, na.rm = T), NA, y)) %>%
                        ungroup()
        }
        
        
        ### tests 1&2 ###
        
        u <- d %>% select(c("component", "to", "from", "y", covariates)) %>%
                split(.$component)
        
        u <- u[which(sapply(u, nrow) > 2)]
        
        s <- u %>% future_map_dfr(mantel, covariates = covariates, nreps = nreps, method = method, 
                                  .progress = TRUE)
        
        s <- d %>%
                select(component, expected_wind) %>%
                distinct() %>%
                right_join(s) %>%
                gather(var, spc, covariates) %>%
                filter(abs(spc) <= 1, 
                       is.finite(spc),
                       is.finite(expected_wind),
                       expected_wind %in% c(0, .5, 1))
        
        s13 <- s %>% mutate(mig_var = mig_var, 
                            method = method, 
                            dist_filter = dist_filter,
                            wind_p = p,
                            nreps = nreps)
        
        
        
        ### test 3: inter-genome ###
        
        u <- d %>%
                filter(is.finite(y)) %>%
                intergenome_data(covariates) %>%
                ungroup() %>%
                mutate(dataset_id = paste(dataset_id, species))
        
        u <- split(u, u$dataset_id)
        
        s <- u %>% future_map_dfr(mantel_intergenome, 
                                  covariates = c(covariates, "cp"),
                                  nreps = nreps, method = method, .progress = TRUE)
        
        s2 <- s %>% mutate(mig_var = mig_var, 
                           method = method, 
                           dist_filter = dist_filter, 
                           wind_p = p,
                           nreps = nreps)
        
        
        ### return randomization data ###
        list(s13 = s13, 
             s2 = s2)
}







######### flow ##########
data_h1 <- function(annotated, stat = "dmig"){
        scale <- function(x) as.vector(base::scale(x))
        d <- annotated %>% mutate(mig = .[[stat]]) %>%
                filter(to != from,
                       is.finite(log(mig)),
                       is.finite(log(wind))) %>%
                group_by(component, wind_p) %>%
                mutate(mig = scale(log(mig))) %>%
                mutate(nmigp = nmigp * .998 + .01,
                       nmigp = log(nmigp / (1 - nmigp)),
                       div_ar = log10(div_ar_to / div_ar_from)) %>%
                ungroup() %>%
                
                mutate(wind = scale(log(wind)),
                       clim = scale(log(clim)),
                       distance = scale(log(distance))) %>%
                filter(is.finite(clim),
                       is.finite(distance),
                       is.finite(wind))
        return(d)
}
mantel_h1 <- function(d, nreps = 1000, covariates = c("wind", "distance", "clim")){
        specs <- expand_grid(mig_var = c("mig", "nmigp")[1],
                             method = c("spearman", "pearson")[2],
                             dist_filter = c("alldist", "shortdist")[1],
                             p = "wind_p1")
        pmap(specs, perms, d = d, covariates = covariates, nreps = nreps)
}



############ asymmetry ###########
data_h2 <- function(annotated, stat = "dmig"){
        ratios <- function(x) c(log10(x[1] / x[2]), log10(x[2] / x[1]))
        d <- annotated %>% mutate(mig = .[[stat]]) %>%
                filter(to != from, is.finite(log(wind))) %>%
                arrange(component, edge, from, to) %>%
                group_by(component, edge, wind_p) %>%
                mutate(wind = ratios(wind),
                       mig = ratios(mig),
                       lat_diff = lat_from - lat_to,
                       div_ar = log10(div_ar_to / div_ar_from),
                       swind = sign(wind),
                       smig = sign(mig),
                       aligned = swind == smig) %>%
                ungroup()
        return(d)
}
mantel_h2 <- function(d, nreps = 1000, covariates = c("wind")){
        specs <- expand_grid(mig_var = c("mig", "smig")[1],
                             method = c("spearman", "pearson")[2],
                             dist_filter = c("alldist", "shortdist")[1],
                             p = "wind_p1")
        pmap(specs, perms, d = d, covariates = covariates, nreps = nreps)
}



######## isolation ##########
data_h3 <- function(annotated, stat = "djost"){
        scale <- function(x) as.vector(base::scale(x))
        d <- annotated %>% mutate(diff = .[[stat]]) %>%
                mutate(diff = ifelse(diff < 0, 0, diff)) %>%
                filter(to != from, is.finite(wind), #is.finite(diff), 
                       is.finite(distance)) %>%
                arrange(component, edge, from, to) %>%
                group_by(component, edge, wind_p) %>%
                mutate(diff = mean(diff, na.rm = T),
                       sim = 1 / diff,
                       distance = mean(distance),
                       clim = mean(clim),
                       wind = mean(wind),
                       div_ar = log10(div_ar_to / div_ar_from),
                       lon = mean(c(lon_from, lon_to)),
                       lat = mean(c(lat_from, lat_to))) %>%
                ungroup()  %>%
                filter(is.finite(sim)) %>%
                mutate(sim = scale(log(sim)),
                       wind = scale(log(wind)),
                       clim = scale(log(clim)),
                       distance = scale(log(distance))) %>%
                filter(is.finite(clim),
                       is.finite(distance),
                       is.finite(sim),
                       is.finite(wind))
        return(d)
}
mantel_h3 <- function(d, nreps = 1000, covariates = c("wind", "distance", "clim")){
        specs <- expand_grid(mig_var = "sim",
                             method = c("spearman", "pearson")[2],
                             dist_filter = c("alldist", "shortdist")[1],
                             p = "wind_p1")
        pmap(specs, perms, d = d, 
             covariates = covariates, nreps = nreps)
}



########## diversity ##############
data_h4 <- function(annotated, stat = "dmig"){
        ratios <- function(x) c(log10(x[1] / x[2]), log10(x[2] / x[1]))
        d <- annotated %>% mutate(mig = .[[stat]]) %>%
                filter(to != from, is.finite(log(wind))) %>%
                arrange(component, edge, from, to) %>%
                group_by(component, edge, wind_p) %>%
                mutate(wind = ratios(wind),
                       mig = ratios(mig),
                       lat_diff = lat_from - lat_to,
                       div_ar = log10(div_ar_to / div_ar_from), # ordered to make predicted corr positive
                       div_an = log10(div_an_to / div_an_from)) %>%
                ungroup()
        return(d)
}
mantel_h4 <- function(d, nreps = 1000, covariates = c("wind", "lat_diff")){
        specs <- expand_grid(mig_var = c("div_ar", "div_an")[1],
                             method = c("pearson"),
                             dist_filter = c("alldist"),
                             p = "wind_p1")
        pmap(specs, perms, d = d, 
             covariates = covariates, nreps = nreps)
}


