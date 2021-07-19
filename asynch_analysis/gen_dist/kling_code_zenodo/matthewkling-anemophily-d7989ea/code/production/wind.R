

unwrap <- function(x, width=20){
        x <- x %>% crop(extent(180-width, 180, -90, 90)) %>% shift(-360) %>% merge(x)
        x <- x %>% crop(extent(-180, -180+width, -90, 90)) %>% shift(360) %>% merge(x)
        x
}


add_coords <- function(windrose){
        rows <- cols <- windrose[[1]]
        rows[] <- rep(1:nrow(rows), each=ncol(rows))
        cols[] <- rep(1:ncol(rows), nrow(rows))
        windrose <- stack(windrose, rows, cols)
        names(windrose) <- c("SW", "W", "NW", "N", "NE", "E", "SE", "S", "row", "col")
        return(windrose)
}


compile_windrose <- function(months = 1:12, outfile,
                             years = 1980:2009, p = 1, overwrite = F,
                             fd = "E:/wind/winds_of_change/tailwinds/data/windrose/monthly_p1_wnd10m"){
        
        if(file.exists(outfile) & !overwrite) return("skipping")
        
        f <- list.files(fd, full.names=T)
        yr <- substr(f, nchar(f)-11, nchar(f)-8) %>% as.integer()
        mo <- substr(f, nchar(f)-7, nchar(f)-6) %>% as.integer()
        f <- f[yr %in% years & mo %in% months]
        if(length(f) %% (length(years) * length(months)) != 0) stop("non-factorial data")
        
        #p <- ifelse(p==0, 1, p)
        y <- f %>% 
                lapply(stack) %>%
                Reduce("+", .) %>%
                "/"(24 * 365/12 * length(months) * length(years)) %>% # number of hours
                #"^"(1/p) %>%
                writeRaster(outfile, overwrite=T)
        return("completed")
}


make_windroses <- function(traits){
        
        m <- traits %>% 
                select(months) %>% distinct() %>% na.omit() %>%
                filter(months != "NA")
        m <- expand_grid(months = m$months,
                         p = 0:3) %>%
                mutate(windrose = paste0("data/wind/roses/p", p, "/", months, ".tif"),
                       vec = str_split(months, "_"),
                       vec = sapply(vec, as.integer),
                       indir = paste0("E:/wind/winds_of_change/tailwinds/data/windrose/monthly_p", p, "_wnd10m"))
        
        y <- list(months = m$vec, outfile = m$windrose, p = m$p, fd = m$indir) %>%
                pmap(compile_windrose)
        
        return(m$windrose)
}




make_wind_graph <- function(infile, outfile = NULL,
                            water = .1, # accessibility multiplier
                            direction = "downwind",
                            logtrans = FALSE # make cost-paths products by log-transforming transitions
){
        
        # create raster of weights for land/water
        weights <- raster("f:/cfsr/land.tif") %>%  
                rotate() %>% 
                unwrap(180) %>%
                reclassify(c(NA, NA, water))
        
        # load windrose data and apply weights
        rose <- stack(infile) %>% 
                rotate() %>% 
                unwrap(width = 180) %>%
                "*"(weights)
        
        # apply log transformation if requested
        # (not currently working -- may need to modify geometery in transition_stack() as well)
        # if(logtrans){
        #         v <- na.omit(values(rose))
        #         mv <- min(v[v > 0])
        #         rose[rose == 0] <- mv
        #         rose <- rose / mv
        #         rose <- calc(rose, fun = log10)
        # }
        
        # make transition graph
        wind <- rose %>% 
                add_coords() %>% 
                transition_stack(windflow, directions = 8, 
                                 symm = F, direction = direction)
        
        if(is.null(outfile)){
                outfile <- infile %>%
                        str_replace("roses", "graphs") %>%
                        str_replace("\\.tif", ".rds")
        }
        
        saveRDS(wind, outfile)
        return(outfile)
}


# pairwise wind hours
wind_times <- function(pops, wind_graphs, traits, annual = FALSE){
        
        mos <- traits$months[traits$component == pops$component[1]]
        if(annual | is.na(mos) | mos == "NA") mos <- "1_2_3_4_5_6_7_8_9_10_11_12"
        
        files <- wind_graphs[grepl(paste0("/", mos, ".rds"), wind_graphs)]
        
        graphs <-  files %>% map(readRDS)
        
        pops <- spTransform(pops, crs(graphs[[1]]))
        
        cd <- graphs %>%
                map(costDistance, fromCoords = pops) %>%
                map(function(x) x / 3600) # seconds to hours
        
        names(cd) <- basename(dirname(files))
        if(annual) names(cd) <- paste0("ann", names(cd))
        
        return(cd)
}
