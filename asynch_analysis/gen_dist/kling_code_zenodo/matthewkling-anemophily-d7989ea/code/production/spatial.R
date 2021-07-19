

spatialize <- function(pops){
        coordinates(pops) <- c("longitude", "latitude")
        crs(pops) <- crs('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
        pops
}


# merge all populations in the same grid cell
# (and remove genepop entries without individuals)
merge_pops <- function(d, grid){
        
        #d <- readd(datasets_raw)[1,]
        message(d$pop)
        # load datasets
        g <- g0 <- read_lines(d$gen)
        p <- read_csv(d$pop)
        if(is.character(p$longitude)) p$longitude <- p$longitude %>% enc2native() %>% as.numeric()
        
        # identify grid cell for each population
        ps <- spatialize(p)
        grid <- crop(grid, extent(ps)*2)
        grid[] <- 1:ncell(grid)
        p$cell <- extract(grid, ps)
        
        # remove individuals that have missing data for all loci
        geno <- g %>% str_sub(13)
        good <- str_detect(geno, paste(1:9, collapse = "|")) | nchar(geno) == 0
        g <- g[good]
        
        # remove empty populations from genepop
        gg <- g
        g <- gg[c(1, which(gg != lag(gg, 1)))]
        for(i in 1:nrow(p)) if(!any(grepl(paste0("p", p$ID[i], "_"), g))) p$ID[i] <- NA
        p <- filter(p, is.finite(ID))
        
        # merge populations
        gm <- g[1:2]
        for(cid in unique(p$cell)){
                gm <- c(gm, "pop")
                ids <- p$ID[p$cell == cid]
                #if(length(ids) > 1) stop()
                for(id in ids){
                        pattern <- paste0("p", id, "_i")
                        inds <- g[grepl(pattern, g)]
                        tag <- letters[match(id, ids) + 8]
                        replacement <- paste0("p", ids[1], "_", tag)
                        replacement <- str_pad(replacement, nchar(pattern), "right", tag)
                        inds <- sub(pattern, replacement, inds)
                        gm <- c(gm, inds)
                }
        }
        
        p$n <- NA
        for(i in 1:nrow(p)) p$n[i] <- sum(grepl(paste0("p", p$ID[i], "_"), gm))
        
        pm <- p %>%
                group_by(cell) %>%
                summarize(ID = ID[1],
                          name = name[1],
                          longitude = mean(longitude),
                          latitude = mean(latitude),
                          n = sum(n)) %>%
                ungroup()
        
        
        g <- gm
        p <- pm
        
        p <- p %>% select(-cell)
        
        # eliminate pops with one individual (they break downstream code)
        if(min(p$n) < 2){
                x <- p$ID[p$n < 2] # pops to cut
                xi <- which(grepl(paste(paste0("p", x, "_"), collapse = "|"), g))
                xi <- unique(c(xi, xi - 1))
                g <- g[setdiff(1:length(g), xi)]
                p <- filter(p, n > 1)
        }
        
        # save results
        d$gen <- sub("\\.txt", "_merged.txt", d$gen)
        sink(d$gen)
        writeLines(g)
        sink()
        
        d$pop <- sub("\\.csv", "_merged.csv", d$pop)
        write_csv(p, d$pop)
        
        return(d)
}

