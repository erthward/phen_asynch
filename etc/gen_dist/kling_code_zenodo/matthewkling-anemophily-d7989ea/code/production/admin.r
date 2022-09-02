

get_metadata <- function(){
        
        drive_auth(email = "mattkling@berkeley.edu")
        sheets_auth(email = "mattkling@berkeley.edu")
        
        folders <- drive_ls("Plant PopGen Project/formatted_datasets")
        
        namefix <- function(x) x %>% tolower() %>% 
                str_replace_all(" |\\?|-", "_") %>%
                str_replace_all("_\\(apa\\)", "") %>%
                str_replace_all("use_", "use")
        
        # google spreadsheet with metadata
        url <- "https://docs.google.com/spreadsheets/d/15B80Kg9pAFaukfDy-oW8MOHefyyKANIG1VvxwyQDaG4/edit#gid=1398094914"
        studies <- read_sheet(url,  sheet = "reformatting") %>%
                rename_all(namefix) %>%
                mutate(rownum = (1:nrow(.)),
                       component = paste0(dataset_id, component_id)) %>%
                select(dataset_id:sequence_type, use, rownum, component, valid_alignment, refresh)
        
        queries <- read_sheet(url, sheet = "results") %>%
                rename_all(namefix) %>% 
                mutate(citation = paste0(str_split(journal_article_citation, ",", n=2, simplify=T)[,1], 
                                         "_", publication_year)) %>%
                select(dataset_id, citation, location)
        queries <- queries %>%
                mutate(dataset_id = dataset_id + 1000) %>%
                bind_rows(queries)
        
        left_join(studies, queries) %>%
                filter(!is.na(use), use=="yes" #, valid_alignment=="pass"
                ) %>%
                mutate(name = ifelse(component %in% folders$name, 
                                     component, dataset_id)) %>%
                left_join(folders) %>% 
                split(f = .$component)
}



import_dataset <- function(md, outdir="data/datasets"){
        
        drive_auth(email = "mattkling@berkeley.edu")
        sheets_auth(email = "mattkling@berkeley.edu")
        
        md <- md[[1]]
        ci <- paste0(md$dataset_id, md$component_id)
        
        dr <- paste0(outdir, "/", ci)
        md$path <- dr
        
        outfiles <- paste0(dr, "/", c("genepop.txt", "populations.csv"))
        names(outfiles) <- c("gen", "pop")
        out <- tibble(component = md$component, gen = outfiles[1], pop = outfiles[2])
        
        if(all(file.exists(c(out$gen, out$pop))) &
           (is.na(md$refresh) | md$refresh != "yes" )) return(out)
        
        message(paste("importing", ci))
        if(dir.exists(dr)) unlink(dr, recursive=T)
        suppressWarnings(dir.create(dr))
        
        f <- drive_ls(as_id(md$id)) %>%
                filter(name %in% c("populations.csv", paste0("populations_", ci, ".csv"), 
                                   "genepop.txt", paste0("genepop_", ci, ".txt"))) %>%
                arrange(name)
        
        for(i in 1:nrow(f)) suppressMessages(
                drive_download(as_id(f$id[i]), 
                               paste0(dr, "/", sub(paste0("_", ci), "", f$name[i])),
                               overwrite=T))
        return(out)
}


# melt a distance matrix
flatten <- function(x, i, varname="dist"){
        # x = distance matrix
        # i = vector of site numbers
        x <- as.data.frame(x)
        names(x) <- paste0("p", i)
        x$from <- names(x)
        x <- gather(x, to, dist, -from)
        x$from <- as.integer(sub("p", "", x$from))
        x$to <- as.integer(sub("p", "", x$to))
        names(x)[3] <- varname
        x
}


# combine variables into a single data frame, for one dataset
collate <- function(datasets, pops, div, pw){
        
        # unnest and flatten 
        pw2 <- list()
        for(i in names(pw)){
                pwi <- pw[[i]]
                
                if(class(pwi) == "list"){
                        names(pwi) <- paste0(i, "_", names(pwi))
                        pw2 <- c(pw2, pwi)
                }else{
                        pw2[[i]] <- pwi
                }
        }
        if(length(unique(as.vector(sapply(pw2, dim)))) > 1) return(NULL)
        pw <- map(names(pw2), function(n) flatten(pw2[[n]], pops$ID, n))
        
        # diversity
        dvt <- dvf <- div[c("Allelic_richness", "Allele_number")] %>% 
                map(function(x) x["overall",]) %>% 
                map(as.data.frame) %>%
                map(rownames_to_column) %>%
                bind_cols() %>%
                select(1, 2, 4) %>%
                separate(rowname, c("pop", "junk"), sep = "_") %>%
                select(-junk) %>%
                mutate(pop = as.integer(str_remove(pop, "p")))
        names(dvf) <- c("from", "div_ar_from", "div_an_from")
        names(dvt) <- c("to", "div_ar_to", "div_an_to")
        
        # lat-long data for each edge
        lon_from <- matrix(rep(pops$longitude, each=nrow(pops)), nrow=nrow(pops), byrow=T) %>%
                flatten(pops$ID, "lon_from")
        lon_to <- matrix(rep(pops$longitude, each=nrow(pops)), nrow=nrow(pops), byrow=F) %>%
                flatten(pops$ID, "lon_to")
        lat_from <- matrix(rep(pops$latitude, each=nrow(pops)), nrow=nrow(pops), byrow=T) %>%
                flatten(pops$ID, "lat_from")
        lat_to <- matrix(rep(pops$latitude, each=nrow(pops)), nrow=nrow(pops), byrow=F) %>%
                flatten(pops$ID, "lat_to")
        
        # combine
        Reduce("full_join", pw) %>%
                full_join(lat_from) %>% full_join(lat_to) %>% full_join(lon_from) %>% full_join(lon_to) %>% 
                full_join(dvf) %>% full_join(dvt) %>%
                mutate(component = datasets$component) %>%
                mutate(edge = paste(pmin(from, to), pmax(from, to), sep="_"))
        
}


# identify pipeline steps at which datasets are lost due to errors
track_attrition <- function(datasets_validated, datasets, collated, unified){
        
        cd <- collated %>%
                filter(is.finite(distance + clim + mig + nmigp + wind + fst)) %>%
                count(component) %>%
                mutate(pairwise = is.finite(n)) %>%
                select(-n)
        
        ud <- unified %>%
                filter(is.finite(distance + clim + mig + nmigp + wind + fst)) %>%
                select(component, dispersal_months:expected_wind) %>%
                distinct() %>%
                mutate_at(vars(dispersal_months:expected_wind), function(x) !is.na(x))
        
        d <- datasets_validated %>%
                select(component, valid) %>%
                left_join(cd) %>%
                left_join(ud) %>%
                mutate_all(function(x) ifelse(is.na(x), F, x)) %>%
                
                select(-pollination_months, -dispersal_months) %>%
                mutate(all = valid & pairwise & expected_wind)
        
        return(d)
        
}


# build dictionary of order, family, genus for every species
get_taxonomy <- function(d){
        classs <- function(genus){
                #require(taxize)
                Sys.sleep(1)
                cls <- get_ids(genus, db = "ncbi", ask = F, messges = F)
                cls <- classification(cls$ncbi, db = "ncbi")[[1]]
                cls %>%
                        filter(rank %in% c("order", "family")) %>%
                        dplyr::select(name, rank) %>%
                        spread(rank, name) %>%
                        mutate(genus = genus)
        }
        #browser()
        d <- mutate(d, genus = as.vector(str_split_fixed(species, pattern = " ", n = 2)[,1]))
        tax <- unique(d$genus) %>% map_df(possibly(classs, NULL))
        d %>% select(species, genus) %>% distinct() %>% left_join(tax)
}


# semipartial correlation
fspcor <- function(x, y, get, ..., method = "spearman"){
        require(ppcor)
        m <- cbind(x, y, ...) %>% na.omit()
        r <- try(spcor(m, method=method))
        if(class(r) == "try-error") return(NA)
        if(get=="est") return(r$estimate[1,2])
        if(get=="sig") return(r$p.value[1,2])
}



annotate <- function(d, taxonomy){
        
        d %>% 
                mutate(syndquant = factor(expected_wind),
                       sperm = case_when(plastid_inheritance == "seed" ~ "angiosperm",
                                         plastid_inheritance == "pollen" ~ "gymnosperm",
                                         TRUE ~ NA_character_),
                       distbin = ifelse(distance > median(distance), "long", "short")) %>%
                
                left_join(taxonomy)
}

