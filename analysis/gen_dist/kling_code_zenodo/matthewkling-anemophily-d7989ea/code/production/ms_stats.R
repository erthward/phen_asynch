

################ metadata stats ##################

pg_stats <- function(infile){
        # infile = datasets$gen[1]
        dat <- haploRGP(infile)
        tibble(component = basename(dirname(infile)),
               npops = length(dat$genos),
               nind = dat$genos %>% map(dim) %>% map_int(1) %>% mean(),
               nloci = length(dat$af))
}


export_global_summaries <- function(annotated, pgs){
        
        a <- annotated %>%
                filter(wind_p == "wind_p1", 
                       is.finite(wind), is.finite(djost)) %>%
                select(component, sequence_type, expected_wind) %>% distinct()
        
        e <- read_csv("data_github/mantel_results.csv") %>%
                filter(syndrome != "0.5~0") %>%
                left_join(a) %>%
                left_join(pgs %>% rename(npops2 = npops))
        
        # number of datasets
        gs <- tibble(stat = "datsets", n = nrow(e))
        
        # number of species
        gs <- rbind(gs, c("species", length(unique(e$species))))
        
        # number of populations
        gs <- rbind(gs, c("populations", sum(e$npops, na.rm = T)))
        gs <- rbind(gs, c("median populations per dataset", median(e$npops, na.rm = T)))
        gs <- rbind(gs, c("mean populations per dataset", mean(e$npops, na.rm = T)))
        gs <- rbind(gs, c("stdev populations per dataset", sd(e$npops, na.rm = T)))
        gs <- rbind(gs, c("q25 populations per dataset", quantile(e$npops, .25, na.rm = T)))
        gs <- rbind(gs, c("q75 populations per dataset", quantile(e$npops, .75, na.rm = T)))
        gs <- rbind(gs, c("population pairs", sum(e$npairs, na.rm = T)))
        gs <- rbind(gs, c("median individuals per population", median(e$nind, na.rm = T)))
        gs <- rbind(gs, c("q25 individuals per population", quantile(e$nind, .25, na.rm = T)))
        gs <- rbind(gs, c("q75 individuals per population", quantile(e$nind, .75, na.rm = T)))
        gs <- rbind(gs, c("median loci: SSR datasets", median(e$nloci[e$sequence_type=="SSR"], na.rm = T)))
        gs <- rbind(gs, c("median loci: SNP datasets", median(e$nloci[e$sequence_type=="SNP"], na.rm = T)))
        
        # breakdowns by data category
        gs <- count(e, expected_wind) %>% rename(stat = expected_wind) %>% 
                mutate(stat = paste("datasets with expected wind =", stat),
                       n = as.character(n)) %>% bind_rows(gs, .)
        gs <- count(e, sequence_type) %>% rename(stat = sequence_type) %>% 
                mutate(stat = paste("datasets with sequence type =", stat),
                       n = as.character(n)) %>% bind_rows(gs, .)
        gs <- count(e, genome) %>% rename(stat = genome) %>% 
                mutate(stat = paste("datasets with genome type =", stat),
                       n = as.character(n)) %>% bind_rows(gs, .)
        
        prop_pos <- e %>%
                select(component, syndrome, asymmetry_p:isolation_r) %>%
                gather(stat, value, -component, -syndrome) %>%
                separate(stat, c("h", "stat")) %>%
                spread(stat, value) %>%
                na.omit() %>%
                group_by(syndrome, h) %>%
                summarize(prop_p = mean(p < .5),
                          prop_r = mean(r > 0))
        prop_pos %>% filter(syndrome == "1") %>% pull(prop_p) %>% mean()
        
        write_csv(gs, "data_github/global_summaries.csv")
}



################ data on every species/publication ###############

export_study_metadata <- function(annotated){
        
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
                select(dataset_id, dataset_url, citation, journal_article_citation)
        queries <- queries %>%
                mutate(dataset_id = dataset_id + 1000) %>%
                bind_rows(queries)
        
        d <- annotated %>%
                select(dataset_id, component, species, genome, sequence_type) %>%
                distinct() %>%
                left_join(queries, .) %>%
                filter(!is.na(component)) %>%
                rename(reference = journal_article_citation) %>%
                group_by(reference) %>%
                select(-dataset_id, -citation, -component) %>%
                arrange(reference) %>%
                select(reference, dryad_url = dataset_url, species, genome, sequence_type) %>%
                group_by(reference, dryad_url) %>%
                summarize_all(function(x) paste(unique(x), collapse = ", "))
        
        write_csv(d, "data_github/dataset_references.csv")
        return(d)
}



export_trait_data <- function(annotated){
        
        url <- "https://docs.google.com/spreadsheets/d/1onSN4dJ0_4A6kFQvL0sR5yWJ3RrNjcso8zaEby3O45Y/edit#gid=1840736075"
        traits <- read_sheet(url,  sheet = "traits", col_types = "c")
        
        d <- traits %>%
                filter(species %in% unique(annotated$species)) %>%
                select(-notes)
        
        write_csv(d, "data_github/species_traits.csv")
        return(d)
}

export_search_genera <- function(annotated){
        
        url <- "https://docs.google.com/spreadsheets/d/15B80Kg9pAFaukfDy-oW8MOHefyyKANIG1VvxwyQDaG4/edit#gid=1398094914"
        q <- read_sheet(url,  sheet = "queries") %>%
                janitor::clean_names() %>%
                select(genus = search_term, n_search_results = number_search_results)
        
        a <- annotated %>% select(genus, species, component) %>% distinct() %>%
                count(genus)
        
        d <- left_join(q, a) %>%
                mutate(n = ifelse(is.na(n), 0, n)) %>%
                rename(n_datasets_analyzed = n)
        
        write_csv(d, "data_github/search_genera.csv")
        return(d)
}

export_combined_xls <- function(etd, esg, emd, dsr = read_csv("data_github/mantel_results.csv")){
        require(openxlsx)
        
        readme = c("This file contains supplementary data for the study 'Isolation by wind.'",
                   "",
                   "The 'search_genera' sheet lists the tree genera used as search terms in Dryad, as well as the number of search results and the number of datasets included in our final study.",
                   "",
                   "The 'genetic_data_sources' sheet lists the published journal articles and Dryad data URLs for the genetic data we analyzed, as well as the species (sometimes more than one) and data type (sometimes more than one) analyzed for each study.",
                   "",
                   "The 'species_traits' sheet contains trait data for each species, as follows:",
                   "wind_pollination: 0 = non-wind-pollinated, 1 = wind-pollinated",
                   "wind_dispersal: 0 = non-wind-dispersed, 1 = wind-dispersed",
                   "pollination_months and dispersal_months: integers from 1-12 indicating months of species phenology; note that these data are not used for non-wind-influenced genomes and so are not listed for some of those taxa.",
                   "plastid_inheritance: pollen or seed",
                   "",
                   "The 'dataset_results' sheet reports, for each of the 120 datasets as well as the 6 genome control analyses, the number of populations, individuals, and loci in the dataset as well as the dataset-level partial correlation and associated Mantel p-value for each of the four hypotheses."
        )
        
        ds <- list("readme" = readme, 
                   "search_genera" = esg, 
                   "genetic_data_sources" = emd, 
                   "species_traits" = etd,
                   "dataset_results" = dsr %>% select(-component) %>% filter(syndrome != "0.5~0"))
        write.xlsx(ds, file = "data_github/supplementary_data.xlsx")
}
