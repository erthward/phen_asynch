
assemble_traits <- function(metadata){
        
        ### import data
        
        url <- "https://docs.google.com/spreadsheets/d/1onSN4dJ0_4A6kFQvL0sR5yWJ3RrNjcso8zaEby3O45Y/edit#gid=1840736075"
        traits <- read_sheet(url,  sheet = "traits", col_types = "c")
        
        metadata <- bind_rows(metadata) %>%
                select(species, genome, name, component)
        
        d <- traits %>%
                select(species, trait, value) %>% distinct() %>%
                spread(trait, value) %>%
                mutate(wind_dispersal = as.numeric(wind_dispersal),
                       wind_pollination = as.numeric(wind_pollination)) %>%
                right_join(metadata) %>%
                mutate(
                        expected_wind = case_when(
                                genome == "cp" & plastid_inheritance == "seed" ~ wind_dispersal,
                                genome == "cp" & plastid_inheritance == "pollen" ~ wind_pollination,
                                genome == "nu" ~ (wind_dispersal + wind_pollination) / 2),
                        genome_wind_poll = case_when(
                                genome == "cp" & plastid_inheritance == "pollen" & wind_pollination == 1 ~ TRUE,
                                genome == "cp" ~ FALSE,
                                genome == "nu" & wind_pollination == 1 ~ TRUE,
                                genome == "nu" ~ FALSE,
                                TRUE ~ NA),
                        genome_wind_disp = case_when(
                                genome == "cp" & plastid_inheritance == "seed" & wind_dispersal == 1 ~ TRUE,
                                genome == "cp" ~ FALSE,
                                genome == "nu" & wind_dispersal == 1 ~ TRUE,
                                genome == "nu" ~ FALSE,
                                TRUE ~ NA),
                        syndrome = case_when(
                                genome_wind_poll & genome_wind_disp ~ "both",
                                genome_wind_poll ~ "wind-pollination",
                                genome_wind_disp ~ "wind-dispersal",
                                !genome_wind_poll & !genome_wind_disp ~ "neither",
                                TRUE ~ NA_character_),
                        syndrome = factor(syndrome, levels=c("both", "wind-pollination", 
                                                             "wind-dispersal", "neither")))
        
        # identify relevant wind months for each dataset
        all_mos <- "1_2_3_4_5_6_7_8_9_10_11_12"
        w <- bind_rows(metadata) %>%
                left_join(d) %>%
                mutate(months = case_when(
                        expected_wind == 0 | is.na(expected_wind) ~ all_mos,
                        genome == "cp" & plastid_inheritance == "seed" & wind_dispersal ~ dispersal_months,
                        genome == "cp" & plastid_inheritance == "pollen" & wind_pollination ~ pollination_months,
                        genome == "nu" & wind_pollination == 1 & wind_dispersal == 1 ~
                                paste(pollination_months, dispersal_months, sep = "_"),
                        genome == "nu" & wind_pollination == 1 & wind_dispersal == 0 ~
                                pollination_months,
                        genome == "nu" & wind_pollination == 0 & wind_dispersal == 1 ~
                                dispersal_months,
                        TRUE ~ all_mos),
                       months = ifelse(months == "NA_NA", NA, months),
                       months = gsub("NA_|_NA", "", months),
                       months = ifelse(is.na(months), all_mos, months)
                       ) %>%
                select(component, months) %>%
                distinct()
        w$months <- sapply(w$months, 
                           function(x){x %>% str_split("_") %>% "[["(1) %>% unique() %>%
                                   as.numeric() %>% sort() %>% paste(collapse = "_")})
        d <- left_join(d, w)
        
        return(d)
}

generate_trait_template <- function(metadata){
        
        sp <-  bind_rows(metadata) %>%
                select(species) %>%
                distinct() %>%
                pull(species)
        
        d <- expand_grid(species = sp,
                         trait = c("wind_pollination", "wind_dispersal",
                                   "pollination_months", "dispersal_months",
                                   "plastid_inheritance"),
                         value = NA,
                         source = NA)
        
        write_csv(d, "data/traits/trait_template.csv")
}


update_trait_template <- function(traits){
        
        stop("run this manually, not for auto use")
        
        loadd(traits)
        sp <- traits %>% 
                filter(is.na(plastid_inheritance)) %>%
                select(species) %>%
                distinct() %>%
                pull(species)
        
        d <- expand_grid(species = sp,
                         trait = c("wind_pollination", "wind_dispersal",
                                   "pollination_months", "dispersal_months",
                                   "plastid_inheritance"),
                         value = NA,
                         source = NA)
        
        write_csv(d, "data/traits/trait_template_new_spp.csv")
 }
