library(rnpn)
library(ggplot2)


# function for converting complex list of phenophase output to data.frame
get_phenophase_df = function(phenophase_list){
  # coerce to data.frame
  pp_df_list = list('family_name' = c(),
                    'family_id' = c(),
                    'phenophase_name' = c(),
                    'phenophase_category' = c(),
                    'phenophase_id' = c(),
                    'phenophase_definition' = c()
                    )
  for (i in seq(1: length(phenophase_list))){
     p = phenophase_list[[i]]
     for (j in seq(1: length(p$phenophases))){
        pp = p$phenophases[[j]]
        for (n in names(pp_df_list)){
           if (length(grep('family', n, value=T)) > 0){
              pp_df_list[[n]] = c(pp_df_list[[n]], p[[n]])
           } else {
              pp_df_list[[n]] = c(pp_df_list[[n]], pp[[n]])
           }
        }
     }
  }
  pp_df = as.data.frame(pp_df_list)
  return(pp_df)
}

# request all years' data
start_yr = 2014
end_yr = 2022
yrs = as.character(seq(start_yr, end_yr))

# get all stations
stations = rnpn::npn_stations()
network_stations = stations[stations['network_id']!='', ]

# get all species and family IDs of dominant tree genera
species <- npn_species()
dom_decid = species[species$genus %in% c('Quercus', 'Populus', 'Acer', 'Robinia', 'Prosopis',
                                         'Betula', 'Fagus', 'Liquidambar', 'Liriodendron',
                                         'Carya', 'Prunus', 'Juglans', 'Alnus'), ]
dom_tree_ids = unique(dom_decid$species_id)
dom_tree_family_ids = unique(as.character(dom_decid$family_id))

# get phenophases
tree_phases = npn_get_phenophases_for_taxon(family_ids=dom_tree_family_ids,
                                            return_all=1)
tp_df = get_phenophase_df(tree_phases)

# get the 'Leaves' phenophase id
leaf_ids = unique(tp_df[tp_df$phenophase_name == 'Leaves', ]$phenophase_id)
# NOTE: there should be two ids for the "Leaves" phase; the majority are 483,
#       but a couple rows show 488 for Fabaceae;
#       might be an error, but doesn't matter because they don't wind up in the data table anyhow
stopifnot("more than one phenophase id found!"=length(leaf_ids) == 2)

# get site leaf phenometrics for all dominant deciduous trees
pms = npn_download_site_phenometrics(request_source = 'Drew Terasaki Hart',
                                     years = yrs,
                                     num_days_quality_filter = '15', # max diff btwn individ's first 'yes' and prev 'no'
                                     species_ids = dom_tree_ids,
                                     phenophase_ids = leaf_ids,
                                     six_leaf_layer = T              # combine with the SI-x spring index dataset?
                                    )

# drop missing data
pms = pms[pms$mean_first_yes_doy >= 0, ]

# drop very small samples
pms = pms[pms$first_yes_sample_size < 5, ]

# change horrible column name for the SI-x data
colnames(pms) = c(colnames(pms)[1:length(colnames(pms))-1], 'six_leaf_val')

# sort by 'mean_first_yes_doy'
pms = pms[order(pms$mean_first_yes_doy), ]

# merge on the species table for some additional info
pms = base::merge(pms, species[, c('species_id', 'family_id', 'family_common_name', 'order_id', 'order_common_name',
                                   'class_id', 'class_common_name', 'functional_type')], by='species_id', all.x=T, all.y=F)

# write to disk
write.csv(pms, 'NPN_leaves_data_dom_tree_spp.csv', row.names=F)

