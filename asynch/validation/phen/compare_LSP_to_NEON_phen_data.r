library(rnpn)

# TODO:
  # figure out a leaf-out-type phenophase
  # get site data for that across a bunch of dominant species
  # output to file
  # use Python to correlate with change-point detection on our LSP curves

start_yr = 2010
end_yr = 2015
yrs = as.character(seq(start_yr, end_yr))

# get all stations
stations = rnpn::npn)stations()
network_stations = stations[stations['network_id']!='', ]

# get all species IDs of dominant tree genera
species <- npn_species()

dom_trees = species[species$genus %in% c('Quercus'), ]#, 'Pinus', 'Populus', 'Acacia'), ]
#dom_trees = species[species$genus %in% c('Quercus', 'Pinus', 'Populus', 'Acacia'), ]
dom_tree_ids = dom_trees$species_id
dom_tree_family_ids = as.character(dom_trees$family_id)

# get oak phenophases
tree_phases = npn_get_phenophases_for_taxon(family_ids=dom_tree_family_ids,
                                            return_all=1)

# get site phenometrics for all oaks
tree_pms = <-npn_download_site_phenometrics(request_source = 'Drew Terasaki Hart',
                                           years = yrs,
                                           num_days_quality_filter = '30',
                                           species_ids = '35',
                                           phenophase_ids = '373')
