library(rtry)

try = rtry_import(input = '/media/deth/SLAB/diss/3-phn/try_db/32735.txt')

# get just the rows with plant growth form data
try = try[try$TraitID == 42, ]

# get just rows for major PGFs
forb_rows = grepl('(forb)|(herb)|(^H$)', try$OrigValueStr, ignore.case=T)
wood_rows = grepl('(tree)|(^T$)|((?<!sub)-?shrub)', try$OrigValueStr, ignore.case=T, perl=T)
vine_rows = grepl('(liana)|(vine)',  try$OrigValueStr, ignore.case=T)
epip_rows = grepl('epiphyte',  try$OrigValueStr, ignore.case=T)
xero_rows = grepl('xerophyte',  try$OrigValueStr, ignore.case=T)
try$forb = forb_rows
try$wood = wood_rows
try$vine = vine_rows
try$epip = epip_rows
try$xero = xero_rows
keep_rows = (forb_rows + wood_rows + vine_rows + epip_rows + xero_rows) > 0
try = try[keep_rows, ]

# subset columns
try = try[ , c('AccSpeciesID', 'SpeciesName', 'AccSpeciesName', 'ObservationID', 'OrigValueStr', 'forb', 'wood', 'vine', 'epip', 'xero')]

# write out
write.csv(try, './TRY_PGF_values_reclassed_filtered.csv', row.names=F)

