library(phenocamapi)

# get all regions of interest across all sites,
# and save to disk for site-level information used in analysis
rois = get_rois()
rois_filepath = "/home/deth/Desktop/CAL/research/projects/seasonality/seasonal_asynchrony/src/phen/eval/phenocam_GCC/ROIs.csv"
write.csv(rois, rois_filepath, sep=',')

