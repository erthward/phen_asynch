# NOTE: to be run from laptop, within the external hard drive's 'final_maps_and_results' directory,
#       which serves as the main map-data directory for all locally-run follow-on analysis

# copy all mosaicked GeoTIFFs back to Google Drive
rclone copy bdrive:LSP_asynch_drivers_analysis_covars/MODIS_IGBP_veg_entropy.tif /global/scratch/users/drewhart/seasonality/asynch_drivers_analysis/
rclone copy bdrive:LSP_asynch_drivers_analysis_covars/vrm_100KMmd_GMTEDmd.tif /global/scratch/users/drewhart/seasonality/asynch_drivers_analysis/

