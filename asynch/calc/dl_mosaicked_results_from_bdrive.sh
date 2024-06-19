# NOTE: to be run from laptop, within the external hard drive's 'final_maps_and_results' directory,
#       which serves as the main map-data directory for all locally-run follow-on analysis

# copy all mosaicked GeoTIFFs back to Google Drive
rclone copy bdrive:LSP_and_asynch_mosaicked_outputs/ /media/deth/SLAB/diss/3-phn/final_maps_and_results/

# ensure that the following prints out '45' (5 layers for each of 9 variables)
ls *.tif | wc -l
