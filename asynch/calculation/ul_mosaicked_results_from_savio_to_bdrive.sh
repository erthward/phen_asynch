# NOTE: to be run from Savio, within the 'GEE_outputs' directory, within the subdirs of which
#       all the mosaicked results files are found
#       (the subdirs being NIRv, NIRv_STRICT, SIF, SIF_STRICT, tmmn, tmmx, pr, def, cloud)

# copy all mosaicked GeoTIFFs back to Google Drive
rclone copy --include *tif ./NIRv/ bdrive:LSP_and_asynch_mosaicked_outputs/
rclone copy --include *tif ./NIRv_STRICT/ bdrive:LSP_and_asynch_mosaicked_outputs/
rclone copy --include *tif ./SIF/ bdrive:LSP_and_asynch_mosaicked_outputs/
rclone copy --include *tif ./SIF_STRICT/ bdrive:LSP_and_asynch_mosaicked_outputs/
rclone copy --include *tif ./tmmn/ bdrive:LSP_and_asynch_mosaicked_outputs/
rclone copy --include *tif ./tmmx/ bdrive:LSP_and_asynch_mosaicked_outputs/
rclone copy --include *tif ./pr/ bdrive:LSP_and_asynch_mosaicked_outputs/
rclone copy --include *tif ./def/ bdrive:LSP_and_asynch_mosaicked_outputs/
rclone copy --include *tif ./cloud/ bdrive:LSP_and_asynch_mosaicked_outputs/

