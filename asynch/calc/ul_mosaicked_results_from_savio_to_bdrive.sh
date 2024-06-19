# NOTE: to be run from Savio, within the 'GEE_outputs' directory, within the subdirs of which
#       all the mosaicked results files are found
#       (the subdirs being NIRv, NIRv_STRICT, SIF, SIF_STRICT, tmmn, tmmx, pr, def, cloud)

# copy all mosaicked GeoTIFFs back to Google Drive
rclone copy --include *tif /global/scratch/users/drewhart/seasonality/GEE_outputs/NIRv/ bdrive:LSP_and_asynch_mosaicked_outputs/
rclone copy --include *tif /global/scratch/users/drewhart/seasonality/GEE_outputs/NIRv_STRICT/ bdrive:LSP_and_asynch_mosaicked_outputs/
rclone copy --include *tif /global/scratch/users/drewhart/seasonality/GEE_outputs/SIF/ bdrive:LSP_and_asynch_mosaicked_outputs/
rclone copy --include *tif /global/scratch/users/drewhart/seasonality/GEE_outputs/SIF_STRICT/ bdrive:LSP_and_asynch_mosaicked_outputs/
rclone copy --include *tif /global/scratch/users/drewhart/seasonality/GEE_outputs/tmmn/ bdrive:LSP_and_asynch_mosaicked_outputs/
rclone copy --include *tif /global/scratch/users/drewhart/seasonality/GEE_outputs/tmmx/ bdrive:LSP_and_asynch_mosaicked_outputs/
rclone copy --include *tif /global/scratch/users/drewhart/seasonality/GEE_outputs/pr/ bdrive:LSP_and_asynch_mosaicked_outputs/
rclone copy --include *tif /global/scratch/users/drewhart/seasonality/GEE_outputs/def/ bdrive:LSP_and_asynch_mosaicked_outputs/
rclone copy --include *tif /global/scratch/users/drewhart/seasonality/GEE_outputs/cloud/ bdrive:LSP_and_asynch_mosaicked_outputs/

