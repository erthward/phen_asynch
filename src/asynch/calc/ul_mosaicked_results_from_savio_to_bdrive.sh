module load rclone

# NOTE: to be run from Savio, within the 'GEE_outputs_and_derivatives' directory, within the subdirs of which
#       all the mosaicked results files are found
#       (the subdirs being NIRv, NIRv_STRICT, SIF, SIF_STRICT, tmmn, tmmx, pr, def, cloud)

# copy all mosaicked GeoTIFFs back to Google Drive
rclone copy --include *tif /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/NIRv/ bdrive:LSP_and_asynch_mosaicked_outputs/
rclone copy --include *tif /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/NIRv_STRICT/ bdrive:LSP_and_asynch_mosaicked_outputs/
rclone copy --include *tif /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/SIF/ bdrive:LSP_and_asynch_mosaicked_outputs/
rclone copy --include *tif /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/SIF_STRICT/ bdrive:LSP_and_asynch_mosaicked_outputs/
rclone copy --include *tif /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/tmmn/ bdrive:LSP_and_asynch_mosaicked_outputs/
rclone copy --include *tif /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/tmmx/ bdrive:LSP_and_asynch_mosaicked_outputs/
rclone copy --include *tif /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/pr/ bdrive:LSP_and_asynch_mosaicked_outputs/
rclone copy --include *tif /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/def/ bdrive:LSP_and_asynch_mosaicked_outputs/
rclone copy --include *tif /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/cloud/ bdrive:LSP_and_asynch_mosaicked_outputs/

