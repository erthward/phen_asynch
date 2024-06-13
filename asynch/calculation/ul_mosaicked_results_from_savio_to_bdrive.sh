# NOTE: to be run from Savio, within the 'GEE_outputs' directory, within the subdirs of which
#       all the mosaicked results files are found
#       (the subdirs being NIRv, NIRv_STRICT, SIF, SIF_STRICT, tmmn, tmmx, pr, def, cloud)

# copy all mosaicked GeoTIFFs back to Google Drive
rclone copy ./NIRv/*tif bdrive:LSP_and_asynch_mosaicked_outputs/
rclone copy ./NIRv_STRICT/*tif bdrive:LSP_and_asynch_mosaicked_outputs/
rclone copy ./SIF/*tif bdrive:LSP_and_asynch_mosaicked_outputs/
rclone copy ./SIF_STRICT/*tif bdrive:LSP_and_asynch_mosaicked_outputs/
rclone copy ./tmmn/*tif bdrive:LSP_and_asynch_mosaicked_outputs/
rclone copy ./tmmx/*tif bdrive:LSP_and_asynch_mosaicked_outputs/
rclone copy ./pr/*tif bdrive:LSP_and_asynch_mosaicked_outputs/
rclone copy ./def/*tif bdrive:LSP_and_asynch_mosaicked_outputs/
rclone copy ./cloud/*tif bdrive:LSP_and_asynch_mosaicked_outputs/

