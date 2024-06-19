# NOTE: to be run within the 'GEE_outputs/LSP_masks' directory, where all GEE LSP-fitting masks should be stored

# copy everything to local from Google Drive (where GEE wrote results)
rclone copy bdrive:LSP_mask_outputs_from_GEE/ /media/deth/SLAB/diss/3-phn/final_seas_and_asynch_maps/

# check that the following line prints '6 files downloaded'
echo "`ls * | wc -l` files downloaded"

