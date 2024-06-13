# NOTE: to be run within the 'GEE_outputs/LSP_masks' directory, where all GEE LSP-fitting masks should be stored

# copy everything to local from Google Drive (where GEE wrote results)
rclone copy bdrive:LSP_mask_outputs_from_GEE/ .

# check that the following line prints '6'
ls * | wc -l

