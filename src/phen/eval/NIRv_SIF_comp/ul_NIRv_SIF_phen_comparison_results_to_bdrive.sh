module load rclone

# NOTE: to be run from Savio, within the 'GEE_outputs_and_derivatives' directory, where the NIRv-SIF phen-comparison results tif is saved
rclone copy /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/NIRv_SIF_phen_R2s.tif bdrive:LSP_and_asynch_mosaicked_outputs/

