module load rclone

# NOTE: to be run from Savio, within the 'GEE_outputs/NIRv' directory

# copy all EOF results to Google Drive
rclone copy /global/scratch/users/drewhart/seasonality/GEE_outputs/NIRv/NIRv_4_EOFs_sqrt_coswts_normts.tif bdrive:LSP_EOF_results/
rclone copy /global/scratch/users/drewhart/seasonality/GEE_outputs/NIRv/NIRv_EOF_PC_table.csv bdrive:LSP_EOF_results/
rclone copy /global/scratch/users/drewhart/seasonality/GEE_outputs/NIRv/NIRv_EOF_scree_plot.png bdrive:LSP_EOF_results/

