module load rclone

# NOTE: to be run from Savio, within the 'GEE_outputs_and_derivatives/NIRv' directory

# copy all EOF results to Google Drive
rclone copy /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/NIRv/NIRv_4_EOFs_sqrt_coswts_standts.tif bdrive:LSP_EOF_results/
rclone copy /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/NIRv/NIRv_EOF_PC_table.csv bdrive:LSP_EOF_results/
rclone copy /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/NIRv/NIRv_EOF_scree_plot.png bdrive:LSP_EOF_results/
rclone copy /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/NIRv/NIRv_LSP_clust_scree_GLOBAL.png bdrive:LSP_EOF_results/
rclone copy /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/NIRv/NIRv_LSP_7_clusts_arr_GLOBAL.txt bdrive:LSP_EOF_results/
rclone copy /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/NIRv/NIRv_LSP_7_clusts_cents_GLOBAL.pkl bdrive:LSP_EOF_results/

