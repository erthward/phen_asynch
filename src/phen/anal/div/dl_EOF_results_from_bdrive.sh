# NOTE: to be run from laptop, within the external hard drive's 'final_maps_and_results' directory,
#       which serves as the main map-data directory for all locally-run follow-on analysis

# copy all EOF results to laptop
rclone copy bdrive:LSP_EOF_results/NIRv_4_EOFs_sqrt_coswts_standts.tif /media/deth/SLAB/diss/3-phn/final_maps_and_results/
rclone copy bdrive:LSP_EOF_results/NIRv_EOF_PC_table.csv /media/deth/SLAB/diss/3-phn/final_maps_and_results/
rclone copy bdrive:LSP_EOF_results/NIRv_EOF_scree_plot.png /media/deth/SLAB/diss/3-phn/final_maps_and_results/

