module load rclone
# NOTE: to be run from Savio

# copy to Google Drive all results required for final analyses and figures
rclone copy --include rf_SHAP_importance_* /global/scratch/users/drewhart/seasonality/rf_data/ bdrive:LSP_asynch_drivers_rf_results
rclone copy --include rf_permut_importance* /global/scratch/users/drewhart/seasonality/rf_data/ bdrive:LSP_asynch_drivers_rf_results
rclone copy --include err_map_* /global/scratch/users/drewhart/seasonality/rf_data/ bdrive:LSP_asynch_drivers_rf_results
rclone copy --include SHAP_map_* /global/scratch/users/drewhart/seasonality/rf_data/ bdrive:LSP_asynch_drivers_rf_results

# OPTIONAL: copy other stuff (tuning results, summary plots, intermediate datasets, and slurm *.out files)
#rclone copy --include asynch_model_all_vars_* /global/scratch/users/drewhart/seasonality/rf_data/ bdrive:LSP_asynch_drivers_rf_results
#rclone copy --include preds_plot_* /global/scratch/users/drewhart/seasonality/rf_data/ bdrive:LSP_asynch_drivers_rf_results
#rclone copy --include var_import_plots_* /global/scratch/users/drewhart/seasonality/rf_data/ bdrive:LSP_asynch_drivers_rf_results
#rclone copy --include rf_SHAP_vals_* /global/scratch/users/drewhart/seasonality/rf_data/ bdrive:LSP_asynch_drivers_rf_results
#rclone copy --include rf_full_preds_* /global/scratch/users/drewhart/seasonality/rf_data/ bdrive:LSP_asynch_drivers_rf_results
#rclone copy --include tuning_results_* /global/scratch/users/drewhart/seasonality/rf_data/ bdrive:LSP_asynch_drivers_rf_results
#rclone copy --include ch3_rf_*Rout /global/scratch/users/drewhart/seasonality/rf_data/ bdrive:LSP_asynch_drivers_rf_results
#rclone copy --include slurm*out /global/scratch/users/drewhart/seasonality/rf_data/ bdrive:LSP_asynch_drivers_rf_results

#rclone copy /global/scratch/users/drewhart/seasonality/rf_data/boruta_boxplot_NIRv_100km.jpg bdrive:LSP_asynch_drivers_rf_results
#rclone copy /global/scratch/users/drewhart/seasonality/rf_data/rast_err.pyout bdrive:LSP_asynch_drivers_rf_results
#rclone copy /global/scratch/users/drewhart/seasonality/rf_data/rast_shap.pyout bdrive:LSP_asynch_drivers_rf_results

