# NOTE: to be run within the 'GEE_outputs_and_derivatives' directory where all GEE output data will be stored,
#       and that directory must have a properly named subdirectory for each variable
#       (NIRv, NIRv_STRICT, SIF, SIF_STRICT, tmmn, tmmx, pr, def, cloud)

# copy everything to local from Google Drive (where GEE wrote results)
rclone copy bdrive:LSP_outputs_from_GEE/ /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/

# move everything into subdirs
mv /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/*_NIRv_STRICT* /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/NIRv_STRICT/
mv /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/*_NIRv* /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/NIRv/
mv /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/*_SIF_STRICT* /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/SIF_STRICT/
mv /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/*_SIF* /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/SIF/
mv /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/*_tmmn* /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/tmmn/
mv /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/*_tmmx* /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/tmmx/
mv /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/*_pr* /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/pr/
mv /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/*_def* /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/def/
mv /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/*_MODISCloud* /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/MODISCloud/

# assert that the following prints '208' nine times
ls /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/NIRv_STRICT/* | wc -l
ls /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/NIRv/* | wc -l
ls /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/SIF_STRICT/* | wc -l
ls /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/SIF/* | wc -l
ls /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/tmmn/* | wc -l
ls /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/tmmx/* | wc -l
ls /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/pr/* | wc -l
ls /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/def/* | wc -l
ls /global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/cloud/* | wc -l

