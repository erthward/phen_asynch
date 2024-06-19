# NOTE: to be run within the 'GEE_outputs' directory where all GEE output data will be stored,
#       and that directory must have a properly named subdirectory for each variable
#       (NIRv, NIRv_STRICT, SIF, SIF_STRICT, tmmn, tmmx, pr, def, cloud)

# copy everything to local from Google Drive (where GEE wrote results)
rclone copy bdrive:LSP_outputs_from_GEE/ /global/scratch/users/drewhart/seasonality/GEE_outputs/

# move everything into subdirs
mv /global/scratch/users/drewhart/seasonality/GEE_outputs/*_NIRv_STRICT* /global/scratch/users/drewhart/seasonality/GEE_outputs/NIRv_STRICT/
mv /global/scratch/users/drewhart/seasonality/GEE_outputs/*_NIRv* /global/scratch/users/drewhart/seasonality/GEE_outputs/NIRv/
mv /global/scratch/users/drewhart/seasonality/GEE_outputs/*_SIF_STRICT* /global/scratch/users/drewhart/seasonality/GEE_outputs/SIF_STRICT/
mv /global/scratch/users/drewhart/seasonality/GEE_outputs/*_SIF* /global/scratch/users/drewhart/seasonality/GEE_outputs/SIF/
mv /global/scratch/users/drewhart/seasonality/GEE_outputs/*_tmmn* /global/scratch/users/drewhart/seasonality/GEE_outputs/tmmn/
mv /global/scratch/users/drewhart/seasonality/GEE_outputs/*_tmmx* /global/scratch/users/drewhart/seasonality/GEE_outputs/tmmx/
mv /global/scratch/users/drewhart/seasonality/GEE_outputs/*_pr* /global/scratch/users/drewhart/seasonality/GEE_outputs/pr/
mv /global/scratch/users/drewhart/seasonality/GEE_outputs/*_def* /global/scratch/users/drewhart/seasonality/GEE_outputs/def/
mv /global/scratch/users/drewhart/seasonality/GEE_outputs/*_MODISCloud* /global/scratch/users/drewhart/seasonality/GEE_outputs/MODISCloud/

# assert that the following prints '208' nine times
ls /global/scratch/users/drewhart/seasonality/GEE_outputs/NIRv_STRICT/* | wc -l
ls /global/scratch/users/drewhart/seasonality/GEE_outputs/NIRv/* | wc -l
ls /global/scratch/users/drewhart/seasonality/GEE_outputs/SIF_STRICT/* | wc -l
ls /global/scratch/users/drewhart/seasonality/GEE_outputs/SIF/* | wc -l
ls /global/scratch/users/drewhart/seasonality/GEE_outputs/tmmn/* | wc -l
ls /global/scratch/users/drewhart/seasonality/GEE_outputs/tmmx/* | wc -l
ls /global/scratch/users/drewhart/seasonality/GEE_outputs/pr/* | wc -l
ls /global/scratch/users/drewhart/seasonality/GEE_outputs/def/* | wc -l
ls /global/scratch/users/drewhart/seasonality/GEE_outputs/cloud/* | wc -l

