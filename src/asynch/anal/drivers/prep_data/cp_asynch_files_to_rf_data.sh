# NOTE: don't copy the non-strict-masked NIRv or SIF files, to avoid any risk
#       of accidentally using them as the random forest response variables
for dir in `echo NIRv_STRICT SIF_STRICT tmmn tmmx pr def cloud`
  do cp /global/scratch/users/drewhart/seasonality/GEE_outputs/$dir/*asynch*.tif /global/scratch/users/drewhart/seasonality/rf_data/
  done
