for dir in `echo NIRv NIRv_STRICT SIF SIF_STRICT tmmn tmmx pr def cloud`
  do cp /global/scratch/users/drewhart/seasonality/GEE_outputs/$dir/*asynch*.tif /global/scratch/users/drewhart/seasonality/rf_data/
  done
