# NOTE: to be run within the 'GEE_outputs' directory where all GEE output data will be stored,
#       and that directory must have a properly named subdirectory for each variable
#       (NIRv, NIRv_STRICT, SIF, SIF_STRICT, tmmn, tmmx, pr, def, cloud)

# copy everything to local from Google Drive (where GEE wrote results)
rclone copy bdrive:LSP_outputs_from_GEE/ .

# move everything into subdirs
mv *_NIRv_STRICT* ./NIRv_STRICT/
mv *_NIRv* ./NIRv/
mv *_SIF_STRICT* ./SIF_STRICT/
mv *_SIF* ./SIF/
mv *_tmmn* ./tmmn/
mv *_tmmx* ./tmmx/
mv *_pr* ./pr/
mv *_def/
mv *_MODISCloud* ./MODISCloud/

# assert that the following prints '208' nine times
ls ./NIRv_STRICT/* | wc -l
ls ./NIRv/* | wc -l
ls ./SIF_STRICT/* | wc -l
ls ./SIF/* | wc -l
ls ./tmmn/* | wc -l
ls ./tmmx/* | wc -l
ls ./pr/* | wc -l
ls ./def/* | wc -l
ls ./cloud/* | wc -l

