cd /run/media/drew/SLAB/seasonality/SIF/OCO-2/gridded/Global_High_Res_SIF_OCO2_1696/data
gsutil -m cp *.tif gs://seasonality_data/OCO2_SIF_ANN/
cd /run/media/drew/SLAB/seasonality/other/cloud
gsutil -m cp *.tif gs://seasonality_data/EarthEnv_MODIS_cloud/
