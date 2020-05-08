for f in `ls /run/media/drew/SLAB/seasonality/SIF/OCO-2/gridded/Global_High_Res_SIF_OCO2_1696/data`;
do 
    gdal_translate -a_srs EPSG:4326 $f -of 'Gtiff' "${f%.*}.tif";
done
