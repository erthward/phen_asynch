cd /media/deth/SLAB/diss/3-phn/SIF/OCO-2/gridded/Global_High_Res_SIF_OCO2_1696/data;
for f in `ls *nc`;
do 
    gdal_translate -a_srs EPSG:4326 $f -of 'Gtiff' "${f%.*}.tif";
done
cd -
