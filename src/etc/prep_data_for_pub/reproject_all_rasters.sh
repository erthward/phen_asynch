for f in `ls *.tif`;
  do gdalwarp -t_srs EPSG:3857 -of GTiff -co "COMPRESS=DEFLATE" -co "PREDICTOR=2" $f ${f/%.tif/_EPSG3857.tif};
  done
