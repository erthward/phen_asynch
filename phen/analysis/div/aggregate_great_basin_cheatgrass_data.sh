# aggregating to same resolution as our phenology data
gdalwarp -tr 0.050000000000004 0.050000000000004 -r cubic ./WGA_Weighted_Mean_Annual_Herbaceous_Cover_2016_2017_2018.tif WGA_Weighted_Mean_Annual_Herbaceous_Cover_2016_2017_2018_AGG5KM.tif 

