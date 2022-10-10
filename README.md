## overview:
This repo contains all the code used to produce the analysis in
"Global phenology maps reveal patterns of regional divergence,
intercontinental convergence, and climate-independent tropical asynchrony"
(Terasaki Hart et al. 2022), chapter 3 of Drew E. Terasaki Hart's PhD dissertation.

# < PUT CITATION HERE >

Code was written by Drew Ellison Terasaki Hart,
Thao-Nguyen Bui, and Lauren Di Maggio.
It is made freely available under the MIT License,
to be distributed and/or modified with proper attribution.


Any questions, concerns, or other interests can be directed
to drew <dot> hart <at> berkeley <dot> edu. 


### contents

Analyses were run in two major stages, so the two major directories of content are organized to reflect this:
  1. `/phen/`: calculation of global maps of characteristic seasonal phenology, and associated analysis
  2. `/asynch/`: calculation of global maps of phenological asynchrony, and associated analysis
Other directories include:
  - `/data_prep`: code used to download the SIF dataset used, convert to Geotiff, and upload to Google Earth Enginge (GEE)
  - `/data`: ancillary data (i.e., not our main datasets) that lives in this repo and is used in analysis
  - `/etc`: a script full of helper functions, and other odds and ends

***NOTE:*** **All GEE Javascript code must be run in GEE. Other code was designed either to be run on UC Berkeley's Savio compute cluster or on a local machine.**


### workflow:

## data prep:

1. Download SIF data (`data_prep/data_dl/dl_orig_OCO2.py`).
2. Translate SIF data from NetCDF to GeoTIFF (`data_prep/gee_data_ul/convert_SIF_OCO2_ANN_NetCDF_to_GTiff.sh`)
3. Create CSV of SIF metada (`data_prep/gee_data_ul/prep_SIF_OCO2_ANN_upload_metadata.py`)
4. Upload SIF data to GEE as individual Images (`data_prep/gee_data_ul/upload_SIF_OCO2_ANN_data_to_GEE_collections.py`)
5. Combine GEE SIF data into an ImageCollection (`data_prep/gee_data_ul/make_image_collection.sh`)


## calculate masking maps for seasonality-fitting procedure:

1. Call `phen/calculation/GEE/masking/calc_MODIS_month_evenness.js` to produce a global map of monthly Pielou's evenness of the MODIS data, as a GEE asset, to be used in the masking procedure for the final phenology maps.
2. Call `phen/calculation/GEE/masking/calc_land_cover_stability.js` to produce a global map indicating how 'stable' (i.e., time-invariant) each pixel's MODIS land cover is, as a GEE asset, to be used in the masking procedure for the final phenology maps.
3. Call `phen/calculation/GEE/masking/calc_permutation_fail_count_maps.js` numerous times, each time varying the values of `metric`, `seedStart`, and `nIts` so as to produce GEE assets representing, in total, 100 permutations of the seasonal phenology-fitting code for the SIF dataset and 25 for the NIRv dataset (because of compute limitations). See explanatory notes in the script, which document this workflow.
4. Call `phen/calculation/GEE/masking/calc_permutation_signif_mask.js` twice (once with `dataset` set to 'NIRv', once set to 'SIF') to combine the outputs of the previous step into the overall significance masks, to be used in the masking procedure for the final phenology maps.


## calculate fitted seasonality maps:

1. Call `phen/calculation/GEE/main.js` to launch a GEE job that will save tiled results to Google Drive.
  NOTE: Results will be formatted as overlapping TFRecord tiles.
  NOTE: In order to produce all results, this script must be manually called once for each combination of NIRv or SIF and default or strict masking, as well as once per climate variable (TerraClimate min temperature, precipitation, and climatic water deficit, as well as MODIS cloud). To do this, the `datasetName`, `climateVar`, and `maskingMode` variables must be manually swapped in this script.
  NOTE: To execute the full phenology-mapping workflow, main.js will in turn call a variety of other scripts located in `phen/calculation`.


## download seasonality results:

1. After all seasonality-fitting jobs successfully finish, use rclone to download all results into a common parent directory on UC Berkeley's Savio compute cluster, then manually move each run's results into a separately named child directory.


## calculate asynchrony maps:

1. Feed `asynch/calculation/asynch_job.sh` to slurm's sbatch command to calculate asynchrony maps for all fitted phenological and climatic seasonality datasets.
  NOTE: Do this once for each of the three neighborhood radii (50 km, 100 km, 150 km) for which we calculate asynchrony, each time concordantly changing the line in `asynch_job.sh` that sets the `neigh_rad` variable and the line that creates the slrum `--job-name` flag's value.
2. Once all 3 jobs are complete, feed `asynch/calculation/mosaic_job.sh` to slurm's sbatch command to mosaic all TFRecord files into a single GeoTIFF file for each phenological and climatic dataset (yielding a map of seasonality coefficients, a map of seasonality-fitting $R^{2}$s, and 3 maps of asynchrony, one per neighborhood radius).


## prepare other physiographic covariates:

1. Run `phen/calculation/GEE/other_datasets/calc_dist_to_ecotone.js` to produce the ecotone-distance map that will be used as a covariate in the phenological asynchrony predictive model.
2. Run `phen/calculation/GEE/other_datasets/calc_dist_to_water.js` to produce the river-distance map that will be used as a covariate in the phenological asynchrony predictive model.
3. Download SRTM-based 50 km median vector ruggedness metric (file 'vrm_50KMmd_SRTM/tif') from the [EarthEnv website](http://www.earthenv.org/topography).
4. Download the CHELSA bio6 (daily min temperature of the coldest month) and bio15 (precipitation seasonality) data from the [CHELSA website](URL: https://chelsa-climate.org/bioclim/) using wget on the files in ***FILENAME HERE***!
5. Run `asynch/analysis/corr/calc_circular_moving_window_chelsa_rasters.r` to calculate the neighborhood mean and standard deviation of the CHELSA bio6 layer and the neighborhood standard deviation of the bio15 layer within 10-cell (i.e., ~55km at the equator) radii.


## run gridded SIF orbital-gap seasonality validation:

1. Manually download TROPOMI data from ***WHERE FROM?***.
2. Run `phen/validation/orbital_gap_seasonality/validate_ANN-gridded_OCO2_orbital_gaps.py` to validate the seasonality of the OCO2-SIF ANN-interpolated data products within OCO2 orbital gaps for three regions across the pantropics.


## run flux-tower validation:

1. Manually download all subset data products (using DownThemAll!) from the Fluxnet network's [download page](https://fluxnet.org/data/download-data/).
2. Call `phen/validation/flux_tower_seasonality_validation/validate_at_all_fluxnet_sites.py` to validate the fitted NIRv seasonality against GPP seasonality at all usable Fluxnet sites (producing **Fig. 1**).
3. 


## compare NIRv-based and SIF-based phenological asynchrony maps:

1. Run `asynch/validation/compare_SIF_and_NIRv_asynch.py` to compare the two datasets' phenological asynchrony maps across all three neighborhood radii (50 km, 100 km, 150 km).


## produce RGB phenology map:

1. Run `phen/analysis/div/FIG2_rotate_eof0_and_map.py` to produce **Fig. 2**'s global and regionally-zoomed RGB land surface phenology maps.


## produce asynchrony map and conceptual figure:


## run phenological asynchrony modeling workflow:


## run climate-distance analysis:


