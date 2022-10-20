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
  - `/data`: ancillary data (not our main datasets) that lives in this repo and is used in analysis
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

1. Run `phen/calculation/GEE/other_datasets/calc_veg_entropy.js` to produce the vegetation cover entropy map that will be used as a covariate in the phenological asynchrony predictive model.
2. Run `phen/calculation/GEE/other_datasets/calc_dist_to_water.js` to produce the river-distance map that will be used as a covariate in the phenological asynchrony predictive model.
3. Download SRTM-based 50 km median vector ruggedness metric (file 'vrm_50KMmd_SRTM/tif') from the [EarthEnv website](http://www.earthenv.org/topography).
4. Download the CHELSA bio6 (daily min temperature of the coldest month) and bio15 (precipitation seasonality) data from the [CHELSA website](URL: https://chelsa-climate.org/bioclim/) using wget on the files in `asynch/analysis/rf/envidatS3paths.txt`!
5. Run `asynch/analysis/rf/calc_circular_moving_window_chelsa_rasters.r` to calculate the neighborhood mean and standard deviation of the CHELSA bio6 layer and the neighborhood standard deviation of the bio15 layer within 10-cell (i.e., ~55km at the equator) radii.
NOTE: Climate asynchrony maps calculated by `asynch/calculation/asynch_job.sh` will also be used as covariates in the phenological asynchrony predictive model.


## run gridded SIF orbital-gap seasonality validation:

1. Manually download TROPOMI data from ***WHERE FROM?***.
2. Run `phen/validation/orbital_gap_seasonality/validate_ANN-gridded_OCO2_orbital_gaps.py` to validate the seasonality of the OCO2-SIF ANN-interpolated data products within OCO2 orbital gaps for three regions across the pantropics.


## run flux-tower validation:

1. Manually download all subset data products (using DownThemAll!) from the Fluxnet network's [download page](https://fluxnet.org/data/download-data/).
2. Call `phen/validation/flux_tower_seasonality_validation/validate_at_all_fluxnet_sites.py <DS>` twice, once where '<DS>' is replaced with 'NIRv' and once with 'SIF', to run validation on both the fitted NIRv and SIF seasonality results against GPP seasonality at all usable Fluxnet sites (producing **Fig. 1**).
3. Run `phen/validation/compare_NIRv_SIF_maps/compare_NIRv_SIF_fitted_phenology.py` (on Savio) to calculate a global map of the $R^2$s between the NIRv and SIF fitted phenology time series.
4. Run `phen/validation/plot_phen_validation_results_FIG_2.py` to combine both of those validations' results to make Fig. 1.


## run asynchrony validation:

1. Run `asynch/validation/compare_SIF_and_NIRv_asynch.py` to compare the two datasets' phenological asynchrony maps across all three neighborhood radii (50 km, 100 km, 150 km).
2. Run `asynch/validation/calc_asynch_r2s_btwn_neighborhood_radii.py` to produce Table S2, containing R2s for all neighborhood radius comparisons and for all variables for which we produced asynchrony maps.


## produce RGB phenology map:

1. Download ancillary cheatgrass data from [Maestas *et. al*](https://www.sciencebase.gov/catalog/item/5ec5159482ce476925eac3b7) (to be used in a statistical test embedded in `phen/analysis/div/plot_eof_maps_and_ts_FIG_2.py`).
2. Run `phen/analysis/div/aggregate_great_basin_cheatgrass_data.sh` to aggregate that dataset to our analysis resolution of $0.05^{circ}$.
3. Run `phen/analysis/div/plot_eof_maps_and_ts_FIG_1.py` to produce **Fig. 2**'s global and regionally-zoomed RGB land surface phenology maps.


## produce asynchrony map and conceptual figure:

1. Run `asynch/viz/make_conceptual_fig_and_asynch_maps_FIG_3_S5-11.py` to create the asynch map figures for the main paper (**Fig. 3**) and the supplements (**Figs. S3, S4**).


## run phenological asynchrony modeling workflow:

1. On Savio, run `asynch/analysis/rf/prep_phen_asynch_rf_data.r NIRv 100` to prep data for random forest analysis of the main phenological asynchrony dataset (i.e., NIRv-based phenological asynchrony using a 100 km radial neighborhood).
2. In an RStudio session on Savio, run `asynch/analysis/rf/run_phen_asynch_rf.r` with var set to 'NIRv' and neigh.rad set to '100' (i.e., uncommenting lines at top), to execute the random forest analysis on the main phenological asynchrony dataset (i.e., NIRv-based phenological asynchrony using a 100 km radial neighborhood). Be sure the execute the code blocks captured by `if (F){ ... }`, to run hyperparameter-tuning, Boruta feature selection, and other interactive analyses.
3. Manually inspect the results of the interactive analysis. Use the results of that to set the hyperparameters (in the code block starting at line 410 in `asynch/analysis/rf/run_phen_asynch_rf.r`) and the feature selection (code block starting at line 486 in the same file) for the main global RF model that will be used for both datasets (NIRv and SIF) and all 3 neighborhood radii (50 km, 100 km, 150 km).
4. Run `asynch/analysis/rf/ch3_rf_job.sh` to loop over vars (NIRv, SIF) and neighborhood radii, each time prepping data layers, running the random forest analysis, and generating identical results.
5. Run `asynch/analysis/rf/ch3_rasterize_SHAP_job.sh` to convert output CSVs of global SHAP values to GeoTIFFs.
6. Run `asynch/analysis/rf/ch3_rasterize_err_job.sh` to convert output CSVs of global RF prediction errors to GeoTIFFs.
7. Run `asynch/analysis/rf/tabulate_model_summaries.py` to combine all permuation-based and SHAP-based importance values and model $R^2$s and MSEs into a single output table, for supplmental materials.
8. Run `python asynch/analysis/rf/make_shap_summary_map.py NIRv 100 y` to produce the SHAP-value interpretation map (for the 100 km-neighborhood NIRv-asynchrony analysis that included the geo-coordinate polynomials as covariates) and save result as a GeoTIFF.
9. Run `asynch/analysis/rf/plot_rf_summary_FIG_4.py` to produce final figure summarizing random forest results.


## run climate-distance analysis:

1. Run `asynch/analysis/clim_dist/compare_phen_clim_geog_dist_FIG_5.py` to run all iterations of the analysis of the latitudinal trend in the phenological distance~climatic distance relationship and produce **Fig. 5**.

