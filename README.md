# < PUT TITLE HERE >
### < PUT CITATION HERE >
This repo contains all the code used to produce the analysis
for this paper and for its precursor, chapter 3 of my PhD dissertation.

Code was written by me (Drew),
with contributions from Lauren Di Maggio and Thao-Nguyen Bui.

Code is made freely available under the MIT License,
to be distributed and/or modified with due attribution.

Data are archived **PUT DATA DOI HERE**

Any questions, concerns, interests, or requests can be directed
to me at: drew *dot* hart *at* berkeley *dot* edu. 



-------------------------------------------



## contents
Analyses were run in two major stages, so the two major directories of content are organized to reflect this:

  1. `/phen/`: calculation of global maps of characteristic seasonal phenology, and associated analysis
  2. `/asynch/`: calculation of global maps of phenological asynchrony, and associated analysis

Other directories include:

  - `/data_prep`: code used to download the SIF dataset used, convert to Geotiff, and upload to Google Earth Enginge (GEE)
  - `/data`: ancillary data (not our main datasets) that lives in this repo and is used in analysis
  - `/etc`: includes a script of commonly used helper functions, a Python port of [MMRR] (https://onlinelibrary.wiley.com/doi/10.1111/evo.12134), and other odds and ends



-------------------------------------------



## workflow:
Each step of the following workflow was executed in the environment indicated *in italics* (either on GEE, on the UC Berkeley Savio cluster, or on my laptop). See **working environments**, below, for the software specifications of each of these environments.


### prep and evaluate SIF dataset:
1. *On laptop*, download [SIF dataset](https://doi.org/10.3334/ORNLDAAC/1696)
2. *On laptop*, translate SIF data from NetCDF to GeoTIFF (`data_prep/sif_gee_ul/convert_SIF_OCO2_ANN_NetCDF_to_GTiff.sh`)
3. *On laptop*, create CSV of SIF metadata (`data_prep/sif_gee_ul/prep_SIF_OCO2_ANN_upload_metadata.py`)
4. *On laptop*, upload SIF data to GEE as individual Images (`data_prep/sif_gee_ul/upload_SIF_OCO2_ANN_data_to_GEE_collections.py`)
5. *On laptop*, combine GEE SIF data into an ImageCollection (`data_prep/sif_gee_ul/make_image_collection.sh`)
6. *On laptop*, manually download [gridded TROPOMI data](ftp://fluo.gps.caltech.edu/data/tropomi/gridded/).
7. *On laptop*, run `phen/evaluation/orbital_gaps/evaluate_ANN-gridded_OCO2_orbital_gaps_FIG_S3_S4.py` to check that the seasonality of the OCO2-SIF ANN-interpolated data products within OCO2 orbital gaps compares favorably to another gridded SIF dataset, for three regions across the tropics.

### calculate masking and preprocessing maps for seasonality-fitting procedure:
1. *On GEE*, establish all desired parameter values for data reading and masking, harmonic regression fitting, significance calculation, and data exports in `phen/calculation/GEE/params.js`. (Note: To produce all needed outputs, this file will be resaved with altered commenting a few times during the following steps.)
2. *On GEE*, run `phen/calculation/GEE/masking/calc_month_data_avail_props.js` to produce a set of 12 global Image assets, each indicating the average proportion of monthly MODIS NBAR data availability at each pixel, to be used for masking based on proportional data availability and to derive the monthly data-availability evenness mask. (**12 tasks, each ~20-30m runtime, but note that a few may fail the first time and need to be rerun until successful**)
3. *On GEE*, run `phen/calculation/GEE/masking/calc_month_evenness.js` to produce a global map of monthly Pielou's evenness of the MODIS NBAR data, as a GEE asset, derived from the outputs of `phen/calculation/GEE/masking/calc_month_data_avail_props.js`. (**2 tasks, both <10m runtime**)
4. *On GEE*, run `phen/calculation/GEE/masking/calc_land_cover_mask.js` to produce a global map indicating all pixels that are valid (i.e., not water, urban, barren, or permanent ice/snow) and that are always 'natural' across our time series (i.e., forest, savanna, shrub, grass, wetland) and/or always agricultural. (All of those pixels will be used in the LSP analyses, but only the 'always natural' pixels will be used in the asynchrony analyses, to try to avoid places subject to anthropogenically created spatial asynchrony.) (**2 tasks, both <10m runtime**)
5. *On GEE*, run `phen/calculation/GEE/masking/calc_ts_pct_data_availability.js` to produce a map of overall data availability for the entire 20-year archive. (**1 task, ~4h runtime**)
6. *On GEE*, run `phen/calculation/GEE/masking/calc_permutation_fail_count_maps.js` to produce GEE assets representing, in total, 20 permutations of the seasonal phenology-fitting code for the NIRv data (unreasonable to produce more because of compute limitations). (**10 tasks, each ~6-20h runtime**)
7. *On GEE*, run `phen/calculation/GEE/masking/calc_permutation_signif_mask.js` to combine the outputs of the previous step into the overall significance masks, to be used in the masking procedure for the final phenology maps. (**1 task, <10m runtime**)
8. *On GEE*, run `phen/calculation/GEE/masking/create_overall_mask.js` twice, once with 'maskingMode' set to 'default' in `phen/calculation/GEE/params.js` and once with it set to 'strict', to combine all previously created masks into a pair of default (used for LSP analyses) and strict (used for asynchrony analyses) mask assets. (**7 tasks total, <10m runtime**)
9. *On GEE*, run `phen/calculation/GEE/io/calc_min_pos_NIRv.js` to calculate a global map asset of each pixel's minimum positive NIRv value (to be used to backfill negative NIRv values that occur, mainly in places and times with solid snow cover. (** 1 task, ~8h runtime**)


### calculate fitted seasonality maps:
1. *On GEE*, run `phen/calculation/GEE/fft/create_NIRv_harmonic_regression_asset.js` to run the main harmonic regression analysis on the MODIS NIRv dataset. (This is by far the most computationally intensive, so it's best to just run it once and save the result as an asset, which is then loaded, masked, and export to Drive by `phen/calculation/GEE/main.js`.) (**1 task, ~20h runtime**)
2. *On GEE*, run `phen/calculation/GEE/main.js` to launch a GEE task that will save tiled results to Google Drive.
  NOTE: Results will be formatted as overlapping TFRecord tiles.
  NOTE: In order to produce all results, this script must be manually called once for each combination of NIRv or SIF and default or strict masking, as well as once per climate variable (TerraClimate min temperature, precipitation, and climatic water deficit, as well as MODIS cloud). To do this, the `datasetName`, `climateVar`, and `maskingMode` variables must be manually swapped in this script.
  NOTE: To execute the full phenology-mapping workflow, main.js will in turn call a variety of other scripts located in `phen/calculation`.
  (**9 tasks total, ranging from 10s of minutes (NIRv, SIF and TerraClimate tasks) to a handful of hours (MODIS cloud task)**)


### download seasonality results:
1. *On both laptop and Savio*, navigate to the 'GEE\_outputs' directory where all fitted LSP and seasonality files from GEE should be stored, then run `phen/calculation/dl_and_organize_GEE_data.sh` to a.) download all results into the parent directory (using `rclone`) (both on UC Berkeley's Savio cluster and locally on external hard drive), then b.) move each run's results into corresponding child directory.


### calculate asynchrony:
1. *On Savio*, run `asynch/calculation/run_all_asynch_jobs.sh` to feed all three job scripts, one per asynchrony neighborhood radius, to slurm's `sbatch` command, calculating each neighborhood radius' set of asynchrony maps for all fitted phenological and climatic seasonality datasets (**3 slurm jobs total, runtimes of ~13h for the 50km neighborhood job, ~30h for 100km, and ~46h for 150km**).


### mosaic and store all results:
1. *On Savio*, run `asynch/calculation/mosaic_job.sh` to mosaic the regression coefficient, regression $R^2$, and asynchrony-result files for all LSP and climate variables and for all three asynchrony neighborhoods (50 km, 100 km, 150 km), producing a set of GeoTIFF outputs for downstream plotting and analysis (**1 job, ~45m runtime**).
2. *On Savio*, run `asynch/calculation/ul_mosaicked_results_from_savio_to_bdrive.sh` to copy all mosaicked results from Savio back up to Google Drive (**1 task, ~15m runtime**)
3. *On laptop*, run `asynch/calculation/dl_mosaicked_results_from_bdrive.sh` to then also copy those results down to their intended location on laptop external driv (**1 task, ~1h runtime**).


### map masks and R2s from harmonic regressions:
1. *On laptop*, navigate to the directory where the mask files should be stored, then run `phen/calculation/masking/dl_GEE_masks.sh` to download all 6 mask GeoTIFFs output by GEE (**1 task, <5m runtime**).
2. *On laptop*, run `phen/calculation/masking/make_masking_maps_supp.py` to produce **Fig. SXXX**, showing all masks used on the LSP datasets (**1 task, <5m runtime**).
3. *On laptop*, run `phen/calculation/R2/make_R2_maps_supp.py` to produce **Fig. SXXX**, showing the $R^2$s of the harmonic regressions fitted to all LSP and climate datasets (**1 task, <5m runtime**).


### produce RGB phenology map:
1. *On Savio*, run `phen/analysis/div/calc_LSP_EOFs.py` on Savio to calculate the global NIRv LSP EOF map and save results (**1 job, <1h total runtime**).
2. *On Savio*, run `phen/analaysis/div/ul_LSP_EOFs.sh` to push the EOF results back to BDrive (**1 job, <5m runtime**).
3. *On laptop*, run `phen/analysis/div/dl_LSP_EOFs.sh` to download the EOF results to the local hard drive (**1 job, <5m runtime**).
4. *On laptop*, download ancillary cheatgrass data from [Maestas *et. al*](https://www.sciencebase.gov/catalog/item/5ec5159482ce476925eac3b7) (to be used in a statistical test embedded in `phen/analysis/div/plot_eof_maps_and_ts_FIG_S3.py`).
5. *On laptop*, run `phen/analysis/div/aggregate_great_basin_cheatgrass_data.sh` to aggregate that dataset to our analysis resolution of $0.05^{\circ}$.
6. *On laptop*, work through the manual steps listed in the notes at the top of `phen/analysis/div/compose_ITCZ_shapefile.py` to digitize the Dec-Jan-Feb and Jun-Jul-Aug mean ITCZ locations delineated by [Zhisheng et al. 2015](annualreviews.org/content/journals/10.1146/annurev-earth-060313-054623), save the output CSVs to the 'data/' subdirectory of the local clone of this repo, then run `phen/analysis/div/compose_ITCZ_shapefile.py` to produce a Shapefile of those digitized lines.
7. *On laptop*, run `phen/analysis/div/make_EOF_and_RGB_map_figs.sh` to produce **Fig. 1**'s global EOF map, the focal-region EOF maps, the raw EOF-map figure, and the unfolded EOF maps (**1 task, <30m runtime**).
8. *On laptop*, manually assemble the focal-region map figure odg files using the outputs of the previous step.


### run NPN and SI-x evaluation:
1. *On laptop*, run `phen/evaluation/NPN_and_SI-x/get_NPN_leaf_data.r` to download, for a wide range of dominant US tree genera and at all NPN sites, both the day of year of first leaf based on NPN ground observations and mean day of year of start of season (SOS) based on MODIS-derived SI-x phenology maps (**1 task, <1h runtime**).
2. *On laptop*, run `phen/evaluation/NPN_and_SI-x/compare_NIRv_LSP_to_NPN_first_leaf.py` to evaluate SOS estimates derived from our NIRv LSP data against both the NPN first-leaf and SI-x SOS datasets (**1 task, <1m runtime**).


### run NIRv-SIF asynchrony comparison:
1. *On Savio*, run `phen/evaluation/compare_NIRv_SIF_maps/ch3_phen_comparison_job.sh` to calculate a global map of $R^2$ values between the fitted annual NIRv-based and SIF-based LSP patterns (**1 task, ~XXXXh runtime**).
2. **ADD UL SCRIPT**
3. **ADD DL SCRIPT**


### run FLUXNET evaluation:
1. *On laptop*, manually download all subset data products (using DownThemAll!) from the FLUXNET network's [download page](https://fluxnet.org/data/download-data/) (**1 task, runs roughly overnight**).
2. *On laptop*, run `phen/evaluation/flux_tower_GPP/run_flux_evaluations.sh` to run the flux-tower GPP comparison, at all usable FLUXNET sites, for both the fitted NIRv and SIF LSP results (producing Fig. XXX) (**1 task, ~XXXXh runtime**).
4. *On laptop*, run `phen/evaluation/plot_phen_evaluation_results_FIG_S5.py` to combine both LSP datasets' FLUXNET evaluations and the NIRv-SIF comparison evaluation to make Fig. XXX (**1 task, ~XXXXh runtime**).
5. *On laptop*, run `phen/evaluation/flux_tower_GPP/combine_fluxnet_val_outputs_TableS1.py` to combine of flux-tower evaluation results, listed by FLUXNET site, into a single supplemental table (**1 task, <1m runtime**).


### run asynchrony neighborhood-comparison evaluation:
1. *On laptop*, run `asynch/evaluation/calc_asynch_r2s_btwn_neighborhood_radii.py` to produce Table S2, containing R2s for all neighborhood radius comparisons and for all variables for which we produced asynchrony maps.


### run asynchrony NIRv-SIF comparison:
1. *On laptop*, run `asynch/evaluation/SIF_comp/compare_SIF_and_NIRv_asynch.py` to compare the two datasets' phenological asynchrony maps across all three neighborhood radii (50 km, 100 km, 150 km).


### produce asynchrony conceptual figure:
1. *On laptop*, run `asynch/viz/make_conceptual_FIG_SXXX.py` to create the asynchrony-calculation conceptual figure (Fig. SXXXX) (**1 task, <1m runtime**).


### prepare physiographic covariates:
1. *On GEE*, run `phen/calculation/GEE/other_datasets/calc_veg_entropy.js` to produce the vegetation cover entropy map that will be used as a covariate in the phenological asynchrony predictive model. (**1 task, <10m runtime**)
2. *On Savio*, download SRTM-based 50 km median vector ruggedness metric (file 'vrm_50KMmd_SRTM/tif') from the [EarthEnv website](http://www.earthenv.org/topography).
NOTE: Climate asynchrony maps calculated by `asynch/calculation/asynch_job.sh` will also be used as covariates in the phenological asynchrony predictive model.


### run phenological asynchrony drivers analysis:
1. On Savio, run `asynch/analysis/drivers/prep_data/prep_phen_asynch_rf_data.r NIRv 100` to prep data for random forest analysis of the main phenological asynchrony dataset (i.e., NIRv-based phenological asynchrony using a 100 km radial neighborhood).
2. In an RStudio session on Savio, run `asynch/analysis/drivers/run_rf/run_phen_asynch_rf.r` with var set to 'NIRv' and neigh.rad set to '100' (i.e., uncommenting lines at top), to execute the random forest analysis on the main phenological asynchrony dataset (i.e., NIRv-based phenological asynchrony using a 100 km radial neighborhood). Be sure the execute the code blocks captured by `if (F){ ... }`, to run hyperparameter-tuning, Boruta feature selection, and other interactive analyses.
3. Manually inspect the results of the interactive analysis. Use the results of that to set the hyperparameters (in the code block starting at line 410 in `asynch/analysis/drivers/run_rf/run_phen_asynch_rf.r`) and the feature selection (code block starting at line 486 in the same file) for the main global RF model that will be used for both datasets (NIRv and SIF) and all 3 neighborhood radii (50 km, 100 km, 150 km).
4. Run `asynch/analysis/drivers/run_rf/ch3_rf_job.sh` to loop over vars (NIRv, SIF) and neighborhood radii, each time prepping data layers, running the random forest analysis, and generating identical results.
5. Run `asynch/analysis/drivers/summ_results/ch3_rasterize_SHAP_job.sh` to convert output CSVs of global SHAP values to GeoTIFFs.
6. Run `asynch/analysis/drivers/summ_results/ch3_rasterize_err_job.sh` to convert output CSVs of global RF prediction errors to GeoTIFFs.
7. Run `asynch/analysis/drivers/summ_results/tabulate_model_summaries.py` to combine all permuation-based and SHAP-based importance values and model $R^2$s and MSEs into a single output table, for supplmental materials.
8. Run `asynch/analysis/drivers/make_figs/make_asynch_viz_and_analysis_maps_FIG_3_S10-18.py` 4 times, once each with a different command line arg ('main', 'asynch_supps', 'error_supp', 'predom_supp'), to produce the main RF-summary fig (Fig. XXX, which maps predominance for the two top-importance covariates, ppt.asy and tmp.min.asy) and the supplementary figs presenting all asynchrony maps (Figs. SXXX-SXXX), the RF error map (Fig. SXXX), and the RF summary map showing all predominant covariates (Fig. SXXX).


### run climate-distance analysis:
1. Run `asynch/analysis/clim_dep/compare_phen_clim_geog_dist_FIG_5.py` to run all iterations of the analysis of the latitudinal trend in the phenological distance~climatic distance relationship and produce the analysis summary in **Fig. XXX**.


### run iNaturalist flowering-phenology analysis:
1. Run `asynch/analysis/phen/get_all_inat_plant_phen_taxa.py` to save a table of the counts of all phenology-annotated and valid (i.e., native, non-captive, research-grade, with positional accuracy ≤ 1000 m) observations for all iNaturalist taxa with at least one observation. (**~5 min runtime; last run 2024-06-05T11:24:00UTC**)
2. Run `asynch/analysis/phen/get_all_inat_phen_obs_and_fit_npeaks.py` to download phenology observation histograms and raw observations for all iNat taxa with at least 50 phenology-annotated and otherwise valid observations, fit alpha hulls to the first ≤ 5000 observations (α = 0.75, the mid-value from the climate-dependence analysis above), and estimate the number of peaks (0, 1, or 2) in their flowering-week histograms (using a KDE with bandwidth = 5, a peak detection algorithm detecting all peaks ≥ 60% of the max histogram height, and significance of the number of peaks being estimated using a 100-permutation test). (**~4 day runtime; last run started about 2024-06-05T23:00:00UTC and ran until about 2024-06-09T00:00:00UTC, with a couple overnight stops to avoid excessive obstacles because of API rate-limiting; script has try/except blocks that seem to handle API rate limiting in an unsophisticated but technically passable way, but it still works best if it is stopped for some hours when 429 errors become too common, which is fine becuase it stashes results along the way and picks up wherever it left off**)
3. Run `asynch/analysis/phen/analyze_inat_flow_phen_results.py` to run an MMRR, predicting flowering observation-time distance as a function of geographic and LSP distances, for all iNat taxa with non-unimodal flowering week histograms (i.e., 0 or 2 peaks).

### run genetic analyses:
1. Download the [supplemental data in Dryad](https://datadryad.org/stash/dataset/doi:10.5061/dryad.pc866t1p4) from [Thomé et al. 2021](https://www.nature.com/articles/s41437-021-00460-7), the only genomic test of the asynchrony of seasons hypothesis (ASH) of which we are aware.
2. Run `asynch/analysis/gen/rhinella/convert_STRUCTURE_to_gendist_mat.r` to produce a genetic distance matrix for the _Rhinella granulosa_ genetic data (i.e., the Thomé et al. supplemental data).
3. Run `asynch/analysis/gen/rhinella/test_rhinella_granulosa.py` to run an MMRR for _Rhinella granulosa_, predicting genetic distance as a function of our NIRv-based phenological distance instead of the Thomé et al. precipitation-seasonality distance variable (and controlling for geographic and environmental distances).
4. Manually compile the sample locations (from the [supplemental data in Zenodo](https://zenodo.org/records/5012226)) and FASTA-format sequences (from NCBI, based on sample voucher IDs in the supplemental data) for all samples of the only speciesin [Quintero et al. 2014](https://www.journals.uchicago.edu/doi/full/10.1086/677261) (a multi-species test of the ASH using archived microsatellite data) that overlaps with the eastern Brazilian Thomé et al. study region.
5. Run `asynch/analysis/gen/xiphorhynchus/align_xiphorhynchus_fuscus.sh` to use [MAFFT v7.520](https://mafft.cbrc.jp/alignment/software/) to align all raw sequence data for the Xiphorhynchus fuscus samples (i.e., to Quintero et al. species that co-occurs in the Thomé et al. study region).
6. Run `asynch/analysis/gen/xiphorhynchus/calc_gen_dist_mat_xiphorhynchus_fuscus.r` to produce a genetic distance matrix from the aligned _Xiphorhynchus fuscus_ sequences.
7. Run `asynch/analysis/gen/xiphorhynchus/test_xiphorhynchus_fuscus.py` to run an MMRR for _Xiphorhynchus fuscus_, predicting genetic distance as a function of our NIRv-based phenological distance instead of the Quintero et al. precipitation asynchrony variable (and controlling for geographic and environmental distances).
8. **COMPILE FIGURE 5 USING WHAT SCRIPTS/TOOLS?**



-------------------------------------------



## working environments
**local**:
  - Linux laptop running Pop!\_OS 22.04 LTS (except Ubuntu 18 used for SIF data prep and evaluation)
  - Python 3.9 (except Python 3.7 used for SIF data prep and evaluation)
    - numpy 1.22.4
    - rasterio 1.2.10
    - xarray 2022.3.0
    - rioxarray 0.11.1
    - pandas 2.2.1
    - shapely 2.0.3
    - geopandas 0.14.3
    - zipfile36 0.1.3
    - tensorflow 2.4.1
    - pyproj 3.6.1
    - cartopy 0.20.2
    - geopy 1.13.0
    - scipy 1.13.0
    - sklearn 1.0.2
    - alphashape 1.3.1
    - rasterstats 0.16.0
    - statsmodels 0.13.2
    - seaborn 0.11.2
    - fuzzywuzzy 0.18.0
    - Bio 1.79
    - pyinaturalist 0.19.0
    - json 2.0.9
    - nlmpy
    - matplotlib 3.7.0
    - h3 3.7.4
    - descartes
    - contextily 1.2.0
    - xyzservices 2022.4.0
    - palettable 3.3.0
    - cmocean 2.0
    - cmcrameri 1.8
    - colormap 1.0.4
    - imageio 2.19.0
    - cv2 4.9.0 
  - R 
    - adegent 2.1.5
    - poppr 2.9.5
    - ape 5.6.2
  - Bash 5 (except Bash 4 used for SIF data prep and evaluation)
  - GDAL 3.3.1
  - MAFFT 7.520

**UC Berkeley Savio Cluster**:
  - Python 3.7
    - numpy 1.21.5
    - rasterio 1.1.5
    - xarray 0.20.2
    - rioxarray 0.9.1
    - pandas 1.3.5
    - geopandas 0.8.1
    - json 2.0.9
    - tensorflow 2.3.1
    - affine 2.3.0
    - haversine
    - eofs 1.4.0
    - scipy 1.4.1
    - sklearn 0.21.3
    - matplotlib 3.1.1
  - R **XXXXX**
    - sp  **XXXX**
    - sf **XXXX**
    - spdep **XXXX**
    - raster **XXXX**
    - terra **USE THIS??**
    - rsample  **XXXX**
    - RRF **XXXX**
    - ranger **XXXX**
    - randomForest  **USE THIS??**
    - spatialRF **USE THIS??**
    - rfUtilities  **USE THIS??**
    - h2o  **USE THIS??**
    - Boruta **XXXX**
    - fastshap **XXXX**
    - SpatialML **USE THIS??**
    - GWmodel **USE THIS??**
    - vip **XXXX**
    - pdp **XXXX**
    - DALEX **XXXX**
    - ggplot2 **XXXX**
    - ggthemes **XXXX**
    - grid **XXXX**
    - cowplot **XXXX**
    - tmap **XXXX**
    - maps **XXXX**
    - RColorBrewer **XXXX**
    - cmocean **XXXX**
    - wesanderson **XXXX**
    - dplyr **XXXX**
    - caret **USE THIS??**
  - Julia 1.4.1
    - Distributed
    - OrderedCollections 1.4.1
    - StaticArrays 1.3.5
    - Glob 1.3.0
    - TFRecord 0.1.0
    - JSON 0.21.3
    - ArchGDAL 0.7.4
    - Distances 0.10.7
    - NearestNeighbors 0.4.9
    - Statistics
    - StatsBase 0.33.16
    - GLM 1.6.1
    - Colors 0.12.8
  - Bash 4

**GEE**:
  - browser-based Javascript IDE (GEE API ≥0.1.404)

