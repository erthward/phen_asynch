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
All analysis code is located in `src`. Analyses were run in two major stages, so the two major subdirectories of `src` are organized to reflect this:

  1. `src/phen/`: calculation of global maps of characteristic seasonal phenology, and associated analysis
  2. `src/asynch/`: calculation of global maps of phenological asynchrony, and associated analysis

Additional `src` subdirectories include:
  - `src/data_prep`: code used to download the SIF dataset used, convert to Geotiff, and upload to Google Earth Enginge (GEE)
  - `src/etc`: includes a script of commonly used helper functions, a Python port of [MMRR] (http://onlinelibrary.wiley.com/doi/10.1111/evo.12134), and other odds and ends

Other directories include:
  - `data`: ancillary data (not our main datasets) that lives in this repo and is used in analysis
  - `res`: results (organized into subdirectories for figures and tables)



-------------------------------------------



## workflow:
Each step of the following workflow was executed in the environment indicated *in italics* (either on GEE, on the UC Berkeley Savio cluster, or on my laptop). See **working environments**, below, for the software specifications of each of these environments.


### prep and evaluate SIF dataset:
1. *On laptop*, download [SIF dataset](http://doi.org/10.3334/ORNLDAAC/1696)
2. *On laptop*, run `bash src/data_prep/sif_gee_ul/convert_SIF_OCO2_ANN_NetCDF_to_GTiff.sh` to translate SIF data from NetCDF to GeoTIFF.
3. *On laptop*, run `python src/data_prep/sif_gee_ul/prep_SIF_OCO2_ANN_upload_metadata.py` to create CSV of SIF metadata.
4. *On laptop*, python `src/data_prep/sif_gee_ul/upload_SIF_OCO2_ANN_data_to_GEE_collections.py` to upload SIF data to GEE as individual Images.
5. *On laptop*, run `bash src/data_prep/sif_gee_ul/make_image_collection.sh` to combine GEE SIF data into an ImageCollection.
6. *On laptop*, manually download [gridded TROPOMI data](ftp://fluo.gps.caltech.edu/data/tropomi/gridded/).
7. *On laptop*, run `python src/phen/eval/orbital_gaps/evaluate_ANN-gridded_OCO2_orbital_gaps.py` to check that the seasonality of the OCO2-SIF ANN-interpolated data products within OCO2 orbital gaps compares favorably to another gridded SIF dataset, for three regions across the tropics.


### prep jurisdictional boundaries datasets for plotting:
1. *On laptop*, manually download and unzip into the 'data' subdirectory the global level-1 administrative boundaries Shapefile at 1:10m scale from [Natural Earth](http://www.naturalearthdata.com/downloads/10m-cultural-vectors/10m-admin-1-states-provinces/).
2. *On laptop*, manuall download and unzip into the 'data' subdirectory the level-1 jurisdictional boundaries GeoJSON datasets for select large-area nations, mostly $>2.5\times10^6\ km^2$ land area (Canada, USA, Mexico, Brazil, Argentina, Kazakhstan, India, Russia, China, Australia) from [GADM](http://gadm.org/data.html).
3. *On laptop*, run `python src/data_prep/juris_bounds/simplify_global_adm1_bounds.py` to transform those files into three files containing global level-0 bounds, global level-1 bounds, and level-1 bounds for just the large-area nations (**1 task, <1m runtime**).


### calculate masking and preprocessing maps for seasonality-fitting procedure:
1. *On GEE*, establish all desired parameter values for data reading and masking, harmonic regression fitting, significance calculation, and data exports in `src/phen/calc/GEE/params.js`. (Note: To produce all needed outputs, this file will be resaved with altered commenting a few times during the following steps.)
2. *On GEE*, run `src/phen/calc/GEE/masking/calc_month_data_avail_props.js` to produce a set of 12 global Image assets, each indicating the average proportion of monthly MODIS NBAR data availability at each pixel, to be used for masking based on proportional data availability and to derive the monthly data-availability evenness mask. (**12 tasks, each ~20-30m runtime, but note that a few may fail the first time and need to be rerun until successful**)
3. *On GEE*, run `src/phen/calc/GEE/masking/calc_month_evenness.js` to produce a global map of monthly Pielou's evenness of the MODIS NBAR data, as a GEE asset, derived from the outputs of `src/phen/calc/GEE/masking/calc_month_data_avail_props.js`. (**2 tasks, both <10m runtime**)
4. *On GEE*, run `src/phen/calc/GEE/masking/calc_land_cover_mask.js` to produce a global map indicating all pixels that are valid (i.e., not water, urban, barren, or permanent ice/snow) and that are always 'natural' across our time series (i.e., forest, savanna, shrub, grass, wetland) and/or always agricultural. (All of those pixels will be used in the LSP analyses, but only the 'always natural' pixels will be used in the asynchrony analyses, to try to avoid places subject to anthropogenically created spatial asynchrony.) (**2 tasks, both <10m runtime**)
5. *On GEE*, run `src/phen/calc/GEE/masking/calc_ts_pct_data_availability.js` to produce a map of overall data availability for the entire 20-year archive. (**1 task, ~4h runtime**)
6. *On GEE*, run `src/phen/calc/GEE/masking/calc_permutation_fail_count_maps.js` to produce GEE assets representing, in total, 20 permutations of the seasonal phenology-fitting code for the NIRv data (unreasonable to produce more because of compute limitations). (**10 tasks, each ~6-20h runtime**)
7. *On GEE*, run `src/phen/calc/GEE/masking/calc_permutation_signif_mask.js` to combine the outputs of the previous step into the overall significance masks, to be used in the masking procedure for the final phenology maps. (**1 task, <10m runtime**)
8. *On GEE*, run `src/phen/calc/GEE/masking/create_overall_mask.js` twice, once with 'maskingMode' set to 'default' in `src/phen/calc/GEE/params.js` and once with it set to 'strict', to combine all previously created masks into a pair of default (used for LSP analyses) and strict (used for asynchrony analyses) mask assets. (**7 tasks total, <10m runtime**)
9. *On GEE*, run `src/phen/calc/GEE/io/calc_min_pos_NIRv.js` to calculate a global map asset of each pixel's minimum positive NIRv value (to be used to backfill negative NIRv values that occur, mainly in places and times with solid snow cover. (** 1 task, ~8h runtime**)


### calculate fitted seasonality maps:
1. *On GEE*, run `src/phen/calc/GEE/fft/create_NIRv_harmonic_regression_asset.js` to run the main harmonic regression analysis on the MODIS NIRv dataset. (This is by far the most computationally intensive, so it's best to just run it once and save the result as an asset, which is then loaded, masked, and export to Drive by `src/phen/calc/GEE/main.js`.) (**1 task, ~20h runtime**)
2. *On GEE*, run `src/phen/calc/GEE/main.js` to launch a GEE task that will save tiled results to Google Drive.
  NOTE: Results will be formatted as overlapping TFRecord tiles.
  NOTE: In order to produce all results, this script must be manually called once for each combination of NIRv or SIF and default or strict masking, as well as once per climate variable (TerraClimate min temperature, precipitation, and climatic water deficit, as well as MODIS cloud). To do this, the `datasetName`, `climateVar`, and `maskingMode` variables must be manually swapped in this script.
  NOTE: To execute the full phenology-mapping workflow, main.js will in turn call a variety of other scripts located in `src/phen/calc`.
  (**9 tasks total, ranging from 10s of minutes for the NIRv, SIF and TerraClimate tasks to a handful of hours for the MODIS cloud task**)


### download seasonality results:
1. *On both laptop and Savio*, navigate to the 'GEE\_outputs' directory where all fitted LSP and seasonality files from GEE should be stored, then run `bash src/phen/calc/dl_and_organize_GEE_data.sh` to a.) download all results into the parent directory (using `rclone`) (both on UC Berkeley's Savio cluster and locally on external hard drive), then b.) move each run's results into corresponding child directory.


### calculate asynchrony:
1. *On Savio*, run `bash src/asynch/calc/run_all_asynch_jobs.sh` to feed all three job scripts, one per asynchrony neighborhood radius, to slurm's `sbatch` command, calculating each neighborhood radius' set of asynchrony maps for all fitted phenological and climatic seasonality datasets (**3 slurm jobs total, runtimes of ~13h for the 50km neighborhood job, ~30h for 100km, and ~46h for 150km**).


### mosaic and store all results:
1. *On Savio*, run `sbatch src/asynch/calc/mosaic_job.sh` to mosaic the regression coefficient, regression $R^2$, and asynchrony-result files for all LSP and climate variables and for all three asynchrony neighborhoods (50 km, 100 km, 150 km), producing a set of GeoTIFF outputs for downstream plotting and analysis (**1 job, ~30m runtime**).
2. *On Savio*, run `bash src/asynch/calc/ul_mosaicked_results_from_savio_to_bdrive.sh` to copy all mosaicked results from Savio back up to Google Drive (**1 task, ~15m runtime**)
3. *On laptop*, run `bash src/asynch/calc/dl_mosaicked_results_from_bdrive.sh` to then also copy those results down to their intended location on laptop external hard drive (**1 task, ~1h runtime**).


### map masks and R2s from harmonic regressions:
1. *On laptop*, run `bash src/phen/calc/masking/dl_GEE_masks.sh` to download all 6 mask GeoTIFFs output by GEE into the appropriate directory on the external hard drive (**1 task, <5m runtime**).
2. *On laptop*, run `python src/phen/calc/masking/make_masking_maps_supp.py` to produce supplemental figure showing all masks used on the LSP datasets (**1 task, <5m runtime**).
3. *On laptop*, run `python src/phen/calc/R2/make_R2_maps_supp.py` to produce supplemental figure showing the $R^2$s of the harmonic regressions fitted to all LSP and climate datasets (**1 task, <5m runtime**).


### produce RGB phenology map:
1. *On Savio*, run `sbatch src/phen/anal/div/LSP_EOF_job.sh` to calculate the global NIRv LSP EOF map and save results (**1 job, <1h total runtime**).
2. *On Savio*, run `bash src/phen/anal/div/ul_EOF_results_to_bdrive.sh` to push the EOF results back to BDrive (**1 job, <5m runtime**).
3. *On laptop*, run `bash src/phen/anal/div/dl_LSP_EOFs.sh` to download the EOF results to the local hard drive (**1 job, <5m runtime**).
4. *On laptop*, download ancillary cheatgrass data from [Maestas *et. al*](http://www.sciencebase.gov/catalog/item/5ec5159482ce476925eac3b7) (to be used in a statistical test embedded in `src/phen/anal/div/plot_EOF_and_RGB_results.py`).
5. *On laptop*, run `bash src/phen/anal/div/aggregate_great_basin_cheatgrass_data.sh` to aggregate that dataset to our analysis resolution of $0.05^{\circ}$.
6. *On laptop*, work through the manual steps listed in the notes at the top of `python src/phen/anal/div/compose_ITCZ_shapefile.py` to digitize the Dec-Jan-Feb and Jun-Jul-Aug mean ITCZ locations delineated by [Zhisheng et al. 2015](annualreviews.org/content/journals/10.1146/annurev-earth-060313-054623), save the output CSVs to the 'data/' subdirectory of the local clone of this repo, then run `python src/phen/anal/div/compose_ITCZ_shapefile.py` to produce a Shapefile of those digitized lines.
7. *On laptop*, run `bash src/phen/anal/div/make_EOF_and_RGB_map_figs.sh` to produce main figures showing the global EOF RGB map and the focal-region EOF RGB maps, and the supplemental figures showing the raw EOF maps and the unfolded EOF RGB maps (**1 task, <30m runtime**).
8. *On laptop*, manually assemble each focal-region map figure's .odg file (for both the main regional figure and the supplementals) in LibreOffice Draw, using the outputs of the previous step, then save a static image (by selecting all, clicking "File" > "Export...", checking the 'Selection' checkbox in the bottom left corner, selecting .png as the file type, clicking "Save", clicking the "Modify resolution" radio-button on the pop-up and setting the resolution to 500 pixels/inch, then clicking the "Modify dimensions" radio button, setting the width to 7 inches and letting the height auto-adjust, then clicking "OK").


### map examples of LSP fitted to raw NIRv data:
1. *On GEE*, run `src/phen/calc/GEE/viz/extract_raw_NIRv_at_CA_sites.js` to extract the raw, 20-year NIRv time series at a transect of 3 California FLUXNET sites with divergent phenologies (**1 task, 1m runtime**).
2. *On laptop*, run `src/phen/viz/plot_NIRv_raw_and_fitted_at_example_sites.py` to produce a supplemental figure demonstrating the LSP time series fitted to the raw NIRv data for those three sites (**1 task, <5m runtime**).


### produce LSP modality map:
1. *On laptop*, run `python src/phen/eval/modality/map_ann_sem_LSP_modality.py` to calculate the global map of LSP modality, by searching for <=2 peaks to each pixel's fitted, minmax-scaled LSP time series, then calculating the absolute difference between the 2 peaks (or assigning a difference of 1 if only a single peak exists), such that 0s refer to places with 2 exactly equal-height fitted LSP peaks per year, 1s refer to places with only a single peak, and values in between have 2 peaks varying in relative height (**1 task, ~1h runtime**).


### produce LSP video:
1. *On laptop*, run `python src/phen/viz/make_LSP_global_video.py` to produce an mp4 video of one year of global, minmax-scaled LSP variability (**1 task, ~3h runtime**).


### run NPN and SI-x evaluation:
1. *On laptop*, run `Rscript --vanilla src/phen/eval/NPN_and_SI-x/get_NPN_leaf_data.r` to download, for a wide range of dominant US tree genera and at all NPN sites, both the day of year of first leaf based on NPN ground observations and mean day of year of start of season (SOS) based on MODIS-derived SI-x phenology maps (**1 task, <1h runtime**).
2. *On laptop*, run `python src/phen/eval/NPN_and_SI-x/compare_NIRv_LSP_to_NPN_first_leaf.py` to evaluate SOS estimates derived from our NIRv LSP data against both the NPN first-leaf and SI-x SOS datasets (**1 task, <1m runtime**).


### run NIRv-SIF LSP comparison:
1. *On Savio*, run `sbatch src/phen/eval/compare_NIRv_and_SIF_maps/ch3_phen_comparison_job.sh` to calculate a global map of $R^2$ values between the fitted annual NIRv-based and SIF-based LSP patterns (**1 task, ~12h runtime**).
2. *On Savio*, run `bash src/phen/eval/compare_NIRv_and_SIF_maps/ul_NIRv_SIF_phen_comparison_results_to_bdrive.sh` to upload the resulting GeoTIFF to Google Drive (**1 task, ~1m runtime**).
3. *On laptop*, run `bash src/phen/eval/compare_NIRv_and_SIF_maps/dl_NIRv_SIF_phen_comparison_results_from_bdrive.sh` to download the GeoTIFF to the necessary local directory (**1 task, ~1m runtime**).
4 *On laptop*, run `python src/phen/eval/compare_NIRv_and_SIF_maps/map_NIRv_SIF_LSP_comparison.py` to map the $R^2$ values between the fitted NIRv and SIF LSP datasets (**1 task, <1m runtime**).


### run FLUXNET evaluation:
1. *On laptop*, manually download all subset data products (using DownThemAll!) from the FLUXNET network's [download page](http://fluxnet.org/data/download-data/) (**1 task, runs roughly overnight**).
2. *On laptop*, run `bash src/phen/eval/flux_tower_GPP/run_flux_evaluations.sh` to run the flux-tower GPP comparison, at all usable FLUXNET2015 sites, for both the fitted NIRv and SIF LSP results (**1 task, ~15m runtime**).
4. *On laptop*, run `python src/phen/eval/plot_phen_evaluation_results.py` to combine both LSP datasets' FLUXNET evaluations and the NIRv-SIF comparison evaluation to make supplemental figure (**1 task, <5m runtime**).
5. *On laptop*, run `python src/phen/eval/flux_tower_GPP/combine_fluxnet_val_outputs.py` to combine of flux-tower evaluation results, listed by FLUXNET site, into a single supplemental table (**1 task, <10s runtime**).
6. *On laptop*, run `python src/phen/eval/flux_tower_GPP/plot_flux_eval_results.py` to generate the flux-tower evaluation supplemental figure (**1 task, <1m runtime**).


### run asynchrony neighborhood-comparison evaluation:
1. *On laptop*, run `python src/asynch/eval/calc_asynch_r2s_btwn_neighborhood_radii.py` to produce Table S2, containing R2s for all neighborhood radius comparisons and for all variables for which we produced asynchrony maps (**1 task, <1m runtime**).


### run NIRv-SIF asynchrony comparison:
1. *On laptop*, run `python src/asynch/eval/NIRv_SIF_comp/compare_SIF_and_NIRv_asynch.py` to compare the two datasets' phenological asynchrony maps across all three neighborhood radii (50 km, 100 km, 150 km) (**1 task, <5m runtime**).


### produce asynchrony supplemental figures:
1. *On laptop*, run `python src/asynch/viz/make_conceptual_asynch_fig.py` to create the asynchrony-calculation conceptual figure (**1 task, <1m runtime**).
2. *On laptop*, run `python src/asynch/viz/plot_all_asynch_maps.py` to create each LSP and climate variable's supplemental figure, displaying the asynchrony maps for all three neighborhood radii (50, 100, and 150 km) (**1 task, <30m runtime**).


### collect all covariates for asynch drivers analysis:
1. *On GEE*, run `src/phen/calc/GEE/other_datasets/calc_veg_entropy.js` to produce the vegetation cover entropy map that will be used as a covariate in the phenological asynchrony predictive model. (**1 task, <10m runtime**)
2. *On laptop*, manually download the GMTED-derived, 100-km median-aggregated vector ruggedness metric (VRM) dataset from [EarthEnv](http://www.earthenv.org/topography) (file 'vrm_100KMmd_GMTEDmd.tif'), then uploade to the 'LSP_ancillary_datasets' folder on Google Drive (**1 task, <10s runtime**)
3. *On Savio*, run `bash src/asynch/anal/drivers/prep_data/dl_asynch_drivers_covariates_from_bdrive.sh` to download both of these covariate maps into the correct directory (**1 task, <1m runtime**).
4. *On Savio*, run `bash src/asynch/anal/drivers/prep_data/cp_asynch_files_to_rf_data.sh` to copy all asynch maps into the random forest-modelling data directory (**1 task, <1m runtime**).


### run phenological asynchrony drivers analysis:
1. *On Savio*, run `bash src/asynch/anal/drivers/prep_data/prep_NIRv_100km_rf_data_for_model_tuning.sh` to prep data for random forest analysis using 100 km-neighborhood NIRv asynchrony dataset (**1 task, <5m runtime**).
2. *On Savio*, in an RStudio session hosted in an Open OnDemand session, running R version 4.0.3 on a standard `savio3` node, interactively run the code in `src/asynch/anal/drivers/run_rf/run_phen_asynch_rf.r`, with `var` set to 'NIRv', `neigh.rad` set to '100', and `coords.as.covars` set to 'y' (by switching the commenting on the corresponding lines at the top of the file), to execute the random forest analysis on the main phenological asynchrony dataset (NIRv-based LSP asynchrony calculated within 100 km radial neighborhoods) and save hyperparameter tuning results. (Be sure to interactively execute the code blocks captured by the `if (F){ ... }` blocks, which will save the hyperparameter-tuning results as CSVs and the Boruta feature selection results as a boxplot PNG.) Manually inspect the results, then use them to set the hyperparameters and variable inclusion to be used in all 12 random forest models (2 LSP datasets × 3 neighborhood radii × 2 options for including/excluding geographic coordinates as covariates) (set these in the code block starting at line 410 and 486) (**1 task, ~2h interactive runtime**).
3. *On Savio*, run `sbatch src/asynch/anal/drivers/run_rf/ch3_rf_job.sh` to loop over vars (NIRv, SIF) and neighborhood radii, each time prepping data layers, running the random forest analysis, and generating identical results (**1 task, ~18h runtime**).
4. *On Savio*, run `sbatch src/asynch/anal/drivers/summ_results/ch3_rasterize_SHAP_job.sh` to convert output CSVs of global SHAP values to GeoTIFFs (**1 task, ~7.5h runtime**).
5. *On Savio*, run `sbatch src/asynch/anal/drivers/summ_results/ch3_rasterize_err_job.sh` to convert output CSVs of global RF prediction errors to GeoTIFFs (**1 task, ~1h runtime**).
6. *On Savio*, run `bash src/asynch/anal/drivers/summ_results/ul_asynch_drivers_results_to_bdrive.sh` to upload all drivers-modeling results to Google Drive (**1 task, <1 h runtime**).
7. *On laptop*, run `bash src/asynch/anal/drivers/summ_results/dl_asynch_drivers_results_from_bdrive.sh` to download all drivers-modeling results into the appropriate directory on the external hard drive (**1 task, <1h runtime**).
8. *On laptop*, run `python src/asynch/anal/drivers/summ_results/tabulate_model_summaries.py` to combine all permuation-based and SHAP-based importance values and model $R^2$s and MSEs into a single output table, for supplmental materials (**1 task, <1m runtime**).
9. *On laptop*, manually copy the contents of each of the three sheets in both of the .xlsx files created by the last step into their corresponding sections in `figs_and_tabs/tabs/TAB_SUPP_model_summary_table_full.xlsx` using LibreOffice Calc, save the file, then save a static image of the table (by manually selecting the entire gray boxed area, clicking "File" > "Export...", checking the 'Selection' checkbox in the bottom left corner, selecting .png as the file type, clicking "Save", clicking the "Modify resolution" radio-button on the pop-up and setting the resolution to 500 pixels/inch, then clicking the "Modify dimensions" radio button, setting the width to 10 inches and letting the height auto-adjust, then clicking "OK").
10. *On laptop*, run `bash src/asynch/anal/drivers/make_figs/make_all_asynch_analysis_figs.sh` to produce the main asynch figure (global map, as well as map summarizing the predominance of the two top-importance covariates), as well as the supplemental figures showing the random forest error map and the map of SHAP-value predominance across all random forest covariates (**1 task, <10m runtime**).


### run analysis of isoclimatic phenological asynchrony:
1. *On laptop*, run `python src/asynch/anal/isoclim/compare_phen_clim_geog_dist.py` to run all iterations of the analysis of the latitudinal trend in the phenological distance~climatic distance relationship and produce the analysis summary figure (NOTE: analysis saves results for each paramater-combination as it runs, so can be stopped and restarted if necessary and will pick up where it left off.) (**1 task, ~26h runtime**).


### run iNaturalist flowering-phenology analysis:
1. *On laptop*, run `python src/asynch/anal/phen/get_all_inat_plant_phen_taxa.py` to save a table of the counts of all phenology-annotated and valid (i.e., native, non-captive, research-grade, with positional accuracy ≤ 1000 m) observations for all iNaturalist taxa with at least one observation. (**~5 min runtime; last run 2024-06-05T11:24:00UTC**)
2. *On laptop*, run `python src/asynch/anal/phen/get_all_inat_phen_obs_and_fit_npeaks.py` to download phenology observation histograms and raw observations for all iNat taxa with at least 50 phenology-annotated and otherwise valid observations, fit alpha hulls to the first ≤ 5000 observations (α = 0.75, the mid-value from the climate-dependence analysis above), and estimate the number of peaks (0, 1, or 2) in their flowering-week histograms (using a KDE with bandwidth = 5, a peak detection algorithm detecting all peaks ≥ 60% of the max histogram height, and significance of the number of peaks being estimated using a 100-permutation test). (**~4 day runtime; run started about 2024-06-05T23:00:00UTC and ran until about 2024-06-09T00:00:00UTC, with a couple overnight stops to avoid excessive obstacles because of API rate-limiting; script has try/except blocks that seem to handle API rate limiting in an unsophisticated but technically passable way, but it still works best if it is stopped for some hours when 429 errors become too common, which is fine becuase it stashes results along the way and picks up wherever it left off**)
3. *On laptop*, run `python src/asynch/anal/phen/fit_inat_flow_phen_LSP_MMRR_models.py` to run an MMRR, predicting flowering observation-time distance as a function of geographic and LSP distances, for all iNat taxa with non-unimodal flowering-week histograms (i.e., 0 or 2 peaks) (**1 task, ~10h runtime**).

### run genetic analyses:
1. *On laptop*, download the [supplemental data in Dryad](http://datadryad.org/stash/dataset/doi:10.5061/dryad.pc866t1p4) from [Thomé et al. 2021](http://www.nature.com/articles/s41437-021-00460-7), the only genomic test of the asynchrony of seasons hypothesis (ASH).
2. *On laptop*, run `Rscript --vanilla src/asynch/anal/gen/rhinella/convert_STRUCTURE_to_gendist_mat.r` to produce a genetic distance matrix for the _Rhinella granulosa_ genetic data (i.e., the Thomé et al. supplemental data).
3. *On laptop*, manually compile the sample locations (from the [supplemental data in Zenodo](http://zenodo.org/records/5012226)) and FASTA-format sequences (from NCBI, based on sample voucher IDs in the supplemental data) for all samples of *Xiphorhychus fuscus*, the only species in [Quintero et al. 2014](http://www.journals.uchicago.edu/doi/full/10.1086/677261) (a multi-species test of the ASH using archived microsatellite data) that is sympatric with the *Rhinella granulosa*, the eastern Brazilian toad studied by Thomé et al.
4. *On laptop*, run `bash src/asynch/anal/gen/xiphorhynchus/align_xiphorhynchus_fuscus.sh` to use [MAFFT v7.520](http://mafft.cbrc.jp/alignment/software/) to align all raw *X. fuscus* sequence data.
5. *On laptop*, run `Rscript --vanilla src/asynch/anal/gen/xiphorhynchus/calc_gen_dist_mat_xiphorhynchus_fuscus.r` to produce a genetic distance matrix from the aligned _Xiphorhynchus fuscus_ sequences.
6. *On laptop*, run `python src/asynch/anal/plot_flowphen_and_landgen_results.py` to execute both landscape genetic analyses (*R. granulosa* and *X. fuscus*) and then visualize results and example taxa from both the iNaturalist flowering phenology peak-fitting and MMRR analyses and the landscape genetic analyses (**1 task, ~15m runtime**).


-------------------------------------------



## working environments
- **laptop**:
  - Pop!\_OS 22.04 LTS (*except Ubuntu 18 was used for SIF data prep and evaluation*)
  - Bash 5.1.16 (*except Bash 4 was used for SIF data prep and evaluation*)
  - GDAL 3.3.1
  - Python 3.9 (*except Python 3.7 was used for SIF data prep and evaluation*)
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
    - ape 5.6.2
  - MAFFT 7.520
  - LibreOffice 7.3.7.2 30(Build:2)
- **UC Berkeley Savio Cluster**:
  - Scientific Linux 7.9
  - Bash 4.2.26
  - GDAL 2.2.3
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
    - scikit-learn 0.21.3
    - matplotlib 3.1.1
  - R 4.0.3
    - dplyr 1.0.8
    - rgdal 1.5.18
    - sp 1.4.6
    - sf 0.9.7
    - raster 3.4.5
    - maps 3.4.0
    - rsample  0.1.1
    - ranger 0.13.1
    - Boruta 7.0.0
    - fastshap 0.0.7
    - vip 0.3.2
    - pdp 0.7.0
    - ggplot2 3.3.5
    - ggthemes 4.2.4
    - grid 4.0.3
    - RColorBrewer 1.1.2
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
- **GEE Javascript API (browser-hosted)**:
  - GEE API ≥0.1.404

