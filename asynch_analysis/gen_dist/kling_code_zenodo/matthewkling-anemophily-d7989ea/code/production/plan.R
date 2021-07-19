
plan <- drake_plan(
        
        ##### import and validate pop gen datasets #####
        
        # triggers
        refresh_metadata = F,
        refresh_datasets = F,
        refresh_traits = F,
        
        # get metadata from google sheet
        metadata = target(get_metadata(),
                          trigger = trigger(condition = refresh_metadata)),
        
        # import datasets from google drive
        datasets_raw = target(import_dataset(metadata),
                              dynamic = map(metadata),
                              trigger = trigger(condition = refresh_datasets)),
        
        # merge populations in the same grid cell
        datasets_merged = target(merge_pops(datasets_raw, 
                                            grid = rotate(raster(list.files("data/wind/roses", 
                                                                            full.names = T)[1]))),
                                 dynamic = map(datasets_raw),
                                 trigger = trigger(condition = F)),
        
        # identify and remove datasets with structural problems
        datasets_validated = target(validate_dataset(datasets_merged),
                                    dynamic = map(datasets_merged),
                                    trigger = trigger(condition = F)),
        datasets = filter(datasets_validated, valid),
        
        
        ##### prep intermediate datasets #####
        
        # format data on population localities
        pops = target(read_csv(datasets$pop) %>% 
                              mutate(component = datasets$component) %>% 
                              arrange(ID), 
                      dynamic = map(datasets)),
        pops_sp = target(spatialize(pops), 
                         dynamic = map(pops)),
        
        # import species traits from google sheet
        traits = target(assemble_traits(metadata),
                        trigger = trigger(condition = refresh_traits)),
        
        # compile windrose rasters for each unique dispersal season
        windroses = target(make_windroses(traits), 
                           trigger = trigger(condition = F)),
        
        # create wind connectivity graphs
        wind_graphs = target(make_wind_graph(windroses),
                             dynamic = map(windroses)),
        
        # prep climate data used for IBE
        clim_files = climate_files(),
        clim_pc = climate_pca(clim_files),
        
        
        ##### compute pairwise metrics for each population pair #####
        
        # genetic stats
        div = target(divBasic(datasets$gen),
                        dynamic = map(datasets)),
        pw_mig = target(haploMigrate(datasets$gen, boots=0),
                        dynamic = map(datasets)),
        pw_nmig = target(nullMigrate(datasets$gen, boots=1000),
                         dynamic = map(datasets)),
        pw_fst = target(fst(datasets$gen, metadata),
                        dynamic = map(datasets)),
        pw_diff = target(haploDiff(datasets$gen, pairwise = TRUE),
                        dynamic = map(datasets)),
        
        # distance and environment
        pw_distance = target(distm(pops_sp), 
                             dynamic = map(pops_sp)),
        pw_clim = target(climate_dist(pops_sp, clim_files, clim_pc), 
                         dynamic = map(pops_sp)),
        
        # wind
        pw_wind = target(wind_times(pops_sp, wind_graphs, traits),
                         dynamic = map(pops_sp)),
        pw_wind_ann = target(wind_times(pops_sp, wind_graphs, traits, annual = TRUE),
                         dynamic = map(pops_sp)),
        
        
        ##### collate and aggregate pairwise datasets #####
        
        collated = target(collate(datasets, pops, div,
                                  pw = list(distance = pw_distance, clim = pw_clim, wind = pw_wind, windann = pw_wind_ann,
                                            gst = pw_diff$pairwise$gst, Gst = pw_diff$pairwise$Gst, djost = pw_diff$pairwise$D, fst = pw_fst,
                                            dmig = pw_mig$dRelMig, gmig = pw_mig$gRelMig, nmigp = pw_nmig$dRelMigP)),
                          dynamic = map(datasets, pops, 
                                        pw_distance, pw_clim, pw_wind, pw_wind_ann, 
                                        pw_diff, pw_mig, pw_nmig, pw_fst)),
        unified = bind_rows(collated) %>% 
                left_join(bind_rows(metadata)) %>%
                left_join(traits) %>%
                gather(wind_p, wind, wind_p0:wind_p3) %>%
                mutate(wind = 1 / wind),
        
        # attrition = track_attrition(datasets_validated, datasets, collated, unified),
        
        
        ##### test hypotheses #####
        
        # add additional metadata, useful for summary analyses
        taxonomy = target(get_taxonomy(unified),
                          trigger = trigger(condition = F)),
        annotated = annotate(unified, taxonomy),
        
        # reformat data to address each hypothesis
        h1_data = data_h1(annotated, stat = "dmig"), # flow
        h2_data = data_h2(annotated, stat = "dmig"), # asymmetry
        h3_data = data_h3(annotated, stat = "fst"), # isolation
        h4_data = data_h4(annotated), # diversity
        
        # mantel randomizations
        nreps = 1000,
        h1_rand = mantel_h1(h1_data, nreps),
        h2_rand = mantel_h2(h2_data, nreps),
        h3_rand = mantel_h3(h3_data, nreps),
        h4_rand = mantel_h4(h4_data, nreps),
        
        # figures
        # must be done outside make() due to self-invalidation error triggered by ggs()
        # build_plots(), 
        
        # export data for manuscript
        pgs = datasets$gen %>% map_df(pg_stats),
        egs = export_global_summaries(annotated, pgs),
        emd = export_study_metadata(annotated),
        etd = export_trait_data(annotated),
        esg = export_search_genera(annotated),
        exp = export_combined_xls(etd, esg, emd)
)