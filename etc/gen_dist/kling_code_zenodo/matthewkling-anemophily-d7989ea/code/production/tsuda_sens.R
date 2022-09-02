
# for the tsuda et al 2017 dataset:
# 1: recalculate flow, isolation, asymmetry hypotheses with alternative statistics
# 2: recalculate all hypotheses with genetic controls

tsuda_sens <- function(){
        
        loadd(annotated)
        
        a <- annotated %>% 
                select(-(sperm:order), -(use:drive_resource)) %>%
                filter(component == "17d")
        
        
        #### analysis 1 ####  covariates = c("wind", "distance", "clim")
        
        s13 <- function(x) x[[1]]$s13
        
        # flow
        h1_dmig <- a %>% data_h1("dmig") %>% mantel_h1() %>% s13() %>% 
                mutate(h = "flow", stat = "djost", control = "-")
        h1_gmig <- a %>% data_h1("gmig") %>% mantel_h1() %>% s13() %>% 
                mutate(h = "flow", stat = "gst", control = "-")
        h1_dmig_c <- a %>% data_h1("dmig") %>% mantel_h1(covariates = c("wind", "distance", "clim", "div_ar")) %>% s13() %>% 
                mutate(h = "flow", stat = "djost", control = "diversity")
        h1_gmig_c <- a %>% data_h1("gmig") %>% mantel_h1(covariates = c("wind", "distance", "clim", "div_ar")) %>% s13() %>% 
                mutate(h = "flow", stat = "gst", control = "diversity")
        
        # asymmetry
        h2_dmig <- a %>% data_h2("dmig") %>% mantel_h2() %>% s13() %>% 
                mutate(h = "asymmetry", stat = "djost", control = "-")
        h2_gmig <- a %>% data_h2("gmig") %>% mantel_h2() %>% s13() %>% 
                mutate(h = "asymmetry", stat = "gst", control = "-")
        h2_dmig_c <- a %>% data_h2("dmig") %>% mantel_h2(covariates = c("wind", "div_ar")) %>% s13() %>% 
                mutate(h = "asymmetry", stat = "djost", control = "diversity")
        h2_gmig_c <- a %>% data_h2("gmig") %>% mantel_h2(covariates = c("wind", "div_ar")) %>% s13() %>% 
                mutate(h = "asymmetry", stat = "gst", control = "diversity")
        
        # isolation
        h3_djost <- a %>% data_h3("djost") %>% mantel_h3() %>% s13() %>% 
                mutate(h = "isolation", stat = "djost", control = "-")
        h3_gst <- a %>% data_h3("gst") %>% mantel_h3() %>% s13() %>% 
                mutate(h = "isolation", stat = "gst", control = "-")
        h3_fst <- a %>% data_h3("fst") %>% mantel_h3() %>% s13() %>% 
                mutate(h = "isolation", stat = "fst", control = "-")
        h3_djost_c <- a %>% data_h3("djost") %>% mantel_h3(covariates = c("wind", "distance", "clim", "div_ar")) %>% s13() %>% 
                mutate(h = "isolation", stat = "djost", control = "diversity")
        h3_gst_c <- a %>% data_h3("gst") %>% mantel_h3(covariates = c("wind", "distance", "clim", "div_ar")) %>% s13() %>% 
                mutate(h = "isolation", stat = "gst", control = "diversity")
        h3_fst_c <- a %>% data_h3("fst") %>% mantel_h3(covariates = c("wind", "distance", "clim", "div_ar")) %>% s13() %>% 
                mutate(h = "isolation", stat = "fst", control = "diversity")
        
        # diversity
        h4_ar <- a %>% data_h4() %>% mantel_h4(covariates = c("wind", "lat_diff")) %>% s13() %>% 
                mutate(h = "diversity", stat = "ar", control = "-")
        h4_ar_c <- a %>% data_h4() %>% mantel_h4(covariates = c("wind", "lat_diff", "mig")) %>% s13() %>% 
                mutate(h = "diversity", stat = "ar", control = "asymmetry")
        
        
        d <- bind_rows(h1_dmig, h1_gmig, h1_dmig_c, h1_gmig_c,
                       h2_dmig, h2_gmig, h2_dmig_c, h2_gmig_c,
                       h3_djost, h3_gst, h3_fst, h3_djost_c, h3_gst_c, h3_fst_c,
                       h4_ar, h4_ar_c)
        y <- d %>%
                filter(var == "wind") %>%
                mutate(stat = case_when(stat == "djost" ~ "D", 
                                        stat == "fst" ~ "Fst",
                                        stat == "gst" ~ "Gst",
                                        stat == "ar" ~ "AR")) %>%
                mutate(h = factor(h, levels = c("flow", "isolation", "asymmetry", "diversity"))) %>%
                group_by(h, stat, control) %>%
                summarize(r = spc[i == 0],
                          p = 1 - ecdf(spc[i != 0])(spc[i == 0])) %>%
                ungroup() %>%
                rename(metric = h) %>%
                mutate(r = round(r, 3))
        
        y %>% write_csv("data_github/Tsuda_2017_sensitivity.csv")
        
}
