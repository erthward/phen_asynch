


climate_files <- function(){
        c("F:/chelsa/bio19/CHELSA_bio10_5.tif",
          "F:/chelsa/bio19/CHELSA_bio10_6.tif",
          "F:/chelsa/derived/CHELSA_annual_CWD_1979_2013.tif",
          "F:/chelsa/derived/CHELSA_annual_AET_1979_2013.tif")
}

climate_pca <- function(files){
        
        files %>%
                stack() %>%
                extract(data.frame(x = runif(30000, -180, 180),
                                   y = runif(30000, -90, 90))) %>%
                na.omit() %>%
                as.data.frame() %>%
                mutate_at(3:4, function(x) log10(x + 1)) %>%
                prcomp(center = TRUE, scale. = TRUE)
}


climate_dist <- function(pops, files, pca){
        
        d <- stack(files) %>%
                extract(pops) %>%
                as.data.frame() %>%
                mutate_at(3:4, function(x) log10(x + 1)) %>%
                predict(pca, .)
        for(i in 1:ncol(d)) d[,i] <- d[,i] / pca$sdev[i]
        
        d %>% dist(diag = TRUE, upper = T) %>% as.matrix()
}






