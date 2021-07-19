

wind_stats_plot <- function(annotated){
        
        annotate <- ggplot2::annotate
        
        flw <- annotated %>%
                filter(wind_p == "wind_p1") %>%
                filter(to != from,
                       is.finite(clim),
                       is.finite(distance),
                       is.finite(wind)) %>%
                mutate(wind = 1/wind,
                       distance = distance / 1000,
                       speed = distance/wind,
                       group = paste(component, edge)) %>%
                group_by(group) %>%
                mutate(mean_wind = mean(wind),
                       mean_speed = mean(speed),
                       speed_ratio = max(speed)/min(speed)) %>%
                ungroup()
        
        mean(flw$speed_ratio)
        flw %>% group_by(component) %>%
                summarize(speed = max(mean_speed)/min(mean_speed)) %>%
                pull(speed) %>% median()
        
        f <- flw %>%
                #filter(component == "17d")
                filter(group %in% sample(flw$group, 500))
        
        diag <- tibble(x = c(10, 10000, 10, 10000, 10, 10000),
                       y = c(3, 3000, 10, 10000, 30, 30000),
                       group = c(1, 1, 2, 2, 3, 3))
        
        p1 <- f %>%
                ggplot(aes(distance, wind, group = group)) +
                geom_smooth(data = diag, aes(x, y), method = lm, se = F, #fullrange = T, 
                            color = "black", size = .5) +
                geom_path(aes(color = speed_ratio)) +
                geom_point(aes(y = mean_wind), size = .5) +
                annotate(geom = "text", x = 10, y = c(3, 10, 30)*1.4, label = c("3", "1", "0.3"),
                         angle = 37, hjust = 0) +
                annotate(geom = "text", x = 10, y = 10, label = "speed (km/h)",
                         angle = 90, hjust = .5, vjust = -.75) +
                scale_fill_gradientn(colors = c("white", "black"), trans = "log10") +
                scale_color_gradientn(colors = c("blue", "red", "orange"),
                                      breaks = c(1, 2, 5, 10), limits = c(1, NA), trans = "log10") +
                guides(color = guide_colorbar(order = 0)) +
                guides(fill = guide_colorbar(order = 1)) +
                scale_x_log10(limits = c(7, NA)) + scale_y_log10() +
                theme_bw() +
                labs(fill = "mean speed (km/h)",
                     color = "speed ratio",
                     x = "distance (km)",
                     y = "wind travel time (h)")
        
        p2 <- f %>% 
                ggplot(aes(speed, speed_ratio, group = group, color = distance)) +
                geom_path() +
                geom_point(aes(x = distance/mean_wind), size = .5, color = "black") +
                scale_color_gradientn(colors = c("blue", "red", "orange"), trans = "log10") +
                scale_x_log10() + scale_y_log10(breaks = c(1, 2, 5, 10), limits = c(1, NA)) +
                theme_bw() +
                labs(x = "speed (km/h)",
                     y = "speed ratio",
                     color = "distance (km)")
        
        p <- p1 + p2 & theme(legend.position = c(.15, .77),
                             legend.background = element_blank())
        
        ggs("figures/manuscript/figure_3.pdf", p, width = 10, height = 5, units = "in",
            add = grid.text(paste0("(", letters[1:2], ")"), 
                            x=c(.04, .53), 
                            y=c(.95),
                            gp=gpar(fontsize=14, fontface="bold", col="black")))
}