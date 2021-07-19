


# this has to be run manually, because ggs returns a weird self-invalidation error I can't solve
build_plots <- function(){
        
        loadd(h3_data, h4_rand, h1_rand, h2_rand, h3_rand, h3_rand_nowind, annotated, nreps, traits, wind_graphs)
        
        fig1 = methods_figure(annotated, traits, wind_graphs, h1_rand, h2_rand, h3_rand, h4_rand)
        
        fig2 = syndrome_diagrams(traits)
        
        fig3 = wind_stats_plot(annotated)
        
        fig4 = ibw_boxplots(h4_rand, h1_rand, h2_rand, h3_rand, annotated, nreps)
}



### FIGURE 4 ###
ibw_boxplots <- function(h4_rand, h1_rand, h2_rand, h3_rand, annotated, nreps, pgs){
        #loadd(h4_rand, h1_rand, h2_rand, h3_rand, annotated, nreps, pgs)
        
        # s13
        rand <- list(h1_rand, h2_rand, h3_rand, h4_rand) %>%
                map(function(x) map_df(x, "s13")) 
        for(i in 1:4) rand[[i]]$h <- c("flow", "asymmetry", "isolation", "diversity")[i] 
        rand <- bind_rows(rand)
        rand <- filter(rand, wind_p == "wind_p1")
        
        
        # s2
        rand2 <- list(h1_rand, h2_rand, h3_rand, h4_rand) %>%
                map(function(x) map_df(x, "s2")) 
        for(i in 1:4) rand2[[i]]$h <- c("flow", "asymmetry", "isolation", "diversity")[i] 
        rand2 <- bind_rows(rand2)
        rand2 <- filter(rand2, wind_p == "wind_p1")
        
        
        n <- annotated %>%
                group_by(component) %>%
                summarize(n = length(unique(from)))
        
        
        specs <- rand %>% select(h, mig_var) %>% distinct()
        specs <- expand_grid(flow = specs$mig_var[specs$h == "flow"],
                             asymmetry = specs$mig_var[specs$h == "asymmetry"],
                             isolation = specs$mig_var[specs$h == "isolation"],
                             diversity = specs$mig_var[specs$h == "diversity"],
                             method = unique(rand$method),
                             dist_filter = unique(rand$dist_filter))
        
        spec <- specs %>%
                filter(flow == "mig",
                       asymmetry == "mig",
                       isolation == "sim",
                       diversity == "div_ar",
                       method == "pearson",
                       dist_filter == "alldist") %>%
                gather(h, mig_var, flow:diversity)
        
        sum_fun <- median
        
        s <- left_join(spec, rand) %>% 
                filter(var %in% c("wind", "swind")) %>%
                mutate(h = factor(h, levels = c("flow", "isolation", "asymmetry", "diversity")),
                       syndrome = expected_wind) %>%
                filter(!is.na(component)) %>%
                group_by(h, syndrome, component) %>%
                mutate(n = n()) %>% filter(n == nreps+1)
        
        ig <- left_join(spec, rand2) %>%
                select(h, dataset_id:cp) %>%
                gather(var, spc, wind:cp) %>%
                mutate(h = factor(h, levels = c("flow", "isolation", "asymmetry", "diversity"))) %>%
                mutate(syndrome = "1~0") %>%
                filter(var == "wind") %>%
                mutate(component = dataset_id) %>%
                mutate(group = "genome control") %>%
                filter(is.finite(spc))
        
        spr <- s %>% 
                ungroup() %>%
                mutate(group = "functional groups") %>% 
                mutate(syndrome = as.character(syndrome)) %>%
                bind_rows(ig) %>% 
                mutate(syndrome = factor(syndrome, levels = c("0", "0.5", "1", "1~0"))) %>%
                group_by(h, component, syndrome) %>%
                summarize(p = 1 - ecdf(spc[i != 0])(spc[i == 0]),
                          spc = spc[i == 0]) %>%
                ungroup() %>%
                gather(stat, value, spc, p) %>%
                mutate(x = case_when(syndrome == "0" ~ 0,
                                     syndrome == "0.5" ~ .5,
                                     syndrome == "1" ~ 1,
                                     syndrome == "1~0" ~ 1.5))
        
        
        
        fun <- median
        
        # significance of correlations, assessed by comparison to Mantel nulls
        sig_ig <- ig %>%
                group_by(h, i) %>%
                summarize(spc = fun(spc)) %>%
                group_by(h) %>%
                summarize(test = 3,
                          p = 1 - ecdf(spc[i != 0])(spc[i == 0]),
                          spc = spc[i == 0])
        sig_single <- s %>%
                group_by(h, syndrome, i) %>%
                summarize(spc = fun(spc)) %>%
                group_by(h, syndrome) %>%
                summarize(p = 1 - ecdf(spc[i != 0])(spc[i == 0]),
                          spc = spc[i == 0]) %>%
                ungroup() %>%
                mutate(syndrome = as.character(as.integer(factor(syndrome))))
        
        sig_diff <- s %>% ungroup() %>%
                mutate(syndrome = syndrome > 0) %>%
                group_by(h, syndrome, i) %>%
                summarize(spc = fun(spc)) %>% ungroup() %>%
                mutate(syndrome = paste0("w", syndrome)) %>%
                spread(syndrome, spc) %>%
                mutate(spc = wTRUE - wFALSE) %>% 
                group_by(h) %>%
                summarize(test = 2,
                          p = 1 - ecdf(spc[i != 0])(spc[i == 0]),
                          spc = spc[i == 0])
        sig_trend <- s %>%
                group_by(h, i) %>%
                summarize(r = cor(spc, syndrome, method = "spearman")) %>%
                group_by(h) %>%
                summarize(test = 1,
                          p = 1 - ecdf(r[i != 0])(r[i == 0]),
                          spc = r[i == 0])
        
        
        
        # significance of p-value distributions, assessed as described below
        alpha <- .5
        # univariate distirbution: assessed via binomial test
        psig_ig <- ig %>%
                mutate(test = 3) %>%
                group_by(component, h, test) %>%
                summarize(p = 1 - ecdf(spc[i != 0])(spc[i == 0]),
                          spc = spc[i == 0]) %>%
                group_by(h, test) %>%
                summarize(x = sum(p < alpha),
                          n = n(),
                          p = binom.test(x = x, n = n, p = alpha, alternative = "greater")$p.value,
                          spc = median(spc)) %>%
                ungroup() %>% select(-x, -n)
        pcomp <- s %>% 
                group_by(h, component, syndrome) %>%
                summarize(p = 1 - ecdf(spc[i != 0])(spc[i == 0]),
                          spc = spc[i == 0]) %>%
                ungroup()
        psig_diff <- pcomp %>% 
                mutate(syndrome = syndrome > 0) %>%
                group_by(h) %>%
                summarize(test = 2,
                          p = fisher.test(table(syndrome, p < .5), alternative="greater")$p.value)
        psig_trend <- pcomp %>%
                group_by(h) %>%
                summarize(test = 1,
                          p = cor.test(syndrome, p, alternative = "less", method = "spearman")$p.value) %>%
                mutate(test = 1)
        
        
        
        
        ####################### plots -- horizontal ########################
        
        egsig <- tibble(test = 1:3,
                        text = c("test 1:\ncorrelation increases\nacross wind\ndispersal levels",
                                 "test 2:\ncorrelation is\npositive for\nwind dispersers",
                                 "test 3:\ncorrelation is\npositive with\ngenome control"),
                        x = c(0, .75, 1.5),
                        y = c(0, 0, 0) - .33)
        
        rsig <- bind_rows(sig_ig, sig_diff, sig_trend) %>% left_join(egsig) %>% mutate(stat = "spc")
        psig <- bind_rows(psig_ig, psig_diff, psig_trend) %>% left_join(egsig) %>% mutate(stat = "p")
        
        # plot for correlations
        sp <- spr %>% filter(stat == "spc") %>% mutate(stat = "correlation")
        means <- sp %>% group_by(h, x, syndrome) %>% summarize(value = mean(value))
        sig <- rsig %>% mutate(stat = "correlation")
        pr <- ggplot() +
                facet_grid(stat ~ h) +
                geom_hline(yintercept = 0, color = "gray") +
                geom_boxplot(data = sp, aes(x, value, group = syndrome, color = syndrome, fill = syndrome), alpha = .15) +
                geom_point(data = sp, aes(x, value, color = syndrome), alpha = .4) +
                
                # means
                geom_point(data = means, aes(x, value), size = 3, color = "black", shape = 17) +
                geom_smooth(data = filter(means, syndrome != "1~0"), aes(x, value), method = lm, se = F, color = "black") +
                
                geom_vline(xintercept = 1.25, size = 1, color = "gray", linetype = "dotted") +
                
                scale_color_manual(values = c("red", "purple", "dodgerblue", "cyan")) +
                scale_fill_manual(values = c("red", "purple", "dodgerblue", "cyan")) +
                scale_shape_manual(values = c(21, 19)) +
                scale_y_continuous(breaks = c(-.2, 0, .2), position = "left") +
                scale_x_continuous(breaks = c(0, .5, 1, 1.5), labels = c("none", "partial", "full", "partial\n~ none")) +
                theme_bw(base_size = 18) +
                theme(legend.position = "none",
                      panel.grid = element_blank(),
                      axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                      strip.background = element_rect(fill = 'black', color = "black"),
                      strip.text = element_text(color = "white")) +
                coord_cartesian(ylim = c(-.3, .3)) +
                labs(x = "wind dispersal                       ",
                     y = "partial correlation between\nwind and genetic pattern")
        
        # plot for p-values
        sp <- spr %>% filter(stat == "p") %>% mutate(stat = "significance")
        means <- sp %>% group_by(h, x, syndrome) %>% summarize(value = mean(value))
        sig <- psig %>% mutate(stat = "significance")
        pp <- ggplot() +
                facet_grid(stat ~ h) +
                geom_hline(yintercept = 0.5, color = "gray") +
                geom_boxplot(data = sp, aes(x, 1-value, group = syndrome, color = syndrome, fill = syndrome), alpha = .15) +
                geom_point(data = sp, aes(x, 1-value, color = syndrome), alpha = .4) +
                
                # means
                geom_point(data = means, aes(x, 1-value), size = 3, color = "black", shape = 17) +
                geom_smooth(data = filter(means, syndrome != "1~0"), aes(x, 1-value), method = lm, se = F, color = "black") +
                
                geom_vline(xintercept = 1.25, size = 1, color = "gray", linetype = "dotted") +
                
                scale_color_manual(values = c("red", "purple", "dodgerblue", "cyan")) +
                scale_fill_manual(values = c("red", "purple", "dodgerblue", "cyan")) +
                scale_shape_manual(values = c(21, 19)) +
                scale_y_continuous(breaks = c(0, .5, 1), position = "left") +
                scale_x_continuous(breaks = c(0, .5, 1, 1.5), labels = c("none", "partial", "full", "partial\n~ none")) +
                theme_bw(base_size = 18) +
                theme(legend.position = "none",
                      panel.grid = element_blank(),
                      strip.background.x = element_blank(),
                      strip.background = element_rect(fill = 'black', color = "black"),
                      strip.text.x = element_blank(),
                      strip.text = element_text(color = "white"),
                      axis.text.x = element_text(angle = 45, hjust = 1, lineheight = .7)) +
                coord_cartesian(ylim = c(0, 1)) +
                labs(x = "wind dispersal level",
                     y = "partial Mantel\nsignificance")
        
        
        
        
        sig <- bind_rows(rsig %>% mutate(stat = "p-value", test = "global"),
                         psig %>% mutate(stat = "p-value", test = "dataset"))
        
        lbls <- c("wind\ndispersers", "syndrome\ncomparison", "genome\ncontrol")
        
        pz <- ggplot() +
                facet_grid(stat ~ h) +
                geom_text(data = sig, aes(x, test, label = round(p, 3), color = p < .1), 
                          fontface = "bold", size=5) +
                scale_x_continuous(breaks = c(0, .75, 1.5), labels = lbls, 
                                   limits = c(-.25, 1.75)) +
                scale_color_manual(values = c("gray50", "black")) +
                theme_bw(base_size = 18) +
                theme(legend.position = "none",
                      panel.grid = element_blank(),
                      panel.background = element_rect(fill = "gray85"),
                      strip.background.x = element_blank(),
                      strip.background = element_rect(fill = 'black', color = "black"),
                      strip.text.x = element_blank(),
                      strip.text = element_text(color = "white"),
                      axis.text = element_text(angle = 45, hjust = 1, lineheight = .7)) +
                labs(y = "null\ntest",
                     x = "prediction")
        
        p <- pr + pp + pz + plot_layout(nrow = 3, heights = c(1, 1, .3))
        
        ggs("figures/manuscript/figure_4.pdf",
            p, width = 12, height = 10, units = "in",
            add = grid.text(paste0("(", letters[1:3], ")"), 
                            x=c(.03), 
                            y=c(.98, .62, .25),
                            gp=gpar(fontsize=20, fontface="bold", col="black")))
        
        
        e <- spr %>%
                select(-x) %>%
                mutate(stat = ifelse(stat == "spc", "r", stat)) %>%
                unite(stat, h, stat) %>%
                spread(stat, value) %>%
                left_join(annotated %>% select(species, component, genome) %>% distinct()) %>%
                mutate(syndrome = as.character(syndrome),
                       syndrome = ifelse(syndrome == "1~0", "0.5~0", syndrome)) %>%
                arrange(syndrome, species) %>%
                left_join(pgs)
        e %>%
                select(component, species:npops, syndrome:isolation_r) %>%
                write_csv("data_github/mantel_results.csv")
        
        
        ############################# SNP-SSR comparison ##############################
        
        st <- annotated %>% select(component, sequence_type) %>% distinct()
        
        pd <- e %>%
                filter(syndrome != "0.5~0") %>%
                gather(stat, value, asymmetry_p:isolation_r) %>%
                separate(stat, c("h", "stat")) %>%
                mutate(fsyndrome = factor(syndrome, levels = c("0", "0.5", "1"))) %>%
                mutate(h = factor(h, levels = c("flow", "isolation", "asymmetry", "diversity"))) %>%
                left_join(st)
        
        # plot for correlations
        sp <- pd %>% filter(fsyndrome != "0") %>% filter(stat == "r") %>% mutate(stat = "correlation")
        means <- sp %>% group_by(h, sequence_type) %>% summarize(value = mean(value, na.rm = T))
        pr <- ggplot() +
                facet_grid(stat ~ h) +
                geom_hline(yintercept = 0, color = "gray") +
                geom_boxplot(data = sp, aes(sequence_type, value, fill = sequence_type, color = sequence_type), alpha = .15) +
                geom_point(data = sp, aes(sequence_type, value, fill = sequence_type, color = sequence_type), alpha = .4) +
                
                # means
                geom_point(data = means, aes(sequence_type, value), size = 3, color = "black", shape = 17) +
                
                geom_vline(xintercept = 1.25, size = 1, color = "gray", linetype = "dotted") +
                
                scale_fill_manual(values = c("forestgreen", "orange")) +
                scale_color_manual(values = c("forestgreen", "orange")) +
                scale_shape_manual(values = c(21, 19)) +
                scale_y_continuous(breaks = c(-.2, 0, .2), position = "left") +
                theme_bw(base_size = 18) +
                theme(#legend.position = "none",
                        panel.grid = element_blank(),
                        #strip.placement = "outside",
                        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                        strip.background = element_rect(fill = 'black', color = "black"),
                        strip.text = element_text(color = "white")) +
                coord_cartesian(ylim = c(-.3, .3)) +
                labs(x = "wind dispersal                       ",
                     y = "partial correlation between\nwind and genetic pattern",
                     fill = "sequence\ntype", color = "sequence\ntype")
        
        # plot for p-values
        sp <- pd %>% filter(fsyndrome != "0") %>% filter(stat == "p") %>% mutate(stat = "significance")
        means <- sp %>% group_by(h, sequence_type) %>% summarize(value = mean(value, na.rm = T))
        pp <- ggplot() +
                facet_grid(stat ~ h) +
                geom_hline(yintercept = 0.5, color = "gray") +
                geom_boxplot(data = sp, aes(sequence_type, 1-value, fill = sequence_type, color = sequence_type), alpha = .15) +
                geom_point(data = sp, aes(sequence_type, 1-value, fill = sequence_type, color = sequence_type), alpha = .4) +
                
                geom_point(data = means, aes(sequence_type, 1-value), size = 3, color = "black", shape = 17) +
                
                geom_vline(xintercept = 1.25, size = 1, color = "gray", linetype = "dotted") +
                
                scale_fill_manual(values = c("forestgreen", "orange")) +
                scale_color_manual(values = c("forestgreen", "orange")) +
                scale_shape_manual(values = c(21, 19)) +
                scale_y_continuous(breaks = c(0, .5, 1), position = "left") +
                #scale_x_continuous(breaks = c(0, .5, 1, 1.5), labels = c("none", "partial", "full", "partial\n~ none")) +
                theme_bw(base_size = 18) +
                theme(#legend.position = "none",
                        panel.grid = element_blank(),
                        strip.background.x = element_blank(),
                        strip.background = element_rect(fill = 'black', color = "black"),
                        strip.text.x = element_blank(),
                        strip.text = element_text(color = "white"),
                        axis.text.x = element_text(angle = 45, hjust = 1, lineheight = .7)) +
                coord_cartesian(ylim = c(0, 1)) +
                labs(x = NULL,
                     y = "partial Mantel\nsignificance",
                     fill = "sequence\ntype", color = "sequence\ntype")
        
        p <- pr + pp + plot_layout(nrow = 2, guides = "collect")
        
        ggs("figures/manuscript/figure_s5.pdf",
            p, width = 12, height = 10, units = "in",
            add = grid.text(paste0("(", letters[1:2], ")"), 
                            x=c(.03), 
                            y=c(.96, .48),
                            gp=gpar(fontsize=20, fontface="bold", col="black")))
        
        
        # significance tests: global and specific null approaches
        
        sig_diff <- s %>% ungroup() %>%
                left_join(st) %>%
                filter(syndrome != 0) %>%
                group_by(h, sequence_type, i) %>%
                summarize(spc = fun(spc)) %>% ungroup() %>%
                spread(sequence_type, spc) %>%
                mutate(spc = SNP - SSR) %>% 
                group_by(h) %>%
                summarize(test = 2,
                          p = 1 - ecdf(spc[i != 0])(spc[i == 0]),
                          spc = spc[i == 0])
        
        pcomp <- s %>%
                left_join(st) %>%
                filter(syndrome != 0) %>% 
                group_by(h, component, sequence_type) %>%
                summarize(p = 1 - ecdf(spc[i != 0])(spc[i == 0]),
                          spc = spc[i == 0]) %>%
                ungroup()
        # difference assessed via fisher exact test
        psig_diff <- pcomp %>% 
                group_by(h) %>%
                summarize(test = 2,
                          p = fisher.test(table(sequence_type, p < .5), alternative="greater")$p.value)
        
        ############################# end SNP-SSR comparison ##############################
        
        
        
        #######
        
        pdr <- spr %>% filter(stat == "spc")  %>%
                filter(syndrome != "1~0") %>%
                filter(syndrome != "0") %>%
                select(-x) %>%
                spread(h, value)
        pdp <- spr %>% filter(stat == "p")  %>%
                filter(syndrome != "1~0") %>%
                filter(syndrome != "0") %>%
                select(-x) %>%
                spread(h, value)
        pd <- bind_cols(pdr, pdp)
        # as.data.frame() %>%
        # ecoclim::pairsData(c("flow", "isolation", "asymmetry", "diversity"), "stat")
        
        # almost none are significant at p<.1
        pvals <- c(cor.test(pd$isolation, pd$asymmetry, method = "spearman")$p.value,
                   cor.test(pd$isolation, pd$diversity, method = "spearman")$p.value,
                   cor.test(pd$diversity, pd$asymmetry, method = "spearman")$p.value,
                   cor.test(pd$isolation1, pd$asymmetry1, method = "spearman")$p.value,
                   cor.test(pd$isolation1, pd$diversity1, method = "spearman")$p.value,
                   cor.test(pd$diversity1, pd$asymmetry1, method = "spearman")$p.value)
        
        
        p1 <- pd %>% ggplot(aes(isolation, diversity)) + geom_point(color = "orange") +
                labs(x = "r: isolation", y = "r: diversity")
        p2 <- pd %>% ggplot(aes(isolation, asymmetry)) + geom_point(color = "orange") +
                labs(x = "r: isolation", y = "r: asymmetry")
        p3 <- pd %>% ggplot(aes(asymmetry, diversity)) + geom_point(color = "orange") +
                labs(x = "r: asymmetry", y = "r: diversity")
        
        p4 <- pd %>% ggplot(aes(asymmetry1, isolation1)) + geom_point(color = "darkred") +
                labs(x = "p: asymmetry", y = "p: isolation") + facet_grid(. ~ "asymmetry")
        p5 <- pd %>% ggplot(aes(diversity1, isolation1)) + geom_point(color = "darkred") +
                labs(x = "p: diversity", y = "p: isolation") + facet_grid("isolation" ~ "diversity")
        p6 <- pd %>% ggplot(aes(diversity1, asymmetry1)) + geom_point(color = "darkred") +
                labs(x = "p: diversity", y = "p: asymmetry") + facet_grid("asymmetry" ~ .)
        
        p7 <- pd %>% ggplot(aes(isolation, isolation1)) + geom_point() +
                labs(x = "r: isolation", y = "p: isolation") + facet_grid(. ~ "isolation")
        p8 <- pd %>% ggplot(aes(asymmetry, asymmetry1)) + geom_point() +
                labs(x = "r: asymmetry", y = "p: asymmetry")
        p9 <- pd %>% ggplot(aes(diversity, diversity1)) + geom_point() +
                labs(x = "r: diversity", y = "p: diversity") + facet_grid("diversity" ~ .)
        
        
        p <- p7 + p4 + p5 +
                p2 + p8 + p6 +
                p1 + p3 + p9 +
                plot_layout(nrow = 3) &
                theme_bw() &
                theme(strip.background = element_rect(fill = "black"),
                      strip.text = element_text(color = "white"))
        #ggsave("figures/manuscript/figure_s3.pdf", p, width = 8.5, height = 8, units = "in")
        
        
        
        hyps <- c("flow", "isolation", "asymmetry", "diversity")
        pd <- spr %>%
                filter(syndrome != "1~0") %>%
                filter(syndrome != "0") %>%
                select(-x) %>%
                spread(h, value) %>%
                as.data.frame() %>%
                ecoclim::pairsData(hyps, "stat", mirror = F) %>%
                mutate(x_var = factor(x_var, levels = hyps),
                       y_var = factor(y_var, levels = hyps))
        rd <- pd %>%
                group_by(stat, x_var, y_var) %>%
                summarize(r2 = cor(x_value, y_value, use = "pairwise.complete.obs")^2,
                          r2 = round(r2, 3))
        
        pr <- ggplot(data = pd %>% filter(stat == "spc"), aes(x_value, y_value)) +
                facet_grid(y_var ~ x_var) +
                geom_point() +
                shadowtext::geom_shadowtext(data = rd %>% filter(stat == "spc"), 
                                            aes(x = -1, y = 1, label = r2), 
                                            hjust = 0, vjust = 1, color = "red", bg.color = "white", 
                                            fontface = "bold", size = 5) +
                scale_x_continuous(breaks = c(-1, 0, 1)) +
                scale_y_continuous(breaks = c(-1, 0, 1)) +
                labs(x = "partial wind-genetic correltion coefficient",
                     y = "partial wind-genetic correltion coefficient") +
                theme_bw() +
                theme(strip.background = element_rect(fill = "black"),
                      strip.text = element_text(color = "white"))
        
        pp <- ggplot(data = pd %>% filter(stat == "p"), aes(x_value, y_value)) +
                facet_grid(y_var ~ x_var) +
                geom_point() +
                shadowtext::geom_shadowtext(data = rd %>% filter(stat == "p"), 
                                            aes(x = 0, y = 1, label = r2), 
                                            hjust = 0, vjust = 1, color = "red", bg.color = "white", 
                                            fontface = "bold", size = 5) +
                scale_x_continuous(breaks = c(0, .5, 1)) +
                scale_y_continuous(breaks = c(0, .5, 1)) +
                labs(x = "partial wind-genetic correltion coefficient",
                     y = "partial wind-genetic correltion coefficient") +
                theme_bw() +
                theme(strip.background = element_rect(fill = "black"),
                      strip.text = element_text(color = "white"))
        
        
        ratios <- function(x) c(log(x[1] / x[2]), log(x[2] / x[1]))
        f <- annotated %>%
                filter(wind_p == "wind_p1") %>%
                filter(to != from,
                       is.finite(clim),
                       is.finite(distance),
                       is.finite(wind)) %>%
                group_by(component, edge) %>%
                mutate(mig = log(dmig),
                       mig_ratio = ratios(mig),
                       div_ratio = log(div_ar_to / div_ar_from),
                       fst = ifelse(fst < 0, 0, fst)) %>%
                select(flow = mig, isolation = fst, asymmetry = mig_ratio, diversity = div_ratio) %>%
                filter(is.finite(flow), is.finite(isolation), is.finite(asymmetry), is.finite(diversity)) %>%
                mutate(group = paste(component, edge)) %>%
                as.data.frame() %>%
                ecoclim::pairsData(hyps, c("component", "group"), mirror = F) %>%
                mutate(x_var = factor(x_var, levels = hyps),
                       y_var = factor(y_var, levels = hyps))
        
        ff <- f %>% filter(group %in% sample(group, 1000)) %>%
                mutate(x_value = ifelse(x_var == "asymmetry" & abs(x_value) > 2, NA, x_value),
                       y_value = ifelse(y_var == "asymmetry" & abs(y_value) > 2, NA, y_value))
        
        rdf <- f %>%
                group_by(component, x_var, y_var) %>%
                summarize(r2 = cor(x_value, y_value, use = "pairwise.complete.obs")^2,
                          r2 = round(r2, 3)) %>%
                group_by(x_var, y_var) %>%
                summarize(r2 = median(r2, na.rm = T)) %>%
                right_join(ff) %>%
                group_by(x_var, y_var) %>%
                summarize(x_value = min(x_value, na.rm = T), 
                          y_value = max(y_value, na.rm = T),
                          r2 = median(r2))
        
        
        
        pf <- ggplot(data = ff, aes(x_value, y_value)) +
                facet_grid(y_var ~ x_var, scales = "free") +
                geom_point(size = .5) +
                shadowtext::geom_shadowtext(data = rdf, 
                                            aes(x = x_value, y = y_value, label = r2), 
                                            hjust = 0, vjust = 1, color = "red", bg.color = "white", 
                                            fontface = "bold", size = 5) +
                scale_x_continuous() +
                scale_y_continuous() +
                labs(x = "            Fst                log gene flow ratio     log diversity ratio",
                     y = "log gene flow ratio                   Fst                     log gene flow     ") +
                theme_bw() +
                theme(strip.background = element_rect(fill = "black"),
                      strip.text = element_text(color = "white"))
        
        p <- pf + pr + pp & theme(strip.text = element_text(size = 12))
        ggs("figures/manuscript/figure_s6.pdf", p, 
            width = 15, height = 5, units = "in",
            add = grid.text(paste0("(", letters[1:3], ")"), 
                            x=c(.02, .35, .68), 
                            y=c(.96),
                            gp=gpar(fontsize=14, fontface="bold", col="black")))
        
        
        
        
        #### compare statistical tests ####
        
        pd <- bind_rows(rsig, psig) %>%
                select(h, test, stat, p) %>%
                spread(stat, p)
        p <- pd %>%
                ggplot(aes(p, spc)) +
                geom_vline(xintercept = .5) +
                geom_hline(yintercept = .5) +
                ggplot2::annotate(geom = "line", x = c(.001, 1), y = c(.001, 1), linetype = "dotted") +
                geom_point() +
                ggplot2::annotate(geom = "text", x = .001, y = .23, label = paste("r =", round(cor(pd$p, pd$spc), 2)),
                                  color = "red", hjust = 0) +
                scale_x_log10(limits = c(.001, 1), breaks = c(.01, .1, .5, 1)) +
                scale_y_log10(limits = c(.001, 1), breaks = c(.01, .1, .5, 1)) +
                coord_fixed() +
                theme_bw() + 
                theme(panel.grid.minor = element_blank()) +
                labs(x = "p-value, dataset null method",
                     y = "p-value, global null method")
        ggs("figures/manuscript/figure_s4.pdf", p, width = 4, height = 4, units = "in",
            add = grid.text(" ", x = .5, y = .5))
}

