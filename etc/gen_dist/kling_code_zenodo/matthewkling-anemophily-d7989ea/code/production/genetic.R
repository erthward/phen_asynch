
fst <- function(infile, metadata){
        
        #infile <- readd(datasets)$gen[1]
        
        # import
        d <- haploRGP(infile)
        d <- d$genos
        genome <- metadata %>%
                bind_rows() %>%
                filter(component == basename(dirname(infile))) %>%
                pull(genome)
        
        # convert to hierfstat format
        d <- lapply(d, function(x) apply(x, 1:2, paste, collapse=""))
        for(i in 1:length(d)){
                d[[i]] <- cbind(rep(i, nrow(d[[i]])), d[[i]])
                d[[i]][d[[i]]=="NANA"] <- NA
        }
        d <- do.call("rbind", d) %>% 
                as.data.frame() %>%
                mutate_each(as.character) %>%
                mutate_each(as.integer)
        
        # compute fst
        npop <- length(unique(d[,1]))
        p <- t(combn(1:npop, 2))
        p <- cbind(p,
                   apply(p, 1, function(x){
                           y <- try(wc(filter(d, V1 %in% x), diploid = genome=="nu")$FST)
                           if(class(y)=="try-error") y <- NA # happens when all genotypes are idendical
                           y
                   }))
        
        fst <- matrix(NA, npop, npop)
        for(i in 1:nrow(p)) fst[p[i,1],p[i,2]] <- fst[p[i,2],p[i,1]] <- p[i,3]
        return(fst)
}



# modifying diveRsity::divMigrate (and diveRsity::RGP) to support haplotype data

haploRGP <- function (infile) {
        
        fastScan <- function(fname) {
                s <- file.info(fname)$size
                buf <- readChar(fname, s, useBytes = TRUE)
                if (length(grep("\r", buf)) != 0L) {
                        buf <- gsub("\r", "\n", buf)
                        buf <- gsub("\n\n", "\n", buf)
                }
                return(strsplit(buf, "\n", fixed = TRUE, useBytes = TRUE)[[1]])
        }
        if(is.list(infile)) {
                infile <- as.matrix(infile)
                dat <- apply(infile, 1, function(x) {
                        x <- x[!is.na(x)]
                        return(paste(x, collapse = "\t"))
                })
        }else{
                dat <- fastScan(infile)
                dat <- sapply(dat, function(x) {
                        sub("^\\s+", "", x)
                })
                dat <- sapply(dat, function(x) {
                        return(sub("\\s+$", "", x))
                })
                names(dat) <- NULL
        }
        popLoc <- grep("^([[:space:]]*)pop([[:space:]]*)$", tolower(dat))
        if(popLoc[1] == 3) {
                if (length(strsplit(dat[4], split = "\\s+")[[1]][-1]) > 
                    1) {
                        locs <- strsplit(dat[2], split = "\\s+")[[1]]
                        if (length(locs) == 1) {
                                locs <- strsplit(dat[2], split = ",")[[1]]
                        }
                        locs <- as.character(sapply(locs, function(x) {
                                x <- strsplit(x, split = "")[[1]]
                                if (is.element(",", x)) {
                                        x <- x[-(which(x == ","))]
                                }
                                return(paste(x, collapse = ""))
                        }))
                        dat <- c(dat[1], locs, dat[-(1:2)])
                }
        }else{
                locs <- as.character(dat[2:(popLoc[1] - 1)])
        }
        locs <- as.character(sapply(locs, function(x) {
                return(strsplit(x, split = "\\s+")[[1]][1])
        }))
        popLoc <- grep("^([[:space:]]*)pop([[:space:]]*)$", tolower(dat))
        npops <- length(popLoc)
        no_col <- length(locs) + 1
        nloci <- length(locs)
        strt <- popLoc + 1
        ends <- c(popLoc[-1] - 1, length(dat))
        genoRet <- function(strt, ends, x) {
                out <- strsplit(x[strt:ends], split = "\\s+")
                x <- do.call("rbind", c(out, deparse.level = 0))
                if (round(mean(nchar(x[, 2]))) == 1L) {
                        x[, 1] <- paste(x[, 1], x[, 2], sep = "")
                        x <- x[, (-2)]
                }
                x[x == "-9"] <- NA
                x[x == "0000"] <- NA
                x[x == "000000"] <- NA
                list(ls = x[, (-1)], nms = as.vector(x[, 1]))
        }
        genos <- mapply(genoRet, strt = strt, ends = ends, MoreArgs = list(x = dat), 
                        SIMPLIFY = FALSE)
        indNames <- lapply(genos, "[[", 2)
        genos <- lapply(genos, "[[", 1)
        badLoc <- apply(genos[[1]], 2, function(x) {
                sum(is.na(x)) == length(x)
        })
        badLoc <- which(!badLoc)
        
        dig <- mean(nchar(na.omit(genos[[badLoc[1]]][, 1])))
        if(dig==3){
                types <- 1
                gp <- dig
        }else{
                types <- 2
                gp <- round(mean(dig/2))
        }
        
        # convert genotypes to arrays
        genos <- lapply(genos, function(x) {
                #x <- genos[[1]]
                if(types==2){
                        al1 <- substr(x, 1, gp)
                        al2 <- substr(x, (gp + 1), (gp * 2))
                        out <- array(NA, dim = c(nrow(x), ncol(x), 2))
                        out[, , 1] <- al1
                        out[, , 2] <- al2
                }
                if(types==1){
                        al1 <- substr(x, 1, gp)
                        out <- array(NA, dim = c(nrow(x), ncol(x), 1))
                        out[, , 1] <- al1
                }
                return(out)
        })
        
        statFun <- function(x, cl = NULL) {
                popSizes <- apply(x, 2, function(y) {
                        length(na.omit(y[, 1])) * 2
                })
                af <- lapply(1:dim(x)[2], function(i) {
                        y <- as.vector(na.omit(x[, i, ])) # this summarizes over the diplotypes, so haplo should work
                        nms <- unique(y)[order(unique(y))]
                        ot <- myTab(y)
                        names(ot) <- nms
                        return(ot)
                })
                popSizes <- popSizes/types
                list(af = af, ps = popSizes)
        }
        
        
        obsAllSize <- lapply(genos, statFun)
        
        af <- lapply(obsAllSize, function(x) {
                out <- x$af
                x$af <- NULL
                return(out)
        })
        ps <- lapply(obsAllSize, function(x) {
                out <- x$ps
                x$ps <- NULL
                return(out)
        })
        
        af <- lapply(1:(nloci), function(i) {
                return(lapply(af, "[[", i))
        })
        ps <- lapply(1:(nloci), function(i) {
                return(sapply(ps, "[", i))
        })
        
        
        # rearrange data by loci
        check <- function(args, gp) {
                #args <- af[[1]]
                npops <- length(args)
                pad <- paste("%0", gp, "g", sep = "")
                rnames <- sprintf(pad, unique(sort(as.numeric(unlist(lapply(args, names))))))
                out <- matrix(0, nrow = length(rnames), ncol = npops)
                rownames(out) <- as.character(rnames)
                for (i in 1:npops) {
                        out[match(names(args[[i]]), rownames(out)), i] <- as.numeric(args[[i]])
                }
                return(out)
        }
        
        af <- lapply(af, check, gp = gp)
        
        gc()
        list(af = af, genos = genos, ps = ps, gp = gp, indnms = indNames, 
             locs = locs)
}



haploMigrate <- function(infile = NULL, outfile = NULL, boots = 0, stat = "all", 
                         filter_threshold = 0, plot_network = FALSE, plot_col = "darkblue", 
                         para = FALSE){
        
        dat <- haploRGP(infile)
        npops <- length(dat$genos)
        nloci <- length(dat$af)
        dat$af <- lapply(dat$af, function(x) {
                cs <- colSums(x)
                x[, cs == 0] <- NA
                return(x)
        })
        if(!is.null(outfile)) {
                dir.create(path = paste(getwd(), "/", outfile, "-[divMigrate]", 
                                        "/", sep = ""))
                of <- paste(getwd(), "/", outfile, "-[divMigrate]", "/", 
                            sep = "")
        }
        pw <- combn(npops, 2)
        hths <- lapply(dat$af, pwHt, pw = pw - 1)
        ht <- lapply(hths, "[[", "ht")
        hs <- lapply(hths, "[[", "hs")
        if (stat == "d" || stat == "all" || stat == "Nm") {
                d <- function(ht, hs) {
                        return(((ht - hs)/(1 - hs)) * 2)
                }
        }
        if (stat == "gst" || stat == "all" || stat == "Nm") {
                g <- function(ht, hs) {
                        ot <- (ht - hs)/ht
                        diag(ot) <- 0
                        return(ot)
                }
        }
        if (stat == "Nm" || stat == "all") {
                Nm <- function(g, d, n) {
                        t1 <- (1 - g)/g
                        t2 <- ((n - 1)/n)^2
                        t3 <- ((1 - d)/(1 - ((n - 1)/n) * d))
                        return(0.25 * t1 * t2 * t3)
                }
        }
        if (stat == "d" || stat == "all" || stat == "Nm") {
                dloc <- mapply(d, ht = ht, hs = hs, SIMPLIFY = "array")
                dloc[is.nan(dloc)] <- 1
                hrmD <- apply(dloc, c(1, 2), function(x) {
                        mn <- mean(x, na.rm = TRUE)
                        vr <- var(x, na.rm = TRUE)
                        return(1/((1/mn) + vr * (1/mn)^3))
                })
                dMig <- (1 - hrmD)/hrmD
                dMig[is.infinite(dMig)] <- NA
                dRel <- dMig/max(dMig, na.rm = TRUE)
                dRel[is.nan(dRel)] <- NA
        }
        if (stat == "gst" || stat == "all" || stat == "Nm") {
                g <- function(ht, hs) {
                        ot <- (ht - hs)/ht
                        diag(ot) <- 0
                        return(ot)
                }
                hsAr <- array(unlist(hs), dim = c(npops, npops, nloci))
                mnHs <- apply(hsAr, c(1, 2), mean, na.rm = TRUE)
                htAr <- array(unlist(ht), dim = c(npops, npops, nloci))
                mnHt <- apply(htAr, c(1, 2), mean, na.rm = TRUE)
                hrmGst <- g(mnHt, mnHs)
                gMig <- ((1/hrmGst) - 1)/4
                gMig[is.infinite(gMig)] <- NA
                gRel <- gMig/max(gMig, na.rm = TRUE)
        }
        
        if (stat == "all" || stat == "Nm") {
                nm <- Nm(hrmGst, hrmD, 2)
                diag(nm) <- NA
                nmRel <- nm/max(nm, na.rm = TRUE)
        }
        
        if (boots != 0L) {
                if (stat == "d" || stat == "all" || stat == "Nm") {
                        dRelPlt <- dRel
                        dRelPlt[dRelPlt < filter_threshold] <- 0
                }
                if (stat == "gst" || stat == "all" || stat == "Nm") {
                        gRelPlt <- gRel
                        gRelPlt[gRelPlt < filter_threshold] <- 0
                }
                if (stat == "all" || stat == "Nm") {
                        nmRelPlt <- nmRel
                        nmRelPlt[nmRelPlt < filter_threshold] <- 0
                }
                
                ps <- sapply(dat$indnms, length)
                idx <- lapply(1:boots, function(i) {
                        lapply(ps, function(x) {
                                return(sample(x, size = x, replace = TRUE))
                        })
                })
                if (para) {
                        cl <- parallel::makeCluster(detectCores())
                        parallel::clusterExport(cl, c("bsFun", "dat", "pw", 
                                                      "stat"), envir = environment())
                        bsStat <- parallel::parLapply(cl, idx, function(x) {
                                return(bsFun(genos = dat$genos, idx = x, af = dat$af, 
                                             pw = pw, stat = stat))
                        })
                        parallel::stopCluster(cl)
                }
                else {
                        bsStat <- lapply(idx, function(x) {
                                return(bsFun(genos = dat$genos, idx = x, af = dat$af, 
                                             pw = pw, stat = stat))
                        })
                }
                if (stat == "d" || stat == "all") {
                        bsD <- sapply(bsStat, "[[", "dRel", simplify = "array")
                }
                if (stat == "gst" || stat == "all") {
                        bsG <- sapply(bsStat, "[[", "gRel", simplify = "array")
                }
                if (stat == "Nm" || stat == "all") {
                        bsNm <- sapply(bsStat, "[[", "nmRel", simplify = "array")
                }
                sigDiff <- function(x, y) {
                        if(any(is.na(c(x, y)))) return(NA)
                        if (x[1] < y[1] && x[2] < y[1]) {
                                return(TRUE)
                        }
                        else {
                                return(FALSE)
                        }
                }
                if (stat == "d" || stat == "all") {
                        sigMatD <- matrix(NA, nrow = ncol(dRel), ncol(dRel))
                        for (i in 1:ncol(pw)) {
                                p1 <- quantile(bsD[pw[1, i], pw[2, i], ], prob = c(0.025, 0.975), na.rm=T)
                                p2 <- quantile(bsD[pw[2, i], pw[1, i], ], prob = c(0.025, 0.975), na.rm=T)
                                sigMatD[pw[2, i], pw[1, i]] <- sigDiff(p1, p2)
                                sigMatD[pw[1, i], pw[2, i]] <- sigDiff(p2, p1)
                        }
                        dRelPlt[!sigMatD] <- 0
                }
                if (stat == "gst" || stat == "all") {
                        sigMatG <- matrix(NA, nrow = ncol(gRel), ncol(gRel))
                        for (i in 1:ncol(pw)) {
                                p1 <- quantile(bsG[pw[1, i], pw[2, i], ], prob = c(0.025, 0.975), na.rm=T)
                                p2 <- quantile(bsG[pw[2, i], pw[1, i], ], prob = c(0.025, 0.975), na.rm=T)
                                sigMatG[pw[2, i], pw[1, i]] <- sigDiff(p1, p2)
                                sigMatG[pw[1, i], pw[2, i]] <- sigDiff(p2, p1)
                        }
                        gRelPlt[!sigMatG] <- 0
                }
                if (stat == "Nm" || stat == "all") {
                        sigMatNm <- matrix(NA, nrow = ncol(nmRel), ncol(nmRel))
                        for (i in 1:ncol(pw)) {
                                p1 <- quantile(bsNm[pw[1, i], pw[2, i], ], prob = c(0.025, 0.975), na.rm=T)
                                p2 <- quantile(bsNm[pw[2, i], pw[1, i], ], prob = c(0.025, 0.975), na.rm=T)
                                sigMatNm[pw[2, i], pw[1, i]] <- sigDiff(p1, p2)
                                sigMatNm[pw[1, i], pw[2, i]] <- sigDiff(p2, p1)
                        }
                        nmRelPlt[!sigMatNm] <- 0
                }
                
                if (stat == "d") {
                        list(dRelMig = dRel, dRelMigSig = dRelPlt)
                }
                else if (stat == "gst") {
                        list(gRelMig = gRel, gRelMigSig = gRelPlt)
                }
                else if (stat == "Nm") {
                        list(nmRelMig = nmRel, nmRelMigSig = nmRelPlt)
                }
                else if (stat == "all") {
                        list(dRelMig = dRel, dRelMigSig = dRelPlt, gRelMig = gRel, 
                             gRelMigSig = gRelPlt, nmRelMig = nmRel, nmRelMigSig = nmRelPlt)
                }
        }else {
                if (stat == "d") {
                        list(dRelMig = dRel)
                }else if (stat == "gst") {
                        list(gRelMig = gRel)
                }else if (stat == "Nm") {
                        list(nmRelMig = nmRel)
                }else if (stat == "all") {
                        list(dRelMig = dRel, gRelMig = gRel, nmRelMig = nmRel)
                }
        }
}


nullMigrate <- function(infile, boots = 1000){
        # infile = s$filename
        # boots = 100
        message(infile)
        
        # function to convert genos to af
        geno2af <- function(genos, nloci){
                lapply(1:nloci, function(l){
                        gl <- lapply(genos, function(x) table(as.vector(x[,l,])))
                        gl <- lapply(gl, function(x) x / sum(x))
                        alleles <- unique(as.vector(unlist(sapply(gl, names))))
                        f <- t(sapply(alleles, function(x) sapply(gl, function(y) y[x])))
                        f[is.na(f)] <- 0
                        return(f)
                })
        }
        
        # function to calculate dRelMig from af, adapted from divMigrate
        af2mig <- function(af, npops){
                pw <- combn(npops, 2)
                hths <- lapply(af, pwHt, pw = pw - 1)
                ht <- lapply(hths, "[[", "ht")
                hs <- lapply(hths, "[[", "hs")
                d <- function(ht, hs) return(((ht - hs)/(1 - hs))    * 2) 
                
                dloc <- mapply(d, ht = ht, hs = hs, SIMPLIFY = "array")
                dloc[is.nan(dloc)] <- 1
                hrmD <- apply(dloc, c(1, 2), function(x) {
                        mn <- mean(x, na.rm = TRUE)
                        vr <- var(x, na.rm = TRUE)
                        return(1/((1/mn) + vr * (1/mn)^3))
                })
                dMig <- (1 - hrmD)/hrmD
                dMig[is.infinite(dMig)] <- NA
                dRel <- dMig/max(dMig, na.rm = TRUE)
                dRel[is.nan(dRel)] <- NA
                #return(dMig)
                return(dRel)
        }
        
        #infile <- filename
        dat <- haploRGP(infile)
        npops <- length(dat$genos)
        genos <- dat$genos 
        
        # remove loci with no data
        allgenos <- abind:::abind(genos, along = 1)
        empty <- apply(allgenos, 2, function(x) all(is.na(x)))
        genos <- lapply(genos, function(x) x[,!empty,])
        genos <- lapply(genos, function(x) array(x, dim = c(dim(x)[1:2], dim(allgenos)[3]))) # necessary for haplo data
        nloci <- sum(!empty)
        
        # calculate allele freqs and migration
        af <- geno2af(genos, nloci)
        mig <- af2mig(af, npops)
        
        # some vars needed for resampling
        allgenos <- abind:::abind(genos, along = 1)
        popind <- sapply(genos, function(x) dim(x)[1])
        totind <- sum(popind)
        cumind <- c(0, cumsum(popind))
        bmig <- array(dim = c(npops, npops, boots))
        
        # resampling
        for(i in 1:boots){
                #setTxtProgressBar(pb, i/boots)
                
                # shuffle samples among popns
                agi <- allgenos[sample(totind, totind, replace = F),,]
                agi <- array(agi, dim(allgenos)) # necessary for haplo data
                gi <- list()
                for(p in 1:length(popind)){
                        agip <- agi[(cumind[p]+1):cumind[p+1],,]
                        
                        gi[[p]] <- array(agip, c(dim(agip)[1:2], dim(allgenos)[3]))
                        
                } 
                
                # calculate allele freqs and migration
                afi <- geno2af(gi, nloci)
                migi <- af2mig(afi, npops)
                bmig[,,i] <- migi
        }
        
        # quantile of true mig in nulls (conservative approach to dealing with equalities)
        # if less than 50% of the nulls are less than or equal to x, then p
        # if more than 50% are greater or equal to x, then 1 - p 
        # otherwise, .5
        migs <- abind::abind(mig, bmig, along = 3)
        ple <- apply(migs, c(1, 2), function(x) mean(x[1] >= x[2:(boots+1)]))
        pge <- apply(migs, c(1, 2), function(x) mean(x[1] <= x[2:(boots+1)]))
        p <- matrix(.5, nrow(pge), ncol(pge))
        p[ple < .5 & is.finite(ple)] <- ple[ple < .5 & is.finite(ple)]
        p[pge < .5 & is.finite(pge)] <- 1 - pge[pge < .5 & is.finite(pge)]
        p[diag(p)] <- NA
        
        return(list(dRelMig = mig,
                    dRelMigP = p))
}


# modified diveRsity::diffCalc to support haplotype data


haploDiff <- function (infile = NULL, outfile = NULL, fst = FALSE, pairwise = FALSE, 
                       bs_locus = FALSE, bs_pairwise = FALSE, boots = NULL, ci_type = "individuals", 
                       alpha = 0.05, para = FALSE) 
{ message(infile)
        dCalc <- function(ht, hs, n = NULL) {
                if (!is.null(n)) {
                        return(((ht - hs)/(1 - hs)) * (n/(n - 1)))
                }
                else {
                        return(2 * ((ht - hs)/(1 - hs)))
                }
        }
        gstCalc <- function(ht, hs) {
                return((ht - hs)/ht)
        }
        GstCalc <- function(ht, hs, n = NULL) {
                if (!is.null(n)) {
                        htmax <- ((n - 1) + hs)/n
                }
                else {
                        htmax <- (1 + hs)/2
                }
                return(((ht - hs)/ht)/((htmax - hs)/htmax))
        }
        GgstCalc <- function(ht, hs, n = NULL) {
                if (!is.null(n)) {
                        return((n * (ht - hs))/(((n * ht) - hs) * (1 - hs)))
                }
                else {
                        return((2 * (ht - hs))/(((2 * ht) - hs) * (1 - hs)))
                }
        }
        thetaCalc <- function(a, b, cdat) {
                return(a/(a + b + cdat))
        }
        diffCalcbCor <- function(act, bs) {
                if (is.matrix(bs)) {
                        mn <- rowMeans(bs, na.rm = TRUE) - act
                        mn[is.nan(mn)] <- NA
                        bc <- mapply(`-`, split(bs, row(bs)), mn)
                        return(bc)
                }
                else {
                        mn <- mean(bs, na.rm = TRUE) - act
                        mn[is.nan(mn)] <- NA
                        return(bs - mn)
                }
        }
        bs <- boots
        ip <- haploRGP(infile) ################ changed from rgp() ##########
        ps = sapply(ip$genos, function(x) {
                dim(x)[1]
        })
        np = ncol(ip$af[[1]])
        # if (ci_type == "loci") {
        #         nl <- dim(ip$genos[[1]])[2]
        # }
        pw <- combn(np, 2)
        # if (para) {
        #         ncor <- parallel::detectCores()
        # }
        # if (ci_type == "individuals") {
        #         if (!is.null(bs)) {
        #                 idx <- lapply(1:bs, function(i) {
        #                         lapply(1:np, function(j) {
        #                                 sample(ps[j], size = ps[j], replace = TRUE)
        #                         })
        #                 })
        #         }
        # }
        # else {
        #         idx <- lapply(1:bs, function(i) {
        #                 sample(nl, nl, replace = TRUE)
        #         })
        # }
        preStats <- statCalc(rsDat = ip$genos, al = ip$af, fst = fst, 
                             bs = FALSE)
        if (fst) {
                hsum <- lapply(preStats$hsum, tabMerge)
                hsum <- lapply(hsum, function(x) {
                        x[!(names(x) == "NA")]
                })
                varC <- mapply("glbWCcpp", hsum = hsum, af = preStats$alOut, 
                               indtyp = preStats$indtyp, SIMPLIFY = FALSE)
                rm(hsum)
                locFst <- sapply(varC, function(x) {
                        return(x$a/(x$a + x$b + x$c))
                })
                locFit <- sapply(varC, function(x) {
                        return(1 - (x$c/(x$a + x$b + x$c)))
                })
                locFis <- sapply(varC, function(x) {
                        return(1 - (x$c/(x$b + x$c)))
                })
                glba <- mean(sapply(varC, "[[", "a"), na.rm = TRUE)
                glbb <- mean(sapply(varC, "[[", "b"), na.rm = TRUE)
                glbc <- mean(sapply(varC, "[[", "c"), na.rm = TRUE)
                glbTheta <- glba/(glba + glbb + glbc)
                glbFit <- 1 - (glbc/(glba + glbb + glbc))
                glbFis <- 1 - (glbc/(glbb + glbc))
        }
        hrm <- function(x) {
                return(length(x)/(sum(1/x)))
        }
        harmLoc <- sapply(preStats$indtyp, hrm)
        hthsLoc <- mapply(varFunc, af = preStats$alOut, sHarm = harmLoc, 
                          SIMPLIFY = FALSE)
        dLoc <- dCalc(ht = sapply(hthsLoc, "[[", "htEst"), 
                      hs = sapply(hthsLoc, "[[", "hsEst"), n = np)
        mnD <- mean(dLoc, na.rm = TRUE)
        varD <- var(dLoc, na.rm = TRUE)
        dGlb <- 1/((1/mnD) + varD * (1/mnD)^3)
        gLoc <- gstCalc(sapply(hthsLoc, "[[", "htEst"), 
                        sapply(hthsLoc, "[[", "hsEst"))
        gGlb <- gstCalc(mean(sapply(hthsLoc, "[[", "htEst"), 
                             na.rm = TRUE), mean(sapply(hthsLoc, "[[", "hsEst"), 
                                                 na.rm = TRUE))
        GLoc <- GstCalc(sapply(hthsLoc, "[[", "htEst"), 
                        sapply(hthsLoc, "[[", "hsEst"), n = np)
        GGlb <- GstCalc(mean(sapply(hthsLoc, "[[", "htEst"), 
                             na.rm = TRUE), mean(sapply(hthsLoc, "[[", "hsEst"), 
                                                 na.rm = TRUE), n = np)
        GGLoc <- GgstCalc(sapply(hthsLoc, "[[", "htEst"), 
                          sapply(hthsLoc, "[[", "hsEst"), n = np)
        GGGlb <- GgstCalc(mean(sapply(hthsLoc, "[[", "htEst"), 
                               na.rm = TRUE), mean(sapply(hthsLoc, "[[", "hsEst"), 
                                                   na.rm = TRUE), n = np)
        if (!fst) {
                est <- data.frame(loci = c(ip$locs, "Global"), 
                                  gst = round(c(gLoc, gGlb), 4), Gst = round(c(GLoc, 
                                                                               GGlb), 4), GGst = round(c(GGLoc, GGGlb), 4), 
                                  D = round(c(dLoc, dGlb), 4))
                rownames(est) <- NULL
        }
        else {
                est <- data.frame(loci = c(ip$locs, "Global"), 
                                  gst = round(c(gLoc, gGlb), 4), Gst = round(c(GLoc, 
                                                                               GGlb), 4), GGst = round(c(GGLoc, GGGlb), 4), 
                                  D = round(c(dLoc, dGlb), 4), Fis = round(c(locFis, 
                                                                             glbFis), 4), Fst = round(c(locFst, glbTheta), 
                                                                                                      4), Fit = round(c(locFit, glbFit), 4))
                rownames(est) <- NULL
        }
        if (bs_locus || bs_pairwise) {
                al <- lapply(ip$af, function(x) {
                        nms <- rownames(x)
                        x[x != 0] <- 0
                        rownames(x) <- nms
                        return(x)
                })
                if (para) {
                        cl <- parallel::makeCluster(ncor)
                        bsDat <- parallel::parLapply(cl, idx, statCalc, rsDat = ip$genos, 
                                                     al = al, fst = fst, ci_type = ci_type)
                        parallel::stopCluster(cl)
                }
                else {
                        system.time({
                                bsDat <- lapply(idx, statCalc, rsDat = ip$genos, 
                                                al = al, fst = fst, ci_type = ci_type)
                        })
                }
                indtyp <- lapply(bsDat, "[[", "indtyp")
                rm(idx)
        }
        if (bs_locus) {
                if (fst) {
                        bsVarC <- lapply(bsDat, function(x) {
                                hsum <- lapply(x$hsum, tabMerge)
                                hsum <- lapply(hsum, function(x) {
                                        x[!(names(x) == "NA")]
                                })
                                stat <- mapply(glbWCcpp, hsum = hsum, af = x$alOut, 
                                               indtyp = x$indtyp, SIMPLIFY = FALSE)
                                bsFst <- sapply(stat, function(y) {
                                        return(y$a/(y$a + y$b + y$c))
                                })
                                bsFit <- sapply(stat, function(y) {
                                        return(1 - (y$c/(y$a + y$b + y$c)))
                                })
                                bsFis <- sapply(stat, function(y) {
                                        return(1 - (y$c/(y$b + y$c)))
                                })
                                glba <- mean(sapply(stat, "[[", "a"), 
                                             na.rm = TRUE)
                                glbb <- mean(sapply(stat, "[[", "b"), 
                                             na.rm = TRUE)
                                glbc <- mean(sapply(stat, "[[", "c"), 
                                             na.rm = TRUE)
                                glbst <- glba/(glba + glbb + glbc)
                                glbit <- 1 - (glbc/(glba + glbb + glbc))
                                glbis <- 1 - (glbc/(glbb + glbc))
                                list(bsFstLoc = bsFst, bsFitLoc = bsFit, bsFisLoc = bsFis, 
                                     bsFstAll = glbst, bsFitAll = glbit, bsFisAll = glbis)
                        })
                        bsFstL <- sapply(bsVarC, "[[", "bsFstLoc")
                        bsFisL <- sapply(bsVarC, "[[", "bsFisLoc")
                        bsFitL <- sapply(bsVarC, "[[", "bsFitLoc")
                        bsFstA <- sapply(bsVarC, "[[", "bsFstAll")
                        bsFisA <- sapply(bsVarC, "[[", "bsFisAll")
                        bsFitA <- sapply(bsVarC, "[[", "bsFitAll")
                        bsFstL <- diffCalcbCor(locFst, bsFstL)
                        bsFisL <- diffCalcbCor(locFis, bsFisL)
                        bsFitL <- diffCalcbCor(locFit, bsFitL)
                        bsFstA <- diffCalcbCor(glbTheta, bsFstA)
                        bsFisA <- diffCalcbCor(glbFis, bsFisA)
                        bsFitA <- diffCalcbCor(glbFit, bsFitA)
                        LCI <- data.frame(Fst = apply(bsFstL, 2, quantile, 
                                                      probs = alpha/2, na.rm = TRUE), Fis = apply(bsFisL, 
                                                                                                  2, quantile, probs = alpha/2, na.rm = TRUE), 
                                          Fit = apply(bsFitL, 2, quantile, probs = alpha/2, 
                                                      na.rm = TRUE))
                        UCI <- data.frame(Fst = apply(bsFstL, 2, quantile, 
                                                      probs = 1 - (alpha/2), na.rm = TRUE), Fis = apply(bsFisL, 
                                                                                                        2, quantile, probs = 1 - (alpha/2), na.rm = TRUE), 
                                          Fit = apply(bsFitL, 2, quantile, probs = 1 - 
                                                              (alpha/2), na.rm = TRUE))
                }
                allHarm <- lapply(bsDat, function(x) {
                        sapply(x$indtyp, function(y) {
                                return(((1/length(y)) * (sum(y^-1)))^-1)
                        })
                })
                listAdd <- function(bsDat, sHarm) {
                        bsDat$nHarm <- sHarm
                        return(bsDat)
                }
                bsDat <- mapply(listAdd, bsDat = bsDat, sHarm = allHarm, 
                                SIMPLIFY = FALSE)
                hetVar <- lapply(bsDat, function(x) {
                        stat <- mapply("varFunc", af = x$alOut, sHarm = x$nHarm, 
                                       SIMPLIFY = FALSE)
                        ht <- sapply(stat, "[[", "htEst")
                        hs <- sapply(stat, "[[", "hsEst")
                        n <- ncol(x$alOut[[1]])
                        bsDLoc <- dCalc(ht = ht, hs = hs, n = n)
                        mnD <- mean(bsDLoc, na.rm = TRUE)
                        varD <- var(bsDLoc, na.rm = TRUE)
                        bsDAll <- 1/((1/mnD) + varD * (1/mnD)^3)
                        bsgLoc <- gstCalc(ht = ht, hs = hs)
                        bsgAll <- gstCalc(ht = mean(ht, na.rm = TRUE), hs = mean(hs, 
                                                                                 na.rm = TRUE))
                        bsGLoc <- GstCalc(ht = ht, hs = hs, n = n)
                        bsGAll <- GstCalc(ht = mean(ht, na.rm = TRUE), hs = mean(hs, 
                                                                                 na.rm = TRUE), n = n)
                        bsGGLoc <- GgstCalc(ht = ht, hs = hs, n = n)
                        bsGGAll <- GgstCalc(ht = mean(ht, na.rm = TRUE), 
                                            hs = mean(hs, na.rm = TRUE), n = n)
                        list(bsDLoc = bsDLoc, bsgLoc = bsgLoc, bsGLoc = bsGLoc, 
                             bsDAll = bsDAll, bsgAll = bsgAll, bsGAll = bsGAll, 
                             bsGGLoc = bsGGLoc, bsGGAll = bsGGAll)
                })
                bsDL <- sapply(hetVar, "[[", "bsDLoc")
                bsgL <- sapply(hetVar, "[[", "bsgLoc")
                bsGL <- sapply(hetVar, "[[", "bsGLoc")
                bsGGL <- sapply(hetVar, "[[", "bsGGLoc")
                bsDA <- sapply(hetVar, "[[", "bsDAll")
                bsgA <- sapply(hetVar, "[[", "bsgAll")
                bsGA <- sapply(hetVar, "[[", "bsGAll")
                bsGGA <- sapply(hetVar, "[[", "bsGGAll")
                bsDL <- diffCalcbCor(dLoc, bsDL)
                bsgL <- diffCalcbCor(gLoc, bsgL)
                bsGL <- diffCalcbCor(GLoc, bsGL)
                bsGGL <- diffCalcbCor(GGLoc, bsGGL)
                bsDA <- diffCalcbCor(dGlb, bsDA)
                bsgA <- diffCalcbCor(gGlb, bsgA)
                bsGA <- diffCalcbCor(GGlb, bsGA)
                bsGGA <- diffCalcbCor(GGGlb, bsGGA)
                if (fst) {
                        LCI$D <- apply(bsDL, 2, quantile, probs = alpha/2, 
                                       na.rm = TRUE)
                        LCI$gst <- apply(bsgL, 2, quantile, probs = alpha/2, 
                                         na.rm = TRUE)
                        LCI$Gst <- apply(bsGL, 2, quantile, probs = alpha/2, 
                                         na.rm = TRUE)
                        LCI$GGst <- apply(bsGGL, 2, quantile, probs = alpha/2, 
                                          na.rm = TRUE)
                        UCI$D <- apply(bsDL, 2, quantile, probs = 1 - (alpha/2), 
                                       na.rm = TRUE)
                        UCI$gst <- apply(bsgL, 2, quantile, probs = 1 - (alpha/2), 
                                         na.rm = TRUE)
                        UCI$Gst <- apply(bsGL, 2, quantile, probs = 1 - (alpha/2), 
                                         na.rm = TRUE)
                        UCI$GGst <- apply(bsGGL, 2, quantile, probs = 1 - 
                                                  (alpha/2), na.rm = TRUE)
                }
                else {
                        LCI <- data.frame(D = apply(bsDL, 2, quantile, probs = alpha/2, 
                                                    na.rm = TRUE))
                        LCI$gst <- apply(bsgL, 2, quantile, probs = alpha/2, 
                                         na.rm = TRUE)
                        LCI$Gst <- apply(bsGL, 2, quantile, probs = alpha/2, 
                                         na.rm = TRUE)
                        LCI$GGst <- apply(bsGGL, 2, quantile, probs = alpha/2, 
                                          na.rm = TRUE)
                        UCI <- data.frame(D = apply(bsDL, 2, quantile, probs = 1 - 
                                                            (alpha/2), na.rm = TRUE))
                        UCI$gst <- apply(bsgL, 2, quantile, probs = 1 - (alpha/2), 
                                         na.rm = TRUE)
                        UCI$Gst <- apply(bsGL, 2, quantile, probs = 1 - (alpha/2), 
                                         na.rm = TRUE)
                        UCI$GGst <- apply(bsGGL, 2, quantile, probs = 1 - 
                                                  (alpha/2), na.rm = TRUE)
                }
                if (fst) {
                        glbBS <- data.frame(gst = bsgA, Gst = bsGA, GGst = bsGGA, 
                                            d = bsDA, fst = bsFstA, fit = bsFitA, fis = bsFisA)
                        rm(bsFstA, bsFisA, bsFitA, bsDA, bsgA, bsGA, bsGGA)
                }
                else {
                        glbBS <- data.frame(gst = bsgA, Gst = bsGA, GGst = bsGGA, 
                                            d = bsDA)
                        rm(bsDA, bsgA, bsGA, bsGGA)
                }
                glbLCI <- apply(glbBS, 2, quantile, probs = alpha/2, 
                                na.rm = TRUE)
                glbUCI <- apply(glbBS, 2, quantile, probs = 1 - (alpha/2), 
                                na.rm = TRUE)
                if (fst) {
                        glbOut <- data.frame(stat = c("gst", "Gst", 
                                                      "GGst", "D", "Fst", "Fit", 
                                                      "Fis"), actual = c(gGlb, GGlb, GGGlb, dGlb, 
                                                                         glbTheta, glbFit, glbFis), lower = glbLCI, upper = glbUCI, 
                                             row.names = NULL)
                }
                else {
                        glbOut <- data.frame(stat = c("gst", "Gst", 
                                                      "GGst", "D"), actual = c(gGlb, GGlb, 
                                                                               GGGlb, dGlb), lower = glbLCI, upper = glbUCI, 
                                             row.names = NULL)
                }
        }
        popnms <- sapply(ip$indnms, "[[", 1)
        pwpops <- paste(popnms[pw[1, ]], " vs ", popnms[pw[2, 
        ]])
        if (pairwise || bs_pairwise) {
                preNharm <- lapply(preStats$indtyp, diffCalcHarm, pw = pw - 
                                           1)
                preStats$sHarm <- preNharm
                nbshths <- mapply(pwHCalc, af = preStats$alOut, sHarm = preStats$sHarm, 
                                  MoreArgs = list(pw = pw - 1), SIMPLIFY = FALSE)
                nbshths <- lapply(c("hsEst", "htEst"), function(x) {
                        lapply(nbshths, "[[", x)
                })
                names(nbshths) <- c("hsEst", "htEst")
                hsMn <- vapply(1:ncol(pw), function(i) {
                        mean(vapply(nbshths$hsEst, "[[", i, FUN.VALUE = numeric(1)), 
                             na.rm = TRUE)
                }, FUN.VALUE = numeric(1))
                htMn <- vapply(1:ncol(pw), function(i) {
                        mean(vapply(nbshths$htEst, "[[", i, FUN.VALUE = numeric(1)), 
                             na.rm = TRUE)
                }, FUN.VALUE = numeric(1))
                pwDLoc <- mapply(dCalc, ht = nbshths$htEst, hs = nbshths$hsEst)
                if(is.null(dim(pwDLoc))) return(NULL) ################# added for bugfix on ds102a, with only 2 pops
                pwDall <- apply(pwDLoc, 1, function(x) {
                        mnD <- mean(x, na.rm = TRUE)
                        vrD <- var(x, na.rm = TRUE)
                        1/((1/mnD) + ((vrD * ((1/mnD)^3))))
                })
                pwgLoc <- mapply(gstCalc, ht = nbshths$htEst, hs = nbshths$hsEst)
                pwgAll <- gstCalc(ht = htMn, hs = hsMn)
                pwGLoc <- mapply(GstCalc, ht = nbshths$htEst, hs = nbshths$hsEst)
                pwGAll <- GstCalc(ht = htMn, hs = hsMn)
                pwGGLoc <- mapply(GgstCalc, ht = nbshths$htEst, hs = nbshths$hsEst)
                pwGGAll <- GgstCalc(ht = htMn, hs = hsMn)
                if (fst) {
                        hsum <- lapply(preStats$hsum, function(x) {
                                pwTabMerge(x, pw - 1)
                        })
                        hsum <- lapply(hsum, function(x) {
                                x <- lapply(x, function(y) {
                                        y <- y[!(names(y) == "NA")]
                                        return(y)
                                })
                        })
                        pwVar <- mapply("pwWCcpp", hsum1 = hsum, af1 = preStats$alOut, 
                                        indtyp1 = preStats$indtyp, MoreArgs = list(pw = pw - 
                                                                                           1), SIMPLIFY = FALSE)
                        pwVar <- lapply(c("a", "b", "c"), 
                                        function(x) {
                                                return(lapply(pwVar, "[[", x))
                                        })
                        names(pwVar) <- c("a", "b", "cdat")
                        pwFstLoc <- mapply(thetaCalc, a = pwVar$a, b = pwVar$b, 
                                           cdat = pwVar$cdat)
                        mnA <- vapply(1:ncol(pw), function(i) {
                                mean(vapply(pwVar$a, "[", i, FUN.VALUE = numeric(1)), 
                                     na.rm = TRUE)
                        }, FUN.VALUE = numeric(1))
                        mnB <- vapply(1:ncol(pw), function(i) {
                                mean(vapply(pwVar$b, "[", i, FUN.VALUE = numeric(1)), 
                                     na.rm = TRUE)
                        }, FUN.VALUE = numeric(1))
                        mnC <- vapply(1:ncol(pw), function(i) {
                                mean(vapply(pwVar$cdat, "[", i, FUN.VALUE = numeric(1)), 
                                     na.rm = TRUE)
                        }, FUN.VALUE = numeric(1))
                        pwFstAll <- mnA/(mnA + mnB + mnC)
                }
        }
        if (bs_pairwise) {
                sHarm <- lapply(indtyp, function(x) {
                        lapply(x, diffCalcHarm, pw = pw - 1)
                })
                rm(indtyp)
                listAdd <- function(bsDat, sHarm) {
                        bsDat$sHarm <- sHarm
                        bsDat$nHarm <- NULL
                        return(bsDat)
                }
                bsDat <- mapply(listAdd, bsDat = bsDat, sHarm = sHarm, 
                                SIMPLIFY = FALSE)
                rm(sHarm)
                hths <- lapply(bsDat, function(x) {
                        lapply(1:length(x$alOut), function(i) {
                                pwHCalc(x$alOut[[i]], sHarm = x$sHarm[[i]], pw = pw - 
                                                1)
                        })
                })
                hths <- lapply(hths, function(x) {
                        hsEst <- lapply(x, "[[", 1)
                        htEst <- lapply(x, "[[", 2)
                        return(list(hsEst = hsEst, htEst = htEst))
                })
                pwDLocbs <- sapply(hths, function(x) {
                        return(mapply(dCalc, ht = x$htEst, hs = x$hsEst, 
                                      SIMPLIFY = TRUE))
                }, simplify = "array")
                pwDAllbs <- apply(pwDLocbs, c(1, 3), function(x) {
                        mn <- mean(x, na.rm = TRUE)
                        vr <- var(x, na.rm = TRUE)
                        1/((1/mn) + ((vr * ((1/mn)^3))))
                })
                rm(pwDLocbs)
                pwgAllbs <- sapply(hths, function(x) {
                        ht <- rowMeans(matrix(unlist(x$htEst), ncol = length(x$htEst)))
                        hs <- rowMeans(matrix(unlist(x$hsEst), ncol = length(x$hsEst)))
                        return(gstCalc(ht, hs))
                }, simplify = "array")
                pwGAllbs <- sapply(hths, function(x) {
                        ht <- rowMeans(matrix(unlist(x$htEst), ncol = length(x$htEst)))
                        hs <- rowMeans(matrix(unlist(x$hsEst), ncol = length(x$hsEst)))
                        return(GstCalc(ht, hs))
                }, simplify = "array")
                pwGGAllbs <- sapply(hths, function(x) {
                        ht <- rowMeans(matrix(unlist(x$htEst), ncol = length(x$htEst)))
                        hs <- rowMeans(matrix(unlist(x$hsEst), ncol = length(x$hsEst)))
                        return(GgstCalc(ht, hs))
                }, simplify = "array")
                if (fst) {
                        wcVar <- lapply(bsDat, function(x) {
                                x$hsum <- lapply(x$hsum, pwTabMerge, pw = pw - 
                                                         1)
                                return(mapply("pwWCcpp", hsum1 = x$hsum, 
                                              indtyp1 = x$indtyp, af1 = x$alOut, MoreArgs = list(pw = pw - 
                                                                                                         1), SIMPLIFY = FALSE))
                        })
                        pwFstAllbs <- sapply(wcVar, function(x) {
                                a <- rowMeans(sapply(x, "[[", "a"))
                                b <- rowMeans(sapply(x, "[[", "b"))
                                cdat <- rowMeans(sapply(x, "[[", "c"))
                                return(thetaCalc(a, b, cdat))
                        }, simplify = "array")
                        rm(wcVar, bsDat)
                        z <- gc()
                        pwFstAllbs <- t(mapply(diffCalcbCor, act = pwFstAll, 
                                               bs = split(pwFstAllbs, row(pwFstAllbs))))
                        pwFstAllLCI <- apply(pwFstAllbs, 1, quantile, prob = 0.025, 
                                             na.rm = TRUE)
                        pwFstAllUCI <- apply(pwFstAllbs, 1, quantile, prob = 0.975, 
                                             na.rm = TRUE)
                }
                pwDAllbs <- t(mapply(diffCalcbCor, act = pwDall, bs = split(pwDAllbs, 
                                                                            row(pwDAllbs))))
                pwDAllLCI <- apply(pwDAllbs, 1, quantile, prob = 0.025, 
                                   na.rm = TRUE)
                pwDAllUCI <- apply(pwDAllbs, 1, quantile, prob = 0.975, 
                                   na.rm = TRUE)
                pwgAllbs <- t(mapply(diffCalcbCor, act = pwgAll, bs = split(pwgAllbs, 
                                                                            row(pwgAllbs))))
                pwgAllLCI <- apply(pwgAllbs, 1, quantile, prob = 0.025, 
                                   na.rm = TRUE)
                pwgAllUCI <- apply(pwgAllbs, 1, quantile, prob = 0.975, 
                                   na.rm = TRUE)
                pwGAllbs <- t(mapply(diffCalcbCor, act = pwGAll, bs = split(pwGAllbs, 
                                                                            row(pwGAllbs))))
                pwGAllLCI <- apply(pwGAllbs, 1, quantile, prob = 0.025, 
                                   na.rm = TRUE)
                pwGAllUCI <- apply(pwGAllbs, 1, quantile, prob = 0.975, 
                                   na.rm = TRUE)
                pwGGAllbs <- t(mapply(diffCalcbCor, act = pwGAll, bs = split(pwGGAllbs, 
                                                                             row(pwGGAllbs))))
                pwGGAllLCI <- apply(pwGGAllbs, 1, quantile, prob = 0.025, 
                                    na.rm = TRUE)
                pwGGAllUCI <- apply(pwGGAllbs, 1, quantile, prob = 0.975, 
                                    na.rm = TRUE)
        }
        op <- list(std_stats = est)
        if (bs_locus) {
                op$global_bs <- data.frame(stat = glbOut[, 1], round(glbOut[, 
                                                                            -1], 4))
        }
        if (bs_locus) {
                if (fst) {
                        statnms <- c("gst", "Gst", "GGst", 
                                     "D", "Fst", "Fis", "Fit")
                }
                else {
                        statnms <- c("gst", "Gst", "GGst", 
                                     "D")
                }
                locCI <- lapply(statnms, function(x) {
                        return(data.frame(locus = ip$locs, actual = est[-nrow(est), 
                                                                        x], lower = round(LCI[, x], 4), upper = round(UCI[, 
                                                                                                                          x], 4), row.names = NULL))
                })
                names(locCI) <- statnms
                op$bs_locus <- locCI
                rm(locCI)
                z <- gc(reset = TRUE)
        }
        if (pairwise) {
                if (fst) {
                        pwgLoc <- data.frame(round(t(pwgLoc), 4))
                        dimnames(pwgLoc) <- list(ip$locs, pwpops)
                        op$pw_locus$gst <- pwgLoc
                        rm(pwgLoc)
                        pwGLoc <- data.frame(round(t(pwGLoc), 4))
                        dimnames(pwGLoc) <- list(ip$locs, pwpops)
                        op$pw_locus$Gst <- pwGLoc
                        rm(pwGLoc)
                        pwGGLoc <- data.frame(round(t(pwGGLoc), 4))
                        dimnames(pwGGLoc) <- list(ip$locs, pwpops)
                        op$pw_locus$GGst <- pwGGLoc
                        rm(pwGGLoc)
                        pwDLoc <- data.frame(round(t(pwDLoc), 4))
                        dimnames(pwDLoc) <- list(ip$locs, pwpops)
                        op$pw_locus$D <- pwDLoc
                        rm(pwDLoc)
                        pwFstLoc <- data.frame(round(t(pwFstLoc), 4))
                        dimnames(pwFstLoc) <- list(ip$locs, pwpops)
                        op$pw_locus$Fst <- pwFstLoc
                        rm(pwFstLoc)
                }
                else {
                        pwgLoc <- data.frame(round(t(pwgLoc), 4))
                        dimnames(pwgLoc) <- list(ip$locs, pwpops)
                        op$pw_locus$gst <- pwgLoc
                        rm(pwgLoc)
                        pwGLoc <- data.frame(round(t(pwGLoc), 4))
                        dimnames(pwGLoc) <- list(ip$locs, pwpops)
                        op$pw_locus$Gst <- pwGLoc
                        rm(pwGLoc)
                        pwGGLoc <- data.frame(round(t(pwGGLoc), 4))
                        dimnames(pwGGLoc) <- list(ip$locs, pwpops)
                        op$pw_locus$GGst <- pwGGLoc
                        rm(pwGGLoc)
                        pwDLoc <- data.frame(round(t(pwDLoc), 4))
                        dimnames(pwDLoc) <- list(ip$locs, pwpops)
                        op$pw_locus$D <- pwDLoc
                        rm(pwDLoc)
                }
                if (fst) {
                        op$pairwise <- list(gst = matrix(NA, nrow = np, ncol = np))
                        op$pairwise$gst[lower.tri(op$pairwise$gst)] <- round(pwgAll, 
                                                                             4)
                        dimnames(op$pairwise$gst) <- list(popnms, popnms)
                        op$pairwise$Gst <- op$pairwise$gst
                        op$pairwise$Gst[lower.tri(op$pairwise$Gst)] <- round(pwGAll, 
                                                                             4)
                        dimnames(op$pairwise$Gst) <- list(popnms, popnms)
                        op$pairwise$GGst <- op$pairwise$gst
                        op$pairwise$GGst[lower.tri(op$pairwise$GGst)] <- round(pwGGAll, 
                                                                               4)
                        dimnames(op$pairwise$GGst) <- list(popnms, popnms)
                        op$pairwise$D <- op$pairwise$gst
                        op$pairwise$D[lower.tri(op$pairwise$D)] <- round(pwDall, 
                                                                         4)
                        dimnames(op$pairwise$D) <- list(popnms, popnms)
                        op$pairwise$Fst <- op$pairwise$gst
                        op$pairwise$Fst[lower.tri(op$pairwise$Fst)] <- round(pwFstAll, 
                                                                             4)
                        dimnames(op$pairwise$Fst) <- list(popnms, popnms)
                }
                else {
                        op$pairwise <- list(gst = matrix(NA, nrow = np, ncol = np))
                        op$pairwise$gst[lower.tri(op$pairwise$gst)] <- round(pwgAll, 
                                                                             4)
                        dimnames(op$pairwise$gst) <- list(popnms, popnms)
                        op$pairwise$Gst <- op$pairwise$gst
                        op$pairwise$Gst[lower.tri(op$pairwise$Gst)] <- round(pwGAll, 
                                                                             4)
                        dimnames(op$pairwise$Gst) <- list(popnms, popnms)
                        op$pairwise$GGst <- op$pairwise$gst
                        op$pairwise$GGst[lower.tri(op$pairwise$GGst)] <- round(pwGGAll, 
                                                                               4)
                        dimnames(op$pairwise$GGst) <- list(popnms, popnms)
                        op$pairwise$D <- op$pairwise$gst
                        op$pairwise$D[lower.tri(op$pairwise$D)] <- round(pwDall, 
                                                                         4)
                        dimnames(op$pairwise$D) <- list(popnms, popnms)
                }
        }
        if (bs_pairwise) {
                if (fst) {
                        op$bs_pairwise <- list(gst = data.frame(populations = pwpops, 
                                                                actual = round(pwgAll, 4), lower = round(pwgAllLCI, 
                                                                                                         4), upper = round(pwgAllUCI, 4), row.names = NULL))
                        rm(pwgAll, pwgAllLCI, pwgAllUCI)
                        op$bs_pairwise$Gst <- data.frame(populations = pwpops, 
                                                         actual = round(pwGAll, 4), lower = round(pwGAllLCI, 
                                                                                                  4), upper = round(pwGAllUCI, 4), row.names = NULL)
                        rm(pwGAll, pwGAllLCI, pwGAllUCI)
                        op$bs_pairwise$GGst <- data.frame(populations = pwpops, 
                                                          actual = round(pwGGAll, 4), lower = round(pwGGAllLCI, 
                                                                                                    4), upper = round(pwGGAllUCI, 4), row.names = NULL)
                        rm(pwGGAll, pwGGAllLCI, pwGGAllUCI)
                        op$bs_pairwise$D <- data.frame(populations = pwpops, 
                                                       actual = round(pwDall, 4), lower = round(pwDAllLCI, 
                                                                                                4), upper = round(pwDAllUCI, 4), row.names = NULL)
                        rm(pwDall, pwDAllLCI, pwDAllUCI)
                        op$bs_pairwise$Fst <- data.frame(populations = pwpops, 
                                                         actual = round(pwFstAll, 4), lower = round(pwFstAllLCI, 
                                                                                                    4), upper = round(pwFstAllUCI, 4), row.names = NULL)
                        rm(pwFstAll, pwFstAllLCI, pwFstAllUCI)
                }
                else {
                        op$bs_pairwise <- list(gst = data.frame(populations = pwpops, 
                                                                actual = round(pwgAll, 4), lower = round(pwgAllLCI, 
                                                                                                         4), upper = round(pwgAllUCI, 4), row.names = NULL))
                        rm(pwgAll, pwgAllLCI, pwgAllUCI)
                        op$bs_pairwise$Gst <- data.frame(populations = pwpops, 
                                                         actual = round(pwGAll, 4), lower = round(pwGAllLCI, 
                                                                                                  4), upper = round(pwGAllUCI, 4), row.names = NULL)
                        rm(pwGAll, pwGAllLCI, pwGAllUCI)
                        op$bs_pairwise$GGst <- data.frame(populations = pwpops, 
                                                          actual = round(pwGGAll, 4), lower = round(pwGGAllLCI, 
                                                                                                    4), upper = round(pwGGAllUCI, 4), row.names = NULL)
                        rm(pwGGAll, pwGGAllLCI, pwGGAllUCI)
                        op$bs_pairwise$D <- data.frame(populations = pwpops, 
                                                       actual = round(pwDall, 4), lower = round(pwDAllLCI, 
                                                                                                4), upper = round(pwDAllUCI, 4), row.names = NULL)
                        rm(pwDall, pwDAllLCI, pwDAllUCI)
                        z <- gc(reset = TRUE)
                }
        }
        if (!is.null(outfile)) {
                opf <- paste(getwd(), "/", outfile, "-[diffCalc]/", 
                             sep = "")
                dir.create(opf, showWarnings = FALSE)
                outnms <- names(op)
                out <- sapply(outnms, function(x) {
                        if (x == "std_stats" || x == "global_bs") {
                                ot <- paste(colnames(op[x][[1]]), collapse = "\t")
                                preot <- apply(op[x][[1]], 1, paste, collapse = "\t")
                                ot <- c(ot, preot)
                                writeLines(paste(ot, collapse = "\n"), 
                                           paste(opf, x, ".txt", sep = ""))
                                ot <- NULL
                        }
                        else if (x == "pairwise") {
                                statnms <- names(op[x][[1]])
                                ot <- lapply(statnms, function(y) {
                                        dat <- op[x][[1]][y][[1]]
                                        dat[is.na(dat)] <- ""
                                        dimnames(dat) <- list(NULL, NULL)
                                        opt <- apply(dat, 1, paste0, collapse = "\t", 
                                                     na.rm = "")
                                        popnmsOut <- paste(popnms, "\t", sep = "")
                                        opt <- mapply(paste, popnmsOut, opt, MoreArgs = list(collapse = "\t"))
                                        opt <- c(y, "", paste("pops", paste(popnms, 
                                                                            collapse = "\t"), sep = "\t"), 
                                                 opt, "")
                                        return(opt)
                                })
                                if (fst) {
                                        ot <- c("Pairwise stats", "gst = Nei & Chesser, (1983)", 
                                                "Gst = Hedrick, (2005)", "GGst = Meirmans & Hedrick, (2011)", 
                                                "D = Jost, (2008)", "Fst = Weir & Cockerham's theta, (1984)", 
                                                "", unlist(ot))
                                }
                                else {
                                        ot <- c("Pairwise stats", "gst = Nei & Chesser, (1983)", 
                                                "Gst = Hedrick, (2005)", "GGst = Meirmans & Hedrick, (2011)", 
                                                "D = Jost, (2008)", "", unlist(ot))
                                }
                                writeLines(paste(ot, collapse = "\n"), 
                                           paste(opf, x, ".txt", sep = ""))
                                ot <- NULL
                        }
                        else if (x == "bs_locus") {
                                statnms <- names(op[x][[1]])
                                ot <- lapply(statnms, function(y) {
                                        ot1 <- c("", y, "", paste(colnames(op[x][[1]][y][[1]]), 
                                                                  collapse = "\t"))
                                        ot2 <- apply(op[x][[1]][y][[1]], 1, paste, 
                                                     collapse = "\t")
                                        return(c(ot1, ot2))
                                })
                                if (fst) {
                                        ot <- c("Locus 95% CIs", "", "gst = Nei & Chesser, 1983", 
                                                "Gst = Hedrick, 2005", "GGst = Meirmans & Hedrick, (2011)", 
                                                "D = Jost, 2008", "Fst = Weir & Cockerham, 1984", 
                                                "Fis = Weir & Cockerham, 1984", "Fit = Weir & Cockerham, 1984", 
                                                "", unlist(ot))
                                }
                                else {
                                        ot <- c("Locus 95% CIs", "", "gst = Nei & Chesser, 1983", 
                                                "Gst = Hedrick, 2005", "GGst = Meirmans & Hedrick, (2011)", 
                                                "D = Jost, 2008", unlist(ot))
                                }
                                writeLines(paste(ot, collapse = "\n"), 
                                           paste(opf, x, ".txt", sep = ""))
                                ot <- NULL
                        }
                        else if (x == "pw_locus") {
                                statnms <- names(op[x][[1]])
                                ot <- lapply(statnms, function(y) {
                                        ot1 <- c("", y, "", paste("Loci", 
                                                                  paste(colnames(op[x][[1]][y][[1]]), collapse = "\t"), 
                                                                  sep = "\t"))
                                        ot2 <- apply(op[x][[1]][y][[1]], 1, paste, 
                                                     collapse = "\t")
                                        ot2 <- mapply(paste, rownames(op[x][[1]][y][[1]]), 
                                                      ot2, MoreArgs = list(sep = "\t"))
                                        return(c(ot1, ot2))
                                })
                                if (fst) {
                                        ot <- c("Locus Pairwise estimates", "", 
                                                "gst = Nei & Chesser, 1983", "Gst = Hedrick, 2005", 
                                                "GGst = Meirmans & Hedrick, (2011)", 
                                                "D = Jost, 2008", "Fst = Weir & Cockerham, 1984", 
                                                unlist(ot))
                                }
                                else {
                                        ot <- c("Locus Pairwise estimates", "", 
                                                "gst = Nei & Chesser, 1983", "Gst = Hedrick, 2005", 
                                                "GGst = Meirmans & Hedrick, (2011)", 
                                                "D = Jost, 2008", unlist(ot))
                                }
                                writeLines(paste(ot, collapse = "\n"), 
                                           paste(opf, x, ".txt", sep = ""))
                                ot <- NULL
                        }
                        else if (x == "bs_pairwise") {
                                statnms <- names(op[x][[1]])
                                ot <- lapply(statnms, function(y) {
                                        ot1 <- c("", y, "", paste(colnames(op[x][[1]][y][[1]]), 
                                                                  collapse = "\t"))
                                        ot2 <- apply(op[x][[1]][y][[1]], 1, paste, 
                                                     collapse = "\t")
                                        return(c(ot1, ot2))
                                })
                                if (fst) {
                                        ot <- c("Pairwise 95% CIs", "", 
                                                "gst = Nei & Chesser, 1983", "Gst = Hedrick, 2005", 
                                                "GGst = Meirmans & Hedrick, (2011)", 
                                                "D = Jost, 2008", "Fst = Weir & Cockerham, 1984", 
                                                unlist(ot))
                                }
                                else {
                                        ot <- c("Pairwise 95% CIs", "", 
                                                "gst = Nei & Chesser, 1983", "Gst = Hedrick, 2005", 
                                                "GGst = Meirmans & Hedrick, (2011)", 
                                                "D = Jost, 2008", unlist(ot))
                                }
                                writeLines(paste(ot, collapse = "\n"), 
                                           paste(opf, x, ".txt", sep = ""))
                                ot <- NULL
                        }
                })
                rm(out)
        }
        return(op)
}
