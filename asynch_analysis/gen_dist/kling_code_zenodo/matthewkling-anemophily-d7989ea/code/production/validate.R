


validate_dataset <- function(x){
        #browser()
        result <- try(validate_alignment(validate_populations(x$pop), 
                                        validate_genepop(x$gen)))
        x$valid <- class(result) != "try-error"
        x
}



# check populations file
validate_populations <- function(infile){
        
        p <- read_csv(infile, col_types = cols())
        
        if(! ncol(p) %in% 4:5) stop("incorrect number of columns in population location file")
        if(all.equal(names(p)[1:4], c("ID", "name", "longitude", "latitude")) != TRUE &
           all.equal(names(p)[1:4], c("ID", "names", "longitude", "latitude")) != TRUE){
                stop("incorrect variable names in population location file")
        }
        if(ncol(p)==5 & names(p)[5] != "n") stop("incorrect variable names in population location file")
        
        if(any(is.na(p))) stop("NA values in location file")
        
        sp <- p
        coordinates(sp) <- c("longitude", "latitude")
        return(p)
}


# check genepop file
validate_genepop <- function(infile){
        gp <- readGenepop(infile)
        hm <- haploMigrate(infile)
        return(gp)
}

# check genepop-populations agreement
validate_alignment <- function(pop, gen){
        #browser()
        if(nrow(pop) != gen$npops) stop("populations.csv and genepop.txt have different number of populations")
        
        gen_ids <- gen$pop_names %>% str_split("_", simplify=T) 
        gen_ids <- gen_ids[,1] %>% sub("p", "", .) %>% as.integer()
        
        if(all.equal(sort(pop$ID), sort(gen_ids)) != TRUE) stop("Population IDs do not match between genepop and populations file.")
        
        if("n" %in% names(pop)){
                if(!all.equal(gen$pop_sizes[order(gen_ids)], pop$n[order(pop$ID)])){
                stop("Population sizes do not match.")}}
        return("AOK")
}
