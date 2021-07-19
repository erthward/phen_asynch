library(adegenet)
library(PopGenReport)
#library(R.oo)
#library(diveRsity)
library(hierfstat)
#library(FinePop)
library(plyr)


##############################################################
# TODO:

    # 1. decide on the genetic distance metric(s) to use
        # for now, using pop-based F_ST (Weir and Cockerham 1984)
        #            and ind-based dist (Smouse and Peakall 1999)
        # But between PopGenReport, hierfstat, and FinePop, there
        # appear to be a whole ton!
        # Do I need to choose based on whether haploid or diploid?
        #                        or on   whether SNPs or STRs?

##############################################################


data_dir = "/home/deth/Desktop/UCB/research/projects/seasonality/seasonal_asynchrony/asynch_analysis/data/kling_gendata/"

# get all the genepop files
genepop_files = grep('gen$', list.files(data_dir), value=T)

#read_diploid_genepop_file = function(fn, ncode=3, output_format="table"){
#    genind_obj = read.genepop(fn, ncode)
#    return(genind_obj)
#}


# function to read either haploid or diploid genetic data
# from a genepop-formatted file into a genind object
read_genepop_file = function(fn, ncode=3, output_format='genind', delete_tmp=T){
    # read in the lines that can be coerced to a table
    lines = readLines(fn)
    newlines = c()
    for (line in lines){
        if (!line == 'title' & !line == 'pop'){
            if (startsWith(line, 'locus1, ')){
                line = gsub(",? ", ",", line)
                newlines = c(newlines, paste0('pop,', line))
            }
            else {
                line = gsub(",? ", ",", line)
                newlines = c(newlines, line)
            }
        }
    }
    # write to a temporary file
    writeLines(newlines, 'genepop.tmp')
    # read the tmp file to a table
    output = read.csv('genepop.tmp',
                      colClasses = c('character'),
                      stringsAsFactors = F)
    output$pop = gsub('_+$', '', output$pop)
    if (output_format == 'genind'){
        # or read the tmp file to a genind object
        pop_ind = ldply(strsplit(output$pop, split='_', fixed=T))
        pop = as.character(pop_ind[,1])
        inds = as.character(pop_ind[,2])
        # use length of string in 1st row, second col, to determine ploidy
        ploidy = nchar(output[1,2])/3
        output = df2genind(output[,2:ncol(output)],
                           ploidy=ploidy,
                           ind.name=output$pop,
                           pop=pop,
                           #sep='',
                           ncode=ncode)
    }
    if (delete_tmp){
        if (file.exists('genepop.tmp')) {
            #Delete file if it exists
            file.remove('genepop.tmp')
        }
    }
    return(output)
}


# return a pairwise genetic distance matrix, given a genind object
# NOTE: default is pairwise F_ST following Weir & Cockerham (1984)
calc_gen_dist = function(genind_obj, method='WC84'){
    if (method == 'smouse'){
        dist = as.matrix(gd.smouse(genind_obj))
    } else {
        dist = as.matrix(genet.dist(genind_obj, method=method))
    }
    return(dist)
}


for (filename in genepop_files){
    tryCatch({
        cat('\nNOW PROCESSING: ', filename, '\n\n')
        # read data
        genind = read_genepop_file(paste0(data_dir, filename), output_format='genind')

        # calculate population-based (Weir and Cockerham 1984) distance matrix
        dist.wc84 = calc_gen_dist(genind, method='WC84')
        # write to file
        write.csv(as.data.frame(dist.wc84),
                  paste0(data_dir, filename, '.wc84.csv'))

        # calculate individual-based (Smouse and Peakall, 1999) distance matrix
        dist.smouse99 = calc_gen_dist(genind, method='smouse')
        # write to file
        write.csv(as.data.frame(dist.smouse99),
                  paste0(data_dir, filename, '.smouse99.csv'))
    }, error = function(err) {
        cat('\tERROR THROWN:\n')
        print(err)
    }, finally = {cat('\n\tFINISHED\n-----------------------------------\n')})
}


#for (filename in genepop_files){
#    tryCatch({
#        cat('\nNOW PROCESSING: ', filename, '\n\n')
#        # try to read the file using length-3 allele codes
#        genind_obj = tryCatch({
#            read_diploid_genepop_file(paste0(data_dir, filename), 3)
#        }, error = function(err) {
#            read_haploid_genepop_file(paste0(data_dir, filename), 3)
#        }, finally = {NA})
#        if (is.na(genind_obj)){
#        }
#        else {
#            cat('\nWRITING: ', filename, '\n\n')
#
#            # TODO: if that fails, must be a diff format, so print out info to inspect
#
#            # calculate the distance matrices and save
#            #dm_kos = as.data.frame(as.matrix(gd.kosman(genind_obj)))
#            #write.csv(dm_kos, paste0(data_dir, filename, '.kosman.csv'))
#
#            dm_smo = as.data.frame(as.matrix(gd.smouse(genind_obj)))
#            write.csv(dm_smo, paste0(data_dir, filename, '.smouse.csv'))
#        }
#    }, error = function(err){
#    })
#}
