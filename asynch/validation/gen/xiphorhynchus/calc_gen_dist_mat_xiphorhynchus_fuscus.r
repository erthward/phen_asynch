library(ape)

# NOTE: Quintero et al. used the GTR+G substitution model
#       when aligning and calclating distances, but I wasn't
#       sure how to do that, and all sequences had non-missing data
#       for everything but the first ~10% of sites, so I just
#       pairwise deleted all those sites with any missing data
#       before distance calculation

# load the aligned sequences
seqs = read.FASTA('./all_xiphorhynchus_fuscus_ALIGNED.fasta')

# calculate genetic distance matrix
dist = dist.dna(x=seqs, model='K80', pairwise.deletion=T, as.matrix=T)

# write out matrix
write.table(dist,
            './all_xiphorhynchus_fuscus_dist.csv',
            sep=',',
            row.names=T,
            col.names=T,
           )

