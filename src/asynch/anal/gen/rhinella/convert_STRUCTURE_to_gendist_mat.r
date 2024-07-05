# NOTE: need to run this simple script interactively and answer all questions:
  # number of genotypes: 80
  # number of markers: 7674
  # col containing genotype labels: 1
  # col containing pop factor: 2
  # other optional cols: <Enter> (none)
  # row containing marker names: 0
  # genotypes coded by a single row? n

library(adegenet)

# read genetic data
# NOTE: appears it is in 2-row STRUCTURE format,
#       with each integer in a column indicating a unique allele for that locus
#       and -9 indicating missing data, based on the STRUCTURE2 documentation
#       (web.stanford.edu/group/pritchardlab/software/readme_structure2.pdf)
# NOTE: column 1 encodes locality (51 serial integers)
# NOTE: because values in the locus columns are all in [-9, 1, 2, 3, 4],
#       I suspect each int encodes a base, and for purposes of their analyses
#       they don't consider any mutation model, so we're just treating each
#       integer value as an alternative genotype for that locus
genind_obj = read.structure('./dryad_archive/2-Mantel_test/Rgranulosa_Mantel.str')

# write out simple Euclidean genetic distance matrix
dist_gen = as.matrix(dist(genind_obj))
write.table(dist_gen,
            './dryad_archive/2-Mantel_test/Rgranulosa_Mantel.csv',
            sep=',',
            row.names=T,
            col.names=T,
            )

