#! /usr/bin/env Rscript

# Author: Andy Saurin (andrew.saurin@univ-amu.fr)
#
# Simple RScript to calculate RPKMs and TPMs
# based on method for RPKM/TPM calculations shown in http://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/
#
# The input file is the output of featureCounts
#

rpkm <- function(counts, lengths) {
  pm <- sum(counts) /1e6
  rpm <- counts/pm
  rpm/(lengths/1000)
}

tpm <- function(counts, lengths) {
  rpk <- counts/(lengths/1000)
  coef <- sum(rpk) / 1e6
  rpk/coef
}


## read table from featureCounts output
args <- commandArgs(T)

tag <- tools::file_path_sans_ext(args[1])


cat('Reading in featureCounts data...')
ftr.cnt <- read.table(args[1], sep="\t", header=T, quote="") #Important to disable default quote behaviour or else genes with apostrophes will be taken as strings
cat(' Done\n')

if ( ncol(ftr.cnt) < 7 ) { 
	cat(' The input file is not the raw output of featureCounts (number of columns > 6) \n')
	quit('no')
}

lengths = ftr.cnt[,6]

counts <- ftr.cnt[,7:ncol(ftr.cnt)]

cat('Performing RPKM calculations...')

rpkms <- apply(counts, 2, function(x) rpkm(x, lengths) )
ftr.rpkm <- cbind(ftr.cnt[,1:6], rpkms)

rpkms <- apply(counts, 2, function(x) rpkm(x, lengths) )
ftr.rpkm <- cbind(ftr.cnt[,1:6], rpkms)
write.table(ftr.rpkm, file=paste0(tag, "_rpkm.txt"), sep="\t", row.names=FALSE, quote=FALSE)
cat(' Done.\n\tSaved as ')
cat ( paste0(tag, "_rpkm.txt", '\n') )

cat('Performing TPM calculations...')

tpms <- apply(counts, 2, function(x) tpm(x, lengths) )

ftr.tpm <- cbind(ftr.cnt[,1:6], tpms)

write.table(ftr.tpm, file=paste0(tag, "_tpm.txt"), sep="\t", row.names=FALSE, quote=FALSE)
cat(' Done.\n\tSaved as ')
cat ( paste0(tag, "_tpm.txt", '\n') )


quit('no')