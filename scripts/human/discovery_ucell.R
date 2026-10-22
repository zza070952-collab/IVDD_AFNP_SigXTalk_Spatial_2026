#!/usr/bin/env Rscript
# Own-code adapter of the executed discovery UCell wrapper, not package source.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2L) stop("Usage: discovery_ucell.R REBUILD_WORKSPACE FROZEN_PROGRAMS_TSV")
root <- normalizePath(args[[1]], mustWork=TRUE)
sigfile <- normalizePath(args[[2]], mustWork=TRUE)
suppressPackageStartupMessages({library(Matrix);library(hdf5r);library(UCell)})
if (as.character(packageVersion("UCell")) != "2.4.0") {
  stop("Historical discovery scoring used UCell2.4.0; restore that isolated environment, or disclose/version any separately authorised alternative")
}
f <- H5File$new(file.path(root,"work/discovery_identity.h5ad"),mode="r")
shape <- h5attr(f[["X"]],"shape")
x <- sparseMatrix(i=as.integer(f[["X/indices"]][])+1L,p=as.integer(f[["X/indptr"]][]),
                  x=as.numeric(f[["X/data"]][]),dims=rev(as.integer(shape)),
                  dimnames=list(f[["var/_index"]][],f[["obs/_index"]][]))
f$close_all()
signatures <- read.delim(sigfile,stringsAsFactors=FALSE)
features <- split(signatures$gene,signatures$program)
names(features) <- paste0(names(features),"_top25")
ranks <- UCell::StoreRankings_UCell(x,maxRank=1500,ncores=1,chunk.size=200,ties.method="average")
scores <- UCell::ScoreSignatures_UCell(features=features,precalc.ranks=ranks,maxRank=1500,
                                       ncores=1,chunk.size=200,ties.method="average")
con <- gzfile(file.path(root,"tables/discovery_programs_UCell.tsv.gz"),"wt")
write.table(data.frame(cell_id=rownames(scores),scores,check.names=FALSE),con,
            sep="\t",quote=FALSE,row.names=FALSE)
close(con)
writeLines(capture.output(sessionInfo()),file.path(root,"envs/UCell_discovery_session.txt"))
cat("Actual discovery UCell cell scores:",nrow(scores),"\n")
