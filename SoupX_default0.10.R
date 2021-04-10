### --- Rscript takes removes the ambient RNA from 10x counts using SoupX  --- ###
#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
require(optparse)
options <- list(
make_option(c("-p", "--path"), action = "store", default = getwd(), type = "character", help=" path for  outs folder from the cellranger .")
)
arguments <- parse_args(OptionParser(option_list = options))
setwd(arguments$p)
cat(arguments$p)

dir.create("SoupX")

message("Loading the required packages ....")
suppressPackageStartupMessages(library(SoupX))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(matrixStats))


message("Creating the soup Channel ...")
sc = load10X(arguments$p, keepDroplets = TRUE)
sc = setContaminationFraction(sc, 0.10)  ### High contamination fraction because of Nulcei Isolation
head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 20)
out = adjustCounts(sc,roundToInt=TRUE)
cntSoggy = rowSums(sc$toc > 0)
cntStrained = rowSums(out > 0)
mostZeroed = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
tail(sort(rowSums(sc$toc > out)/rowSums(sc$toc > 0)), n = 20)
DropletUtils:::write10xCounts("./soupX-contamination-fraction-0.10", out)
rm(list=ls())


