### --- Rscript takes removes the ambient RNA from 10x counts using SoupX  --- ###
#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
require(optparse)
options <- list(
make_option(c("-p", "--path"), action = "store", default = getwd(), type = "character", help=" path for  outs folder from the cellranger ."),
make_option(c("-c", "--fraction"), type="integer", default=0.10, help=" contamination fraction percentage, default is 0.10", metavar="number")
)
arguments <- parse_args(OptionParser(option_list = options))
setwd(arguments$p)
cat(arguments$p, arguments$c)

message("Loading the required packages ....")
suppressPackageStartupMessages(library(SoupX))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(matrixStats))


message("Creating the soup Channel ...")
sc = load10X(arguments$p, keepDroplets = TRUE)
sc = setContaminationFraction(sc, arguments$c)
head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 20)
out = adjustCounts(sc,roundToInt=TRUE)
cntSoggy = rowSums(sc$toc > 0)
cntStrained = rowSums(out > 0)
mostZeroed = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
tail(sort(rowSums(sc$toc > out)/rowSums(sc$toc > 0)), n = 20)
message("writting the mtx files to soupX-contamination-fraction-0.10 directory ...")
DropletUtils:::write10xCounts("./soupX-contamination-fraction-0.10", out)
rm(list=ls())


