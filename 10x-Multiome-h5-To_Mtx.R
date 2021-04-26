### Extract RNA from Multiome h5 file produced from cellranger and convert to sparse matrix format to be used by SoupX
#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
require(optparse)
options <- list(
make_option(c("-p", "--path"), action = "store", default = getwd(), type = "character", help=" path for  folder containing h5 file ."),
make_option(c("-f", "--file"), action = "store", default = getwd(), type = "character", help="h5 file containg rna and adt counts produced from cellranger   ."),
make_option(c("-n", "--name"), action = "store", default = getwd(), type = "character", help=" folder name to store mtx files .")
)
arguments <- parse_args(OptionParser(option_list = options))
cat(arguments$p, arguments$f, arguments$n)
setwd(arguments$p)

### --- loading the required pacakges --- ###
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(DropletUtils))

message("Reading h5 file ...")
inputdata.10x <- Read10X_h5(arguments$f)
message("writting rna counts as mtx in", sep = "" ,  arguments$n, sep = "" , "folder")
write10xCounts(x = inputdata.10x[["Gene Expression"]],  path = arguments$n)


