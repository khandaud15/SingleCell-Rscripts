### Extract RNA and ADT counts from Cite-seq h5 file produced from cellranger and convert to sparse matrix format
#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
require(optparse)
options <- list(
make_option(c("-f", "--file"), action = "store", default = getwd(), type = "character", help="h5 file containg rna and adt counts produced from cellranger   .")
)
arguments <- parse_args(OptionParser(option_list = options))
setwd(getwd())

cat(arguments$f )

### --- loading the required pacakges --- ###
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(DropletUtils))

message("Reading h5 file ...")
inputdata.10x <- Read10X_h5(arguments$f)
message("writting adt counts as mtx in ADT folder ...")
write10xCounts(x = inputdata.10x[["Antibody Capture"]], path = "ADT")
message("writting rna counts as mtx in RNA folder ...")
write10xCounts(x = inputdata.10x[["Gene Expression"]],  path = "RNA")


