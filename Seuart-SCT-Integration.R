#!/usr/bin/env Rscript
### --- Seuart Integration on txt file in the directory rpca mode --- ####
suppressPackageStartupMessages(library(optparse))
require(optparse)
options <- list(
  make_option(c("-p", "--path"), action = "store", default = getwd(), type = "character", help=" path for  folder containing all the text files ."),
  make_option(c("-o", "--outname"), action = "store", default = "output", type = "character", help="Output file(s) base name of rds.")
)
arguments <- parse_args(OptionParser(option_list = options))
setwd(arguments$p)
cat(arguments$p,arguments$o)

message("creating an directory for seurat plots...")
# dir.create("SeuratIntegration")


message(" Loading the required packages .... ")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(hdf5r))
suppressPackageStartupMessages(library(ggplot2))


myFiles <- list.files(path=arguments$p, pattern="*.txt")
message("Reading all the text file and making list....")
myList <- lapply(myFiles, function(x) fread(x))
myList <- lapply(myList , column_to_rownames, "uid")
samples <- gsub(pattern = "_counts.txt", replacement = "", x = basename(myFiles))

### --- Number of cells in each samples --- ###
message("finding number of cells barcode in each sample ...")
names(myList) <- samples
for (i in 1:length(myList)){
  print(paste("no of cell barcodes in ", names(myList[i]), "are ", length(myList[[i]])))
}

### ---  create an seurat object list --- ###
scrna.list = list();
for (i in 1:length(myList)) {
  scrna.list[[i]] = CreateSeuratObject(counts = myList[[i]],  project=samples[i]);
  scrna.list[[i]][["DataSet"]] = samples[i];
}
names(scrna.list) <- samples

scrna.list <- lapply(X = scrna.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = scrna.list, nfeatures = 3000)
scrna.list <- PrepSCTIntegration(object.list = scrna.list, anchor.features = features)


immune.anchors <- FindIntegrationAnchors(object.list = scrna.list, normalization.method = "SCT", 
                                         anchor.features = features)
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")

immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE)
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:40)

saveRDS(immune.combined.sct, file = paste(arguments$o, "_SCT-Integrated.rds"))
