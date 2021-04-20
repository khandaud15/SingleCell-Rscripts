#!/usr/bin/env Rscript
### --- Seuart Integration on txt file in the directory rpca mode --- ####
suppressPackageStartupMessages(library(optparse))
require(optparse)
options <- list(
make_option(c("-p", "--path"), action = "store", default = getwd(), type = "character", help=" path for  folder containing all the text files ."),
make_option(c("-r", "--res"), type="integer", default=1, help="number of resolution for finding the clusters ", metavar="number"),
make_option(c("-n", "--npcs"), type="integer", default=30, help="number of pcs for PCA function ", metavar="number"),
make_option(c("-o", "--outname"), action = "store", default = "output", type = "character", help="Output file(s) base name of rds.")
)
arguments <- parse_args(OptionParser(option_list = options))
setwd(arguments$p)
cat(arguments$p,arguments$r, arguments$n, arguments$o)

message("creating an directory for seurat plots...")
dir.create("SeuratIntegration")


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

### ---  normalize and identify variable features for each dataset independently --- ###
ifnb.list <- lapply(X = scrna.list, FUN = function(x) {
       x <- NormalizeData(x)
       x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

#### ---  select features that are repeatedly variable across datasets for integration run PCA on each ### ---  
features <- SelectIntegrationFeatures(object.list = ifnb.list)
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
          x <- ScaleData(x, features = features, verbose = FALSE)
          x <- RunPCA(x, features = features, verbose = FALSE)
})

M.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features, reduction = "rpca")
#### --- this command creates an 'integrated' data assay #### --- #### --- 
M.combined <- IntegrateData(anchorset = M.anchors)



DefaultAssay(M.combined) <- "integrated"

#### ---  Run the standard workflow for visualization and clustering #### --- 
M.combined <- ScaleData(M.combined, verbose = FALSE)
M.combined <- RunPCA(M.combined, npcs = arguments$n, verbose = FALSE)
M.combined <- RunUMAP(M.combined, reduction = "pca", dims = 1:40)
M.combined <- FindNeighbors(M.combined, reduction = "pca", dims = 1:40)
M.combined <- FindClusters(M.combined, resolution = arguments$r)


### --- save RDS file for further analysis ---###
message("saving the rds file for downstream analysis ...")
saveRDS(M.combined, paste0(arguments$o, "_integrated.rds"))

### --- batch corrected umap plot --- ###
pdf("SeuratIntegration/Integrated-Batch_umap.pdf", width = 10.19, height = 7.61)
DimPlot(M.combined, group.by = "DataSet") + ggtitle("Integrated UMAP Plot")
dev.off()

pdf("SeuratIntegration/Integrated-Clusters_umap.pdf", width = 8.5, height = 7.6)
DimPlot(M.combined, group.by = "seurat_clusters", label = T) + NoLegend()
dev.off()


### --- findMarkers for each cluster --- ###
DefaultAssay(M.combined) <- "integrated"
M.markers <- FindAllMarkers(M.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(M.markers, "SeuratIntegration/Integration_Markers.txt", col.names = T, sep = "\t", row.names = F, quote = F)

sessionInfo()
quit()

