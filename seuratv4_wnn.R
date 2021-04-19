#!/usr/bin/env Rscript
### --- Seuart WNN based analysis of MultiOme dataset --- ####
suppressPackageStartupMessages(library(optparse))
require(optparse)
options <- list(
make_option(c("-p", "--path"), action = "store", default = getwd(), type = "character", help=" path for  folder containing all the text files ."),
make_option(c("-r", "--res"), type="integer", default=1, help="number of resolution for finding the clusters ", metavar="number"),
make_option(c("-o", "--outname"), action = "store", default = "output", type = "character", help="Output file(s) base name of rds.")
)
arguments <- parse_args(OptionParser(option_list = options))
setwd(arguments$p)
cat(arguments$p,arguments$r,  arguments$o)

#### SeuratV4  WNN based basic analysis for QC on mutiome H5 file
message(" Loading the required packages .... ")

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))


message("creating an directory for seurat plots...")
dir.create("Seurat_v4_wnn")
filtered_h5 <- "filtered_feature_bc_matrix.h5"
# the 10x hdf5 file contains both data types. 
inputdata.10x <- Read10X_h5(filtered_h5)

# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

# Create Seurat object
M <- CreateSeuratObject(counts = rna_counts, project = arguments$o)
M[["percent.mt"]] <- PercentageFeatureSet(M, pattern = "^MT-")

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

frag.file <- "atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
M[["ATAC"]] <- chrom_assay

### We perform basic QC based on the number of detected molecules for each modality as well as mitochondrial percentage
pdf("Seurat_v4_wnn/VoiLin.pdf", width = 14, height = 8)
VlnPlot(M, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0) + NoLegend()
dev.off()
### --- VoiLin plot with dots --- ###
pdf("Seurat_v4_wnn/VoiLin_with_dots.pdf", width = 15, height = 10)
VlnPlot(M, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0.2) + NoLegend()
dev.off()

### --- subset the data --- ###
M.subset <- subset(
                  x = M,
                  subset = nCount_ATAC < 7e4 &
                  nCount_ATAC > 5e3 &
                  nCount_RNA < 25000 &
                  nCount_RNA > 1000 &
                  percent.mt < 20
)


### --- We next perform pre-processing and dimensional reduction on both assays independently, using standard approaches for RNA and ATAC-seq data --- ###
# RNA analysis
DefaultAssay(M.subset) <- "RNA"
M.subset<- SCTransform(M.subset, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(M.subset) <- "ATAC"
M.subset <- RunTFIDF(M.subset)
M.subset <- FindTopFeatures(M.subset, min.cutoff = 'q0')
M.subset <- RunSVD(M.subset)
M.subset <- RunUMAP(M.subset, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")


### --- We calculate a WNN graph, representing a weighted combination of RNA and ATAC-seq modalities. We use this graph for UMAP visualization and clustering --- ###
M.subset <- FindMultiModalNeighbors(M.subset, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
M.subset <- RunUMAP(M.subset, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
M.subset <- FindClusters(M.subset, graph.name = "wsnn", algorithm = 3, verbose = FALSE, resolution = arguments$r)

### --- lets make the umap plot --- ###

p1 <- DimPlot(M.subset, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(M.subset, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(M.subset, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("WNN")
pdf("Seurat_v4_wnn/umap.pdf", width = 15.75, height = 7)
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()


### --- findMarkers for each cluster --- ###
DefaultAssay(M.subset) <- "SCT"
M.markers <- FindAllMarkers(M.subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(M.markers, file ="Seurat_v4_wnn/markers.txt", col.names = T, sep = "\t", row.names = F, quote = F)

saveRDS(M.subset, paste0(arguments$o, "_wnn.rds"))
### --done -- ###


### usage Rscript Seurat_wnn.R "data/sal//"  "WM40"










