### --- This scripts takes the h5 files from multiome h5 and makes individual text file as well as Merged text  File --- ###
#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
require(optparse)
#Parse arguments from command line
options <- list(
  make_option(c("-p", "--path"), action = "store", default = getwd(), type = "character", help=" path for MultiModal h5 (GEX+SCATAC files) .")
  
)
arguments <- parse_args(OptionParser(option_list = options))
setwd(arguments$p) 

cat(arguments$p)


### --- Rscript takes the h5 files from gex H5 and return the text file for each and also the merged file --- ###

message(" Loading the required packages .... ")
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(hdf5r))


myFiles <- list.files(path=arguments$p, pattern="*.h5")
message("Reading all the text file and making list....")
myList <- lapply(myFiles, function(x) Read10X_h5(x))
samples <- gsub(pattern = "*.h5", replacement = "", x = basename(myFiles))

message(" creating an list for gene expression counts ...")
counts.list <- list();
for (i in 1:length(myList)) {
  counts.list[[i]] = myList[[i]][["Gene Expression"]]
}

message("apending the samplename to the respective cellbarcodes ....")
names(counts.list) <- samples
for (i in 1:length(myList)) {
  colnames(counts.list[[i]]) = paste(colnames(counts.list[[i]]), names(counts.list[i]) , sep = "_" );
  print(paste("no of cell barcodes in ", names(counts.list[i]), "are ",  ncol(counts.list[[i]])))
}

### --- convert list to an dataframe --- ###
message("converting list to an dataframes ....")
dfs <- lapply(counts.list, data.frame, stringsAsFactors = FALSE)
dfs <- lapply(dfs, rownames_to_column, "uid")

for (i in 1:length(dfs)){
  message("writting the text files ....")
  write_tsv(dfs[[i]], file = paste0(names(dfs[i]), "_counts.txt"), quote_escape = "none", col_names = T)
  print(paste("text file for" , names(dfs[i]) , "completed" ))
}

message("individual file done ...")

###  ---as we alreeady have dataframe which has all different samples saved as list, lets use that --- ###
message("merging all the text file ...")
mymergeddata  <- Reduce(function(x,y) {merge(x,y)}, dfs)
write_tsv(mymergeddata , file = "MergedFiles.txt", quote_escape = "none", col_names = T)
message("Merging completed ...")


### --- usage Rscript 10x-MultiH5-To-Text.R --path "path to the h5 files"









