#### This scripts converts 10x mtx files into text file compatible with altanalyze
### Make sure the files in the directory are gziped 
#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
require(optparse)
options <- list(
make_option(c("-p", "--path"), action = "store", default = getwd(), type = "character", help=" path for  gene expression h5 file produced from normal cellranger   ."),
make_option(c("-o", "--outname"), action = "store", default = "output", type = "character", help="Output file(s) base name of counts.")
)


arguments <- parse_args(OptionParser(option_list = options))
setwd(arguments$p, arguments$o)

cat(arguments$p)

message(" Loading the required packages .... ")
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(data.table))

Load10x <- function(
  mtx.dir = NULL,
  gene.column = 2,
  unique.features = TRUE,
  strip.suffix = FALSE
) {
  full.data <- list()
  for (i in seq_along(along.with = mtx.dir)) {
    run <- mtx.dir[i]
    if (!dir.exists(paths = run)) {
      stop("Directory provided does not exist")
    }
    barcode.loc <- file.path(run, 'barcodes.tsv')
    gene.loc <- file.path(run, 'genes.tsv')
    features.loc <- file.path(run, 'features.tsv.gz')
    matrix.loc <- file.path(run, 'matrix.mtx')
    # Flag to indicate if this data is from CellRanger >= 3.0
    pre_ver_3 <- file.exists(gene.loc)
    if (!pre_ver_3) {
      addgz <- function(s) {
        return(paste0(s, ".gz"))
      }
      barcode.loc <- addgz(s = barcode.loc)
      matrix.loc <- addgz(s = matrix.loc)
    }
    if (!file.exists(barcode.loc)) {
      stop("Barcode file missing. Expecting ", basename(path = barcode.loc))
    }
    if (!pre_ver_3 && !file.exists(features.loc) ) {
      stop("Gene name or features file missing. Expecting ", basename(path = features.loc))
    }
    if (!file.exists(matrix.loc)) {
      stop("Expression matrix file missing. Expecting ", basename(path = matrix.loc))
    }
    data <- readMM(file = matrix.loc)
    cell.names <- readLines(barcode.loc)
    if (all(grepl(pattern = "\\-1$", x = cell.names)) & strip.suffix) {
      cell.names <- as.vector(x = as.character(x = sapply(
        X = cell.names,
        FUN = ExtractField,
        field = 1,
        delim = "-"
      )))
    }
    if (is.null(x = names(x = mtx.dir))) {
      if (i < 2) {
        colnames(x = data) <- cell.names
      } else {
        colnames(x = data) <- paste0(i, "_", cell.names)
      }
    } else {
      colnames(x = data) <- paste0(names(x = mtx.dir)[i], "_", cell.names)
    }
    feature.names <- read.delim(
      file = ifelse(test = pre_ver_3, yes = gene.loc, no = features.loc),
      header = FALSE,
      stringsAsFactors = FALSE
    )
    if (any(is.na(x = feature.names[, gene.column]))) {
      warning(
        'Some features names are NA. Replacing NA names with ID from the opposite column requested',
        call. = FALSE,
        immediate. = TRUE
      )
      na.features <- which(x = is.na(x = feature.names[, gene.column]))
      replacement.column <- ifelse(test = gene.column == 2, yes = 1, no = 2)
      feature.names[na.features, gene.column] <- feature.names[na.features, replacement.column]
    }
    if (unique.features) {
      fcols = ncol(x = feature.names)
      if (fcols < gene.column) {
        stop(paste0("gene.column was set to ", gene.column,
                    " but feature.tsv.gz (or genes.tsv) only has ", fcols, " columns.",
                    " Try setting the gene.column argument to a value <= to ", fcols, "."))
      }
      rownames(x = data) <- make.unique(names = feature.names[, gene.column])
    }
    # In cell ranger 3.0, a third column specifying the type of data was added
    # and we will return each type of data as a separate matrix
    if (ncol(x = feature.names) > 2) {
      data_types <- factor(x = feature.names$V3)
      lvls <- levels(x = data_types)
      if (length(x = lvls) > 1 && length(x = full.data) == 0) {
        message("10X data contains more than one type and is being returned as a list containing matrices of each type.")
      }
      expr_name <- "Gene Expression"
      if (expr_name %in% lvls) { # Return Gene Expression first
        lvls <- c(expr_name, lvls[-which(x = lvls == expr_name)])
      }
      data <- lapply(
        X = lvls,
        FUN = function(l) {
          return(data[data_types == l, , drop = FALSE])
        }
      )
      names(x = data) <- lvls
    } else{
      data <- list(data)
    }
    full.data[[length(x = full.data) + 1]] <- data
  }
  # Combine all the data from different directories into one big matrix, note this
  # assumes that all data directories essentially have the same features files
  list_of_data <- list()
  for (j in 1:length(x = full.data[[1]])) {
    list_of_data[[j]] <- do.call(cbind, lapply(X = full.data, FUN = `[[`, j))
    # Fix for Issue #913
    list_of_data[[j]] <- as(object = list_of_data[[j]], Class = "dgCMatrix")
  }
  names(x = list_of_data) <- names(x = full.data[[1]])
  # If multiple features, will return a list, otherwise
  # a matrix.
  if (length(x = list_of_data) == 1) {
    return(list_of_data[[1]])
  } else {
    return(list_of_data)
  }
}

data = Load10x(mtx.dir = arguments$p)   #### dir containing the mtx file
raw_counts = as.data.frame(as.matrix(data))
setDT(raw_counts, keep.rownames = "UID")
raw_counts = raw_counts[!duplicated(raw_counts$UID),]
write.table(raw_counts, paste0(arguments$o, "_raw_counts.txt"), col.names = T, sep = "\t", row.names = F, quote = F)


#### Usage Rscript 10x-Mtx-To-Txt.R --path /data/salomonis2/outs/bc_filtered  --outname "whatever"
