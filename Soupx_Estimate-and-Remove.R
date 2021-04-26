### --- calcualte the contamination fraction and remove the ambient RNA from the 10x mtx --- ###
suppressPackageStartupMessages(library(optparse))
require(optparse)
options <- list(
  make_option(c("-p", "--path"), action = "store", default = getwd(), type = "character", help=" path for outs folder produced from cellranger."),
  make_option(c("-f", "--file"), type="character", default=NULL, help="MarkerGenes.txt  file produced by ICGS run ...", metavar="filename")
)

arguments <- parse_args(OptionParser(option_list = options))
setwd(arguments$p)
cat(arguments$p,arguments$f)

message("loading all the required packages ....")
suppressPackageStartupMessages(library(SoupX))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(DropletUtils))


message("Reading the markers file ....")
markers <-read.delim(arguments$f, h=T)
markers <- markers[order(-markers$Pearson.rho),] ## sorting by decreasing order
markers <- markers%>%filter(Pearson.rho >= 0.70)
##markers2 = markers%>%group_by(Cell.State)%>% top_n(n=2, wt = Pearson.rho)
nonExpressedGeneList <- list(ICGS = markers$UID); ICGS <- markers$UID


message("Reading the mtx files from the outs folder ...")
sc = load10X(dataDir = arguments$p, keepDroplets = TRUE)
head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 20)
useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = list(IG = ICGS))

#####Calculating the contamination fraction#######
message("estimating contamination from the soup channel ...")
fraction <-capture.output(c_fraction = calculateContaminationFraction(sc, list(IG = ICGS), useToEst = useToEst), type = "message")
fraction <- gsub(pattern = c("Estimated global contamination fraction of ", "%"),  replacement = "", x = fraction[1])
print(paste("total percentage of conatmination estimetedd ", "is", fraction))
fraction$percent <- gsub(pattern = "%", replacement = "", x = fraction)

### ---  Removing Contamination Now --- ###
message("removing the contamination fraction ...")
sc = setContaminationFraction(sc, as.numeric(fraction$percent)/100)
head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 20)
out = adjustCounts(sc,roundToInt=TRUE)
cntSoggy = rowSums(sc$toc > 0)
cntStrained = rowSums(out > 0)
mostZeroed = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
mostZeroed
tail(sort(rowSums(sc$toc > out)/rowSums(sc$toc > 0)), n = 20)
message("writting new mtx files to soupX-contamination-fraction-0.0453 folder ...")
DropletUtils:::write10xCounts(paste("./soupX-contamination-fraction_",fraction$percent, sep = ""), out)
quit()



