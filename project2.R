library(DropSeq.util)
library("affycoretools")
library("DESeq2")
library(SingleCellExperiment)
library(scran)
library(GENIE3)
library(edgeR)
library(scater)
library(M3Drop)
source("project2functions.R")

# start_time <- Sys.time()
# dge1.path <- "C:\\Users\\BackUp Account\\Documents\\H_1stRound_CrossTissue_FibroblastLike_5-3-17.raw.dge.txt.gz"
# dge1.path <- "C:\\Users\\BackUp Account\\Documents\\H_1stRound_CrossTissue_Microglia_Macrophage_5-3-17.raw.dge.txt.gz"
 dge1.path <- "C:\\Users\\BackUp Account\\Documents\\H_1stRound_CrossTissue_Mural_5-3-17.raw.dge.txt.gz"
# dge1.path <- "C:\\Users\\BackUp Account\\Documents\\H_1stRound_CrossTissue_Polydendrocytes_5-3-17.raw.dge.txt.gz"

# imports as compact matrix
d1 <- loadSparseDge(dge1.path) 
d1 <- as.matrix(d1)

options(stringsAsFactors = FALSE)
umi <- SingleCellExperiment(assays = list(counts = d1))
logcounts(umi) <- log2(counts(umi) + 1)
plotPCA(umi)
plotTSNE(umi)
rm(umi)
# TASK 1 - Cell Quality Control

# removes zero count genes
print(paste("number of genes before zero count removal: ", nrow(d1)))
formattedMatrix <- d1[base::rowSums(d1!=0)>0,]
print(paste("number of genes after zero count removal: ", nrow(formattedMatrix)))
# for memory:
rm(d1)
# gets column totals
sampleTotals <- base::colSums(formattedMatrix)
hist(sampleTotals, breaks = 100)

# for clocking
# end_time <- Sys.time()
# total_time <- end_time - start_time
# after abline is set:
# restart_time <- Sys.time()
# parameter specific to cell type
abline(v = 800, col = "red")
rm(sampleTotals)
# prints number of cells before read count QC:
print(paste("original cell count: ",base::ncol(formattedMatrix)))
# removes likely failed reads:
QCmatrix <- formattedMatrix[,base::colSums(formattedMatrix)>1000]
rm(formattedMatrix)

# prints number of cells before read count QC:
print(paste("cell # following read count QC: ",base::ncol(QCmatrix)))

# makes a vector where each element is the number of unique genes per column
originalTotals <- colSums(QCmatrix != 0)
hist(originalTotals, breaks = 100)

# for clocking
# total_time <- (Sys.time() - restart_time)
# after abline is set:
# restart_time <- Sys.time()
# parameter specific to cell type
abline(v = 600, col = "red")
rm(originalTotals)
# removes cells with relatively poor distribution:
QCmatrix <- QCmatrix[,which(600 <= colSums(QCmatrix != 0))]
# prints number of cells before read count QC:
print(paste("cell # following gene detection QC: ",base::ncol(QCmatrix)))

ERCCindeces <- which(TRUE == grepl("Ercc",row.names(QCmatrix), fixed = TRUE))
MTindeces <- which(TRUE == grepl("mt-",row.names(QCmatrix), fixed = TRUE))

QCmatrix <- filterMatrix(QCmatrix, ERCCindeces)
print(paste("cell # after Ercc gene spike QC: ",base::ncol(QCmatrix)))
QCmatrix <- filterMatrix(QCmatrix, MTindeces)
print(paste("cell # after Mt- gene spike QC: ",base::ncol(QCmatrix)))

umi2 <- SingleCellExperiment(assays = list(counts = QCmatrix))
logcounts(umi2) <- log2(counts(umi2) + 1)
plotPCA(umi2)
plotTSNE(umi2)
rm(umi2)
# TASK 2 - Normalization by Library Size

umi3cpm <- SingleCellExperiment(assays = list(counts = QCmatrix))
normcounts(umi3cpm) <- calculateCPM(umi3cpm, use_size_factors = FALSE)
logcounts(umi3cpm) <- log2(normcounts(umi3cpm) + 1)
# logcounts(umi3cpm) <- log2(calculateCPM(umi3cpm, use_size_factors = FALSE)+1)
plotPCA(umi3cpm)
plotTSNE(umi3cpm)

#method 1 for pool dropout, manual normalization
'
library(scater)
umi3pool <- SingleCellExperiment(assays = list(counts = QCmatrix))
factors <- getSizeFactors(QCmatrix)
poolNormalized <- scater::normalizeCounts(QCmatrix, size_factors = factors, return_log = TRUE, log_exprs_offset = 1)
logcounts(umi3pool) <- poolNormalized
summary(sizeFactors(umi3pool))
'

#method 2 for pool dropout
umi3pool <- SingleCellExperiment(assays = list(counts = QCmatrix))
qclust <- quickCluster(umi3pool, min.size=30)
# sizes factor must be removed for fibroblast dataset
umi3pool <- computeSumFactors(umi3pool, sizes = 15, clusters=qclust)
summary(sizeFactors(umi3pool))
# scran uses scater!
# https://bioconductor.org/packages/3.7/bioc/vignettes/scran/inst/doc/scran.html
# saves to normcounts w/ second argument, instead of logcounts
umi3pool <- scater::normalize(umi3pool, return_log = FALSE)
# plots require logcounts assay
logcounts(umi3pool) <- log2(normcounts(umi3pool) + 1)
plotPCA(umi3pool)
plotTSNE(umi3pool)

rm(QCmatrix)
# total_time <- (Sys.time() - restart_time)
# TASK 3 - Gene Selection
restart_time <- Sys.time()
#dgeVersion <- convertTo(umi3pool, type="edgeR")
#matrixVersion <- as.matrix(dgeVersion)

BrenneckGenes <-
  BrenneckeGetVariableGenes(
    normcounts(umi3pool),
    fdr = .01,
    minBiolDisp = 0.5,
    suppress.plot = FALSE
  )

# start_time <- Sys.time()
# myCellLabels <- colnames(QCmatrix)
# my_list <- M3DropCleanData(normcounts(umi3pool), labels = myCellLabels, min_detected_genes = 100, is.counts = TRUE)
# my_matrix <- my_list$data

M3Drop_genes <- M3DropDifferentialExpression(
  normcounts(umi3pool),
  mt_method = "fdr",
  mt_threshold = .01,
  suppress.plot = FALSE
)
# just grab the names
M3Drop_genes <- M3Drop_genes[,1]
# total_time <- (Sys.time() - restart_time)
# TASK 4 - Network Inference
# restart_time <- Sys.time()
# because it's returning dublicates for some reason:
BrenneckGenes <- unique(BrenneckGenes)
# get indeces for selected genes:
# whichOnes <- which(TRUE == BrenneckGenes%in%row.names(normcounts(umi3pool)))
whichOnes <- which(TRUE == row.names(normcounts(umi3pool))%in%BrenneckGenes)
whichOthers <- which(TRUE == row.names(normcounts(umi3pool))%in%M3Drop_genes)
selected <- union(whichOnes,whichOthers)
# makes new matrix of only selected genes
genieMatrix <- normcounts(umi3pool)[selected,]
# reduce number of genes for inference by removing genes populated mainly by zeros (variable adjustable):
cleanGenie <- filterZeros(genieMatrix,0.75)
rm(genieMatrix)
# for clocking
total_time <- (Sys.time() - restart_time)
# to resume after modifying threshold and checking genieMatrix rows
# restart_time <- Sys.time()
# to reduce number of genes by row count. Commented out in favor of selection argument adjustment
# genieMatrix <- preGenieMatrix[rowSums(preGenieMatrix!=0)>2500,]

# packages for multithreading
# library(doParallel)
# library(doRNG)
# parallel processing causes overheat! nCores should not exceed 1

wam <- GENIE3(cleanGenie, regulators = NULL, targets = NULL, treeMethod = "RF",
              K = "all", nTrees = 7, nCores = 1, verbose = TRUE)

shuffledMatrix <- cleanGenie
for (i in 1:nrow(cleanGenie)) shuffledMatrix[i,] <- cleanGenie[i,order(runif(length(cleanGenie[i,])))] 
# reshuffle:
# for (i in 1:nrow(shuffledMatrix)) shuffledMatrix[i,] <- shuffledMatrix[i,order(runif(length(shuffledMatrix[i,])))] 

swam <- GENIE3(shuffledMatrix, regulators = NULL, targets = NULL, treeMethod = "RF",
               K = "all", nTrees = 7, nCores = 1, verbose = TRUE)
rm(shuffledMatrix)
Gresults <- getLinkList(wam)
rm(wam)
shuffledGresults <- getLinkList(swam)
rm(swam)
Goutput <- Gresults[1:50,]

library(igraph)
net <- graph_from_data_frame(Goutput, directed = T)
# combines duplicate nodes
net <- simplify(net, remove.multiple = T, remove.loops = T) 
plot(net,  edge.arrow.size=.02, vertex.label.dist = 1.5, 
     edge.color = "darkgrey", vertex.size = 5, vertex.frame.color = "NA", 
     vertex.label.color= "black", margin = 0, vertex.label = Gresults$V1)

# false positive calculation
index <- round(nrow(Gresults)*0.05)
threshold <- shuffledGresults[index,3]
sigIndex <- min(which(threshold > Gresults[,3]))*2

library(ggplot2)

for (i in 1:5){
  xIndex <- match(Goutput[i,1], row.names(logcounts(umi3pool)))
  vectorX <- logcounts(umi3pool)[xIndex,]
  yIndex <- match(Goutput[i,2], row.names(logcounts(umi3pool)))
  vectorY <- logcounts(umi3pool)[yIndex,]
  plot(x = vectorX, y = vectorY,  xlab = row.names(logcounts(umi3pool))[xIndex],
       ylab = row.names(logcounts(umi3pool))[yIndex], main = paste("Interaction # ",i))
} # end for

# graph for worst result:
xIndex <- match(Gresults[nrow(Gresults),1], row.names(logcounts(umi3pool)))
vectorX <- logcounts(umi3pool)[xIndex,]
yIndex <- match(Gresults[nrow(Gresults),2], row.names(logcounts(umi3pool)))
vectorY <- logcounts(umi3pool)[yIndex,]
plot(x = vectorX, y = vectorY,  xlab = row.names(logcounts(umi3pool))[xIndex],
     ylab = row.names(logcounts(umi3pool))[yIndex], main = "Bottom-ranked Interaction")

# Inference method 2, resets all variables ------------------------
library(netbenchmark)
# transposes
cleanGenie <- t(cleanGenie)

wam <- clr.wrap(cleanGenie)
Gresults <- getLinkList(wam)
# to remove doubles for printing:
evens <-seq(2,length(Gresults[,1]),2)
top50 <-seq(1,50,1)
Gresults <- Gresults[evens,]
Goutput <- (Gresults[1:50,])
row.names(Goutput) <- top50

shuffledMatrix <- cleanGenie
# adjusted for transposition:
for (i in 1:ncol(cleanGenie)) shuffledMatrix[,i] <- cleanGenie[order(runif(length(cleanGenie[,i]))),i] 
# reshuffle:
# for (i in 1:ncol(shuffledMatrix)) shuffledMatrix[,i] <- shuffledMatrix[order(runif(length(cleanGenie[,i]))),i] 
swam <- clr.wrap(shuffledMatrix)
shuffledGresults <- getLinkList(swam)
# to remove doubles for printing:
evens <-seq(2,length(shuffledGresults[,1]),2)
top50 <-seq(1,50,1)
shuffledGresults <- shuffledGresults[evens,]
shuffledGoutput <- (shuffledGresults[1:50,])
row.names(shuffledGoutput) <- top50

# false positive calculation
index <- round(nrow(Gresults)*0.05)
threshold <- shuffledGresults[index,3]
sigIndex <- min(which(threshold > Gresults[,3]))*2

net <- graph_from_data_frame(Goutput, directed = T)
# combines duplicate nodes
net <- simplify(net, remove.multiple = T, remove.loops = T) 
plot(net,  edge.arrow.size=.02, vertex.label.dist = 1.5, 
     edge.color = "darkgrey", vertex.size = 5, vertex.frame.color = "NA", 
     vertex.label.color= "black", margin = -0.1, vertex.label = Gresults$V1)

library(ggplot2)
for (i in 1:5){
  xIndex <- match(Goutput[i,1], row.names(logcounts(umi3pool)))
  vectorX <- logcounts(umi3pool)[xIndex,]
  yIndex <- match(Goutput[i,2], row.names(logcounts(umi3pool)))
  vectorY <- logcounts(umi3pool)[yIndex,]
  plot(x = vectorX, y = vectorY,  xlab = row.names(logcounts(umi3pool))[xIndex],
       ylab = row.names(logcounts(umi3pool))[yIndex], main = paste("Interaction # ",i))
} # end for

# graph for worst result:
xIndex <- match(Gresults[nrow(Gresults),1], row.names(logcounts(umi3pool)))
vectorX <- logcounts(umi3pool)[xIndex,]
yIndex <- match(Gresults[nrow(Gresults),2], row.names(logcounts(umi3pool)))
vectorY <- logcounts(umi3pool)[yIndex,]
plot(x = vectorX, y = vectorY,  xlab = row.names(logcounts(umi3pool))[xIndex],
     ylab = row.names(logcounts(umi3pool))[yIndex], main = "Bottom-ranked Interaction")

#total_time <- (Sys.time() - restart_time)

