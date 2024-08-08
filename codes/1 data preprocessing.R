
library(dplyr)
library(limma)
library(GEOquery)
library(umap)

#-------------------------------------------------------------------------------
#- GSE136324: Generate from RAW data
#-------------------------------------------------------------------------------

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("affy", quietly = TRUE))
  BiocManager::install("affy")

library(affy)
url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE136324&format=file"
options(timeout = max(3600, getOption("timeout")))
download.file(url, destfile = "./gse136324.tar", method = "auto", mode = "wb")

untar("./gse136324.tar", files = NULL, list = F, exdir = "./gse136324/")
files = list.files("./gse136324", full.names = TRUE)
affy.data = ReadAffy(filenames = files)
gse136324 = mas5(affy.data[1:10])
save(list = "gse136324", file= "../data/GSE136324_RAW.RData")

#- merge the gse136324 SOFT format that provides clinical information ----------

s   <- "GSE136324" #- n = 867
gse <- getGEO(s, GSEMatrix = T, AnnotGPL = T)[[1]]
gse$data_processing[1]

#- summary the clinical data ----
f   <- gse@phenoData@data
f$data_processing[1]
f$data_processing.1[1]
f$data_processing.2[1]

id <- f$geo_accession
group  <- gsub("pointgroup: ", "", f[,"characteristics_ch1.8"])
gender <- gsub("gender: ", "", f[,"characteristics_ch1.58"])
race <- gsub("race: ", "", f[,"characteristics_ch1.59"])
age  <- as.numeric(gsub("ageatsampledate: ", "", f[,"characteristics_ch1.60"]))
pts  <- f$characteristics_ch1.1

time <- as.numeric(gsub("monthsos.ar: ", "", f[,"characteristics_ch1.145"]))
time <- as.numeric(time)
time <- ifelse(time == 0, 0.001, time)

status <- gsub("censos.ar: ", "", f[,"characteristics_ch1.143"])
status <- as.numeric(!as.logical(status))

clin <- data.frame(id = id, group = group, pts = pts,
                   gender = gender, race = race,
                   age = age, time = time,
                   status = status)

#- patients selection, retain the sample of the first test

s <- order(clin$group, decreasing = F)
clin <- clin[s, ]

dup <- !duplicated(clin$pts)
clin <- clin[dup, ]

time <- ifelse(clin$time == 0, 0.001, as.numeric(clin$time))

#- add micro_array with probe data
#- load raw data with probes ----

load("../data/GSE136324_RAW.RData")
X <- exprs(gse136324)
samples <- gsub("_.*", "", colnames(X))
colnames(X) <- samples
X <- X[,clin$id]

#- preprocessing X data

# log2 transformation
qx <- as.numeric(quantile(X, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||  (qx[6]-qx[1] > 50 && qx[2] > 0)

if (LogC) {
  X[which(X <= 0)] <- NaN
  X <- log2(X)}

# check 1, box-and-whisker plot

dev.new(width = 3 + ncol(X) / 6, height = 5)
par(mar =c (7, 4, 2, 1))
boxplot(X, boxwex = 0.7, notch = T, outline = FALSE, las = 2)
dev.off()

#- check 2,
par(mar=c(4, 4, 2, 1))
plotDensities(X, legend = "right")

#- check 3, mean-variance trend
X <- na.omit(X)
plotSA(lmFit(X))

# UMAP plot (multi-dimensional scaling)
X <- X[!duplicated(X), ]
ump <- umap(t(X), n_neighbors = 15, random_state = 123)
plot(ump$layout, xlab = "", ylab = "", pch = 20, cex = 1.5)

# Point labels without overlaps
library("maptools")
pointLabel(ump$layout, labels = rownames(ump$layout),
           method = "SANN", cex = 0.6)

save(list = c("clin", "X"), file = "../data/GSE136324.RData")

#-------------------------------------------------------------------------------
#- GSE2658
#-------------------------------------------------------------------------------

s = "GSE2658"
gse <- getGEO(s, GSEMatrix = T, getGPL = T, AnnotGPL = T)[[1]]
X   <- exprs(gse)
f   <- gse@phenoData@data

#- remove isoform
# g   <- gse@featureData@data

f$data_processing[1]
f$data_processing.1[1]
f$data_processing.2[1]

#- GCOS ------

library(stringr)
id <- f$geo_accession
status <- as.numeric(grepl("SUIND=1", f$characteristics_ch1))
time <- str_extract(f$characteristics_ch1.2, "[0-9.]+")
time <- as.numeric(time)
time <- ifelse(time == 0, 0.001, time)

grp  <- ifelse(grepl("TT2", f$title), "TT2", "TT3")

clin <- data.frame(id = id, grp = grp, status = status, time = time)
#- log2 transformation
qx <- as.numeric(quantile(X, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)

if (LogC){
  X[which(X <= 0)] <- NaN
  X <- log2(X)
}


# check 1, box-and-whisker plot

dev.new(width = 3 + ncol(X) / 6, height = 5)
par(mar =c (7, 4, 2, 1))
boxplot(X, boxwex = 0.7, notch = T, outline = FALSE, las = 2)
dev.off()

#- check 2,
par(mar=c(4, 4, 2, 1))
plotDensities(X)

#- check 3, mean-variance trend
X <- na.omit(X)
plotSA(lmFit(X))

# UMAP plot (multi-dimensional scaling)
X <- X[!duplicated(X), ]
ump <- umap(t(X), n_neighbors = 15, random_state = 123)
plot(ump$layout, xlab = "", ylab = "", pch = 20, cex = 1.5)

# Point labels without overlaps
library("maptools")
pointLabel(ump$layout, labels = rownames(ump$layout),
           method = "SANN", cex = 0.6)

save(list = c("X", "clin"), file = "../data/GSE2658.RData")

#-------------------------------------------------------------------------------

s = "GSE9782"
gse <- getGEO(s, GSEMatrix = T, getGPL = T, AnnotGPL = T)[[1]]
X   <- exprs(gse)
f   <- gse@phenoData@data

#- byproduct, save probe annotation profile

# anno <- gse@featureData@data
# save(list = "anno", file = "../data/prob_annotation.RData")

#- byproduct end

id <- f$geo_accession
treat  <- gsub("treatment = ", "", f$characteristics_ch1.1)
gender <- gsub("sex = ", "", f$characteristics_ch1.2)
race <- gsub("race = ", "", f$characteristics_ch1.3)
age  <- as.numeric(gsub("Age_at_Randomization = ", "", f$characteristics_ch1.4))

#-

find <- function(x){
  re <- apply(f, 1, function(y) y[grepl(x, y)])
  gsub(x, "", re)
}

time <- find("Days_Survived_From_Randomization = ")
#-adjust to month
time <- as.numeric(time) / 30
time <- ifelse(time == 0, 0.001, time)
status <- find("Did_Patient_Die")
status <- as.numeric(status == "(0=No,1=Yes) = 1")

clin <- data.frame(id = id, treat = treat, gender = gender,
                   race = race, age = age, time = time, status = status)

nrow(clin) == ncol(X)

# log2 transformation
qx <- as.numeric(quantile(X, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) {
  X[which(X <= 0)] <- NaN
  X <- log2(X)
}

# check 1, box-and-whisker plot

dev.new(width = 3 + ncol(X) / 6, height = 5)
par(mar =c (7, 4, 2, 1))
boxplot(X, boxwex = 0.7, notch = T, outline = FALSE, las = 2)
dev.off()

#- check 2,
par(mar=c(4, 4, 2, 1))
plotDensities(X, legend = "right")

#- check 3, mean-variance trend
X <- na.omit(X)
plotSA(lmFit(X))

# UMAP plot (multi-dimensional scaling)
X <- X[!duplicated(X), ]
ump <- umap(t(X), n_neighbors = 15, random_state = 123)
plot(ump$layout, xlab = "", ylab = "", pch = 20, cex = 1.5)

# Point labels without overlaps
library("maptools")
pointLabel(ump$layout, labels = rownames(ump$layout),
           method = "SANN", cex = 0.6)

save(list = c("clin", "X"), file = "../data/GSE9782.RData")

#-------------------------------------------------------------------------------
#- keep the datasets of shared genes
#-------------------------------------------------------------------------------

load("../data/GSE136324.RData")
g1 <- rownames(X)

load("../data/GSE9782.RData")
g2 <- rownames(X)

load("../data/GSE2658.RData")
g3 <- rownames(X)

#- intersect
g <- intersect(intersect(g1, g2), g3)

#-------------------------------------------------------------------------------

load("../data/GSE2658.RData")
g1 <- rownames(X)

load("../data/GSE9782.RData")
g2 <- rownames(X)

load("../data/GSE136324.RData")
g3 <- rownames(X)

g <- intersect(intersect(g1, g2), g3)

#- load annotation

library(hgu133a.db)
probs <- mapIds(hgu133a.db, keys = g, column = "SYMBOL",
                keytype = "PROBEID", multiVals = "first")

#- For the Probes to Multiple Symbols
length(probs) == length(g)

#- Average the expression for multiple

load("../data/GSE2658.RData")
X <- aggregate(X[g, ] ~ probs, FUN = "mean")
rownames(X) <- X[,1]
X <- X[, -1]
save(list = c("clin", "X"), file = "../data/GSE2658.RData")

load("../data/GSE9782.RData")
X <- aggregate(X[g, ] ~ probs, FUN = "mean")
rownames(X) <- X[,1]
X <- X[, -1]
save(list = c("clin", "X"), file = "../data/GSE9782.RData")

load("../data/GSE136324.RData")
X <- aggregate(X[g, ] ~ probs, FUN = "mean")
rownames(X) <- X[,1]
X <- X[, -1]
save(list = c("clin", "X"), file = "../data/GSE136324.RData")

#-------------------------------------------------------------------------------
#- Batch effects correction
#-------------------------------------------------------------------------------

load("../data/GSE136324.RData")
X136324 <- X

load("../data/GSE2658.RData")
X2658 <- X

load("../data/GSE9782.RData")
X9782 <- X

library(sva)
load("../data/batchCorrectedData.RData")

#-------------------------------------------------------------------------------
#- mean variance standardization -----------------------------------------------
#-------------------------------------------------------------------------------

# batch <- c(rep(1, ncol(X136324)), rep(2, ncol(X2658)), rep(3, ncol(X9782)))
# cb <- ComBat(dat = cbind(X136324, X2658, X9782),
#              batch = batch,
#              mod = NULL,
#              par.prior = FALSE,
#              prior.plots = FALSE,
#              mean.only = FALSE,
#              ref.batch = 1,
#              BPPARAM = bpparam("SerialParam"))

#-------------------------------------------------------------------------------
#- Data compare with other data
#-------------------------------------------------------------------------------

library(sva)

geoget <- function(series){

  nam <- stringr::str_extract(series, "[0-9]+")
  gse <- getGEO(series, GSEMatrix = T, getGPL = T, AnnotGPL = T)[[1]]

  .GlobalEnv[[paste0("X", nam)]] <- exprs(gse)
  .GlobalEnv[[paste0("clin", nam)]] <-gse@phenoData@data
  .GlobalEnv[[paste0("f", nam)]] <-gse@featureData@data
}

geoget("GSE2658") #- 54675, 559
geoget("GSE9782") # 22283, 264

# impute missing values

miss.impute <- function(NewX, training) {

  #- if individual patients in vector, transform to Matrix
  if (is.vector(NewX))
    NewX <- matrix(NewX, ncol = 1)

  #- For valid, select all features including in training features
  sel <- intersect(rownames(NewX), rownames(training))
  NewX <- NewX[sel, ]

  #- tranform valid to training data structure
  ind <- match(rownames(NewX), rownames(training))
  mat <- matrix(NA, nrow = nrow(training), ncol = ncol(NewX),
                dimnames = list(rownames(training),
                                colnames(NewX)))
  mat[ind, ] <- NewX

  for(i in seq_len(ncol(NewX))){

    cat("Processing:", colnames(NewX)[i], "\n")
    sink("nul")
    out <- impute::impute.knn(cbind(training, ref = mat[,i]), k = 10)
    sink()
    mat[,i] <- out$data[,"ref"]
  }
  mat
}

New9782 <- miss.impute(X9782, X2658)

meanvarstandardize <- function(testdata, traindata = ""){

  inters <- intersect(rownames(traindata), rownames(testdata))
  trainmean <- rowMeans(traindata[inters,])
  trainsd <- apply(traindata[inters,], 1, sd)

  testdata <- testdata[inters,] - trainmean
  testdata <- testdata / trainsd
  return(testdata)
}

#- batch correction ------------------------------------------------------------

load("../data/GSE136324.RData")
X136324 <- X

load("../data/GSE9782.RData")
X9782 <- X

load("../data/GSE2658.RData")
X2658 <- X

batch <- c(rep(1, ncol(X136324)), rep(2, ncol(X9782)), rep(3, ncol(X2658)))
cb <- ComBat(dat = cbind(X136324, X9782, X2658),
             batch = batch,
             mod = NULL,
             par.prior = T,
             prior.plots = FALSE,
             mean.only = F,
             ref.batch = 1,
             BPPARAM = bpparam("SerialParam"))

load("../data/GSE9782.RData") #- get clin data
X <- cb[,batch == 2]
save(list = c("X", "clin"), file = "../data/GSE9782B.RData")

load("../data/GSE2658.RData")
X <- cb[,batch == 3]
save(list = c("X", "clin"), file = "../data/GSE2658B.RData")

#-------------------------------------------------------------------------------



















