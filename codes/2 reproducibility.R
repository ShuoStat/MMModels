
pkgs <- c("survival", "dplyr", "glmnet", "survminer", "msigdb",
          "doParallel", "doSNOW", "reshape2", "ggplot2", "tidyr",
          "tidyverse", "ggpubr", "risksetROC", "impute", "sva", "GSVA")

for (p in pkgs) {
  tryCatch({
    test <- require(p, character.only = TRUE)
  }, error = function(e) {
    test <<- FALSE
  })
  
  if (!test) {
    print(paste("Package", p, "not found. Installing Package!"))
    install.packages(p)
    if (requireNamespace("BiocManager", quietly = TRUE)) {
      BiocManager::install(p)
    } else {
      install.packages("BiocManager")
      BiocManager::install(p)
    }
    require(p, character.only = TRUE)
  }
}

source("./fun.R")

#-------------------------------------------------------------------------------
#- load training data
#-------------------------------------------------------------------------------

load("../data/GSE136324.RData")
X <- as.matrix(X)
n = ncol(X)
p = nrow(X)
y = survival::Surv(clin$time, clin$status)

#- Parameters ------------------------------------------------------------------

# {
#   ncv, number of duplicated runs
#   nfolds, folds to generate foldid for glmnet
#   setN, number of pathways used in grplasso
#   allGenes, all genes in X
#   cutoff, cutoff to filter some redundant features
#   minSize, lower limits
#   ncores, cores used for parallel running
# }

nsim = 100
nfolds = 10
setN = 100
allGenes = rownames(X)
cutoff = 1
minSize = 10
ncores = 25
#- get random splits

getsplits(1584812, n = n, nsim = nsim, nfolds = nfolds, percentTraining = 0.7)

#-------------------------------------------------------------------------------
#- get pathways
#-------------------------------------------------------------------------------

#- GO pathways:"GO:BP";GO:CC; GO:MF
#- Canonical pathways: CP:REACTOME; CP:WIKIPATHWAYS; CP:PID; CP:BIOCARTA
#- Immune
#- C6: OncoSig
#- CM: Cancer modules:CM

filter.fun <- function(X, surv, cutoff = 0.05){
  
  if (cutoff == 1) {
    return(rownames(X))
  }
  
  if (cutoff != 1) {
    pvalues <- apply(X, 1, function(x) 
      summary(survival::coxph(surv ~ x))$coefficients[5])
    # adj.p <- p.adjust(pvalues, method = "BH")
    rownames(X)[pvalues <= cutoff]
  }
}

# gsc <- getMsigdb('hs', 'SYM')
# getMsigdbVersions()
# listCollections(gsc)
# listSubCollections(gsc)

getpathways <- function(){
  
  hset <- function(collection = c(), subcollection  = c()){
    h <- subsetCollection(gsc, 
                          collection = collection,
                          subcollection = subcollection)
    hset <- list()
    for(i in seq_along(h)){
      tmp <- h@.Data[[i]]
      hset[[tmp@setName]] <- tmp@geneIds
    }
    return(hset)
  }
  
  gsc = getMsigdb('hs', version = "7.5")
  
  #- C1 collection: positional gene sets. not included
  #- C2: Canonical pathways
  #- C3 collection: regulatory target gene sets, not included
  h.hall <- hset(collection = "h")
  
  #- C2
  h.biocarta <- hset(subcollection = "CP:BIOCARTA")
  # h.reactome <- hset(subcollection = "CP:REACTOME")
  h.pid      <- hset(subcollection = "CP:PID")
  h.cgp      <- hset(subcollection = "CGP")
  h.wiki     <- hset(subcollection = "CP:WIKIPATHWAYS")
  
  #- C3
  h.tft.legency  <- hset(subcollection = "TFT:TFT_Legacy")
  h.tft.gtrd     <- hset(subcollection = "TFT:GTRD")
  
  #- C4, computational gene sets
  # h.cgn <- hset(subcollection = "CGN")
  h.cm  <- hset(subcollection = "CM")
  
  #- C5, ontology gene sets
  h.gobp <- hset(subcollection = "GO:BP")
  h.gocc <- hset(subcollection = "GO:CC")
  h.gomf <- hset(subcollection = "GO:MF")
  # h.hpo  <- hset(subcollection = "HPO")
  
  #- C6, oncogene signature gene
  h.oncogene <- hset(collection = "c6")
  
  #- C7 
  h.immune <- hset(subcollection = "IMMUNESIGDB")
  h.vax    <- hset(subcollection = "VAX")
  
  hs <- list(hall = h.hall,
             biocarta = h.biocarta,
             # reactome = h.reactome,
             pid      = h.pid,
             cgp      = h.cgp,
             wiki     = h.wiki,
             tft.legency = h.tft.legency,
             tft.gtrd    = h.tft.gtrd,
             # cgn = h.cgn,
             cm  = h.cm,
             gobp = h.gobp,
             gocc = h.gocc,
             gomf = h.gomf,
             # hpo  = h.hpo, 
             oncogene = h.oncogene, 
             immune = h.immune,
             vax = h.vax)
  
  hs <- lapply(hs, function(sets){
    genes <- rownames(X)
    getpaths(sets, genes, N = length(sets), min.genes = 10)
  })
  return(hs)
}

hs <- getpathways()
# save("hs", file = "../output/hs.RData")
load("../output/hs.RData")
lapply(hs, length)

#- Table 1 
#- Pathway Information

A <- t(cbind.data.frame(lapply(hs, length)))
write.csv(A, "../results/pathwaysinf.csv")

#- pathways for group lasso ----------------------------------------------------

reducePaths <- function(paths, allGenes, minSize, setN){
  
  #- allGenes, all genes available in X
  #- minSize, remove pathways without much  
  #- setN, number of pathways to retain
  
  ind <- unlist(lapply(paths, function(x) 
    length(intersect(x, allGenes)) >= minSize))
  
  #- select the most distinct pathways
  paths <- paths[ind]
  allGenes <- rownames(X)
  
  #- output pathways
  getpaths(paths, allGenes, N = setN, min.genes = minSize, rm.duplicated = T)
  
}

#- pathway information for group lasso
hs.grp <- lapply(hs, reducePaths, allGenes = allGenes, minSize = minSize,
                 setN = setN)

#- save(list = "hs.grp", file = "../output/hs.grp.RData")
load("../output/hs.grp.RData")
lapply(hs.grp, function(x) length(unlist(x)))

#-------------------------------------------------------------------------------
# Get pathway-based models
#-------------------------------------------------------------------------------

fitMod.scores <- function(X, y, paths, cutoff, foldid,
                          method = c("gsva", "ssgsea", "zscore", "plage"),
                          minSize = 10) {
  
  submod <- function(X, y, paths, method, foldid, cutoff, minSize){
    
    #- trans to pathway scores
    Z <- GSVA::gsva(X, paths, method = method, min.sz = minSize,
                    ssgsea.norm = FALSE)
    
    #- filter pathways based on Uni-cox
    s <- unlist(filter.fun(Z, surv = y, cutoff = cutoff))
    
    #- model fitting
    cv.las <- glmnet::cv.glmnet(t(Z[s,]), y, family = "cox", foldid = foldid)
    las <- glmnet::glmnet(t(Z[s,]), y, family = "cox", lambda = cv.las$lambda.min)
    las <- las$beta[,1]
    #- return coefficients
    las[las != 0]
    
  }
  
  #- 
  mapply(submod, method = method, 
         MoreArgs = list(X = X, y = y, paths = paths, minSize = minSize, 
                         foldid = foldid, cutoff = cutoff),
         SIMPLIFY = F)
  
}

ncores <- pmin(detectCores() - 1, ncores)
cl <- makeCluster(ncores)
registerDoSNOW(cl)

pb <- txtProgressBar(min = 1, max = nsim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

foreach(i = 1:nsim, .packages = c("glmnet"), 
        .options.snow = opts) %dopar% {
          
          sam <- split.sample[[i]]
          f <- foldids[[i]]#- foldid for i.nsim splits and j fold
          
          #- 
          modsPath <- list()
          for(k in names(hs)){
            
            mod <- fitMod.scores(X[,sam], y[sam], 
                                 paths = hs[[k]],
                                 cutoff = cutoff, 
                                 foldid = f,
                                 method = c("gsva", "ssgsea", "zscore", "plage"),
                                 minSize = minSize)
            
            mod[["grp"]] <- grplas.surv(t(X[,sam]), y[sam], 
                                        pathways = hs.grp[[k]], 
                                        foldid = f)
            modsPath[[k]] <- mod
          }
          
          modsGenes <- fitfun(t(X[,sam]), y[sam], 
                              foldid = f, 
                              method = "las")$las
          
          save(list = c("modsPath", "modsGenes"), 
               file = paste0("../output/nsim_", i, ".RData"))
        }

stopCluster(cl)

#- Prediction performance in internal validation -------------------------------
#- get validation

# function to derive IBS (integrated Brier Score)
get.ibs <- function(X.t, X.v, y.t, y.v, beta) {
  
  # X.t, y.t training data
  # X.v, y.v validation data
  
  require(pec)
  
  if (length(beta) == 0){
    lp.t <- rep(1, ncol(X.t))
    lp.v <- rep(1, ncol(X.v))
  } else {
    lp.t <- t(X.t[names(beta),,drop = F]) %*% beta
    lp.v <- t(X.v[names(beta),, drop = F]) %*% beta
  }
  
  train <- data.frame(lp = lp.t, time = y.t[,1], status = y.t[,2])
  valid <- data.frame(lp = lp.v, time = y.v[,1], status = y.v[,2])
  mod <- survival::coxph(Surv(time, status) ~ lp, data = train, init = 1, iter.max = 0, x = T)
  f   <- as.formula("Surv(time, status) ~ 1")
  p   <- pec::pec(mod, f, data = valid, reference = F, verbose = F)
  as.numeric(pec::crps(p, times = p$minmaxtime))
}
  
# modPred, get model prediction used in parallel runing
modPred <- function(modsPath, modsGenes, X.t, X.v, y.t, y.v){
  
  # modsPath, including models for all pathways
  funPaths <- function(X.t,X.v, y.t, y.v, mods, pathBases){
    
    #- mods, mods for different methods
    #- pathBases, one of hs
    pred.c <- pred.ibs <- matrix(NA, 
                                 nrow = 1, 
                                 ncol = length(mods),
                                 dimnames = list(c(),
                                                 names(mods)))
    methods <- names(mods)
    
    if ("grp" %in% methods) {
      beta <- mods$grp$beta
      if (length(beta) > 0) {
        lp <- t(X.v[names(beta),,drop = F]) %*% beta
      } else {
        lp <- rep(1, ncol(X.v))
      }
      
      pred.c[, "grp"]  <- intsurv::cIndex(y.v[,1], y.v[,2], lp)["index"]
      pred.ibs[,"grp"] <- get.ibs(X.t, X.v, y.t, y.v, beta)
    }
    
    for(i in setdiff(methods, "grp")){
      
      beta <- mods[[i]]
      if (length(beta) > 0) {
        Z.v <- GSVA::gsva(X.v, pathBases[names(beta)], method = i, ssgsea.norm = FALSE)
        Z.t <- GSVA::gsva(X.t, pathBases[names(beta)], method = i, ssgsea.norm = FALSE)
        lp  <- t(Z.v) %*% beta
      } else {
        lp  <- rep(1, ncol(X.v))
        Z.t <- t(rep(1, ncol(X.t)))
        Z.v <- t(rep(1, ncol(X.v)))
      }
      pred.c[,i]   <- intsurv::cIndex(y.v[,1], y.v[,2], lp)["index"]
      pred.ibs[,i] <- get.ibs(Z.t, Z.v, y.t, y.v, beta)
    }
    
    pred.c <- as.data.frame(pred.c)
    pred.ibs <- as.data.frame(pred.ibs)
    
    pred.c[,"measure"] <- "cindex"
    pred.ibs[,"measure"] <- "ibs"
    
    out <- rbind(pred.c, pred.ibs) %>% 
      tidyr::pivot_longer(-measure) %>%
      dplyr::rename(method = name)
    #- stop here
    return(out)
  }
  
  #- prediction of genes
  funGenes <- function(X.t,X.v, y.t, y.v, modsGenes){
    
    beta = modsGenes
    if (length(beta) > 0) {
      lp <- t(X.v[names(beta),,drop = F]) %*% beta
    } else {
      lp <- rep(1, ncol(X.v))
    }
    
    pred.c   <- intsurv::cIndex(y.v[,1], y.v[,2], lp)["index"]
    pred.ibs <- get.ibs(X.t, X.v, y.t, y.v, beta)
    
    return(data.frame(measure = c("cindex", "ibs"),
                      method  = "genes",
                      value   = c(pred.c, pred.ibs)))
  }
  
  #- combine into a data.frame
  output <- mapply(funPaths, pathBases = hs, mods = modsPath, 
                   MoreArgs = list(X.t = X.t, X.v = X.v, y.t = y.t, y.v = y.v),
                   SIMPLIFY = F)
  
  output <- output %>% dplyr::bind_rows(.id = "paths")
  
  #- add genes prediction
  outputGenes <- funGenes(X.t, X.v, y.t, y.v, modsGenes)
  output <- rbind(output, cbind(paths = "genes", outputGenes))
  
  return(output)
}

#-
cl <- makeCluster(ncores)
registerDoSNOW(cl)

pb <- txtProgressBar(min = 1, max = nsim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

pred <- foreach(i = 1:nsim, 
                .combine = "rbind",
                .packages = c("glmnet", "survival", "pec", "intsurv", "dplyr"), 
                .options.snow = opts) %dopar% {
                  
                  sam <- split.sample[[i]]
                  load(paste0("../output/nsim_", i, ".RData"))
                  
                  X.v <- X[ ,!sam]
                  y.v <- y[!sam]
                  
                  X.t <- X[ ,sam]
                  y.t <- y[sam]
                  
                  #- validate procedure
                  out <- modPred(modsPath, modsGenes, X.t, X.v, y.t, y.v)
                  out$index <- i
                  return(out)
                }

stopCluster(cl)
# save(list = "pred", file = "../output/pred.nsim.RData")

#-------------------------------------------------------------------------------
# Visualization Internal CV, Figure 2
# In the revised version, we used the mean and 1.96 *se, instead of boxplits
#-------------------------------------------------------------------------------

load("../output/pred.nsim.RData")
library(ggplot2)

pred$paths <- stringr::str_to_title(pred$paths)
pred$paths <- recode(pred$paths, "Hall" = "Hallmark")

ord <- c("Hallmark", "Biocarta", "Wiki", "Pid", "Tft.gtrd", "Tft.legency", "Cgp", 
         "Gobp", "Gocc", "Gomf", "Oncogene", "Cm", "Immune", "Vax", "Genes")

# plage is removed because plage is not suitable for transformation 
ibs <- dplyr::filter(pred, measure == "ibs", method != "plage")
ibs$paths <- stringr::str_to_title(ibs$paths)

ibs.means <- ibs %>% 
  group_by(paths, method) %>%
  summarise(mean = format(round(mean(value), 3), nsmall = 3),
            sd   = format(round(sd(value), 3), nsmall = 3)) %>%
  mutate(ibs = paste0(mean, "(", sd, ")")) %>%
  dplyr::select(!c("mean", "sd")) %>%
  tidyr::pivot_wider(names_from = c("method"),
                     values_from = "ibs")

medianGene <- median(unlist(ibs[ibs$method == "genes", "value"]))

p1 <- ggplot(ibs, aes(x = factor(paths, 
                                 levels = ord, 
                                 order = T), 
                      y = value, fill = method)) +
  geom_boxplot(color = "black", position = position_dodge2(preserve = "single")) +
  xlab(NULL) + 
  ylab("IBS") + 
  theme_bw() + 
  ggtitle("B, IBS") + 
  geom_hline(yintercept = medianGene) +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))

# C-index
cindex <- filter(pred, measure == "cindex", method != "plage")

cindex.means <- cindex %>%
  group_by(paths, method) %>%
  summarise(mean = format(round(mean(value), 3), nsmall = 3),
            sd   = format(round(sd(value), 3), nsmall = 3)) %>%
  mutate(cindex = paste0(mean, "(", sd, ")")) %>%
  dplyr::select(!c("mean", "sd")) %>%
  tidyr::pivot_wider(names_from = c("method"),
                     values_from = "cindex")

medianGene <- median(unlist(cindex[cindex$method == "genes", "value"]))

p2 <- ggplot(cindex, aes(x = factor(paths, 
                                    levels = ord, 
                                    order = T), 
                         y = value, fill = method)) +
  geom_boxplot(color = "black", position = position_dodge2(preserve = "single")) +
  xlab(NULL) + 
  ylab("C-index") + 
  ggtitle("A, C-index") + 
  theme_bw() + 
  geom_hline(yintercept = medianGene) +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) 

p <- ggpubr::ggarrange(p2, p1, nrow = 2, common.legend = TRUE, legend = "bottom")
ggsave("../results/Figure 2.pdf", p, width = 8, height = 6)
ggsave("../results/Figure 2.tiff", p, width = 8, height = 6, units = "in", 
       dpi = 300, compression = "lzw")

# In the revised version, we use mean and se -----------------------------------
# C-index, mean and se

CMeanSe <- cindex %>%
  group_by(paths, method) %>%
  summarise(mean = mean(value),
            se   = sd(value) / sqrt(100)) %>%
  mutate(paths = factor(paths, levels = ord))
           
pd <- position_dodge(width = 0.5)
p1 <- ggplot(CMeanSe, aes(x = paths, y = mean, color = method)) + 
  geom_line(position = pd) +
  geom_pointrange(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se), 
                  position = pd,
                  size = 0.35) +
  theme_bw() + 
  geom_hline(yintercept = CMeanSe$mean[CMeanSe$paths == "Genes"]) +
  ylab("C-index") +
  xlab(NULL) + 
  ggtitle("A, C-index") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank())

# IBS, mean and se

IBSMeanSe <- ibs %>%
  group_by(paths, method) %>%
  summarise(mean = mean(value),
            se   = sd(value) / sqrt(100)) %>%
  mutate(paths = factor(paths, levels = ord))

pd <- position_dodge(width = 0.5)
p2 <- ggplot(IBSMeanSe, aes(x = paths, y = mean, color = method)) + 
  geom_line(position = pd) +
  geom_pointrange(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se), 
                  position = pd,
                  size = 0.35) +
  geom_hline(yintercept = IBSMeanSe$mean[IBSMeanSe$paths == "Genes"]) +
  theme_bw() + 
  ylab("IBS") +
  xlab(NULL) + 
  ggtitle("B, IBS") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank())


p <- ggpubr::ggarrange(p1, p2, nrow = 2, common.legend = TRUE, legend = "bottom")
ggsave("../results/Figure 2R.pdf", p, width = 8, height = 7)
ggsave("../results/Figure 2R.tiff", p, width = 8, height = 7, units = "in", 
       dpi = 300, compression = "lzw")

# in the revision, we conduct statistical test to compare method 

adj.t.test <- function(x, y, n1, n2) {
  
  # Let x and y be the accuracy of algorithms A and B respectively
  # n1 instances are used for training, and the remaining n2 instances for testing
  
  # Check if the lengths of x and y are the same
  if(length(x) != length(y))
    stop("x and y should be the same length")
  
  # Calculate the differences between x and y
  z = as.numeric(x) - as.numeric(y)
  
  # Determine the number of differences
  n = length(z)
  
  # Calculate the variance of the differences
  sigma = var(z)
  
  # Calculate the t-statistic
  t = mean(z) / sqrt((1 / n + n2 / n1) * sigma)
  
  # Calculate the p-value based on the t-statistic and degrees of freedom
  p.value <- pt(q = t, df = n - 1)
  
  # Return the mean of x, mean of y, the t-statistic, and p-value as a named list
  return(list(mean.x = mean(x),
              mean.y = mean(y),
              statistic = t,
              p.value = p.value))
}

#-------------------------------------------------------------------------------

load( "../output/pred.nsim.RData")

# cindex
cMat <- pred %>%
  filter(measure == "cindex" & method  != "plage") %>%
  mutate(path_method = str_to_title(paste0(paths, "(", method, ")"))) %>%
  mutate(ifelse(path_method == "Genes(Genes)", "Gene", path_method)) %>%
  pivot_wider(id_cols = index, id_expand = TRUE, 
              values_from = value, names_from = path_method) 

iMat <- pred %>%
  filter(measure == "ibs" & method  != "plage") %>%
  mutate(path_method = str_to_title(paste0(paths, "(", method, ")"))) %>%
  mutate(ifelse(path_method == "Genes(Genes)", "Gene", path_method)) %>%
  pivot_wider(id_cols = index, id_expand = TRUE, 
              values_from = value, names_from = path_method)
n1 = 306
n2 = 130
cMat <- cMat[, -1]
allMethods <- names(cMat)
nMethod <- ncol(cMat)

pmat <- matrix(NA, nrow = nMethod, ncol = nMethod,
               dimnames = list(allMethods, allMethods))

for(i in seq_len(nMethod)) {
  if(i < (nMethod - 1)) {
    for(j in (i + 1):nMethod){
      pmat[i, j] <- adj.t.test(pull(cMat[,i]), pull(cMat[,j]), 
                               n1 = n1, n2 = n2)$p.value
      pmat[j, i] <- adj.t.test(pull(iMat[,i]), pull(iMat[,j]),
                               n1 = n1, n2 = n2)$p.value
    }
  }
}

library(RColorBrewer)
pmat[pmat > 0.10] <- NA
diag(pmat) <- 0
pdf("../results/stat.pdf", width = 10, height = 10)
jpeg("../results/stat.jpeg", width = 10, height = 10, units = "in", res = 300)

m <- pheatmap::pheatmap(pmat, 
                        color = colorRampPalette(rev(brewer.pal(n = 7, name = "Greens")))(100),
                        cluster_rows = FALSE, 
                        cluster_cols = FALSE,
                        legend_breaks = c(0, 0.01, 0.05, 0.10),
                        legend_labels = c("0", "0.01", "0.05", "0.10"),
                        na_col = "white",
                        fontsize = 8)
print(m)
dev.off()

#-------------------------------------------------------------------------------
# External validation
#-------------------------------------------------------------------------------

#- training models -------------------------------------------------------------

set.seed(6852348)
foldid <- sample(rep(1:10, length = n), n)

modsPath <- list()
for(k in names(hs)){
  
  mod <- fitMod.scores(X, y, 
                       paths = hs[[k]],
                       cutoff = cutoff, 
                       foldid = foldid,
                       method = c("gsva", "ssgsea", "zscore"),
                       minSize = minSize)
  
  mod[["grp"]] <- grplas.surv(t(X), y, 
                              pathways = hs.grp[[k]], 
                              foldid = foldid)
  modsPath[[k]] <- mod
}

modsGenes <- fitfun(t(X), y, 
                    foldid = foldid, 
                    method = "las")$las

# save(list = c("modsPath", "modsGenes"), file = "../output/training.RData")
#-------------------------------------------------------------------------------

load("../output/training.RData")
load("../data/GSE136324.RData")
X.t <- as.matrix(X)
y.t <- y

#- test data, consider batch effects
# GSE9782
load("../data/GSE9782B.RData")
X.v = as.matrix(X)
y.v = survival::Surv(clin$time, clin$status)
re9782 <- modPred(modsPath, modsGenes, X.t, X.v, y.t, y.v)

# GSE2658
load("../data/GSE2658B.RData")
X.v = as.matrix(X)
y.v = Surv(clin$time, clin$status)
re2658 <- modPred(modsPath, modsGenes, X.t, X.v, y.t, y.v)

# save(list = c("re9782", "re2658"), file = "../output/outvalidB.RData")

predict.mods <- function(mods, newX, paths = "hall", method = "ssgsea"){
  
  # mods, either modsPath, or modsGenes
  # method, gsva, ssgsea, zscre grp
  # paths, biocarta, cgn, cpg, cm, genes, gobp, gocc, gomf, hall, hpo, immune, oncogene
  # pid, reactome, tft.gtrd, tft.legency, vax, wiki
  
  if (method == "genes"){
    
    beta <- mods
    if (length(beta) > 0) {
      lp <- t(newX[names(beta),,drop = F]) %*% beta
    } else {
      lp <- rep(1, ncol(newX))
    }
  }
  
  if (method == "grp") {
    
    mod  <- mods[[paths]]
    beta <- mod$grp$beta
    if (length(beta) > 0) {
      lp <- t(newX[names(beta),,drop = F]) %*% beta
    } else {
      lp <- rep(1, ncol(newX))
    }
  }
  
  if(method %in% c("ssgsea", "zscroe", "gsva")){
    
    mod  <- mods[[paths]]
    beta <- mod[[method]]
    
    if (length(beta) > 0) {
      newZ <- GSVA::gsva(newX, hs[[paths]][names(beta)], method = method, 
                        ssgsea.norm = FALSE)
      lp  <- t(newZ) %*% beta
    } else {
      lp  <- rep(1, ncol(newZ))
    }
  }
  
  return(lp)
}
 
load("../output/training.RData")

#- in original data
# predFun <- function(modsPath, modsGenes, newX, newy, path_method, X.t, y.t){
# 
#   newX = as.matrix(newX)
#   X.t  = as.matrix(X.t)
#   # to GSVA method
# 
#   tranX.new <- list()
#   tranX.t   <- list()
# 
#   for(i in seq_along(path_method)) {
# 
#     path <- path_method[[i]][1]
#     method <- path_method[[i]][2]
# 
#     if (method %in% c("ssgsea", "zscore", "gsva")){
# 
#       mod  <- modsPath[[path]]
#       beta <- mod[[method]]
# 
#       tranX.new[[paste0("X_", path, "_", method)]] <-
#         GSVA::gsva(newX, hs[[path]][names(beta)], method = method,
#                    ssgsea.norm = FALSE)
# 
#       tranX.t[[paste0("X.t_", path, "_", method)]] <-
#         GSVA::gsva(X.t, hs[[path]][names(beta)], method = method,
#                    ssgsea.norm = FALSE)
# 
#     }
#   }
# }


#-------------------------------------------------------------------------------
# Get CI using bootstrap
#-------------------------------------------------------------------------------
# 
# set.seed(123)
# nCores = 10
# nRun = 1000
# 
# CIList <- list()
# methods <- c("EMC92Classifier", "UAMS70Classifier", "IFM15Classifier", "MILLENNIUM100Classifier")
# 
# for(method in methods) {
#   require(doParallel)
#   require(doSNOW)
#   nCores  <- pmin(detectCores() - 1, nCores)
#   cluster <- makeCluster(nCores)
#   registerDoParallel(cluster)
# 
# 
#   # adjDat, function in fun.R, used for batch correction by methods
#   dat <- adjDat(method)
#   pb <- txtProgressBar(min = 1, max = nRun, style = 3)
#   progress <- function(n) setTxtProgressBar(pb, n)
#   opts <- list(progress = progress)
# 
#   re <- foreach(i = seq_len(nRun),
#                 .export = c("dat"),
#                 .combine = "rbind",
#                 .options.snow = opts) %dopar% {
#                   source("previous.models.R")
#                   out <- parFun(dat, method = method)
#                   return(out)
#                 }
# 
#   CIList[[method]] <- re
#   stopCluster(cluster)
# }
# 

# using bootstrap to obtain CI, also includes external validation 
# mods, c("modsPath", "modsGenes")
  
BootMyModsFun <- function(modsPath, modsGenes, 
                          newX, newy, 
                          path_method, 
                          samIndList,
                          nCores = 1, 
                          X.t, y.t) {
  
  # path_method, list. The first element is pathway and the second is method
  # example
  # path_method <- list(c("vax", "grp"),
  #                     c("vax", "ssgsea"),
  #                     c("", "genes"))
  
  newX = as.matrix(newX)
  X.t  = as.matrix(X.t)
  # to GSVA method
  
  tranX.new <- list()
  tranX.t   <- list()
  
  for(i in seq_along(path_method)) {
    
    path <- path_method[[i]][1]
    method <- path_method[[i]][2]
    
    if (method %in% c("ssgsea", "zscore", "gsva")){
      
      mod  <- modsPath[[path]]
      beta <- mod[[method]]
      
      tranX.new[[paste0("X_", path, "_", method)]] <- 
        GSVA::gsva(newX, hs[[path]][names(beta)], method = method, 
                   ssgsea.norm = FALSE)
      
      tranX.t[[paste0("X.t_", path, "_", method)]] <- 
        GSVA::gsva(X.t, hs[[path]][names(beta)], method = method, 
                   ssgsea.norm = FALSE)
      
    }
  }
  
  # getMethod
  doParFun <- function(samInd){
    
    out <- list()
    for(j in seq_along(path_method)) {
      
      path   <- path_method[[j]][1]
      method <- path_method[[j]][2]
      
      # gene model
      if (method == "genes"){
        beta <- modsGenes
        if (length(beta) > 0) {
          newZ <- newX[names(beta), samInd, drop = F]
          lp <- t(newZ) %*% beta
        } else {
          lp <- rep(1, ncol(newZ))
        }
      }
      
      # grp
      if (method == "grp") {
        mod  <- modsPath[[path]]
        beta <- mod$grp$beta
        if (length(beta) > 0) {
          newZ <- newX[names(beta), samInd, drop = F]
          lp <- t(newZ) %*% beta
        } else {
          lp <- rep(1, ncol(newZ))
        }
      }
      
      # other
      if(method %in% c("ssgsea", "zscore", "gsva")){
        
        mod  <- modsPath[[path]]
        beta <- mod[[method]]
        
        if (length(beta) > 0) {
          newZ <- tranX.new[[paste0("X_", path, "_", method)]]
          lp   <- t(newZ[ ,samInd, drop = FALSE]) %*% beta
        } else {
          lp  <- rep(1, ncol(newZ))
        }
      }
      
      # get C-index
      cindex <- intsurv::cIndex(newy[samInd, "time"], 
                                newy[samInd, "status"], lp)["index"]
      
      # for gene, 
      if (method %in% c("ssgsea", "zscore", "gsva")) {
        Z.t <- tranX.t[[paste0("X.t_", path, "_", method)]]
      } else {
        Z.t <- X.t
      }
      
      ibs <- get.ibs(Z.t, newZ, y.t, newy[samInd,], beta)
      
      out[[paste0(path, "_", method)]] <- 
        setNames(c(cindex, ibs), c("c-index", "ibs"))
    }
    
    return(out)
  }
  
  # get Cindex and IBS for whole data
  allPred <- doParFun(seq_len(ncol(newX)))

  require(doParallel)
  require(doSNOW)
  nCores  <- pmin(detectCores() - 1, nCores)
  cluster <- makeCluster(nCores)
  registerDoParallel(cluster)  
  
  nRun <- length(samIDList)
  pb <- txtProgressBar(min = 1, max = nRun, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  re <- foreach(i = seq_len(length(samIndList)),
                .packages = "survival",
                .export = c("tranX.t", "tranX.new", "get.ibs"),
                .options.snow = opts) %dopar% {
                  
                  samInd <- samIndList[[i]]
                  return(doParFun(samInd))
                }
  
  stopCluster(cluster)
  
  # extract group
  extract_group <- function(input_list, group_name) {
    do.call(rbind, lapply(input_list, `[[`, group_name))
  }
  
  outList <- list()
  for(i in seq_along(path_method)) {
    path <- path_method[[i]][1]
    method <- path_method[[i]][2]
    outList[[paste0(path, "_", method)]] <- extract_group(re, paste0(path, "_", method))
  }
  
  CI <- lapply(outList, function(x) 
    apply(x, 2, function(xx) quantile(xx, c(0.025, 0.975))))
  
  out <- list()
  for(k in names(allPred)) {
    
    mean <- round(allPred[[k]], 3)
    ci   <- round(CI[[k]], 3)
    ci   <- apply(ci, 2, function(x) paste0("(", x[1], ",", x[2], ")"))
    out[[k]] <- paste0(mean, ci)
  }
  
  return(list(forTab = as.data.frame(out),
              forPlot = list(CI = CI, maen = allPred),
              forTest = outList))
}



# in GSE136324

path_method <- list(c("vax", "grp"),
                    c("vax", "ssgsea"),
                    c("", "genes"))
nRun <- 1000

load("../output/training.RData")
load("../output/hs.RData")
load("../data/GSE136324.RData")
X136324 <- as.matrix(X)
y136324 <- survival::Surv(clin$time, clin$status)

n <- ncol(X136324)
set.seed(123)
samIndList <- sapply(seq_len(nRun), function(x) sample(seq_len(n), n, replace = TRUE),
                     simplify = FALSE)

Bootpred136324 <- BootMyModsFun(modsPath, modsGenes, 
                                newX = X136324, newy = y136324,
                                path_method = path_method, 
                                samIndList = samIndList,
                                nCores = 4, 
                                X.t = X136324, y.t = y136324)

# GSE9782
load("../data/GSE9782B.RData")
X9782 = as.matrix(X)
y9782 = survival::Surv(clin$time, clin$status)
n <- ncol(X9782)
set.seed(123)
samIndList <- sapply(seq_len(nRun), function(x) sample(seq_len(n), n, replace = TRUE),
                     simplify = FALSE)

Bootpred9782 <- BootMyModsFun(modsPath, modsGenes, 
                            newX = X9782, newy = y9782,
                            path_method = path_method, 
                            samIndList = samIndList,
                            nCores = 10, 
                            X.t = X136324, y.t = y136324)


# GSE2658
load("../data/GSE2658B.RData")
X2658 = as.matrix(X)
y2658 = survival::Surv(clin$time, clin$status)

set.seed(123)
n <- ncol(X2658)
samIndList <- sapply(seq_len(nRun), function(x) sample(seq_len(n), n, replace = TRUE),
                     simplify = FALSE)

Bootpred2658 <- BootMyModsFun(modsPath, modsGenes, 
                            newX = X2658, newy = y2658,
                            path_method = path_method, 
                            nCores = 10, 
                            samIndList = samIndList,
                            X.t = X136324, y.t = y136324)

# generate table

tab3 <- cbind(t(Bootpred136324$forTab), t(Bootpred9782$forTab), t(Bootpred2658$forTab))
colnames(tab3) <- paste0(
  rep(c("GSE136324", "GSE9782", "GSE2658"), each = 2), "-",
  rep(c("C-index", "IBS"), time = 3))

write.csv(tab3, file = "../results/Tab3.csv")

# generate plot
# 
# meanCiDat <- c()
# 
# for(i in c("GSE136324", "GSE9782", "GSE2658")) {
#   
#   obj <- get(paste0("Bootpred", str_extract(i, "[0-9]+")))
#   
#   ci <- obj$forPlot$CI
#   mean <- obj$forPlot$maen
#   
#   for(j in names(ci)) {
#     
#     tmpci <- t(ci[[j]]) %>%
#       as.data.frame() %>%
#       tibble::rownames_to_column("metric")
#     
#     tmpmean <- mean[[j]] %>%
#       data.frame(mean = .) %>%
#       tibble::rownames_to_column("metric") %>%
#       left_join(y = tmpci, by = "metric") %>%
#       mutate(model = j, data = i)
#     
#     meanCiDat <- rbind(meanCiDat, tmpmean)
#   }
# }
# 
# meanCiDat <- meanCiDat %>%
#   mutate(model = case_when(model == "vax_grp" ~ "Vax(grp)",
#                            model == "vax_ssgsea" ~ "Vax(ssgsea)",
#                            model == "_genes" ~ "Genes"))
#          
# 
# # plot function                  
# plotInternal <- function(metric, data) {
#   
#   dat <- meanCiDat[meanCiDat$metric == metric & meanCiDat$data == data,]
#   pd <- position_dodge(width = 0.5)
#   p1 <- ggplot(dat, aes(x = model, y = mean, fill = "black")) + 
#     geom_line(position = pd) +
#     geom_pointrange(aes(ymin =`2.5%`, ymax = `97.5%`), 
#                     position = pd,
#                     size = 0.35) +
#     theme_bw() + 
#     ylab("C-index") +
#     xlab(NULL) + 
#     ggtitle("A, C-index") + 
#     theme(legend.title = element_blank(),
#           legend.position = "none")
#   p1
# }
# 
# p1 <- plotInternal(metric = "c-index", data = "GSE9782")

#-------------------------------------------------------------------------------
#- statistical test based on bootstrep resampling
#-------------------------------------------------------------------------------

test.fun <- function(obj) {
  
  obj <- obj$forTest
  dat <- as.data.frame(obj)
  
  n <- ncol(dat)
  pmat <- matrix(NA, n, n, dimnames = list(colnames(dat), colnames(dat)))
  for(i in seq_len(n-1)){
    for(j in (i+1):n){
      pmat[i,j] <- wilcox.test(dat[,i], dat[,j], paired = TRUE)$p.value
    }
  }
  
  return(pmat)
}

test.fun(Bootpred136324)
test.fun(Bootpred9782)
test.fun(Bootpred2658)

#-Statistical test using likelihood ratio --------------------------------------

# FitCoxMod <- function(path, method, newX, newy) {
#   
#   # to matrix     
#   newX = as.matrix(newX)
#   # transDat 
#   if (method %in% c("ssgsea", "zscore", "gsva")){
#     
#     mod  <- modsPath[[path]]
#     beta <- mod[[method]]
#     
#     newX <- GSVA::gsva(newX, hs[[path]][names(beta)], method = method, 
#                        ssgsea.norm = FALSE)
#     
#     if (length(beta) > 0) {
#       lp   <- t(newX[ , , drop = FALSE]) %*% beta
#     } else {
#       lp  <- rep(1, ncol(newX))
#     }
#   }
#   
#   # gene model
#   if (method == "genes"){
#     beta <- modsGenes
#     if (length(beta) > 0) {
#       newX <- newX[names(beta), , drop = F]
#       lp <- t(newX) %*% beta
#     } else {
#       lp <- rep(1, ncol(newX))
#     }
#   }
#   
#   # grp
#   if (method == "grp") {
#     mod  <- modsPath[[path]]
#     beta <- mod$grp$beta
#     if (length(beta) > 0) {
#       newX <- newX[names(beta), , drop = F]
#       lp <- t(newX) %*% beta
#     } else {
#       lp <- rep(1, ncol(newX))
#     }
#   }
#   
#   # fit model
#   require(survival)
#   mod <- survival::coxph(newy ~ lp, init = 1, iter.max = 0, x = T)
#   return(mod)
# }
# 
# load("../output/training.RData")
# load("../data/GSE9782B.RData")
# X9782 = as.matrix(X)
# y9782 = survival::Surv(clin$time, clin$status)
# 
# mod1 <- FitCoxMod(path = "vax", method = "grp", newX = X9782, newy = y9782)
# mod2 <- FitCoxMod(path = "vax", method = "ssgsea", newX = X9782, newy = y9782)
# mod3 <- FitCoxMod(path = "", method = "genes", newX = X9782, newy = y9782)
# 
# nonnestcox::plrtest(mod1, mod2, nested = FALSE, adjusted = "none")
# nonnestcox::plrtest(mod1, mod3, nested = FALSE, adjusted = "none")
# nonnestcox::plrtest(mod2, mod3, nested = FALSE, adjusted = "none")
# 
# p1 <- anova(mod1, mod2)$`Pr(>|Chi|)`[2]
# p2 <- anova(mod1, mod3)$`Pr(>|Chi|)`[2]
# p3 <- anova(mod2, mod3)$`Pr(>|Chi|)`[2]
# 
# p.adjust(c(p1, p2, p3), "bonferroni")
# #
# 
# load("../data/GSE2658B.RData")
# X2658 = as.matrix(X)
# y2658 = survival::Surv(clin$time, clin$status)
# 
# mod1 <- FitCoxMod(path = "vax", method = "grp", newX = X2658, newy = y2658)
# mod2 <- FitCoxMod(path = "vax", method = "ssgsea", newX = X2658, newy = y2658)
# mod3 <- FitCoxMod(path = "", method = "genes", newX = X2658, newy = y2658)
# 
# nonnestcox::plrtest(mod1, mod2, nested = FALSE, adjusted = "none")
# nonnestcox::plrtest(mod1, mod3, nested = FALSE, adjusted = "none")
# nonnestcox::plrtest(mod2, mod3, nested = FALSE, adjusted = "none")
# 
# 
# p1 <- anova(mod1, mod2)$`Pr(>|Chi|)`[2]
# p2 <- anova(mod1, mod3)$`Pr(>|Chi|)`[2]
# p3 <- anova(mod2, mod3)$`Pr(>|Chi|)`[2]
# p.adjust(c(p1, p2, p3), "bonferroni")

# summarize results
#-------------------------------------------------------------------------------
#- compare with previous models
#-------------------------------------------------------------------------------
# data prepare
# parallel function
# get the confidence interval using bootstrap resampling

parFun <- function(dat, method, samInds) {
  
  gFun <- function(dat, datName, method, samInds) {
    # datName, e.g., "GSE9782"
    datID <- gsub("GSE", "", datName)
    n <- ncol(dat[[paste0("X", datID)]])
    
    samInd <- samInds[[datName]]
    parList <- list(trainlogdata = dat$Ref,
                    testlogdata  = dat[[paste0("X", datID)]][,samInd])
    
    lp <- do.call(method, parList)[ ,"score"]
    # get cindex
    cindex <-  intsurv::cIndex(dat[[paste0("y", datID)]][samInd, "time"],
                               dat[[paste0("y", datID)]][samInd, "status"],
                               lp)["index"]
    return(cindex)
  }
  
  # reList <- list()
  # for(method in methods) {
  #   reList[[method]] <- gFun(dat, datName, method, samInd)
  # }
  
  out <- c("GSE136324" = gFun(dat, datName = "GSE136324", method = method, samInds = samInds),
           "GSE9782" =   gFun(dat, datName = "GSE9782", method = method, samInds = samInds),
           "GSE2658" =   gFun(dat, datName = "GSE2658", method = method, samInds = samInds))
  return(out)
}


nCores = 6
nRun = 1000
methods <- c("EMC92Classifier", "UAMS70Classifier", "IFM15Classifier", "MILLENNIUM100Classifier")

CIList <- list()
getSamIndList <- function(datName, nRun) {
  
  # ensure use same samInd with bootstrap for comparing proposed model and gene model
  set.seed(123)
  load(paste0("../data/", datName, ".RData"))
  n <- ncol(X)
  sapply(seq_len(nRun), function(x) sample(seq_len(n), n, replace = TRUE),
         simplify = FALSE)
}

# get sam list
samIndList <- sapply(c("GSE136324", "GSE9782", "GSE2658"),
                     getSamIndList, 
                     nRun = 1000,
                     simplify = FALSE)

for(method in methods) {
  
  require(doParallel)
  require(doSNOW)
  nCores  <- pmin(detectCores() - 1, nCores)
  cluster <- makeCluster(nCores)
  registerDoParallel(cluster)  
  
  # adjDat, function in fun.R, used for batch correction by methods
  dat <- adjDat(method)
  pb <- txtProgressBar(min = 1, max = nRun, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  re <- foreach(i = seq_len(nRun),
                .export = c("dat", "samIndList"),
                .combine = "rbind",
                .options.snow = opts) %dopar% {
                  
                  source("previous.models.R")
                  
                  samInds <- list()
                  for(datNam in c("GSE136324", "GSE9782", "GSE2658")) {
                    samInds[[datNam]] <- samIndList[[datNam]][[i]]
                  }
                  
                  out <- parFun(dat, method = method, samInds) 
                  return(out)
                }
  
  CIList[[method]] <- re
  stopCluster(cluster)
}

re <- lapply(CIList, function(x) {
  
  tidyr::pivot_longer(as.data.frame(x), 
               everything(), 
               names_to = "data") %>% 
    dplyr::group_by(data) %>%
    dplyr::summarise(lower = round(quantile(value, 0.025), 3),
              upper = round(quantile(value, 0.975), 3)) %>%
    mutate(ci = paste0("(", lower, ",", upper, ")"),
           data = gsub(".index", "", data)) %>%
    dplyr::select(data, ci)
})


# summarize ci
ci <- re %>%
  purrr::imap_dfr( ~ .x %>% 
                     mutate(Classifier = .y)) %>%
  tidyr::pivot_wider(names_from = Classifier, values_from = ci)

# summarize mean
meanC <- c()
source("previous.models.R")

# sample ID for all samples
samInds <- list() 

for(i in c("GSE136324", "GSE9782", "GSE2658")){
  # ensure use same samInd with bootstrap for comparing proposed model and gene model
  set.seed(123)
  load(paste0("../data/", i, ".RData"))
  n <- ncol(X)
  samInds[[i]] <- seq_len(n)
}

for(method in methods){
  dat <- adjDat(method)
  meanC <- cbind(meanC, parFun(dat, method = method, samInds = samInds))
}

colnames(meanC) <- methods
# rownames remove index
rownames(meanC) <- gsub(".index", "", rownames(meanC))
# 3 digits
meanC <- round(meanC, 3)

# merge and summarize results
otherC <- meanC %>% 
  as.data.frame() %>%
  tibble::rownames_to_column("data") %>%
  left_join(ci, by = "data") %>%
  mutate(EMC92Classifier = paste0(EMC92Classifier.x, EMC92Classifier.y),
         UAMS70Classifier = paste0(UAMS70Classifier.x, UAMS70Classifier.y),
         IFM15Classifier = paste0(IFM15Classifier.x, IFM15Classifier.y),
         MILLENNIUM100Classifier = paste0(MILLENNIUM100Classifier.x, 
                                          MILLENNIUM100Classifier.y)) %>%
  dplyr::select(data, IFM15Classifier, UAMS70Classifier, EMC92Classifier, MILLENNIUM100Classifier) %>%
  t()

write.csv(otherC, "../results/comparemodels.csv")

# statistical test

pvalues <- c()

for(i in c("GSE136324", "GSE9782", "GSE2658")) {

  myMod <- get(paste0("Bootpred", gsub("GSE", "", i)))$forTest$vax_grp[,1]
  ind <- grep(i, colnames(CIList[[1]]))
  otherMods <- lapply(CIList, function(x) x[,ind])
  
  pvalue <- c()
  for(j in seq_along(otherMods)) {
    pvalue <- c(pvalue,  wilcox.test(myMod, otherMods[[j]], paired = TRUE)$p.value)
  }
  names(pvalue) <- names(otherMods)
  pvalues[[i]] <- pvalue
}

#-------------------------------------------------------------------------------
#- Missing data imputation
#-------------------------------------------------------------------------------

#- add missing values ----------------------------------------------------------
#- C-index

#- load traning data (for missing imputation), models, and hs 

mistest <- function(missing, foldid, mods, cutoff = 1, minSize = 10) {
  
  #- nam, data name
  #- missing, missing variables
  #- mods, two vectors, with the fist being the pathways, the second, method
  
  load("../data/GSE136324.RData")
  X136324 <- as.matrix(X)
  y136324 <- survival::Surv(clin$time, clin$status)
  
  getd <- function(nam, missing){
    load(paste0("../data/", nam, ".RData"))
    New.X <- X[!rownames(X) %in% missing,]
    New.X <- miss.impute(New.X, X136324)
    new.y <- survival::Surv(clin$time, clin$status)
    list(X = New.X, y = new.y)
  }
  
  GSE9782 <- getd("GSE9782B", missing)
  GSE2658 <- getd("GSE2658B", missing)
  
  cindex9782 <- list()
  cindex2658 <- list()
  
  ibs9782 <- list()
  ibs2658 <- list()
  
  for(k in seq_along(mods)){
    
    ind <- mods[[k]]   
    if (ind[2] == "genes") {
      models = modsGenes
    } else {
      models = modsPath
    }
    
    lp  <- predict.mods(models, GSE9782$X, paths = ind[1], method = ind[2])
    #- c-index
    cindex9782[[k]] <- intsurv::cIndex(GSE9782$y[,"time"], GSE9782$y[,"status"], lp)["index"]
    
    #- ibs
    lp.t <- predict.mods(models, X136324, paths = ind[1], method = ind[2])
    ibs9782[[k]] <- get.ibs(t(data.frame(x = lp.t)), t(data.frame(x = lp)), y136324, 
                          GSE9782$y, setNames(1, "x"))
    
    #- c-index
    lp <- predict.mods(models, GSE2658$X, paths = ind[1], method = ind[2])
    cindex2658[[k]] <- intsurv::cIndex(GSE2658$y[,"time"], GSE2658$y[,"status"], lp)["index"]
    
    #- ibs
    ibs2658[[k]] <- get.ibs(t(data.frame(x = lp.t)), t(data.frame(x = lp)), 
                            y136324, GSE2658$y, setNames(1, "x"))
    
  }
  
  return(list(cindex9782 = unlist(cindex9782), 
              cindex2658 = unlist(cindex2658),
              ibs9782    = unlist(ibs9782),
              ibs2658    = unlist(ibs2658)))
}

#-------------------------------------------------------------------------------
#- pts, percentage missing

pts <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.3)
rep = 10
ncores = 25
set.seed(16584358)
load("../output/training.RData")

#- import training model

miss <- list()
for (i in seq_len(rep)) {
  re <- sapply(pts, function(x) sample(allGenes, size = ceiling(p * x)))
  names(re) <- pts
  miss[[i]] <- re
}


ncores <- pmin(detectCores() - 1, ncores)
cl <- makeCluster(ncores)
registerDoParallel(cl)

subruns <- (rep * length(pts))

pb <- txtProgressBar(min = 1, max = subruns, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

foldid <- sample(rep(1:10, length = n), n)

misspred <- foreach(i = 1:subruns, .packages = c("impute", "survival")) %dopar% {

  rep.ind <- ceiling(i / 6)
  pts.ind <- i - 6 * (rep.ind - 1)
  
  out <- list()
  out[["misspred"]] <- mistest(missing = miss[[rep.ind]][[pts.ind]],
                               foldid = foldid, 
                               mods = list(c("vax", "grp"),
                                           c("vax", "ssgsea"),
                                           c("genes", "genes")))
  
  out[["rep"]] <- rep.ind
  out[["pts"]] <- pts[pts.ind]
  return(out)
}

stopCluster(cl)
save(list = "misspred", file = "../output/missing.RData")

#- output the Cindex chagnes of the two important variables

extractmissing <- function(obj, method = c("cindex", "ibs"), data = c("9782", "2658")){
  
  ind <- paste0(method, data)
  out <- lapply(obj, function(x) c(x$misspred[[ind]], x$rep, x$pts))
  out <- Reduce(rbind, out) %>% as.data.frame()
  colnames(out) <- c("Vax.grp", "Vax.ssGSEA", "Genes", "rep", "missing")

  out %>% tidyr::pivot_longer(cols = 1:3)
} 

load("../output/missing.RData")
#- cindex in 9782
cindex <- extractmissing(misspred, method = "cindex", data = "9782")
cmean <- aggregate(value ~ name + missing, FUN = "mean", data = cindex)
csd   <- aggregate(value ~ name + missing, FUN = "sd", data = cindex)
cindex <- merge(cmean, csd, by = c("name", "missing"), suffixes = c(".mean", ".sd"))

pd <- position_dodge(0.15)
cindex <- cindex %>% mutate(missing = as.factor(missing),
                            name = recode(name, "Vax.grp" = "Vax(grp)",
                                          "Vax.ssGSEA" = "Vax(ssGSEA)"))

p1 <- ggplot2::ggplot(data = cindex, aes(x = missing, y = value.mean, color = name, group = name)) +
  geom_errorbar(aes(ymin = value.mean - value.sd, ymax = value.mean + value.sd), 
                width = .15, position = pd) +
  geom_line(position = pd) +
  geom_point(position = pd, size = 1.5) + 
  labs(x = "% of missing", y = "C-index") + 
  ggtitle("A, C-index(GSE9782)") + 
  theme_bw() + 
  theme(legend.title = element_blank())

#- cindex in 2658
cindex <- extractmissing(misspred, method = "cindex", data = "2658")
cmean <- aggregate(value ~ name + missing, FUN = "mean", data = cindex)
csd   <- aggregate(value ~ name + missing, FUN = "sd", data = cindex)
cindex <- merge(cmean, csd, by = c("name", "missing"), suffixes = c(".mean", ".sd"))

pd <- position_dodge(0.15)
cindex <- cindex %>% mutate(missing = as.factor(missing),
                            name = recode(name, "Vax.grp" = "Vax(grp)",
                                          "Vax.ssGSEA" = "Vax(ssGSEA)"))

p2 <- ggplot2::ggplot(data = cindex, aes(x = missing, y = value.mean, color = name, group = name)) +
  geom_errorbar(aes(ymin = value.mean - value.sd, ymax = value.mean + value.sd), 
                width = .15, position = pd) +
  geom_line(position = pd) +
  geom_point(position = pd, size = 1.5) + 
  labs(x = "% of missing", y = "C-index") + 
  ggtitle("B, C-index(GSE2658)") + 
  theme_bw() + 
  theme(legend.title = element_blank())

p <- ggarrange(p1, p2, nrow = 1, common.legend = TRUE, legend = "bottom")
# Figure 3
ggsave("../results/Figure 3.pdf", p, width = 9, height = 3.5)

#- IBS for missing
#- IBS in 9782
ibs <- extractmissing(misspred, method = "ibs", data = "9782")
imean <- aggregate(value ~ name + missing, FUN = "mean", data = ibs)
isd   <- aggregate(value ~ name + missing, FUN = "sd", data = ibs)
ibs <- merge(imean, isd, by = c("name", "missing"), suffixes = c(".mean", ".sd"))

pd <- position_dodge(0.15)
ibs <- ibs %>% mutate(missing = as.factor(missing),
                      name = recode(name, "Vax.grp" = "Vax(grp)",
                                    "Vax.ssGSEA" = "Vax(ssGSEA)"))

p1 <- ggplot2::ggplot(data = ibs, aes(x = missing, y = value.mean, color = name, group = name)) +
  geom_errorbar(aes(ymin = value.mean - value.sd, ymax = value.mean + value.sd), 
                width = .15, position = pd) +
  geom_line(position = pd) +
  geom_point(position = pd, size = 1.5) + 
  labs(x = "% of missing", y = "IBS") + 
  ggtitle("A, IBS(GSE9782)") + 
  theme_bw() + 
  theme(legend.title = element_blank())

#- cindex in 2658
ibs <- extractmissing(misspred, method = "ibs", data = "2658")
imean <- aggregate(value ~ name + missing, FUN = "mean", data = ibs)
isd   <- aggregate(value ~ name + missing, FUN = "sd", data = ibs)
ibs <- merge(imean, isd, by = c("name", "missing"), suffixes = c(".mean", ".sd"))

pd <- position_dodge(0.15)
ibs <- ibs %>% mutate(missing = as.factor(missing),
                      name = recode(name, "Vax.grp" = "Vax(grp)",
                                    "Vax.ssGSEA" = "Vax(ssGSEA)"))

p2 <- ggplot2::ggplot(data = ibs, aes(x = missing, y = value.mean, color = name, group = name)) +
  geom_errorbar(aes(ymin = value.mean - value.sd, ymax = value.mean + value.sd), 
                width = .15, position = pd) +
  geom_line(position = pd) +
  geom_point(position = pd, size = 1.5) + 
  labs(x = "% of missing", y = "IBS") + 
  ggtitle("B, IBS(GSE2658)") + 
  theme_bw() + 
  theme(legend.title = element_blank())

p <- ggarrange(p1, p2, nrow = 1, common.legend = TRUE, legend = "bottom")
ggsave("../results/Figure A1.pdf", p, width = 9, height = 3.5)

#-------------------------------------------------------------------------------



