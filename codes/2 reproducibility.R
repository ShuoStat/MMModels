
pkgs <- c("survival", "dplyr", "glmnet", "survminer", "msigdb",
          "doParallel", "doSNOW", "reshape2", "ggplot2", "tidyr",
          "ggpubr", "risksetROC", "impute", "sva", "GSVA")

test = T
for (p in pkgs) {
  tryCatch(test <- require(p, character.only = T), 
           warning = function(w) return())
  
  if(!test) {
    print(paste("Package", p, "not found. Installing Package!"))
    install.packages(p)
    BiocManager::install(p)
    require(p)
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
#- Pathway-based models
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

get.ibs <- function(X.t, X.v, y.t, y.v, beta) {
  
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
  
#- function to get C-index and ibs

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
                  
                  #- Validate procedure
                  modPred(modsPath, modsGenes, X.t, X.v, y.t, y.v)
                }

stopCluster(cl)
# save(list = "pred", file = "../output/pred.nsim.R")

#-------------------------------------------------------------------------------
#- Visualization Internal CV, Figure 2
#-------------------------------------------------------------------------------

load("../output/pred.nsim.R")
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
  select(!c("mean", "sd")) %>%
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


cindex <- filter(pred, measure == "cindex", method != "plage")

cindex.means <- cindex %>%
  group_by(paths, method) %>%
  summarise(mean = format(round(mean(value), 3), nsmall = 3),
            sd   = format(round(sd(value), 3), nsmall = 3)) %>%
  mutate(cindex = paste0(mean, "(", sd, ")")) %>%
  select(!c("mean", "sd")) %>%
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
load("../data/GSE136324.RData")
X136324 = as.matrix(X)
y136324 = survival::Surv(clin$time, clin$status)
lp136324_vax_grp <- predict.mods(modsPath, X136324, paths = "vax", method = "grp")
# lp136324_cgn_grp <- predict.mods(modsPath, X136324, paths = "cgn", method = "grp")
lp136324_vax_ssgsea <- predict.mods(modsPath, X136324, paths = "vax", method = "ssgsea")
lp136324_genes <- predict.mods(modsGenes, X136324, paths = "", method = "genes")

pred136324 <- list(lp136324_vax_grp = lp136324_vax_grp, 
                   # lp136324_cgn_grp = lp136324_cgn_grp,
                   lp136324_vax_ssgsea = lp136324_vax_ssgsea,
                   lp136324_genes = lp136324_genes)
#-
load("../data/GSE9782B.RData")
X9782 = as.matrix(X)
y9782 = survival::Surv(clin$time, clin$status)
lp9782_vax_grp <- predict.mods(modsPath, X9782, paths = "vax", method = "grp")
lp9782_vax_ssgsea <- predict.mods(modsPath, X9782, paths = "vax", method = "ssgsea")
lp9782_genes <- predict.mods(modsGenes, X9782, paths = "", method = "genes")

pred9782 <- list(lp9782_vax_grp = lp9782_vax_grp, 
                 # lp9782_cgn_grp = lp9782_cgn_grp,
                 lp9782_vax_ssgsea = lp9782_vax_ssgsea,
                 lp9782_genes = lp9782_genes)

load("../data/GSE2658B.RData")
X2658 = as.matrix(X)
y2658 = survival::Surv(clin$time, clin$status)
lp2658_vax_grp <- predict.mods(modsPath, X2658, paths = "vax", method = "grp")
lp2658_vax_ssgsea <- predict.mods(modsPath, X2658, paths = "vax", method = "ssgsea")
lp2658_genes <- predict.mods(modsGenes, X2658, paths = "", method = "genes")

pred2658 <- list(lp2658_vax_grp = lp2658_vax_grp, 
                 # lp2658_cgn_grp = lp2658_cgn_grp,
                 lp2658_vax_ssgsea = lp2658_vax_ssgsea,
                 lp2658_genes = lp2658_genes)

# save(list = c("pred136324", "pred9782", "pred2658"), file = "../output/pred.RData")
#- get C-index and IBS table 

m1 <- re9782 %>% filter((paths == "vax" & method == "grp") |
                          (paths == "vax" & method == "ssgsea") |
                          method == "genes") %>% 
  mutate(data = "GSE9782",
         method = paste0(paths, "(", method, ")")) %>%
  dplyr::select(method, measure, value) %>%
  pivot_wider(names_from = measure)

m2 <- re2658 %>% filter((paths == "vax" & method == "grp") |
                          (paths == "vax" & method == "ssgsea") |
                          method == "genes") %>% 
  mutate(data = "GSE2658",
         method = paste0(paths, "(", method, ")")) %>%
  dplyr::select(method, measure, value) %>% 
  pivot_wider(names_from = measure)

# Table 3
re <- merge(m1, m2, by = "method", suffixes = c("(GSE9782)", "(GSE2658)"))
write.csv(re, file = "../results/outvalid.csv")

#-------------------------------------------------------------------------------
#- compare with previous models
#-------------------------------------------------------------------------------
source("previous.models.R")
#- compare with previous models
#- adjust batch effects

adjDat <- function(method) {
  
  #- load previous models
  load("../data/dataPreviousModel.RData")
  load("../data/batchCorrectedData.RData")
  
  X9782 <- miss.impute(X9782, cbind(X136324, X2658))
  
  if (method == "MILLENNIUM100") {
    Ref <- batchCorrectedData$gepall
    Ref <- Ref[,batchCorrectedData$batchesall["APEX",]]
  } else if(method == "EMC92"){
    Ref <- batchCorrectedData$gepall
    Ref <- Ref[,batchCorrectedData$batchesall["H65",]]
  }  else {
    Ref <- X136324
  }
  
  #- impute for X9782
  probes.in <- rownames(Ref)
  dat <- cbind(Ref[probes.in,], X9782[probes.in,], X2658[probes.in,], 
               X136324[probes.in,])
  
  batch <- c(rep(1, ncol(Ref)), rep(2, ncol(X9782)), rep(3, ncol(X2658)), 
             rep(4, ncol(X136324)))
  
  cb <- ComBat(dat = dat,
               batch = batch,
               mod = NULL,
               par.prior = T,
               prior.plots = FALSE,
               mean.only = F,
               ref.batch = 1,
               BPPARAM = bpparam("SerialParam"))
  
  #- X9782
  X9782 <- cb[,batch == 2]
  
  #- X2658
  X2658 <- cb[,batch == 3]
  
  #- X2658
  X136324 <- cb[,batch == 4]
  
  .GlobalEnv$Ref <- Ref
  .GlobalEnv$X9782 <- X9782
  .GlobalEnv$X2658 <- X2658
  .GlobalEnv$X136324 <- X136324
}

#-
adjDat(method = "EMC92")

#- m92
m92.2658 <- EMC92Classifier(trainlogdata = Ref, testlogdata = X2658)[,"score"]
m92.9782 <- EMC92Classifier(trainlogdata = Ref, testlogdata = X9782)[,"score"]
m92.136324 <- EMC92Classifier(trainlogdata = Ref, testlogdata = X136324)[,"score"]

#- m70
adjDat(method = "UAMS70")
m70.2658 <- UAMS70Classifier(testlogdata = X2658)[,"score"]
m70.9782 <- UAMS70Classifier(testlogdata = X9782)[,"score"]
m70.136324 <- UAMS70Classifier(testlogdata = X136324)[,"score"]

#- m15
adjDat(method = "FM15")
m15.2658 <- IFM15Classifier(testlogdata = X2658, trainlogdata = Ref)[,"score"]
m15.9782 <- IFM15Classifier(testlogdata = X9782, trainlogdata = Ref)[,"score"]
m15.136324 <- IFM15Classifier(testlogdata = X136324, trainlogdata = Ref)[,"score"]

#- m100
adjDat(method = "MILLENNIUM100")
m100.2658 <- MILLENNIUM100Classifier(trainlogdata = Ref, testlogdata = X2658)[,"score"]
m100.9782 <- MILLENNIUM100Classifier(trainlogdata = Ref, testlogdata = X9782)[,"score"]
m100.136324 <- MILLENNIUM100Classifier(trainlogdata = Ref, testlogdata = X136324)[,"score"]

#- grp1
#- import our models
load("../output/pred.RData")

#- C-index
#- GSE136324
load("../data/GSE136324.RData")
c92 <- intsurv::cIndex(clin[,"time"], clin[,"status"], m92.136324)["index"]
c70 <- intsurv::cIndex(clin[,"time"], clin[,"status"], m70.136324)["index"]
c15 <- intsurv::cIndex(clin[,"time"], clin[,"status"], m15.136324)["index"]
c100 <- intsurv::cIndex(clin[,"time"], clin[,"status"], m100.136324)["index"]
#- 1, grp.vax; 2, genes; 3, ssgsea.vax
vaxgrp    <- intsurv::cIndex(clin[,"time"], clin[,"status"], pred136324$lp136324_vax_grp)["index"]
vaxssgsea <- intsurv::cIndex(clin[,"time"], clin[,"status"], pred136324$lp136324_vax_ssgsea)["index"]
genes     <- intsurv::cIndex(clin[,"time"], clin[,"status"], pred136324$lp136324_genes)["index"]

c136324 <- c("IFM15" = c15, "UAMS-70" = c70, "EMC-92" = c92, "MILLENNIUM-100" = c100,
             "vax(grp)" = vaxgrp, "vax(ssgsea)" = vaxssgsea, "genes" = genes)

#- GSE2658
load("../data/GSE2658.RData")
c92 <- intsurv::cIndex(clin[,"time"], clin[,"status"], m92.2658)["index"]
c70 <- intsurv::cIndex(clin[,"time"], clin[,"status"], m70.2658)["index"]
c15 <- intsurv::cIndex(clin[,"time"], clin[,"status"], m15.2658)["index"]
c100 <- intsurv::cIndex(clin[,"time"], clin[,"status"], m100.2658)["index"]
#- 1, grp.vax; 2, genes; 3, ssgsea.vax
vaxgrp    <- intsurv::cIndex(clin[,"time"], clin[,"status"], pred2658$lp2658_vax_grp)["index"]
vaxssgsea <- intsurv::cIndex(clin[,"time"], clin[,"status"], pred2658$lp2658_vax_ssgsea)["index"]
genes   <- intsurv::cIndex(clin[,"time"], clin[,"status"], pred2658$lp2658_genes)["index"]
c2658   <- c("IFM15" = c15, "UAMS-70" = c70, "EMC-92" = c92, "MILLENNIUM-100" = c100,
             "vax(grp)" = vaxgrp, "vax(ssgsea)" = vaxssgsea, "genes" = genes)

#- GSE9782
load("../data/GSE9782.RData")
c92 <- intsurv::cIndex(clin[,"time"], clin[,"status"], m92.9782)["index"]
c70 <- intsurv::cIndex(clin[,"time"], clin[,"status"], m70.9782)["index"]
c15 <- intsurv::cIndex(clin[,"time"], clin[,"status"], m15.9782)["index"]
c100 <- intsurv::cIndex(clin[,"time"], clin[,"status"], m100.9782)["index"]
#- 1, grp.vax; 2, genes; 3, ssgsea.vax
vaxgrp  <- intsurv::cIndex(clin[,"time"], clin[,"status"], pred9782$lp9782_vax_grp)["index"]
vaxssgsea <- intsurv::cIndex(clin[,"time"], clin[,"status"], pred9782$lp9782_vax_ssgsea)["index"]
genes <- intsurv::cIndex(clin[,"time"], clin[,"status"], pred9782$lp9782_genes)["index"]
c9782 <- c("IFM15" = c15, "UAMS-70" = c70, "EMC-92" = c92, "MILLENNIUM-100" = c100,
           "vax(grp)" = vaxgrp, "vax(ssgsea)" = vaxssgsea, "genes" = genes)

cmodels <- round(cbind(c136324 = c136324, c2658 = c2658, c9782 = c9782), 3)
# Table 4
write.csv(cmodels, file = "../results/comparemodels.csv")

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



