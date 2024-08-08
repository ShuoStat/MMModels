

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
source("previous.models.R")

#-------------------------------------------------------------------------------
# computational time
#-------------------------------------------------------------------------------

# time on two data
# missing data imputation
# pathway score
# prediction

# load X136324
# load hs

timeFun <- function(newX, modsPath, modsGenes,
                    path, method) {
    
    # batch time
    t1 <- Sys.time()
    batch <- c(rep(1, ncol(X136324)), rep(2, ncol(newX)))
    cb <- ComBat(dat = cbind(X136324, newX),
                 batch = batch,
                 mod = NULL,
                 par.prior = T,
                 prior.plots = FALSE,
                 mean.only = F,
                 ref.batch = 1,
                 BPPARAM = bpparam("SerialParam"))
    
    t2 <- Sys.time()
    tBatch <- as.numeric(t2 - t1, units = "secs")
    Z <- cb[, batch == 2]
    
    # pathway score time
    if (method %in% c("ssgsea", "zscore", "gsva")){
      
      mod  <- modsPath[[path]]
      beta <- mod[[method]]  
      
      t1 <- Sys.time()
      Z <- GSVA::gsva(Z, hs[[path]][names(beta)], 
                      method = method, 
                      ssgsea.norm = FALSE)
      t2 <- Sys.time()
      tPathScore <- as.numeric(t2 - t1, units = "secs")
    } else {
      tPathScore <- 0
    }
    
    # prediction 
    t1 <- Sys.time()
    if (method == "genes"){
        beta <- modsGenes
        lp <- t(Z[names(beta),]) %*% beta
    }
    
    # grp
    if (method == "grp") {
        mod  <- modsPath[[path]]
        beta <- mod$grp$beta
        lp   <- t(Z[names(beta),]) %*% beta
    }
    
    # pathway score methods
    if(method %in% c("ssgsea", "zscore", "gsva")){
        mod  <- modsPath[[path]]
        beta <- mod[[method]]
        lp   <- t(Z[names(beta),]) %*% beta
    }
    
    t2 <- Sys.time()
    tPred <- as.numeric(t2 - t1, units = "secs")
    
    # missing data imputation
    t1 <- Sys.time()
    n <- nrow(newX)
    # generate 1% missing
    set.seed(1)
    indMiss <- sample(seq_len(n), n * 0.01)
    newX <- newX[-indMiss, ]
    miss.impute(newX, X136324)
    t2 <- Sys.time()
    tMiss <- as.numeric(t2 - t1, units = "secs")
    
    out <- list(tBatch = tBatch,
                tPathScore = tPathScore,
                tPred = tPred, 
                tMiss = tMiss)
    return(out)
}
    

path_method <- list(c("vax", "grp"),
                    c("vax", "ssgsea"),
                    c("", "genes"))

# load data
load("../output/training.RData")
load("../data/GSE136324.RData")
X136324 = as.matrix(X)
load("../output/hs.RData")

# time X9782
load("../data/GSE9782B.RData")
X9782 = as.matrix(X)
# y9782 = survival::Surv(clin$time, clin$status)
t9872 <- mapply(timeFun, 
                path = c("vax", "vax", ""), 
                method = c("grp", "ssgsea", "genes"),
                MoreArgs = list(newX = X9782, 
                                modsPath = modsPath, 
                                modsGenes = modsGenes),
                SIMPLIFY = T)

colnames(t9872) <- c("Vax(grp)", "Vax(ssgsea)", "Genes")

 
# time 2
load("../data/GSE2658B.RData")
X2658 = as.matrix(X)
t2658 <- mapply(timeFun, 
                path = c("vax", "vax", ""), 
                method = c("grp", "ssgsea", "genes"),
                MoreArgs = list(newX = X2658, 
                                modsPath = modsPath, 
                                modsGenes = modsGenes))


colnames(t2658) <- c("Vax(grp)", "Vax(ssgsea)", "Genes")
time <- cbind(t9872, t2658)
write.csv(time, "../results/TableTime.csv")

save(list = c("t9872", "t2658"), file = "../output/time.RData")


#-------------------------------------------------------------------------------
# K-M plot for competitor models
#------------------------------------------------------------------------------- 

lpFun <- function(dat, method, sample = TRUE) {
  
  gFun <- function (dat, datName, method) {
    # datName, e.g., "GSE9782"
    datID <- gsub("GSE", "", datName)
    n <- ncol(dat[[paste0("X", datID)]])
    
    if (sample)
      samInd <- sample(seq_len(n), n, replace = TRUE)
    else
      samInd <- seq_len(n)
    
    parList <- list(trainlogdata = dat$Ref,
                    testlogdata  = dat[[paste0("X", datID)]][,samInd])
    
    lp <- do.call(method, parList)[ ,"score"]
    return(lp)
  }
  
  out <- list("GSE136324" = gFun(dat, "GSE136324", method = method),
              "GSE9782" = gFun(dat, "GSE9782", method = method), 
              "GSE2658" = gFun(dat, "GSE2658", method = method))
  return(out)
}


#- 

LpList <- list()
methods <- c("IFM15Classifier",  "UAMS70Classifier", "EMC92Classifier", "MILLENNIUM100Classifier")

for(method in methods){
  
  dat <- adjDat(method)
  lp  <- lpFun(dat, method = method, sample = FALSE)
  lp  <- c(lp, dat[grepl("^y", names(dat))])
  LpList[[method]] <- lp
  # combine y list
}

# get learner prediction for proposed models, path = vax, method = grp.
predTrain <- function(newX, path, method) {
  
  load("../output/training.RData")
  mod  <- modsPath[[path]]
  beta <- mod[[method]]$beta
  
  newZ <- newX[names(beta), , drop = F]
  lp <- t(newZ) %*% beta
  return(lp)
}

load("../data/GSE136324.RData")
X136324  <- as.matrix(X)
y136324  <- survival::Surv(clin$time, clin$status)
lp136324 <- predTrain(X136324, path = "vax", method = "grp")

load("../data/GSE9782B.RData")
X9782 <- as.matrix(X)
y9782 <- survival::Surv(clin$time, clin$status)
lp9782 <- predTrain(newX = X9782, path = "vax", method = "grp")

load("../data/GSE2658B.RData")
X2658 <- as.matrix(X)
y2658 <- survival::Surv(clin$time, clin$status)
lp2658 <- predTrain(newX = X2658, path = "vax", method = "grp")

# K-M plot by data
# GSE136324

library(survival)
pList <- list()
pList[["g1"]]  <- plot.os.fit(lp136324, y136324, title = "GSE136324 \n Vax(grp)*")
pList[["g2"]]  <- plot.os.fit(LpList$IFM15Classifier$GSE136324, y136324, title = "\n IFM15")
pList[["g3"]]  <- plot.os.fit(LpList$UAMS70Classifier$GSE136324, y136324, title = "\n UAMS-70")
pList[["g4"]]  <- plot.os.fit(LpList$EMC92Classifier$GSE136324, y136324, title = "\n EMC-92")
pList[["g5"]]  <- plot.os.fit(LpList$MILLENNIUM100Classifier$GSE136324, y136324, title = "\n MILLENNIUM-100")

# GSE2658
pList[["g6"]] <- plot.os.fit(lp2658, y2658, title = "GSE2658 \n Vax(grp)")
pList[["g7"]] <- plot.os.fit(LpList$IFM15Classifier$GSE2658, y2658, title = "\n IFM15")
pList[["g8"]] <- plot.os.fit(LpList$UAMS70Classifier$GSE2658, y2658, title = "\n UAMS-70*")
pList[["g9"]] <- plot.os.fit(LpList$EMC92Classifier$GSE2658, y2658, title = "\n EMC-92")
pList[["g10"]] <- plot.os.fit(LpList$MILLENNIUM100Classifier$GSE2658, y2658, title = "\n MILLENNIUM-100")

# GSE9782
pList[["g11"]] <- plot.os.fit(lp9782, y9782, title = "GSE9782 \n Vax(grp)")
pList[["g12"]] <- plot.os.fit(LpList$IFM15Classifier$GSE9782, y9782, title = "\n IFM15")
pList[["g13"]] <- plot.os.fit(LpList$UAMS70Classifier$GSE9782, y9782, title = "\n UAMS-70")
pList[["g14"]] <- plot.os.fit(LpList$EMC92Classifier$GSE9782, y9782, title = "\n EMC-92")
pList[["g15"]] <- plot.os.fit(LpList$MILLENNIUM100Classifier$GSE9782, y9782, title = "\n MILLENNIUM-100*")

pList <- lapply(pList, `[[`, "plot")

for(i in seq_along(pList)){
  print(i)
  pList[[i]] <- pList[[i]] + 
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
    theme(plot.title = element_text(size = 12))
}

# plot(grid_layout)
pout <- ggarrange(plotlist = pList,
                  nrow = 3, ncol = 5,
                  common.legend = TRUE, legend = "bottom", align = "h")

ggsave("../results/KM.pdf", plot = pout, width = 12, height = 7.5)
ggsave("../results/KM.tiff", plot = pout, width = 12, height = 7.5, units = "in", 
       dpi = 300, compression = "lzw")



################################################################################
# Likelihood ratio test
################################################################################

# use LpList obtained from K-M plot 
# use lp136324, lp9782, and lp2658

# data, GSE136324
cox136324  <- coxph(y136324 ~ lp136324, init = 1, iter.max = 0)
cox136324s <- coxph(y136324 ~ scale(lp136324), init = 1, iter.max = 0)

lp136324List  <- lapply(LpList, function(x) {x[[ "GSE136324"]]})
# anova test
lpOthers136324 <- lapply(lp136324List[1:3], function(x){
  m <- coxph(y136324 ~ x, init = 1, iter.max = 0)
  ano <- anova(m, cox136324)
  c(ano$loglik, round(ano$`Pr(>|Chi|)`, 5))
}
)

# data, GSEGSE9782
cox9782  <- coxph(y9782 ~ lp9782, init = 1, iter.max = 0)

lp9782List  <- lapply(LpList, function(x) {x[[ "GSE9782"]]})
lpOthers9782  <- lapply(lp9782List, function(x){
  m <- coxph(y9782 ~ x, init = 1, iter.max = 0)
  ano <- anova(m, cox9782)
  c(ano$loglik, round(ano$`Pr(>|Chi|)`, 5))
  }
)

# data, GSE2658

cox2658  <- coxph(y2658 ~ lp2658, init = 1, iter.max = 0)
lp2658List  <- lapply(LpList, function(x) {x[[ "GSE2658"]]})
lpOthers2658  <- lapply(lp2658List, function(x){
  m <- coxph(y2658 ~ x, init = 1, iter.max = 0)
  ano <- anova(m, cox2658)
  c(ano$loglik, round(ano$`Pr(>|Chi|)`, 5))
}
)

#-------------------------------------------------------------------------------
# Prediction of Single patients
#-------------------------------------------------------------------------------


load("../output/training.RData")
load("../data/GSE136324.RData")
X136324 <- t(as.matrix(X))
y136324 <- survival::Surv(clin$time, clin$status)
timeMax <- max(y136324[y136324[,2] == 1,1])
times <- seq(0, timeMax, length = 20)

# proposed model
beta <- modsPath$vax$grp$beta

# validation data
#- GSE2658
load("../data/GSE2658B.RData")
X2658 <- t(as.matrix(X))

# prediction
probs <- survP(beta, trainX = X136324, trainY = y136324, validX = X2658, times = times)
colnames(probs) <- times

plotDat <- probs %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>%
  tidyr::pivot_longer(cols = c(-1), names_to = "time", values_to = "survP") %>%
  mutate(time = as.numeric(time))


p1 <- ggplot(plotDat, aes(x = time, y = survP, group = sample)) +
  geom_line(linewidth = 0.1) + 
  theme_bw() + 
  labs(x = "Times", y = "Survival Probabilities", title = "A, GSE2658") 

p1


# GSE9782

load("../data/GSE9782B.RData")
X9782 <- t(as.matrix(X))

# prediction
probs <- survP(beta, trainX = X136324, trainY = y136324, validX = X9782, times = times)
colnames(probs) <- times

plotDat <- probs %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>%
  tidyr::pivot_longer(cols = c(-1), names_to = "time", values_to = "survP") %>%
  mutate(time = as.numeric(time))


p2 <- ggplot(plotDat, aes(x = time, y = survP, group = sample)) +
  geom_line(linewidth = 0.1) + 
  theme_bw() + 
  labs(x = "Times", y = "Survival Probabilities", title = "B, GSE9782") 
  
p2

pout <- ggpubr::ggarrange(p1, p2, nrow = 1, ncol = 2)

ggsave("../results/survP.pdf", plot = pout, width = 8, height = 3.5)
ggsave("../results/survP.tiff", plot = pout, width = 8, height = 3.5, units = "in", 
       dpi = 300, compression = "lzw")










