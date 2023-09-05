

getsplits <- function(seed, n, nsim, nfolds = 10, percentTraining = 0.7){
  
  set.seed(seed)
  ntrain = ceiling(percentTraining * n)
  split.sample <- sapply(1:nsim, function(x) seq(n) %in% sample(1:n, ntrain),
                         simplify = F)
  
  foldids      <- sapply(1:nsim, function(x) sample(rep(1:nfolds, length = ntrain), ntrain),
                         simplify = F)
  
  .GlobalEnv$foldids <- foldids
  .GlobalEnv$split.sample <- split.sample
}


#- functions to fit Lasso, ENet, and Ridge -------------------------------------

fitfun <- function(X, y, foldid = NULL, method = c("las", "enet", "rdg")){
  
  require(glmnet)
  
  n <- nrow(X)
  if(is.null(foldid)) 
    foldid <- sample(rep(1:10, length = n), n)
  
  #- lasso
  if ("las" %in% method){
    cv.las <- cv.glmnet(X, y, family = "cox", foldid = foldid)
    las <- glmnet(X, y, family = "cox", lambda = cv.las$lambda.min)
    las <- las$beta[,1]
    las <- las[las != 0]
  }
  
  #- ridge
  if ("rdg" %in% method){
    cv.rdg <- cv.glmnet(X, y, family = "cox", foldid = foldid, alpha = 0)
    rdg <- glmnet(X, y, family = "cox", lambda = cv.rdg$lambda.min, alpha = 0)
    rdg <- rdg$beta[,1]
  }
  
  #- elastic net
  if ("enet" %in% method){
    cv.enet <- cv.elastic(X, y, family = "cox", foldid = foldid,
                          type.measure = "default",
                          alpha = seq(0, 1, 0.25))
    enet <- glmnet(X, y, family = "cox", alpha = cv.enet$alpha.min,
                   lambda = cv.enet$lambda.min)
    enet <- enet$beta[,1]
    enet <- enet[enet != 0]
  }
  
  out <- list()
  for(i in method){
    out[[i]] <- get(i)
  }
  
  return(out)
}

#- prediction performance to get IBS or C-index with IPCW ----------------------

performance <- function(betas, X.t, y.t, X.v, y.v,
                        method = c("ibs", "cindex")){
  
  require(survival)
  require(pec)
  compute_coxph_mod <- function(coeffs, X, y) {
    
    dat <- cbind(y, X)
    if (length(coeffs) == 0){
      coxph(Surv(time, status) ~ 1, data = as.data.frame(dat))
    } else {
      sub_data <- dat[, c("time", "status", names(coeffs)), drop = FALSE]
      coxph(Surv(time, status) ~ .,
            data = as.data.frame(sub_data),
            init = coeffs,
            iter.max = 0,
            x = TRUE,
            y = TRUE)
    }
  }
  
  #- pfm: performance
  
  mod <- lapply(betas, function(x) compute_coxph_mod(x, X = X.t, y = y.t))
  names(mod) <- names(betas)
  
  newdata <- as.data.frame(cbind(y.v, X.v))
  f <- as.formula("Surv(time, status) ~ 1")
  pfm <- c()
  if ("ibs" %in% method){
    
    p <- pec(mod, f, data = newdata, reference = F)
    #- time, the maximum events
    tmp <- as.numeric(crps(p, times = p$minmaxtime))
    names(tmp) <- paste0("ibs.", names(betas))
    pfm <- c(pfm, tmp)
  }
  
  if ("cindex" %in% method){
    tmp <- unlist(cindex(mod, f, data = newdata)$AppCindex)
    names(tmp) <- paste0("cindx.", names(betas))
    pfm <- c(pfm, tmp)
  }
  return(pfm)
}

#- loglik for cox models
loglik.coxph <- function(pred, surv) {	
  
  surv.time <- surv[,1]
  surv.event <- surv[,2]
  
  n <- length(pred)
  r <- rank(surv.time)
  ita <- pred
  epita <- exp(ita)
  d <- rep(0, n)
  dono <- rep(0, n)
  for(i in 1:n) {
    d[i] <- sum(surv.event[r == r[i]])
    dono[i] <- sum(epita[r >= r[i]])
  }
  
  lik <- (ita - log(dono)) * surv.event
  return(lik)
}

#- calculate C-index -----------------------------------------------------------

performance.cindex <- function(betas, X.v, y.v){
  
  mapply(function(beta){
    s <- names(beta)
    pred <- X.v[,s, drop = F] %*% beta
    intsurv::cIndex(y.v[,1], y.v[,2], pred)["index"]
  }, beta = betas)
}

#- calcualte Deviance difference -----------------------------------------------

difdev <- function(pred.v, surv.v, glmnet.deviance = T){
  
  require(survival)
  
  n <- length(pred.v)
  if (glmnet.deviance){
    ll.null <- -glmnet::coxnet.deviance(rep(1, n), surv.v) / 2
    ll <- -glmnet::coxnet.deviance(pred.v, surv.v) / 2
  } else {
    ll.null <- sum(loglik.coxph(rep(1, n), surv.v))
    ll <- sum(loglik.coxph(pred.v, surv.v))
  }
  devdif <- - 2 * (ll - ll.null)
  return(devdif)
}

performance.dev <- function(betas, X.v, y.v){
  
  mapply(function(beta){
    s <- names(beta)
    pred <- X.v[,s, drop = F] %*% beta
    difdev(pred, y.v, T)
  }, beta = betas)
}

#-------------------------------------------------------------------------------
#- Kepaln Meier plot (using median split)
#-------------------------------------------------------------------------------

plot.os.fit <- function(x, os, title = ""){
  
  require(survminer)
  warning("x should be the same order of os")
  x <- as.factor(ifelse(x > median(x), "High", "Low"))
  os.2 <- data.frame(time = os[,1], 
                     status = os[,2],
                     x = x)
  fit <- survfit(Surv(time, status) ~ x, data = os.2)
  
  p <- ggsurvplot(fit, data = os.2, pval = T, legend = c(0.8, 0.8), 
                  conf.int = F, surv.median.line = "hv", 
                  title = title, risk.table = F, 
                  ggtheme = theme_bw(),
                  color = "x",
                  palette = c("firebrick1", "cyan4"))
  
  p$plot <- p$plot + theme(plot.margin = unit(c(10, 10, 10, 20), "points"),
                           plot.title = element_text(size = 15))
  return(p)
}

#- Log-rank test 

os.fit <- function(x, os){
  
  x <- as.factor(ifelse(x > median(x), "High", "Low"))
  os.2 <- cbind(os, x = x)
  fit <- survfit(Surv(time, status) ~ x, data = os.2)
  p <- surv_pvalue(fit, data = os.2)$pval
  m <- setNames(surv_median(fit)$median, surv_median(fit)$strata)
  list(median = m, p = p)
}

#- ENet methods ----------------------------------------------------------------

get.lambda.cox <- function(X, y, 
                           nlam = 100, 
                           lambda.min.ratio = 0.01,
                           alpha, 
                           weights = rep(1, nrow(X)),
                           offset = rep(0, nrow(X)), 
                           exclude = c(),
                           pfs = rep(1, ncol(X))){
  
  mysd <- function(x) sqrt(sum((x-mean(x))^2) / length(x))
  sx <- scale(X, scale = apply(X, 2, mysd))
  sx <- as.matrix(sx)
  
  nobs <- nrow(X); nvars <- ncol(X)
  # extract strata (if any)
  if ("strata" %in% names(attributes(y))) {
    strata <- attr(y, "strata")
  } else {
    strata <- rep(1, nobs)
  }
  if (length(strata) != nobs) stop("length of strata != nobs")
  
  # if some penalty factors are zero, we need to compute eta
  
  if (is.list(pfs)){
    pfs <- lapply(pfs, function(x){
      x = pmax(0, x)
      if (max(x) <= 0) stop("All penalty factors are <= 0")
      as.double(x * nvars / sum(x))
    })}else{
      pfs = pmax(0, pfs)
      if (max(pfs) <= 0) stop("All penalty factors are <= 0")
      pfs <- list(as.double(pfs * nvars / sum(pfs)))
    }
  
  lambda_max <- list()
  
  for(i in seq_along(pfs)){
    
    pf <- pfs[[i]]
    pf_zero <- setdiff(which(pf == 0), exclude)
    
    if (length(pf_zero) > 0) {
      
      tempx <- sx[, pf_zero, drop = FALSE]
      eps <- glmnet.control()$epsnr
      if (length(unique(strata)) == 1) {
        fit <- survival::coxph(y ~ offset(offset) + tempx, 
                               weights = weights, 
                               eps = eps)
      } else {
        fit <- survival::coxph(y ~ offset(offset) + tempx + strata(strata),
                               weights = weights, 
                               eps = eps)
      }
      eta <- predict(fit, reference = "sample") ## Coxph can do strata-specific centering
    } else {
      eta <- offset
    }
    
    # keep numbers small; partial likelihood independent of centering
    eta <- eta - mean(eta)  
    ju  <- rep(1, nvars)
    # we have already included constant variables in exclude
    ju[exclude] <- 0 
    
    # Get cox gradient at "null" point
    # Note that coxgrad already includes weights, so no need to include them
    # In subsequent computations
    null_grad <- glmnet::coxgrad(eta = eta, y = y, w = weights)
    g <- abs(drop(t(null_grad) %*% sx))
    g <- g * ju / ifelse(pf > 0, pf, 1)
    lambda_max[[i]] <- max(g) / max(alpha, 1e-3)
  }
  
  lambda_max <- max(unlist(lambda_max))
  exp(seq(log(lambda_max), log(lambda_max * lambda.min.ratio),
          length.out = nlam))
}

get.lambda <- function(X, y, family, alpha, 
                       intercept = T, 
                       is.offset = F, 
                       offset = rep(0, nrow(X)), 
                       weights = rep(1, nrow(X)), exclude = c(), 
                       pfs = rep(1, ncol(X)), 
                       nlam = 100, 
                       lambda.min.ratio = 0.01) {
  
  family = do.call(family, list())
  nobs <- nrow(X); nvars <- ncol(X)
  
  mysd <- function(x) sqrt(sum((x-mean(x))^2) / length(x))
  sx <- scale(X, scale = apply(X, 2, mysd))
  sx <- as.matrix(sx)
  
  # compute mu and null deviance
  # family = binomial() gives us warnings due to non-integer weights
  # to avoid, suppress warnings
  
  if (intercept) {
    if (is.offset) {
      suppressWarnings(tempfit = glm(y ~ 1, family = family,
                                     weights = weights, offset = offset))
      mu <- tempfit$fitted.values
    } else {
      mu <- rep(weighted.mean(y, weights), times = nobs)
    }
  } else {
    mu <- family$linkinv(offset)
  }
  
  dev_function <- function(y, mu, weights, family) {
    sum(family$dev.resids(y, mu, weights))
  }
  
  nulldev <- dev_function(y, mu, weights, family)
  
  #- standardize pfs, pfs can be a list
  if (is.list(pfs)){
    pfs <- lapply(pfs, function(x){
      x = pmax(0, x)
      if (max(x) <= 0) stop("All penalty factors are <= 0")
      as.double(x * nvars / sum(x))
    })} else {
      pfs = pmax(0, pfs)
      if (max(pfs) <= 0) stop("All penalty factors are <= 0")
      pfs <- list(as.double(pfs * nvars / sum(pfs)))
    }
  
  lambda_max <- list()
  for(i in seq_along(pfs)){
    
    pf <- pfs[[i]]
    # if some penalty factors are zero, we have to recompute mu
    pf_zero <- setdiff(which(pf == 0), exclude)
    if (length(pf_zero) > 0) {
      tempx <- sx[, pf_zero, drop = FALSE]
      
      if (intercept) {
        tempx <- cbind(1,tempx)
      }
      tempfit <- glm.fit(tempx, y, family = family, weights = weights, 
                         offset = offset)
      mu <- tempfit$fitted.values
    }
    
    # compute lambda max
    ju <- rep(1, nvars)
    ju[exclude] <- 0 # we have already included constant variables in exclude
    r <- y - mu
    eta <- family$linkfun(mu)
    v <- family$variance(mu)
    m.e <- family$mu.eta(eta)
    weights <- weights / sum(weights)
    rv <- r / v * m.e * weights
    g <- abs(drop(t(rv) %*% sx))
    
    g <- g * ju / ifelse(pf > 0, pf, 1)
    lambda_max[[i]] <- max(g) / max(alpha, 1e-3)
  }
  
  lambda_max <- max(unlist(lambda_max))
  exp(seq(log(lambda_max), log(lambda_max * lambda.min.ratio),
          length.out = nlam))
}

lambda.fun <- function(X, y, family, alpha, pf, lambda, pf.in.beta){
  
  fit <- glmnet(x = X, y = y,
                family = family,
                lambda = lambda,
                alpha = alpha,
                penalty.factor = pf)
  
  if (!pf.in.beta) {
    
    for(i in seq_along(lambda)){
      sel <- as.matrix(fit$beta)[,i] != 0
      if(sum(sel) < 2)
        next
      obj <- glmnet(x = X[,sel], y = y,
                    family = family,
                    lambda = lambda[i],
                    alpha = alpha)
      fit$a0[i] <- obj$a0
      fit$beta[sel, i, drop = F] <- obj$beta[,1, drop = F]
      fit$beta[,i] = fit$beta[,i] * sel
      fit$dev.ratio[i] <- obj$dev.ratio
    }
  }
  
  #- pmat, selected variables by lambda
  pmat = fit$beta != 0
  
  #- initial beta, a0, dev, nulldev
  beta = fit$beta
  a0   = fit$a0
  dev  = fit$dev.ratio
  
  
  if (family != "cox"){
    #- get ols estimates
    for(i in 1:ncol(pmat)){
      if (any(pmat[,i])) {
        exclude = which(!pmat[,i])
        r.fit = update(fit, x = X, lambda = 0, exclude = exclude,
                       relax = FALSE, penalty.factor = rep(1, ncol(X)))
        
        a0[i] = r.fit$a0
        beta[,i] = r.fit$beta[,1]
        dev[i] = r.fit$dev.ratio
      }
    }
  } else {
    
    for(i in 1:ncol(pmat)){
      
      if (any(pmat[,i])) {
        include = which(pmat[,i])
        #- if number of var > events; just consider 0 for coefficients
        if (length(include) > sum(y[,2])){
          a0[i] = 0
          beta[include, i] = 0
          dev[i] = 0
        }else{
          D <- as.data.frame(X[,include, drop = F])
          mod <- survival::coxph(y ~ ., data = D)
          a0[i] = 0
          beta[include, i] = mod$coefficients
          dev[i] = 0
        }
      }
    }
  }
  
  fit.glm <- fit
  fit$beta <- beta
  fit$a0 <- a0
  fit$dev.ratio <- dev
  fit$alpha <- alpha
  
  fit.glm$relaxed = fit
  class(fit.glm) <- class(fit)
  return(fit.glm)
}

#- used in elastic function 

cvstats <- function(cvstuff,foldid,nfolds,lambda,grouped,...){
  
  if (grouped){
    nlams=rep(dim(cvstuff$cvraw)[2],nfolds) ## All the same - this should go
    cvstuff= cvcompute(cvstuff, foldid, nlams)
  }
  cvm=with(cvstuff,apply(cvraw, 2, weighted.mean, w = weights, na.rm = TRUE))
  cvsd=with(cvstuff, sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, weighted.mean,
                                w = weights, na.rm = TRUE)/(N - 1)))
  nas=is.na(cvsd)
  if(any(nas)){
    lambda=lambda[!nas]
    cvm=cvm[!nas]
    cvsd=cvsd[!nas]
  }
  list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvup = cvm + cvsd, 
       cvlo = cvm - cvsd)
}

#- used in elastic function 
getstat <- function(statlist, lambda, type.measure){
  
  alpha.nam <- as.numeric(gsub("alpha:", "", names(statlist)))
  alp  <- rep(alpha.nam, each = length(lambda))
  lam    <- rep(lambda, times = length(alpha.nam))

  cvm  <- unlist(lapply(statlist, `[[`, "cvm"))
  cvsd <- unlist(lapply(statlist, `[[`, "cvsd"))
  
  if (type.measure %in% c("auc", "C"))
    cvm <- -cvm
  
  #- get min and order them 
  id.min <- which(cvm <= min(cvm), arr.ind = T)
  id.min <- id.min[order(alp[id.min], lam[id.min], decreasing = T)[1]]
  lambda.min <- lam[id.min]
  alpha.min  <- alp[id.min]
  
  #- get 1sd and order them
  cvm.1se <- cvm[id.min] + cvsd[id.min]
  id.1se  <- which(cvm <= cvm.1se, arr.ind = T)
  id.1se  <- id.1se[order(alp[id.1se], lam[id.1se], decreasing = T)[1]]
  lambda.1se <- lam[id.1se]
  alpha.1se  <- alp[id.1se]
  
  list(statlist = statlist,
       alpha = alpha.nam,
       lambda = lambda,
       lambda.min = lambda.min,
       lambda.1se = lambda.1se,
       alpha.min = alpha.min,
       alpha.1se = alpha.1se)
}


#- used in elastic function 

alpha.fun <- function(X, y, family, lambda, alpha, pf, nfolds, 
                      foldid, type.measure, grouped, weights){
  
  alphalist = list()
  for(a in alpha){
    
    lamlist = list()
    for (i in 1:nfolds){
      
      ind = foldid != i
      lamlist[[i]]<- glmnet(X[ind,,drop = F], y[ind], family = family, 
                            alpha = a, penalty.factor = pf, 
                            lambda = lambda)
    }
    
    subclass = class(lamlist[[1]])
    class(lamlist) <- paste0(subclass, "list")
    
    predmat = glmnet:::buildPredmat(lamlist, lambda, X, offset = NULL, foldid,
                                    alignment = "lambda", y = y, weights = weights,
                                    grouped = grouped, type.measure = type.measure,
                                    family = family(lamlist[[1]]))
    
    cvstuff <- switch(family,
                      "gaussian" = glmnet:::cv.elnet(predmat, y, type.measure, weights, foldid, grouped),
                      "cox"      = glmnet:::cv.coxnet(predmat, y, type.measure, weights, foldid, grouped),
                      "binomial" = glmnet:::cv.lognet(predmat, y, type.measure, weights, foldid, grouped))
    
    nam <- paste0("alpha:", a)
    
    alphalist[[nam]] <- 
      cvstats(cvstuff, foldid, nfolds, lambda, grouped = cvstuff$grouped)    
  }
  alphalist
}

#- used in elastic function 

cv.elastic <- function(X, y, family, foldid, type.measure = "default", 
                       lambda = NULL, nlam = 100, alpha = seq(0, 1, 0.25), 
                       lambda.min.ratio = 0.01, weights = rep(1, nrow(X)),
                       penalty.factor = rep(1, ncol(X)), grouped = F){
  
  if (type.measure == "default")
    type.measure <- switch(family,
                           "gaussian" = "mse",
                           "binomial" = "deviance",
                           "cox"      = "deviance")
  
  
  if (is.null(lambda)) {
    if (family == "cox"){
      lambda = get.lambda.cox(X, y, 
                              nlam = nlam, 
                              alpha = alpha,
                              lambda.min.ratio = lambda.min.ratio,
                              weights = weights,
                              pfs = penalty.factor)
    }else{
      lambda = get.lambda(X, y, 
                          family = family,
                          alpha = alpha,
                          nlam = nlam, 
                          lambda.min.ratio = lambda.min.ratio,
                          weights = weights,
                          pfs = penalty.factor)
    }
    
    require(glmnet)
    p = ncol(X); n = nrow(X)
    
    if (is.null(foldid)) 
      foldid <- sample(rep(1:nfolds, length = n), n)
    
    nfolds = max(foldid)
    cv.alpha <- alpha.fun(X, y, family, lambda, alpha, penalty.factor, nfolds, 
                          foldid, type.measure, grouped, weights)
    
    getstat(cv.alpha, lambda, type.measure)
  }
}

#- get pathways ----------------------------------------------------------------

dist.kul <- function(obj){
  
  kul <- function(v1, v2){
    s1 <- sum(v1 %in% v2)
    s2 <- length(v1)
    s3 <- length(v2)
    1 - (s1 / (2 * s2) + s1 / (2 * s3))
  }
  
  n <- length(obj)
  mat <- matrix(NA, n, n, dimnames = list(names(obj), names(obj)))
  for(i in 1:n){
    for(j in i:n){
      v1 <- obj[[i]]
      v2 <- obj[[j]]
      mat[i, j] <- mat[j, i] <- kul(v1, v2)
    }
  }
  
  mat <- mat[lower.tri(mat)]
  
  attr(mat, "Size") <- n
  attr(mat, "Diag") <- F
  attr(mat, "Upper") <- F
  attr(mat, "method") <- "kul"
  attr(mat, "Labels") <-  names(obj)
  return(mat)
}

cluster.go <- function(obj, k = 3, plot = F){
  
  d <- dist.kul(obj)
  h <- hclust(d, method = "average")
  if (plot)
    plot(h, labels = names(obj))
  g <- cutree(h, k = k)
  setNames(g, names(obj))
}

getpaths <- function(paths, genes, N, min.genes = 5, rm.duplicated = T){
  
  paths <- lapply(paths, intersect, y = genes)
  small <- unlist(lapply(paths, function(x) length(x) >= min.genes))
  if (sum(small) > 1)
    cat(paste0("A total of ", sum(!small), 
               " pathways were removed because of less than ", 
               min.genes, " genes in pathways \n\n"))
  
  paths <- paths[small]
  
  #- remove duplicated pathways: from bottom to top
  n1 <- length(paths)
  
  if (rm.duplicated) {
    
    #- detect and remove duplicated pathways
    tmp <- lapply(paths, function(x) paste0(sort(x), collapse = ""))
    paths <- paths[!duplicated(tmp)]
    
    if (any(duplicated(tmp)))
      cat("Duplicated pathways: ", 
          sum(tmp), "duplicated pathways were removed \n\n")
  }
  
  #- clustering
  if (N >= length(paths)) {
    
    cat("N > # of pathways: ", length(paths), "pathways returned \n\n")
    
    return(paths)} else {
      
      paths.c <- cluster.go(paths, k = N, plot = F)
      paths.s <- c()
      
      for(i in 1:N){
        s <- dist.kul(paths[paths.c == i])
        n <- attr(s, "Size")
        mat <- matrix(NA, n, n)
        mat[lower.tri(mat)] <- s
        
        mat <- Matrix::forceSymmetric(mat, uplo = "L")
        diag(mat) <- 0
        paths.s <- c(paths.s, attr(s,"Labels")[which.min(colSums(mat))])
      }
      return(paths[paths.s])
    }
}

#-------------------------------------------------------------------------------
#- missing data imputation
#-------------------------------------------------------------------------------

miss.impute <- function(valid, training) {
  
  require(impute)
  #- For valid, select all features including in training features
  sel <- intersect(rownames(valid), rownames(training))
  valid <- valid[sel, ]
  
  training <- as.matrix(training)
  
  #- tranforma valid to training data structure
  ind <- match(rownames(valid), rownames(training))
  mat <- matrix(NA, nrow = nrow(training), ncol = ncol(valid), 
                dimnames = list(rownames(training),
                                colnames(valid)))
  mat[ind, ] <- valid
  
  for(i in seq_len(ncol(valid))){
    
    cat("Processing:", colnames(valid)[i], "\n")
    # sink("nul")
    out <- impute.knn(cbind(training, ref = mat[,i]), k = 10,
                      colmax = 0.9)
    # sink()
    mat[,i] <- out$data[,"ref"]
  }
  mat
}

#-------------------------------------------------------------------------------
#- Group Lasso
#-------------------------------------------------------------------------------

grplas.surv <- function(X, y, pathways, foldid){
  
  require(dplyr)
  require(grpreg)
  
  colNam <- colnames(X)
  
  if (!is.matrix(X))
    X <- as.matrix(X)
  
  #- check pathways
  ps     <- unique(unlist(pathways))
  if (any(!ps %in% colNam))
    stop(paste0("These variables cannot be found:",
                paste0(setdiff(ps, colNam), collapse = ",")))
  
  #- remove duplicated 
  vs <- c()
  grps <- c()
  for(i in seq_along(pathways)){
    vs   <- c(vs, match(pathways[[i]], colNam))
    grps <- c(grps, rep(i, length(pathways[[i]])))
  }
  
  cv <- cv.grpsurv(X = X[,vs],  
                   y = y,           
                   group = grps,  
                   penalty = "grLasso",         
                   fold = foldid)
  
  fit <- grpsurv(X = X[,vs], 
                 y = y,
                 group  = grps,
                 lambda = cv$lambda.min)

  beta <- fit$beta[,1]
  grpBetas <- aggregate(beta ~ grps, FUN = "sum") %>% filter(beta != 0)
  
  if (is.null(names(pathways)))
    nonZeroGrps <- seq_along(pathways)[grpBetas$grps]
  else
    nonZeroGrps <- names(pathways)[grpBetas$grps]
  
  Betas <- aggregate(beta ~ names(beta), FUN = "sum") %>% filter(beta != 0)
  beta  <- setNames(Betas$beta, Betas$`names(beta)`)
  
  out <- list(beta = beta,
              nonZeroGrps = nonZeroGrps)
  return(out)
}

#-------------------------------------------------------------------------------
#- Multiple ssgesa 
#-------------------------------------------------------------------------------

biMeans <- function(X, sets, minGenes = 10){
  
  require(dplyr)
  
  if (!is.list(sets))
    sets = list(sets)
  
  allGenes  <- rownames(X)
  geneRanks <- X %>% 
    as.data.frame() %>%
    mutate_all( ~ rev(rank(.)))
  
  #- filter sets
  ind <- unlist(lapply(sets, function(x) 
    length(intersect(x, rownames(X))) >= minGenes))
  sets <- sets[ind]
  
  if (sum(!ind) < 0)
    message(paste("A total of", sum(!ind), "pathways are removed with minGenes = 10"))
  
  #- output
  meanDat <- c()
  setsCluster <- list()
  
  for(i in seq_along(sets)){
    
    set <- sets[[i]]
    difGenes  <- setdiff(set, allGenes)
    bothGenes <- intersect(allGenes, set)
    
    #- trans Genes to Ranks
    selGeneRanks <- geneRanks[bothGenes, ]
    setCluster   <- kmeans(selGeneRanks, centers = 2)$cluster
    
    sampleMeans <- selGeneRanks %>% 
      bind_cols(grps = setCluster) %>% 
      group_by(grps) %>%
      summarise(across(everything(), mean)) %>%
      as.data.frame() %>% 
      t()
    
    sampleMeans <- sampleMeans[-1,]
    
    if (is.unsorted(colMeans(sampleMeans)))
      colnames(sampleMeans)  <- c("Low", "High")
    else
      colnames(sampleMeans)  <- c("High", "Low")
    
    setsCluster[[i]] <- setCluster
    meanDat <- cbind(meanDat, sampleMeans)
  }
  
  colnames(meanDat) <- paste0(rep(names(sets), each = 2), "_", colnames(meanDat))
  names(setsCluster) <- names(sets)
  list(Z = t(meanDat),
       clusters = setsCluster)
  #- output, row features
}

#-------------------------------------------------------------------------------




















