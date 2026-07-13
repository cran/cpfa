pficpfa <- 
  function(object, nshuffles = 10, type = c("marginal", "conditional"), 
           conditional.model = c("ridge", "rf"), ridge.lambda = 1e-4, 
           ntree = 500, nodesize = 5, safealign = FALSE, 
           safealign.stat = c("min", "mean", "median"), 
           safealign.threshold = 0.9, parallel = FALSE, cl = NULL) 
{
  if (!(inherits(object, "wrapcpfa"))) {
    stop("Input 'object' must be of class 'wrapcpfa'.")
  }
  opt.model.all <- object$opt.model 
  opt.param <- object$opt.tune
  omethods <- c("PLR", "SVM", "RF", "NN", "RDA", "GBM")
  method <- as.character(which(omethods %in% toupper(object$method) == TRUE))
  nfac <- object$nfac
  family <- object$family
  predstor <- object$predweights
  trweights <- object$train.weights$Classtrain.weights
  y <- object$y
  lnfac <- length(nfac)
  didalign <- object$aligned
  if (!(isTRUE(didalign))) {
    stop("Input 'object' from function 'cpfa' contained 'aligned = FALSE'. \n
         Must set 'cpfa' argument 'align' to 'TRUE' to implement permutation \n
         feature importance calculations.")
  }
  nrep <- length(object$testIDs)
  storrows <- length(object$testIDs[[1]])
  threshold <- 0.5
  if ((!(ceiling(nshuffles) == nshuffles)) || (nshuffles < 1)) {
    stop("Input 'nshuffles' must be an integer greater than 0.")
  }
  numcheck(nshuffles)
  type <- match.arg(type, c("marginal", "conditional"))
  if (type == "conditional") {
    conditional.model <- match.arg(conditional.model, c("ridge", "rf"))
    if (conditional.model == "ridge") {
      numcheck(ridge.lambda)
      if (ridge.lambda <= 0) {
        stop("Input 'ridge.lambda' must be a positive number.")
      }
    } else {
      numcheck(ntree)
      numcheck(nodesize)
      if ((ntree < 1) || (ntree != as.integer(ntree))) {
        stop("Input 'ntree' must be a single positive integer.")
      }
      if ((nodesize < 1) || (nodesize != as.integer(nodesize))) {
        stop("Input 'nodesize' must be a single positive integer.")
      }
    }
  }
  logicheck(safealign)
  safealign.stat <- match.arg(safealign.stat, c("min", "mean", "median"))
  numcheck(safealign.threshold)
  if ((safealign.threshold <= 0) || (safealign.threshold >= 1)) {
    stop("Input 'safealign.threshold' must be between 0 and 1, exclusive.")
  }
  if (safealign == TRUE) {
    statf <- match.fun(safealign.stat)
    goodreps <- lapply(seq_along(object$tccb), function(w) {
                  x <- object$tccb[[w]]
                  if (is.null(x)) return(seq_len(nrep))
                  rowresults <- apply(x[, 2:ncol(x), drop = FALSE], 1, statf)
                  sort(x[which(rowresults > safealign.threshold), 1])})
    goodlengths <- calcflags <- numeric(lnfac)
    for (hh in 1:lnfac) {
       cglength <- length(goodreps[[hh]])
       goodlengths[hh] <- cglength
       if (cglength == 0) {
         calcflags[hh] <- 1
       } else {
         calcflags[hh] <- 0
       }
    }
    remover <- as.logical(calcflags)
    aggreps <- lapply(seq_along(goodreps), function(w) {
                  gr <- goodreps[[w]]
                  tmod <- object$targetmod[w]
                  if (!(is.na(tmod))) gr <- sort(union(gr, tmod)); gr})
  } else {
    remover <- as.logical(numeric(lnfac))
    aggreps <- rep(list(seq_len(nrep)), lnfac)
  }
  logicheck(parallel)
  nfeat <- NULL
  for (ww in 1:lnfac) {
     nfeat <- c(nfeat, ncol(as.matrix(predstor[[1]][[ww]])))
  }
  if (all(remover) == TRUE) {
    stop("Input 'safealign' was TRUE. For the current values of \n 
         'safealign.stat' and 'safealign.threshold', there is no useable \n 
         replication for any component model. As such, feature importance \n 
         could not be calculated for any component model.")
  }
  for (hh in 1:lnfac) {
     if (remover[hh] == TRUE) {
       warning(paste0("Input 'safealign' was TRUE. For the current values \n
                      of 'safealign.stat' and 'safealign.threshold', there \n 
                      is no useable replication for the component model \n
                      with ", nfac[hh], " components. As such, feature \n 
                      importance was not calculated for this component model."))
     }
  }
  if ((type == "conditional") && (any(nfeat == 1L))) {
    warning("Input 'type' was set to 'conditional', but at least one \n
            component model had only one feature. As such, conditional \n
            feature importance could not be calculated for this model.")
  }
  nfacseq <- 1:lnfac
  pnfacseq <- nfacseq[!(remover)]
  repstor <- vector(mode = "list", length = nrep)
  ccluster <- FALSE
  if (parallel == TRUE) {
    if (is.null(cl)) {
      cl <- makeCluster(max(1L, detectCores() - 1L))
      ccluster <- TRUE
    }
    clusterSetRNGStream(cl, iseed = sample.int(.Machine$integer.max, 1)) 
    registerDoParallel(cl)
  }
  ppac <- c("cpfa", "glmnet", "e1071", "randomForest", "nnet", "rda", "xgboost")
  for (i in 1:nrep) {
     opt.model <- opt.model.all[[i]]
     opt.param.i <- opt.param[[i]]
     ytrain <- object$y[object$trainIDs[[i]]]
     ylab <- object$y[object$testIDs[[i]]]
     facstor <- vector(mode = "list", length = length(nfacseq))
     for (w in pnfacseq) {
        if ((safealign == TRUE) && (!(i %in% aggreps[[w]]))) next
        trainweight <- trweights[[i]][[w]]
        C.pred <- as.matrix(predstor[[i]][[w]])
        cnfac <- nfac[w]
        component.order <- seq_len(cnfac)
        if (cnfac > 1L) {
          change <- object$changeorders[[w]]
          change.row <- match(i, change[, 1L])
          if (!(is.na(change.row))) {
            component.order <- as.integer(change[change.row, -1L])
          }
        }
        inverse.order <- order(component.order)
        model.order <- c(inverse.order,
                         if (ncol(C.pred) > cnfac) {
                           seq.int(cnfac + 1L, ncol(C.pred))
                         } else {
                           integer(0L)
                         })
        C.pred.model <- C.pred[, model.order, drop = FALSE]
        trainweight.model <- trainweight[, model.order, drop = FALSE]
        preds <- imphelper(w = w, method = method, opt.model = opt.model,
                           C.pred = C.pred.model, nfac = nfac,
                           threshold = threshold, storrows = storrows,
                           family = family, y = y, ytrain = ytrain, ylab = ylab,
                           opt.param = opt.param.i, 
                           train.weights = trainweight.model)
        ground0 <- cpm.all(x = as.data.frame(preds), y = as.numeric(ylab) - 1)
        ground <- ground0$cpms
        flit <- vector(mode = "list", length = ncol(C.pred))
        for (ss in 1:ncol(C.pred)) {
           cfit <- NULL             
           if (type == "conditional") {
             Xo <- C.pred[, -ss, drop = FALSE]           
             yt <- C.pred[, ss]                          
             if (ncol(Xo) > 0L) {                
               if (conditional.model == "ridge") {
                 Xc   <- sweep(Xo, 2L, colMeans(Xo), "-")
                 ybar <- mean(yt)
                 bhat <- solve(crossprod(Xc) + ridge.lambda * diag(ncol(Xc)),
                               crossprod(Xc, yt - ybar))
                 cfit <- as.numeric(ybar + Xc %*% bhat)
               } else {                                  
                 rf.cond  <- randomForest::randomForest(x = Xo, y = yt, 
                                                        ntree = ntree,
                                                        nodesize = nodesize)
                 cfit <- as.numeric(rf.cond$predicted)
                 if (anyNA(cfit)) cfit[is.na(cfit)] <- mean(yt)
               }
             }
           }
           if (parallel == TRUE) {
             shuflist <- foreach(jj = 1:nshuffles, 
                                .packages = ppac) %dorng% {
                                C.pred.shuf <- C.pred
                                if (type == "conditional") {
                                  if (is.null(cfit)) {                
                                    C.pred.shuf[, ss] <- sample(C.pred[, ss], 
                                                                replace = FALSE)
                                  } else {
                                    C.pred.shuf[, ss] <- cfit +  
                                      sample(C.pred[, ss] - cfit, 
                                             replace = FALSE)
                                  }
                                } else {
                                  C.pred.shuf[, ss] <- sample(C.pred.shuf[, ss], 
                                                              replace = FALSE)
                                }
                                C.pred.shuf.model <- C.pred.shuf[, model.order, 
                                                                 drop = FALSE]
                                preds <- imphelper(w = w, method = method, 
                                                   opt.model = opt.model,
                                                   C.pred = C.pred.shuf.model, 
                                                   nfac = nfac, y = y,
                                                   threshold = threshold, 
                                                   storrows = storrows, 
                                                   family = family, 
                                                   ytrain = ytrain, ylab = ylab, 
                                                   opt.param = opt.param.i, 
                                                   train.weights = 
                                                     trainweight.model)
                                shufper0 <- cpm.all(x = as.data.frame(preds), 
                                                    y = as.numeric(ylab) - 1)
                                shufper <- shufper0$cpms
                                ground - shufper
                              }
           } else {
             shuflist <- vector(mode = "list", length = nshuffles)
             for (jj in 1:nshuffles) {
                C.pred.shuf <- C.pred
                if (type == "conditional") {
                  if (is.null(cfit)) {                
                    C.pred.shuf[, ss] <- sample(C.pred[, ss], replace = FALSE)
                  } else {
                    C.pred.shuf[, ss] <- cfit + sample(C.pred[, ss] - cfit, 
                                                       replace = FALSE)
                  }
                } else {
                  C.pred.shuf[, ss] <- sample(C.pred.shuf[, ss], 
                                              replace = FALSE)
                }
                C.pred.shuf.model <- C.pred.shuf[, model.order, drop = FALSE]
                preds <- imphelper(w = w, method = method, 
                                   opt.model = opt.model,
                                   C.pred = C.pred.shuf.model, 
                                   nfac = nfac, y = y, threshold = threshold, 
                                   storrows = storrows, family = family, 
                                   ytrain = ytrain, ylab = ylab, 
                                   opt.param = opt.param.i, 
                                   train.weights = trainweight.model)
                shufper0 <- cpm.all(x = as.data.frame(preds), 
                                    y = as.numeric(ylab) - 1)
                shufper <- shufper0$cpms
                shuflist[[jj]] <- ground - shufper
             }
           }
           flit[[ss]] <- Reduce("+", shuflist) / length(shuflist)
        }
        facstor[[w]] <- flit
     }
     repstor[[i]] <- facstor
  }
  if ((parallel == TRUE) && (ccluster == TRUE)) {stopCluster(cl)}
  repbymod <- aggreps
  impstats <- vector("list", lnfac)
  for (w in pnfacseq) {
     greps <- intersect(repbymod[[w]], seq_len(nrep))
     greps <- greps[!(vapply(greps, 
                             function(i) is.null(repstor[[i]][[w]]), 
                             logical(1)))]
     if (!(length(greps))) next
     nfeat <- length(repstor[[greps[1L]]][[w]])
     featstats <- vector("list", nfeat)
     for (ss in seq_len(nfeat)) {
        mats <- lapply(greps, function(i) repstor[[i]][[w]][[ss]])
        featstats[[ss]] <- repsums(mats)
     }
     impstats[[w]] <- featstats
  }
  methodlab <- omethods[as.integer(method)]
  flat <- do.call(rbind, lapply(pnfacseq, function(w) {
    fs <- impstats[[w]]; if (is.null(fs)) return(NULL) 
    do.call(rbind, lapply(seq_along(fs), function(ss) {
      m <- fs[[ss]]; nr <- nrow(m$mean); nc <- ncol(m$mean)
      data.frame(nfac = nfac[w], feature = ss, 
                 method = rep(methodlab, times = nc), 
                 metric = rep(colnames(m$mean), each = nr), 
                 mean = as.vector(m$mean), median = as.vector(m$median),
                 sd = as.vector(m$sd), n = m$n, nvalid = as.vector(m$nvalid),
                 row.names = NULL, stringsAsFactors = FALSE)}))}))
  return(flat)
}