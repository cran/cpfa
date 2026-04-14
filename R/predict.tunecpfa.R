predict.tunecpfa <-
  function(object, newdata = NULL, newdata.z = NULL, method = NULL,
           type = c("response", "prob", "classify.weights"), 
           threshold = NULL, ...)                                               
{   
    if (!(inherits(object, "tunecpfa"))) {
      stop("Input 'object' must be of class 'tunecpfa'.")
    }
    model <- object$model
    xold.dim <- object$xdim
    cmode <- object$cmode
    oyold <- object$y
    zout <- object$z
    if (is.null(newdata)) {
      newdata <- object$x
      xdim <- object$xdim
      lxdim <- object$lxdim
      skippy <- TRUE
    } else {
      skippy <- FALSE
    }
    if (!(skippy)) {
      if (is.array(newdata) && (!is.matrix(newdata)) && 
          (model %in% c("parafac", "pca"))) {
        xdim <- dim(newdata)                                                            
        lxdim <- length(xdim)
        if (!((lxdim == 3L) || (lxdim == 4L))) {
          stop("Input 'newdata' must be a 3-way or 4-way array.")
        }
        if (any(is.nan(newdata)) || any(is.infinite(newdata))) {
          stop("Input 'newdata' cannot contain NaN or Inf values.")
        }
        if (any(is.na(newdata))) {
          stop("Input 'newdata' cannot contain missing values.")
        }
        if (model == "parafac") {
          if (cmode != length(xold.dim)) {
            modeval <- 1:lxdim
            mode.re <- c(modeval[-cmode], cmode)
            newdata <- aperm(newdata, mode.re)
            xdim <- dim(newdata)
          }
        } else {
          modeval <- 1:lxdim
          mode.re <- c(cmode, modeval[-cmode])
          newdata <- matrix(aperm(newdata, mode.re), nrow = xdim[cmode])
          xdim <- dim(newdata)
          lxdim <- length(xdim)
        }
      } else if (is.array(newdata) && (!is.matrix(newdata)) 
                 && (model == "parafac2")) {
        xdim <- dim(newdata)                                                            
        lxdim <- length(xdim)
        if (!((lxdim == 3L) || (lxdim == 4L))) {
          stop("Input 'newdata' must be a 3-way or 4-way array.")
        }
        if (any(is.nan(newdata)) || any(is.infinite(newdata))) {
          stop("Input 'newdata' cannot contain NaN or Inf values.")
        }
        if (any(is.na(newdata))) {
          stop("Input 'newdata' cannot contain missing values.")
        }
        if (lxdim == 3L) {
          storlist <- vector("list", xdim[3])
          for (k in 1:xdim[3]) {storlist[[k]] <- newdata[, , k]}
        } else {
          storlist <- vector("list", xdim[4])
          for (k in 1:xdim[4]) {storlist[[k]] <- newdata[, , , k]}
        }
        newdata <- storlist
        rm(storlist)
      } else if (is.list(newdata) && (model == "parafac2")) {
        xdim1 <- dim(newdata[[1]])
        lxdim <- length(xdim1) + 1L
        if (!((lxdim == 3L) || (lxdim == 4L))) {
          stop("Input 'newdata' must be a list of matrices or 3-way arrays.")
        }
        if (any(as.logical(lapply(newdata, 
                                  function(a){return(any(is.nan(a)))})))) {
          stop("Input 'newdata' cannot contain NaN values")
        }
        if (any(as.logical(lapply(newdata,
                                  function(a){return(any(is.infinite(a)))})))) {
          stop("Input 'newdata' cannot contain Inf values")
        }
        if (any(as.logical(lapply(newdata, 
                                  function(a){return(any(is.na(a)))})))) {
          stop("Input 'newdata' cannot contain missing values")
        }
        if (lxdim == 3L) {
          xdim <- rep(NA, 3)
          xdim[2] <- xdim1[2]
          xdim[3] <- length(newdata)
          if (any(unlist(lapply(newdata, ncol)) != xdim[2])) {
            stop("Input 'newdata' must be list of matrices with same number \n
                 of columns.")
          }
        } else {
          xdim <- rep(NA, 4)
          xdim[2] <- xdim1[2]
          xdim[3] <- xdim1[3]
          xdim[4] <- length(newdata)
          index2 <- seq(2, (3 * length(newdata) - 1), by = 3)
          index3 <- seq(3, (3 * length(newdata)), by = 3)
          if (any(unlist(lapply(newdata, dim))[index2] != xdim[2])) {
            stop("Input 'newdata' must be list of arrays with same \n
                 number of columns.")
          }
          if (any(unlist(lapply(newdata, dim))[index3] != xdim[3])) {
            stop("Input 'newdata' must be list of arrays with same \n
                 number of slabs.")
          }
        }
      } else if (is.list(newdata) && (model != "parafac2")) {
        stop("Input 'newdata' cannot be of class 'list' unless \n
             'model = parafac2'.")
      } else if (is.matrix(newdata) && (model != "pca")) {
        stop("Input 'newdata' cannot be of class 'matrix' if model is \n
             'parafac' or 'parafac2'.")
      } else if (is.matrix(newdata) && (model == "pca")) {  
        xdim <- origdim <- dim(newdata)                                                            
        lxdim <- length(xdim)
        if (!((lxdim == 2L))) {
          stop("If 'newdata' has class 'matrix' and model is 'pca', 'newdata' \n 
               must be a 2-way matrix.")
        }
        if (any(is.nan(newdata)) || any(is.infinite(newdata))) {
          stop("Input 'newdata' cannot contain NaN or Inf values.")
        }
        if (any(is.na(newdata))) {
          stop("Input 'newdata' cannot contain missing values.")
        }
        if (!(is.null(cmode))) {
          if (!(cmode %in% (1:lxdim))) {
            stop("Input 'cmode' must be 1 or 2 when model is 'pca' and when \n
                 'newdata' is of class 'matrix'.")
          }
          if (cmode == 1) {
            newdata <- newdata
          } else {
            modeval <- 1:lxdim
            mode.re <- c(modeval[-cmode], cmode)
            newdata <- aperm(newdata, mode.re)
            xdim <- dim(newdata)
          }
        } else {
          cmode <- 1
        }
      } else {
        stop("Input 'newdata' can be of class 'array' for any model, can be \n
             of class 'list' for 'parafac2', and can be of class 'matrix' for \n
             'pca'. Else, input 'newdata' cannot be of a different class.")
      }
    }
    if (model == "parafac") {
      if (xdim[1] != xold.dim[1]) {
        stop("Number of levels for A mode of input 'newdata' must match \n
             number of levels for A mode used in 'object'.")
      }
      if (xdim[2] != xold.dim[2]) {
        stop("Number of levels for B mode of input 'newdata' must match \n
             number of levels for B mode used in 'object'.")
      }
      if (lxdim == 4L) {
        if (xdim[3] != xold.dim[3]) {
          stop("Number of levels for C mode of input 'newdata' must match \n 
               number of levels for C mode used in 'object'.")
        }
      }
    } else if (model == "parafac2") {
      if (xdim[2] != xold.dim[2]) {
        stop("Number of levels for B mode of input 'newdata' must match \n 
             number of levels for B mode used in 'object'.")
      }
      if (lxdim == 4L) {
        if (xdim[3] != xold.dim[3]) {
          stop("Number of levels for C mode of input 'newdata' must match \n 
               number of levels for C mode used in 'object'.")
        }
      }
    } else {
      if (xdim[2] != xold.dim[2]) {
        stop("Number of variables of input 'newdata' must match number of \n 
             variables used in 'object'.")
      }
    }
    if ((!(is.null(newdata.z))) && ((is.null(zout)))) {
      stop("Input 'newdata.z' was provided, but input 'z' from the \n 
           'tunecpfa' object is NULL. When 'newdata.z' is provided, \n
           input 'z' from 'tunecpfa' cannot be NULL.")
    } else if ((is.null(newdata.z)) && (!(is.null(zout)))) {
      stop("Input 'newdata.z' was not provided, but input 'z' from the \n
           'tunecpfa' object is not NULL. When 'z' from 'tunecpfa' is not \n
           NULL, 'newdata.z' must be provided.")
    } else if ((is.null(newdata.z)) && (is.null(zout))) {
      newdata.z <- NULL
    } else {
      newzdim <- dim(newdata.z)                                                            
      nlzdim <- length(newzdim)
      if (!((nlzdim == 2L))) {stop("Input 'newdata.z' must have two modes.")}
      if (!(is.matrix(newdata.z))) {
        stop("Input 'newdata.z', when provided, must be of class 'matrix'.")
      }
      if (any(is.nan(newdata.z)) || any(is.infinite(newdata.z))) {
        stop("Input 'newdata.z' cannot contain NaN or Inf values.")
      }
      if (any(is.na(newdata.z))) {
        stop("Input 'newdata.z' cannot contain missing values.")
      }
      if (newzdim[2] != dim(zout)[2]) {
        stop("Input 'newdata.z', when provided, must have the same number of \n
             number of features as the number of features used in the \n
             original input 'z' from 'tunecpfa'.")
      }
    }
    nfac <- object$opt.param$nfac
    opt.param <- object$opt.param
    lnfac <- length(nfac)
    nfac.names <- paste("fac.", nfac, sep = "")
    if (is.null(method)) {
      method <- object$method
      lmethod <- length(method)
    } else {
      omethods <- c("PLR", "SVM", "RF", "NN", "RDA", "GBM")
      checkmethod <- sum(toupper(method) %in% omethods)
      if (checkmethod != length(method)) {
        stop("Input 'method' contains at least one value that is not valid.")
      }
      method <- which(omethods %in% toupper(method) == TRUE)
      lmethod <- length(method)
      if (sum(method %in% object$method) != lmethod) {
        stop("Input 'method' contains at least one method that was not tuned \n
             in input 'object'.")
      }
    }
    meth.names <- NULL
    if (1 %in% method) {meth.names <- c(meth.names, "PLR")}
    if (2 %in% method) {meth.names <- c(meth.names, "SVM")}
    if (3 %in% method) {meth.names <- c(meth.names, "RF")}
    if (4 %in% method) {meth.names <- c(meth.names, "NN")}
    if (5 %in% method) {
      meth.names <- c(meth.names, "RDA")
      train.weights <- object$train.weights
    }
    if (6 %in% method) {meth.names <- c(meth.names, "GBM")}
    meth.names <- tolower(meth.names)
    opt.model <- object$opt.model
    Aweights <- object$Aweights
    Bweights <- object$Bweights   
    Cweights <- object$Cweights
    Phi <- object$Phi
    scenters <- object$scenters
    sscales <- object$sscales
    const <- object$const
    classify.weights <- vector("list", lnfac)
    names(classify.weights) <- nfac.names
    if (length(type) != 1) {
      stop("Input 'type' must be a character of length equal to 1. Make sure \n
           to specify 'type'.")
    }
    types <- c("response", "prob", "classify.weights")
    if (!(tolower(type) %in% types)) {
      stop("Input 'type' must be 'response', 'prob', or 'classify.weights'.")
    }
    type <- tolower(type)
    family <- object$family
    if (family == "binomial") {
      if (!(is.null(threshold))) {
        if ((threshold < 0) || (threshold > 1) || (length(threshold) > 1)) {
          stop("Input 'threshold' must be a single real number from 0 to 1, \n
              inclusive, for binary classification.")
        }
      }
      if ((is.null(threshold)) && (type != "classify.weights")) {
        threshold <- 0.5
      }
    }
    if (family == "multinomial") {
      if (!(is.null(threshold))) {
        if (any(threshold < 0) || any(threshold > 1)) {
          stop("Input 'threshold' must contain real numbers from 0 to 1 for
             multiclass classification.")
        }
        if (sum(threshold) != 1) {
          stop("Input 'threshold' must sum to 1 for multiclass classification.")
        }
        warning("Argument 'threshold' is not currently implemented for 
                multiclass classification.")
      }
      if ((is.null(threshold)) && (type == "response")) {
        yt <- object$y
        fraction <- table(yt) / length(yt)
        threshold <- as.numeric(fraction)
      }
      if ((is.null(threshold)) && (type == "prob")) {
        threshold <- 0.5
      }
    }
    stor.name <- NULL
    for (i in 1:lnfac) {
       stor.name <- c(stor.name, paste0(nfac.names[i], meth.names)) 
    }
    if (model == "parafac") {
      storrows <- dim(newdata)[lxdim]
    } else if (model == "parafac2") {
      storrows <- length(newdata)
    } else {
      storrows <- xdim[1]
    }
    storfac <- matrix(NA, nrow = storrows, ncol = lmethod * lnfac)
    storprob <- vector("list", lmethod * lnfac)
    for (w in 1:lnfac) {
       colcount <- lmethod * (w - 1) + 1
       if (model == "parafac") {
         Afixed <- Aweights[[w]]
         Bfixed <- Bweights[[w]]
         if (lxdim == 3L) {
           ppfac <- parafac(X = newdata, nfac = nfac[w], nstart = 1, 
                            ctol = sqrt(.Machine$double.eps), verbose = FALSE, 
                            const = const, Afixed = Afixed, Bfixed = Bfixed)                
           C.pred0 <- ppfac$C
           if (!(is.null(newdata.z))) {
             if (nrow(newdata.z) != nrow(C.pred0)) {
               stop("Input 'newdata.z', when provided, must have number of \n
                    observations equal to the number of levels of the \n
                    classification mode of input 'newdata'.")
             } else {
               C.pred0 <- cbind(C.pred0, newdata.z)
             }
           }
           C.pred <- scale(C.pred0, center = scenters[[w]], 
                           scale = sscales[[w]])
           classify.weights[[w]] <- C.pred
         }
         if (lxdim == 4L) {
           Cfixed <- Cweights[[w]]
           ppfac <- parafac(X = newdata, nfac = nfac[w], nstart = 1, 
                            ctol = sqrt(.Machine$double.eps), verbose = FALSE, 
                            const = const, Afixed = Afixed, Bfixed = Bfixed,
                            Cfixed = Cfixed)                                
           C.pred0 <- ppfac$D
           if (!(is.null(newdata.z))) {
             if (nrow(newdata.z) != nrow(C.pred0)) {
               stop("Input 'newdata.z', when provided, must have number of \n
                    observations equal to the number of levels of the \n
                    classification mode of input 'newdata'.")
             } else {
               C.pred0 <- cbind(C.pred0, newdata.z)
             }
           }
           C.pred <- scale(C.pred0, center = scenters[[w]], 
                           scale = sscales[[w]])
           classify.weights[[w]] <- C.pred
         }
       } else if (model == "parafac2") {
         Phifixed <- Phi[[w]]
         phi.pfac2 <- eigen(Phifixed, symmetric = TRUE)
         if (nfac[w] == 1) {
           Gfixed <- t(phi.pfac2$vectors * sqrt(phi.pfac2$values))
         } else {
           Gfixed <- t(phi.pfac2$vectors %*% diag(sqrt(phi.pfac2$values)))
         }
         Bfixed <- Bweights[[w]]
         if (lxdim == 3L) {
           ppfac <- parafac2(X = newdata, nfac = nfac[w], nstart = 1, 
                             ctol = sqrt(.Machine$double.eps), verbose = FALSE, 
                             const = const, Gfixed = Gfixed, Bfixed = Bfixed)                
           C.pred0 <- ppfac$C
           if (!(is.null(newdata.z))) {
             if (nrow(newdata.z) != nrow(C.pred0)) {
               stop("Input 'newdata.z', when provided, must have number of \n
                    observations equal to the number of levels of the \n
                    classification mode of input 'newdata'.")
             } else {
               C.pred0 <- cbind(C.pred0, newdata.z)
             }
           }
           C.pred <- scale(C.pred0, center = scenters[[w]], 
                           scale = sscales[[w]])
           classify.weights[[w]] <- C.pred
         }
         if (lxdim == 4L) {
           Cfixed <- Cweights[[w]]
           ppfac <- parafac2(X = newdata, nfac = nfac[w], nstart = 1, 
                             ctol = sqrt(.Machine$double.eps), verbose = FALSE, 
                             const = const, Gfixed = Gfixed, Bfixed = Bfixed,
                             Cfixed = Cfixed)                                
           C.pred0 <- ppfac$D
           if (!(is.null(newdata.z))) {
             if (nrow(newdata.z) != nrow(C.pred0)) {
               stop("Input 'newdata.z', when provided, must have number of \n
                    observations equal to the number of levels of the \n
                    classification mode of input 'newdata'.")
             } else {
               C.pred0 <- cbind(C.pred0, newdata.z)
             }
           }
           C.pred <- scale(C.pred0, center = scenters[[w]], 
                           scale = sscales[[w]])
           classify.weights[[w]] <- C.pred
         }
       } else {
         pcacenter <- object$pcacenter
         Afixed <- Aweights[[w]]
         ndcent <- scale(newdata, center = pcacenter, scale = FALSE)
         C.pred0 <- ndcent %*% Afixed
         if (!(is.null(newdata.z))) {
           if (nrow(newdata.z) != nrow(C.pred0)) {
             stop("Input 'newdata.z', when provided, must have number of \n
                    observations equal to the number of levels of the \n
                    classification mode of input 'newdata'.")
           } else {
             C.pred0 <- cbind(C.pred0, newdata.z)
           }
         }
         C.pred <- scale(C.pred0, center = scenters[[w]], 
                         scale = sscales[[w]])
         classify.weights[[w]] <- C.pred
       }
       if (type != "classify.weights") {
         C.pred <- as.matrix(C.pred)
           if ('1' %in% method) {
             plr.fit <- opt.model[[w]][[1]]
             if (!(is.null(plr.fit$glmnet.fit))) {
               plr.fit.class <- plr.fit$glmnet.fit$classnames
               lambda.min <- plr.fit$lambda.min
             }
             if (is.null(plr.fit$glmnet.fit)) {
               plr.fit.class <- plr.fit$classnames
               lambda.min <- opt.param[which(
                                       opt.param$nfac == nfac[w]), ]$lambda
             }
             if ((nfac[w] == 1) || (nfac[w] == 1L)) {
               C.pred.plr <- cbind(C.pred, 0)
             } else {
               C.pred.plr <- C.pred
             }
             if (type == "response") {
               type.plr <- "response" 
               if (family == "binomial") {
                 vals <- as.numeric(predict(plr.fit, newx = C.pred.plr, 
                                            type = type.plr, s = lambda.min,
                                            levels = plr.fit.class))
                 storfac[, colcount] <- as.numeric(vals > threshold)
                 colcount <- colcount + 1
               }
               if (family == "multinomial") {
                 vals <- predict(plr.fit, newx = C.pred.plr, type = type.plr, 
                                 s = lambda.min, levels = plr.fit.class)
                 storfac[, colcount] <- as.numeric((apply(vals, 1, 
                                                          which.max))) - 1
                 colcount <- colcount + 1
               }
             }
             if (type == "prob") {
               type.plr <- "response"
               if (family == "binomial") {
                 storprob[[colcount]] <- as.numeric(predict(plr.fit, 
                                                    newx = C.pred.plr, 
                                                    type = type.plr, 
                                                    s = lambda.min,
                                                    levels = plr.fit.class))
                 colcount <- colcount + 1
               }
               if (family == "multinomial") {
                 plr.vals <- as.numeric(predict(plr.fit, 
                                                newx = C.pred.plr, 
                                                type = type.plr, 
                                                s = lambda.min,
                                                levels = plr.fit.class))
                 plr.dim <- dim(C.pred.plr)
                 lplr.class <- length(plr.fit.class)
                 storprob[[colcount]] <- matrix(plr.vals, nrow = plr.dim[1], 
                                                ncol = lplr.class)
                 colcount <- colcount + 1
               }
             }
           }
           if ('2' %in% method) {
             svm.fit <- opt.model[[w]][[2]]
             if (type == "response") {
               if (family == "binomial") {
                 svm.prob <- attr(predict(svm.fit, C.pred, type = type,
                                  probability = TRUE), "probabilities")
                 svm.prob <- svm.prob[, which(colnames(svm.prob) == "1")]
                 storfac[, colcount] <- as.numeric(svm.prob > threshold)
                 colcount <- colcount + 1
               }
               if (family == "multinomial") {
                 svm.prob <- attr(predict(svm.fit, C.pred, type = type,
                                  probability = TRUE), "probabilities")
                 sord <- cbind(1:ncol(svm.prob), 
                               as.numeric(colnames(svm.prob)) + 1)
                 svmord <- sord[order(sord[,2]), ]
                 svm.prob <- svm.prob[, svmord[, 1]]
                 storfac[, colcount] <- as.numeric((apply(svm.prob, 1, 
                                                          which.max))) - 1
                 colcount <- colcount + 1
               }
             }
             if (type == "prob") {
               type.svm <- "response"
               if (family == "binomial") {
                 svm.prob <- attr(predict(svm.fit, C.pred, type = type.svm,
                                  probability = TRUE), "probabilities")
                 svm.prob <- svm.prob[, which(colnames(svm.prob) == "1")]
                 storprob[[colcount]] <- as.numeric(svm.prob)
                 colcount <- colcount + 1
               }
               if (family == "multinomial") {
                 storprob[[colcount]] <- attr(predict(svm.fit, C.pred, 
                                              type = type.svm,
                                              probability = TRUE),
                                              "probabilities")
                 colcount <- colcount + 1
               }
             }
           }
           if ('3' %in% method) {
             rf.fit <- opt.model[[w]][[3]]
             if (type == "response") {
               if (family == "binomial") {
                 rf.prob <- predict(rf.fit, C.pred, type = "prob")
                 rf.prob <- rf.prob[, which(colnames(rf.prob) == "1")]
                 storfac[, colcount] <- as.numeric(rf.prob > threshold)
                 colcount <- colcount + 1
               }
               if (family == "multinomial") {
                 rf.prob <- predict(rf.fit, C.pred, type = "prob")
                 storfac[, colcount] <- as.numeric((apply(rf.prob, 1, 
                                                          which.max))) - 1
                 colcount <- colcount + 1
               }
             }
             if (type == "prob") {
               if (family == "binomial") {
                 rf.prob <- predict(rf.fit, C.pred, type = "prob")
                 rf.prob <- rf.prob[, which(colnames(rf.prob) =="1")]
                 storprob[[colcount]] <- as.numeric(rf.prob)
                 colcount <- colcount + 1
               }
               if (family == "multinomial") {
                 storprob[[colcount]]  <- predict(rf.fit, C.pred, type = "prob")
                 colcount <- colcount + 1
               }
             }
           }
           if ('4' %in% method) {
             nn.fit <- opt.model[[w]][[4]]
             if (type == "response") {
               if (family == "binomial") {
                 nn.prob <- predict(nn.fit, newdata = C.pred, type = "raw")
                 nn.prob <- nn.prob[, which(colnames(nn.prob) == "y1")]
                 storfac[, colcount] <- as.numeric(nn.prob > threshold)
                 colcount <- colcount + 1
               }
               if (family == "multinomial") {
                 nn.prob <- predict(nn.fit, newdata = C.pred, type = "raw")
                 storfac[, colcount] <- as.numeric((apply(nn.prob, 1, 
                                                          which.max))) - 1
                 colcount <- colcount + 1
               }
             }
             if (type == "prob") {
               if (family == "binomial") {
                 nn.prob <- predict(nn.fit, newdata = C.pred, type = "raw")
                 nn.prob <- nn.prob[, which(colnames(nn.prob) =="y1")]
                 storprob[[colcount]] <- as.numeric(nn.prob)
                 colcount <- colcount + 1
               }
               if (family == "multinomial") {
                 storprob[[colcount]]  <- predict(nn.fit, C.pred, type = "raw")
                 colcount <- colcount + 1
               }
             }
           }
           if ('5' %in% method) {
             rda.fit <- opt.model[[w]][[5]]
             if (type == "response") {
               if (family == "binomial") {
                 rda.prob <- predict(rda.fit, x = t(train.weights[[w]]), 
                                     y = as.numeric(oyold) - 1, 
                                     xnew = t(C.pred), type = "posterior")
                 colnames(rda.prob) <- c(0, 1)
                 rda.prob <- rda.prob[, which(colnames(rda.prob) == "1")]
                 storfac[, colcount] <- as.numeric(rda.prob > threshold)
                 colcount <- colcount + 1
               }
               if (family == "multinomial") {
                 rda.prob <- predict(rda.fit, x = t(train.weights[[w]]), 
                                     y = as.numeric(oyold) - 1, 
                                     xnew = t(C.pred), type = "posterior")
                 storfac[, colcount] <- as.numeric((apply(rda.prob, 1, 
                                                          which.max))) - 1
                 colcount <- colcount + 1
               }
             }
             if (type == "prob") {
               if (family == "binomial") {
                 rda.prob <- predict(rda.fit, x = t(train.weights[[w]]), 
                                     y = as.numeric(oyold) - 1, 
                                     xnew = t(C.pred), type = "posterior")
                 colnames(rda.prob) <- c(0, 1)
                 rda.prob <- rda.prob[, which(colnames(rda.prob) == "1")]
                 storprob[[colcount]] <- as.numeric(rda.prob)
                 colcount <- colcount + 1
               }
               if (family == "multinomial") {
                 storprob[[colcount]]  <- predict(rda.fit, 
                                                  x = t(train.weights[[w]]), 
                                                  y = as.numeric(oyold) - 1, 
                                                  xnew = t(C.pred), 
                                                  type = "posterior")
                 colcount <- colcount + 1
               }
             }
           }
           if ('6' %in% method) {
             gbm.fit <- opt.model[[w]][[6]]
             xgdata <- xgb.DMatrix(data = C.pred)
             num_classes <- length(unique(oyold))
             if (type == "response") {
               if (family == "binomial") {
                 gbm.pred <- predict(gbm.fit, xgdata)
                 storfac[, colcount] <- as.numeric(gbm.pred > threshold)       
               }
               if (family == "multinomial") {
                 gbm.pred <- matrix(predict(gbm.fit, xgdata), 
                                    ncol = num_classes)
                 storfac[, colcount] <- as.numeric((apply(gbm.pred, 1, 
                                                          which.max))) - 1
               }
             }
             if (type == "prob") {
               if (family == "binomial") {
                 storprob[[colcount]] <- predict(gbm.fit, xgdata)
               }
               if (family == "multinomial") {
                 storprob[[colcount]] <- matrix(predict(gbm.fit, xgdata), 
                                                ncol = num_classes)
               }
             }
           } 
       }
    }
    storfac <- as.data.frame(storfac)
    colnames(storfac) <- stor.name
    names(storprob) <- stor.name
    if (type == "classify.weights") {
      classify.weight.names <- paste(nfac, "-component(s)", sep ="")
      names(classify.weights) <- classify.weight.names
      return(classify.weights)
    } 
    if (type == "response") {return(storfac)} 
    if (type == "prob") {return(storprob)}
}