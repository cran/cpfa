predict.cpfa <-
  function(object, newdata = NULL, nfac = NULL, method = NULL,
           type = c("response", "prob", "classify.weights"),
           threshold = NULL, ...) 
{   
    if (class(object) != "cpfa") {
      stop("Input 'object' must be of class 'cpfa'.")
    }
    xold <- x <- object$x
    if (is.null(newdata)) {
      newdata <- object$x
    }
    xold.dim <- dim(xold)
    cmode <- object$cmode
    xdim <- dim(newdata)
    lxdim <- length(xdim)
    if (cmode == lxdim) {
      noperm <- TRUE
    }
    if (cmode != lxdim) {
      noperm <- FALSE
    }
    if (!((lxdim == 3L) | (lxdim == 4L)))
      stop("Input 'newdata' must be a 3-way or 4-way array.")
    if (noperm == TRUE) {
      if (xdim[1] != xold.dim[1])
        stop("Number of levels for A mode of input 'newdata' must 
             match number of levels for A mode used in 'object'.")
      if (xdim[2] != xold.dim[2])
        stop("Number of levels for B mode of input 'newdata' must 
             match number of levels for B mode used in 'object'.")
      if (lxdim == 4L) {
        if (xdim[3] != xold.dim[3])
          stop("Number of levels for C mode of input 'newdata' must 
             match number of levels for C mode used in 'object'.")
      }
    }
    if (noperm == FALSE) {
      modeval <- 1:lxdim
      mode.re <- modeval[-cmode]
      for (w in 1:length(mode.re)) {
         if (xdim[mode.re[w]] != xold.dim[mode.re[w]])
           stop("Number of levels for each mode of input 'newdata' must 
                match number of levels for each mode used in 'object', excluding
                the classification mode.")
      }
      newmode <- c(mode.re, cmode)
      newdata <- aperm(newdata, newmode)
    }
    if (any(is.nan(newdata)) | any(is.infinite(newdata))) 
      stop("Input 'newdata' cannot contain NaN or Inf values.")
    if (any(is.na(newdata)))
      stop("Input 'newdata' cannot contain missing values.")
    if (is.null(nfac)) {
      nfac <- object$opt.param$nfac
    }
    opt.param <- object$opt.param
    lnfac <- length(nfac)
    nfac.names <- paste("fac.", nfac, sep ="")
    methods <- c("PLR", "SVM", "RF", "NN")
    method <- pmatch(toupper(method), methods)
    b.lmethod <- length(method)
    if (is.null(method) | b.lmethod == 0) {
      method <- object$method
    }
    lmethod <- length(method)
    meth.names <- NULL
    if (1 %in% method) {meth.names <- c(meth.names, "PLR")}
    if (2 %in% method) {meth.names <- c(meth.names, "SVM")}
    if (3 %in% method) {meth.names <- c(meth.names, "RF")}
    if (4 %in% method) {meth.names <- c(meth.names, "NN")}
    meth.names <- tolower(meth.names)
    opt.model <- object$opt.model
    Aweights <- object$Aweights
    Bweights <- object$Bweights
    Cweights <- object$Cweights
    const <- object$const
    storfac <- vector("list", lnfac)
    names(storfac) <- nfac.names
    classify.weights <- vector("list", lnfac)
    names(classify.weights) <- nfac.names
    types <- c("response", "prob", "classify.weights")
    type.num <- pmatch(toupper(type), types)
    if (length(type.num) == 0) {
      type <- c("response")
    }
    if(length(type.num) >= 2) {
      warning("Input 'type' either not specified or contains more than one type,
              defaulting to 'type = response'.")
      type <- c("response")
    }
    family <- object$family
    if (family == "binomial") {
      if (!(is.null(threshold))) {
        if ((threshold < 0) | (threshold > 1) | (length(threshold) > 1))
          stop("Input 'threshold' must be a single real number from 0 to 1, 
              inclusive, for binary classification.")
      }
      if ((is.null(threshold)) && (type == "response")) {
        yt <- object$y
        fraction <- table(yt)/length(yt)
        threshold <- fraction[which(names(fraction) == "0")]
      }
      if ((is.null(threshold)) && (type == "prob")) {
        threshold <- 0.5
      }
    }
    if (family == "multinomial") {
      if (!(is.null(threshold))) {
        if (any(threshold < 0) | any(threshold > 1))
          stop("Input 'threshold' must contain real numbers from 0 to 1 for
             multiclass classification.")
        if (sum(threshold) != 1)
          stop("Input 'threshold' must sum to 1 for multiclass classification.")
        warning("Argument 'threshold' is not currently implemented for 
                multiclass classification.")
      }
      if ((is.null(threshold)) & (type == "response")) {
        yt <- object$y
        fraction <- table(yt)/length(yt)
        threshold <- as.numeric(fraction)
      }
      if ((is.null(threshold)) & (type == "prob")) {
        threshold <- 0.5
      }
    }
    stor.name <- NULL
    for (i in 1:lnfac) {
       stor.name <- c(stor.name, paste0(nfac.names[i], meth.names)) 
    }
    storfac <- matrix(NA, nrow = dim(newdata)[lxdim], ncol = lmethod*lnfac)
    storprob <- vector("list", lmethod*lnfac)
    for (w in 1:lnfac) {
       colcount <- lmethod*(w - 1) + 1
       Afixed <- Aweights[[w]]
       Bfixed <- Bweights[[w]]
       if (lxdim == 3L) {
         ppfac <- parafac(X = newdata, nfac = nfac[w], nstart = 1, 
                          ctol = sqrt(.Machine$double.eps), 
                          verbose = FALSE, const = const,
                          Afixed = Afixed, Bfixed = Bfixed)
         classify.weights[[w]] <- C.pred <- ppfac$C
       }
       if (lxdim == 4L) {
         Cfixed <- Cweights[[w]]
         ppfac <- parafac(X = newdata, nfac = nfac[w], nstart = 1, 
                          ctol = sqrt(.Machine$double.eps), 
                          verbose = FALSE, const = const,
                          Afixed = Afixed, Bfixed = Bfixed,
                          Cfixed = Cfixed)
         classify.weights[[w]] <- C.pred <- ppfac$D
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
             if (dim(C.pred)[2] == 1) {
               C.pred.plr <- cbind(0, C.pred)
             } 
             if (dim(C.pred)[2] > 1) {
               C.pred.plr <- C.pred
             }
             if (type == "response") {
               type.plr <- "response" 
               if (family == "binomial") {
                 vals <- as.numeric(predict(plr.fit, newx = C.pred.plr, 
                                            type = type.plr, 
                                            s = lambda.min,
                                            levels = plr.fit.class))
                 storfac[, colcount] <- as.numeric(vals > threshold)
                 colcount <- colcount + 1
               }
               if (family == "multinomial") {
                 vals <- predict(plr.fit, newx = C.pred.plr, 
                                            type = type.plr, 
                                            s = lambda.min,
                                            levels = plr.fit.class)
                 storfac[, colcount] <- as.numeric((apply(vals, 
                                                    1, which.max))) - 1
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
                                  probability = TRUE),"probabilities")
                 svm.prob <- svm.prob[, which(colnames(svm.prob) == "1")]
                 storfac[, colcount] <- as.numeric(svm.prob > threshold)
                 colcount <- colcount + 1
               }
               if (family == "multinomial") {
                 svm.prob <- attr(predict(svm.fit, C.pred, type = type,
                                  probability = TRUE),"probabilities")
                 storfac[, colcount] <- as.numeric((apply(svm.prob, 
                                                          1, which.max))) - 1
                 colcount <- colcount + 1
               }
             }
             if (type == "prob") {
               type.svm <- "response"
               if (family == "binomial") {
                 svm.prob <- attr(predict(svm.fit, C.pred, type = type.svm,
                                  probability = TRUE),"probabilities")
                 svm.prob <- svm.prob[, which(colnames(svm.prob) == "1")]
                 storprob[[colcount]] <- as.numeric(svm.prob)
                 colcount <- colcount + 1
               }
               if (family == "multinomial") {
                 storprob[[colcount]] <- attr(predict(svm.fit, C.pred, 
                                            type = type.svm,
                                            probability = TRUE),"probabilities")
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
                 storfac[, colcount] <- as.numeric((apply(rf.prob, 
                                                          1, which.max))) - 1
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
             }
             if (family == "multinomial") {
               nn.prob <- predict(nn.fit, newdata = C.pred, type = "raw")
               storfac[, colcount] <- as.numeric((apply(nn.prob, 
                                                        1, which.max))) - 1
             }
           }
           if (type == "prob") {
             if (family == "binomial") {
               nn.prob <- predict(nn.fit, newdata = C.pred, type = "raw")
               nn.prob <- nn.prob[, which(colnames(nn.prob) =="y1")]
               storprob[[colcount]] <- as.numeric(nn.prob)
             }
             if (family == "multinomial") {
               storprob[[colcount]]  <- predict(nn.fit, C.pred, type = "raw")
             }
           }
         }
       }
    }
    storfac <- as.data.frame(storfac)
    colnames(storfac) <- stor.name
    names(storprob) <- stor.name
    if (type == "classify.weights") {
      classify.weight.names <- paste(nfac, "-factor(s)", sep ="")
      names(classify.weights) <- classify.weight.names
      return(classify.weights)
    }
    if (type == "response") {
      return(storfac)
    }
    if (type == "prob") {
      return(storprob)
    }
}