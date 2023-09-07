cpfa <- 
  function(x, y, nrep = 5, ratio = 0.8, seeds = NULL,
           type.out = c("measures", "descriptives"), nfac = 1, nfolds = 10, 
           foldid = NULL, prior = NULL, model = c("parafac", "parafac2"),
           method = c("PLR", "SVM", "RF", "NN"), 
           family = c("binomial", "multinomial"), alpha = NULL, lambda = NULL, 
           cost = NULL, gamma = NULL, ntree = NULL, nodesize = NULL, 
           size = NULL, decay = NULL, parallel = FALSE, cl = NULL, 
           verbose = TRUE, cmode = NULL, ...) 
{   
    models <- c("parafac", "parafac2")
    model0 <- pmatch(tolower(model), models)
    lmodel <- length(model0)
    if ((lmodel == 0) || (lmodel > 1)) {
      stop("Input 'model' not specified correctly. Input must be one of \n
            either 'parafac' or 'parafac2'.")
    }
    modellow <- tolower(model)
    if (is.array(x) && (model == "parafac")) {
      xdim <- dim(x)                                                            
      lxdim <- length(xdim)
      if (!((lxdim == 3L) || (lxdim == 4L))) 
        stop("Input 'x' must be a 3-way or 4-way array.")
      if (any(is.nan(x)) || any(is.infinite(x))) 
        stop("Input 'x' cannot contain NaN or Inf values.")
      if (any(is.na(x))) 
        stop("Input 'x' cannot contain missing values.")
      if (!(is.null(cmode))) {
        if (!(cmode %in% (1:lxdim))) 
          stop("Input 'cmode' must be 1, 2, or 3 \n 
               (or 4 if 'x' is a 4-way array).")
        modeval <- 1:lxdim
        mode.re <- c(modeval[-cmode], cmode)
        x <- aperm(x, mode.re)
      } else {
        cmode <- lxdim
      }
    } else if (is.array(x) && (model == "parafac2")) {
      xdim <- dim(x)                                                            
      lxdim <- length(xdim)
      if (!((lxdim == 3L) || (lxdim == 4L))) 
        stop("Input 'x' must be a 3-way or 4-way array.")
      if (any(is.nan(x)) || any(is.infinite(x))) 
        stop("Input 'x' cannot contain NaN or Inf values.")
      if (any(is.na(x))) 
        stop("Input 'x' cannot contain missing values.")
      if (!is.null(cmode)) {
        cmode <- lxdim
        warning("Input 'cmode' is ignored when 'model = parafac2'. Last mode \n
                is classification mode by default.")
      } else {
        cmode <- lxdim
      }
      if (lxdim == 3L) {
        storlist <- vector("list", xdim[cmode])
        for (k in 1:xdim[cmode]) {
          storlist[[k]] <- x[, , k]
        }
      } else {
        storlist <- vector("list", xdim[cmode])
        for (k in 1:xdim[cmode]) {
           storlist[[k]] <- x[, , , k]
        }
      }
      x <- storlist
      rm(storlist)
    } else if (is.list(x) && (model == "parafac2")) {
      xdim1 <- dim(x[[1]])
      lxdim <- length(xdim1) + 1L
      if (!((lxdim == 3L) || (lxdim == 4L))) 
        stop("Input 'x' must be a list of matrices or 3-way arrays.")
      if (!(is.null(cmode))) {
        cmode <- lxdim
        warning("Input 'cmode' is ignored if 'model = parafac2'. Last mode \n
                is classification mode by default. First mode is nested \n
                within last mode (i.e., number of levels for first mode can \n
                vary for each level of the last mode).")
      } else {
        cmode <- lxdim
      }
      if (any(as.logical(lapply(x, function(a){return(any(is.nan(a)))}))))
        stop("Input 'x' cannot contain NaN values")
      if (any(as.logical(lapply(x,function(a){return(any(is.infinite(a)))}))))
        stop("Input 'x' cannot contain Inf values")
      if (any(as.logical(lapply(x,function(a){return(any(is.na(a)))}))))
        stop("Input 'x' cannot contain missing values")
      if (lxdim == 3L) {
        xdim <- rep(NA, 3)
        xdim[2] <- xdim1[2]
        xdim[3] <- length(x)
        if (any(unlist(lapply(x, ncol)) != xdim[2])) {
          stop("Input 'x' must be list of matrices with same \n 
                  number of columns.")
        }
      } else {
        xdim <- rep(NA, 4)
        xdim[2] <- xdim1[2]
        xdim[3] <- xdim1[3]
        xdim[4] <- length(x)
        index2 <- seq(2, (3*length(x) - 1), by = 3)
        index3 <- seq(3, (3*length(x)), by = 3)
        if (any(unlist(lapply(x, dim))[index2] != xdim[2]))
          stop("Input 'x' must be list of arrays with same number of columns.")
        if (any(unlist(lapply(x, dim))[index3] != xdim[3]))
          stop("Input 'x' must be list of arrays with same number of slabs.")
      }
    } else if (is.list(x) && (model == "parafac")) {
      stop("Input 'x' must be of class 'array' if 'model = parafac'.")
    } else {
      stop("Input 'x' must be of class 'array' or 'list'.")
    }
    if (!is.factor(y)) 
      stop("Input 'y' must be of class 'factor'.")
    if (!(length(y) == xdim[cmode])) 
      stop("Length of 'y' must match number of levels in \n 
           classification mode/dimension of 'x'.")
    if (!(ceiling(nrep) == nrep) || (nrep < 1)) {
      stop("Input 'nrep' must be an integer greater than 0.")
    }
    if (!(is.numeric(ratio))) {
      stop("Input 'ratio' must be of class numeric.")
    }
    if ((ratio > 1) || (ratio < 0)) {
      stop("Input 'ratio' must be a number between 0 and 1.")
    }
    if (is.null(seeds)) {
      seeds <- 1:nrep
    } else {
      if (!(is.numeric(seeds))) {
        stop("Input 'seeds' must be of class 'numeric'.")
      }
      if ((sum(seeds-floor(seeds) == 0)) != 0) {
        seeds <- round(seeds)
        warning("At least one seed was not an integer. Non-integer seeds were \n
                rounded using 'round()'.")
      }
      if (length(seeds) != nrep) {
        stop("Input 'seeds' must have length equal to input 'nrep'.")
      }
      if (length(unique(seeds)) != length(seeds)) {
        warning("Not all seeds are unique.")
      } 
    }
    types <- c("measures", "descriptives")
    type <- pmatch(tolower(type.out), types)
    ltype <- length(type)
    if ((ltype == 0) || (ltype > 1)) {
      warning("Input 'type.out' not specified correctly. Input must be either \n
              'measures' or 'descriptives'. Defaulting to 'descriptives'.")
      type.out <- "descriptives"
    }
    typelow <- tolower(type.out)
    if (model == "parafac") {
      nobs <- dim(x)[lxdim]
    } else {
      nobs <- length(x)
    }
    ntrain <- ceiling(nobs*ratio)
    if (!(is.null(foldid))) {
      if (length(foldid) != ntrain) {
        stop("Input 'foldid' must have length equal to ceiling(nobs*ratio) \n
             where 'nobs' is the number of observations in the classification \n
             mode.")
      }
    }
    stor <- array(0, dim = c(length(nfac)*length(method), 11, nrep))
    if (cmode <- lxdim) {
      cmode <- NULL
    }
    for (i in 1:nrep) {
       if (verbose == TRUE) {
         cat("nrep =", i, " \n")
       }
       set.seed(seed = seeds[i])
       train.id <- sample.int(nobs, size = ntrain)
       y.train <- y[train.id]
       y.test <- as.numeric(y[-train.id]) - 1 
       if (model == "parafac") {
         if (lxdim == 3L) { 
           X.train <- x[,,train.id]                                               
           X.test <- x[,,-train.id]
         } else {
           X.train <- x[,,,train.id]
           X.test <- x[,,, -train.id]
         }
       } else {
         X.train <- x[train.id]
         X.test <- x[-train.id]
       }
       cpfalist <- tune.cpfa(x = X.train, y = y.train, nfac = nfac,            
                             nfolds = nfolds, foldid = foldid, prior = prior, 
                             model = model, method = method, family = family, 
                             alpha = alpha, lambda = lambda, cost = cost, 
                             gamma = gamma, ntree = ntree, nodesize = nodesize, 
                             size = size, decay = decay, parallel = parallel, 
                             cl = cl, verbose = verbose, cmode = cmode, ...)    
       yhat <- predict(object = cpfalist, newdata = X.test, type = "response")                        
       out <- cpm.all(x = yhat, y = y.test)
       stor[,,i] <- as.matrix(out$cpms)
    }
    rnam <- rownames(out$cpms)
    cnam <- colnames(out$cpms)
    dimnames(stor)[[1]] <- rnam
    dimnames(stor)[[2]] <- cnam
    if (typelow == "measures") {
      return(stor)                                                              
    } else {
      dfun <- c("mean", "median", "sd")
      output <- vector(mode = "list", length = length(dfun))
      for (j in seq_along(dfun)) {
         output[[j]] <- apply(stor, 1:2, 
                      FUN = function(x){return(get(dfun[j])(x, na.rm = TRUE))})
         rownames(output[[j]]) <- rnam
         colnames(output[[j]]) <- cnam
      }
      names(output) <- dfun                                                     
      return(output)
    }
}