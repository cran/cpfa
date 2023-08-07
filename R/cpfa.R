cpfa <- 
  function(x, y, nrep = 5, ratio = 0.8, seeds = NULL,
           type.out = c("measures", "descriptives"), nfac = 1, nfolds = 10, 
           foldid = NULL, prior = NULL, method = c("PLR", "SVM", "RF", "NN"), 
           family = c("binomial", "multinomial"), alpha = NULL, lambda = NULL, 
           cost = NULL, gamma = NULL, ntree = NULL, nodesize = NULL, 
           size = NULL, decay = NULL, parallel = FALSE, cl = NULL, 
           verbose = TRUE, cmode = NULL, ...) 
{
    xdim <- dim(x)
    lxdim <- length(xdim)
    if (!((lxdim == 3L) | (lxdim == 4L))) 
      stop("Input 'x' must be a 3-way or 4-way array.")
    if (any(is.nan(x)) | any(is.infinite(x))) 
      stop("Input 'x' cannot contain NaN or Inf values.")
    if (any(is.na(x))) 
      stop("Input 'x' cannot contain missing values.")
    if (!(is.null(cmode))) {
      if (!(cmode %in% (1:lxdim))) 
        stop("Input 'cmode' must be 1, 2, or 3 (or 4 if 'x' is a 4-way array).")
      modeval <- 1:lxdim
      mode.re <- c(modeval[-cmode], cmode)
      x <- aperm(x, mode.re)
    }
    if (is.null(cmode)) {
      cmode <- lxdim
    }
    if (!is.factor(y)) 
      stop("Input 'y' must be of class 'factor'.")
    if (!(length(y) == xdim[cmode])) 
      stop("Length of 'y' must match number of levels in \n 
           classification mode/dimension of 'x'.")
    if (!(ceiling(nrep) == nrep) | (nrep < 1)) {
      stop("Input 'nrep' must be an integer greater than 0.")
    }
    if (!(is.numeric(ratio))) {
      stop("Input 'ratio' must be of class numeric.")
    }
    if ((ratio > 1) | (ratio < 0)) {
      stop("Input 'ratio' must be a number between 0 and 1.")
    }
    if (is.null(seeds)) {
      seeds <- 1:nrep
    } else {
      if (!(is.numeric(seeds))) {
        stop("Input 'seeds' must be of class 'numeric'.")
      }
      if ((sum(seeds-floor(seeds)==0)) != 0) {
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
    if ((ltype == 0) | (ltype > 1)) {
      warning("Input 'type.out' not specified correctly. Input must be either \n
              'measures' or 'descriptives'. Defaulting to 'descriptives'.")
      type.out <- "descriptives"
    }
    typelow <- tolower(type.out)
    nobs <- dim(x)[lxdim]
    ntrain <- ceiling(nobs*ratio)
    if (!(is.null(foldid))) {
      if (length(foldid) != ntrain) {
        stop("Input 'foldid' must have length equal to ceiling(nobs*ratio) \n
             where 'nobs' is the number of observations in the classification \n
             mode.")
      }
    }
    stor <- array(0, dim = c(length(nfac)*length(method), 11, nrep))
    for (i in 1:nrep) {
       if (verbose == TRUE) {
         cat("nrep =", i, " \n")
       }
       set.seed(seed = seeds[i])
       train.id <- sample.int(nobs, size = ntrain)
       y.train <- y[train.id]
       y.test <- as.numeric(y[-train.id]) - 1
       if (lxdim == 3L) { 
         X.train <- x[,,train.id]
         X.test <- x[,,-train.id]
       } else {
         X.train <- x[,,,train.id]
         X.test <- x[,,,-train.id]
       }
       cpfalist <- tune.cpfa(x = X.train, y = y.train, nfac = nfac, 
                             nfolds = nfolds, foldid = foldid, prior = prior, 
                             method = method, family = family, alpha = alpha, 
                             lambda = lambda, cost = cost, gamma = gamma, 
                             ntree = ntree, nodesize = nodesize, size = size, 
                             decay = decay, parallel = parallel, cl = cl, 
                             verbose = verbose, cmode = cmode)
       yhat <- predict.cpfa(object = cpfalist, newdata = X.test, 
                            type = "response")
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
      for (j in 1:length(dfun)) {
         output[[j]] <- apply(stor, 1:2, 
                      FUN = function(x){return(get(dfun[j])(x, na.rm = TRUE))})
         rownames(output[[j]]) <- rnam
         colnames(output[[j]]) <- cnam
      }
      names(output) <- dfun
      return(output)
    }
}