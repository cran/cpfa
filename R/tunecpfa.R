tunecpfa <- 
  function(x, y, model = c("parafac", "parafac2"), nfac = 1, nfolds = 10,
           method = c("PLR", "SVM", "RF", "NN", "RDA", "GBM"),
           family = c("binomial", "multinomial"), parameters = list(), 
           foldid = NULL, prior = NULL, cmode = NULL, parallel = FALSE, 
           cl = NULL, verbose = TRUE, ...)
{   
    if (!(is.list(parameters))) {
      stop("Input 'parameters' must be of class 'list'.")
    } 
    if (length(parameters) == 0) {
      size <- nodesize <- ntree <- gamma <- cost <- lambda <- alpha <- NULL 
      subsample <-  max.depth <- eta <- delta <- decay <- rda.alpha <- size
      nrounds <- subsample
    } else {
      names(parameters) <- tolower(names(parameters))
      paranames <- c("alpha", "lambda", "cost", "gamma", "ntree", "nodesize", 
                     "size", "decay", "rda.alpha", "delta", "eta", "max.depth", 
                     "subsample", "nrounds")
      paralogical <- paranames %in% names(parameters)
      if (sum(paralogical) == 0) {
        stop("Input 'parameters' was provided but does not contain any \n
             acceptable values. Input 'parameters' must contain only \n 
             acceptable values. See help file for function 'cpfa' under \n
             argument 'parameters' for a list of acceptable values.")
      }
      if (length(unique(names(parameters))) != length(names(parameters))) {
        stop("Input 'parameters' was provided but contains duplicate \n
             parameter names. Input 'parameters' cannot contain duplicates.")
      }
      logicalpara <- names(parameters) %in% paranames
      if ((sum(logicalpara) < length(parameters)) && (sum(logicalpara) > 0)) {
        stop("Input 'parameters' was provided and contains some acceptable \n
             values. However, 'parameters' also contains one or more values \n
             that are not acceptable. Input 'parameters' must contain only \n 
             acceptable values. See help file for function 'cpfa' under \n
             argument 'parameters' for a list of acceptable values.")
      }
      parainput <- paranames[which(paralogical)]
      paranull <- paranames[-which(paralogical)]
      for (j in 1:length(parainput)) {
         inputparam <- paste0(parainput[j], " = ", "parameters$", parainput[j])
         eval(parse(text = inputparam))
      }
      iparam <- paranull
      if (length(iparam) != 0) {
        for (j in 1:length(iparam)) {
           inputparam <- paste0(iparam[j], " = ", "NULL")
           eval(parse(text = inputparam))
        }
      } 
    }
    models <- c("parafac", "parafac2")
    model0 <- sum(tolower(model) %in% models)
    if ((model0 == 0L) || (model0 > 1L)) {
      stop("Input 'model' not specified correctly. Input must be one of \n
           either 'parafac' or 'parafac2'.")
    }
    model <- tolower(model)
    if (is.array(x) && (model == "parafac")) {
      xdim <- origdim <- dim(x)                                                            
      lxdim <- length(xdim)
      if (!((lxdim == 3L) || (lxdim == 4L))) {
        stop("Input 'x' must be a 3-way or 4-way array.")
      }
      if (any(is.nan(x)) || any(is.infinite(x))) {
        stop("Input 'x' cannot contain NaN or Inf values.")
      }
      if (any(is.na(x))) {stop("Input 'x' cannot contain missing values.")}
      if (!(is.null(cmode))) {
        if (!(cmode %in% (1:lxdim))) {
          stop("Input 'cmode' must be 1, 2, or 3 \n 
               (or 4 if 'x' is a 4-way array).")
        }
        modeval <- 1:lxdim
        mode.re <- c(modeval[-cmode], cmode)
        x <- aperm(x, mode.re)
        xdim <- dim(x)
      } else {
        cmode <- lxdim
      }
    } else if (is.array(x) && (model == "parafac2")) {
      xdim <- dim(x)                                                            
      lxdim <- length(xdim)
      if (!((lxdim == 3L) || (lxdim == 4L))) {
        stop("Input 'x' must be a 3-way or 4-way array.")
      }
      if (any(is.nan(x)) || any(is.infinite(x))) {
        stop("Input 'x' cannot contain NaN or Inf values.")
      }
      if (any(is.na(x))) {stop("Input 'x' cannot contain missing values.")}
      if (!is.null(cmode)) {
        cmode <- lxdim
        warning("Input 'cmode' is ignored when 'model = parafac2'. Last mode \n
                is classification mode by default.")
      } else {
        cmode <- lxdim
      }
      if (lxdim == 3L) {
        storlist <- vector("list", xdim[cmode])
        for (k in 1:xdim[cmode]) {storlist[[k]] <- x[, , k]}
      } else {
        storlist <- vector("list", xdim[cmode])
        for (k in 1:xdim[cmode]) {storlist[[k]] <- x[, , , k]}
      }
      x <- storlist
      rm(storlist)
    } else if (is.list(x) && (model == "parafac2")) {
      xdim1 <- dim(x[[1]])
      lxdim <- length(xdim1) + 1L
      if (!((lxdim == 3L) || (lxdim == 4L))) {
        stop("Input 'x' must be a list of matrices or 3-way arrays.")
      }
      if (!(is.null(cmode))) {
        cmode <- lxdim
        warning("Input 'cmode' is ignored if 'model = parafac2'. Last mode \n
                is classification mode by default. First mode is nested \n
                within last mode (i.e., number of levels for first mode can \n
                vary for each level of the last mode).")
      } else {
        cmode <- lxdim
      }
      if (any(as.logical(lapply(x, function(a){return(any(is.nan(a)))})))) {
        stop("Input 'x' cannot contain NaN values")
      }
      if (any(as.logical(lapply(x,function(a){return(any(is.infinite(a)))})))) {
        stop("Input 'x' cannot contain Inf values")
      }
      if (any(as.logical(lapply(x,function(a){return(any(is.na(a)))})))) {
        stop("Input 'x' cannot contain missing values")
      }
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
        index2 <- seq(2, (3 * length(x) - 1), by = 3)
        index3 <- seq(3, (3 * length(x)), by = 3)
        if (any(unlist(lapply(x, dim))[index2] != xdim[2])) {
          stop("Input 'x' must be list of arrays with same number of columns.")
        }
        if (any(unlist(lapply(x, dim))[index3] != xdim[3])) {
          stop("Input 'x' must be list of arrays with same number of slabs.")
        }
      }
    } else if (is.list(x) && (model == "parafac")) {
      stop("Input 'x' must be of class 'array' if 'model = parafac'.")
    } else {
      stop("Input 'x' must be of class 'array' or 'list'.")
    }
    if (!(is.factor(y))) {stop("Input 'y' must be of class 'factor'.")}
    if (model == "parafac") {
      if (!(length(y) == origdim[cmode])) {
        stop("Length of 'y' must match number of levels in classification \n 
             mode of 'x'.")
      }
    } else {
      if (!(length(y) == xdim[cmode])) {
        stop("Length of 'y' must match number of levels in classification \n 
             mode of 'x'.")
      }
    }
    if (length(unique(y)) == 0L || length(unique(y)) == 1L) {
      stop("Input 'y' contains less than two unique labels. Need to provide \n
           at least two unique labels.")
    }
    ylev <- length(levels(y))
    if (!(sum(sort(as.numeric(levels(y))) == 0:(ylev - 1)) == ylev)) {
      stop("Input 'y' must contain labels of 0 and 1 for binary problems, or \n
           0, 1, 2, ..., for multiclass problems.")
    }
    nfac <- sort(nfac)
    lnfac <- length(nfac)
    if (!((is.numeric(nfac)) || (is.integer(nfac)))) {
      stop("Input 'nfac' must be of class integer or numeric.")
    }
    if (!((is.integer(nfolds)) || (is.numeric(nfolds)))) {
      stop("Input 'nfolds' must be of class integer or numeric.")
    }
    if (length(nfolds) != 1L || nfolds < 2L || nfolds != floor(nfolds)) {
      stop("Input 'nfolds' must be a single integer equal to or greater \n
           than 2.")
    }
    if (nfolds > length(y)) {
      stop("Input 'nfolds' must be a single whole number equal to or less \n 
           than the number of labels in 'y'.")
    }
    if (is.null(foldid)) {
      foldid <- sample(rep(1:nfolds, length.out = length(y)))
    } else {
      if (length(foldid) != xdim[cmode]) {
        stop("Input 'foldid' must match number of levels in classification \n
             mode.")
      }
      if (!(is.integer(foldid))) {
        stop("Input 'foldid' must be of class integer.")
      }
      if (length(unique(foldid)) < 2L) {
        stop("Input 'foldid' must contain IDs for two or more folds.")
      }
      if (!(all(foldid == floor(foldid)))) {
        stop("Input 'foldid' must contain integers.")
      }
      if (length(unique(foldid)) != as.integer(nfolds)) {
        stop("Input 'foldid' must contain the number of unique values equal \n
             to nfolds.")
      }
    }
    omethods <- c("PLR", "SVM", "RF", "NN", "RDA", "GBM") 
    checkmethod <- sum(toupper(method) %in% omethods)
    if (checkmethod != length(method)) {
      stop("Input 'method' contains at least one value that is not valid.")
    }
    method <- which(omethods %in% toupper(method) == T)
    if (length(method) == 0) {method <- 1:6}
    if (!(is.logical(verbose))) {
      stop("Input 'verbose' must be of class logical.")
    }
    if (length(verbose) != 1L) {stop("Input 'verbose' must be a single value.")}
    if ('1' %in% method) {
      if (is.null(alpha)) {alpha <- seq(0, 1, length = 6)}
      if (any(alpha < 0L) || any(alpha > 1L)) {
        stop("Input 'alpha' must contain real numbers 
                  between zero and one (inclusive).")
      }
      if (!(is.numeric(alpha))) {stop("Input 'alpha' must be numeric.")}
      alpha <- sort(alpha)
    }
    if ('2' %in% method) {
      if (is.null(gamma)) {gamma <- c(0, 0.01, 0.1, 1, 10, 100, 1000)}
      if (any(gamma < 0L)) {
        stop("Input 'gamma' must contain real numbers equal to or greater \n 
             than zero.")
      }
      if (!(is.numeric(gamma))) {stop("Input 'gamma' must be numeric.")}
      gamma <- sort(gamma)
      if (is.null(cost)) {cost <- c(1, 2, 4, 8, 16, 32, 64)}
      if (any(cost <= 0L)) {
        stop("Input 'cost' must be a real number greater than zero.")
      }
      if (!(is.numeric(cost))) {stop("Input 'cost' must be numeric.")}
      cost <- sort(cost)
      svm.grid <- expand.grid(gamma, cost)
    }
    if ('3' %in% method) {
      if (is.null(ntree)) {ntree <- c(100, 200, 400, 600, 800, 1600, 3200)} 
      if (any(ntree < 1L)) {
        stop("Input 'ntree' must be greater than or equal to one.")
      }
      if (!(all(ntree == floor(ntree)))) {
        stop("Input 'ntree' must contain only integers.")
      }
      if (!(is.numeric(ntree))) {stop("Input 'ntree' must be numeric.")}
      ntree <- sort(ntree)
      if (is.null(nodesize)) {nodesize <- c(1, 2, 4, 8, 16, 32, 64)}
      if (any(nodesize < 1L)) {
        stop("Input 'nodesize' must be greater than or equal to one.") 
      }
      if (!(all(nodesize == floor(nodesize)))) {
        stop("Input 'nodesize' must be integer.")
      }
      if (!(is.numeric(nodesize))) {stop("Input 'nodesize' must be numeric.")}
      nodesize <- sort(nodesize)
      rf.grid <- expand.grid(ntree, nodesize)
    }
    if ('4' %in% method) {
      if (is.null(size)) {size <- c(1, 2, 4, 8, 16, 32, 64)} 
      if (any(size < 0L)) {
        stop("Input 'size' must be greater than or equal to zero.")
      }
      if (!(all(size == floor(size)))) {
        stop("Input 'size' must contain only integers.")
      }
      if (!(is.numeric(size))) {stop("Input 'size' must be numeric.")}
      size <- sort(size)
      if (is.null(decay)) {decay <- c(0.001, 0.01, 0.1, 1, 2, 4, 8, 16)}
      if (!(is.numeric(decay))) {stop("Input 'decay' must be numeric.")}
      decay <- sort(decay)
      nn.grid <- expand.grid(size, decay)
    }
    if ('5' %in% method) {
      if (is.null(rda.alpha)) {rda.alpha <- seq(0, 0.999, length = 6)}
      if (any(rda.alpha < 0L) || any(rda.alpha >= 1L)) {
        stop("Input 'rda.alpha' must contain real numbers equal to or greater \n
             than zero, and less than one.")
      }
      if (!(is.numeric(rda.alpha))) {stop("Input 'rda.alpha' must be numeric.")}
      rda.alpha <- sort(rda.alpha)
      if (is.null(delta)) {delta <- c(0, 0.1, 1, 2, 3, 4)}
      if (any(delta < 0L)) {
        stop("Input 'delta' must be a real number greater than or equal \n 
             to zero.")
      }
      if (!(is.numeric(delta))) {stop("Input 'delta' must be numeric.")}
      delta <- sort(delta)
      rda.grid <- expand.grid(rda.alpha, delta)
    }      
    if ('6' %in% method) {
      if (is.null(eta)) {eta <- c(0.1, 0.3, 0.5, 0.7, 0.9)}
      if (any(eta <= 0L) || any(eta >= 1L)) {
        stop("Input 'eta' must contain real numbers greater than zero and \n 
             less than one.")
      }
      if (!(is.numeric(eta))) {stop("Input 'eta' must be numeric.")}
      eta <- sort(eta)
      if (is.null(max.depth)) {max.depth <- c(1, 2, 3, 4)} 
      if (any(max.depth < 1L)) {
        stop("Input 'max.depth' must be greater than or equal to one.")
      }
      if (!(all(max.depth == floor(max.depth)))) {
        stop("Input 'max.depth' must contain only integers.")
      }
      if (!(is.numeric(max.depth))) {stop("Input 'max.depth' must be numeric.")}
      max.depth <- sort(max.depth)
      if (is.null(subsample)) {subsample <- c(0.6, 0.7, 0.8, 0.9)}
      if (any(subsample <= 0L) || any(subsample > 1L)) {
        stop("Input 'subsample' must be greater than zero and less than or \n
             equal to one.")
      }
      if (!(is.numeric(subsample))) {stop("Input 'subsample' must be numeric.")}
      if (is.null(nrounds)) {nrounds <- c(100, 200, 300, 500)} 
      if (any(nrounds < 1L)) {
        stop("Input 'nrounds' must be greater than or equal to one.")
      }
      if (!(all(nrounds == floor(nrounds)))) {
        stop("Input 'nrounds' must contain only integers.")
      }
      if (!(is.numeric(nrounds))) {stop("Input 'nrounds' must be numeric.")}
      gbm.grid <- expand.grid(eta, max.depth, subsample, nrounds)
    }
    families <- c("binomial", "multinomial")
    family <- pmatch(tolower(family), families)
    lfam <- length(family)
    if (any(is.na(family))) {
      stop("Input 'family' must be a character value, either 'binomial' \n
           or 'multinomial'.")
    }
    if (lfam == 0L) {
      stop("Input 'family' must not be NULL.")
    } else if(lfam > 2L) {
      stop("Input 'family' must be 'binomial' or 'multinomial'.")
    } else if (lfam == 2L && (length(levels(y)) == 2L)) {
      family <- "binomial"
      warning("Input 'family' was not specified. Two classes detected: \n
              defaulting to family = 'binomial'.")
    } else if (lfam == 2L && (length(levels(y)) > 2L)) {
      family <- "multinomial"
      warning("Input 'family' was not specified. Three or more classes \n 
              detected: defaulting to family = 'multinomial'.")
    } else if (lfam == 1L && family == 1L) {
      family <- "binomial"
    } else if (lfam == 1L && family == 2L) {
      family <- "multinomial"
    }
    if (is.null(prior)) {
      frac <- table(y) / length(y)
      prior <- as.numeric(frac)
      if (family == "binomial") {
        names(frac) <- c(0, 1)
        weight <- as.numeric(frac[y])
      }
      if (family == "multinomial") {
        names(frac) <- seq_along(levels(y)) - 1
        weight <- as.numeric(frac[y])
      }
    } else {
      if (!(abs(sum(prior) - 1) < .Machine$double.eps^0.5)) {
        stop("Values within input 'prior' must sum to one.")
      }
      if (!(is.numeric(prior))) {stop("Input 'prior' must be numeric.")}
      if (family == "binomial") {
        if (!(length(prior) == 2L))
          stop("Input 'prior' must contain two values for family of \n
               'binomial'.")
      }
      if (family == "multinomial") {
        if (!(length(prior) >= 3L))
          stop("Input 'prior' must contain three or more values for family \n
               of 'multinomial'.")
      }
      frac <- as.table(prior)                                                  
      if (family == "binomial") {
        names(frac) <- c(0, 1)
        weight <- as.numeric(frac[y])
      }
      if (family == "multinomial") {
        names(frac) <- seq_along(levels(y)) - 1
        weight <- as.numeric(frac[y])
      }
    }
    if ('1' %in% method) {plr.weights <- as.numeric((table(y)[y] / length(y)))}
    Aweights <- Bweights <- Cweights <- Phi <- vector("list", lnfac)
    train.weights <- opt.model <- Aweights                                         
    opt.param <- est.time <- kcv.error <- NULL
    if (!(is.logical(parallel))) {
      stop("Input 'parallel' must be of class logical.")
    }
    if (length(parallel) != 1L) {
      stop("Input 'parallel' must be a single value.")
    }
    if (parallel) {
      if (is.null(cl)) {cl <- makeCluster(detectCores())}
        ce <- clusterEvalQ(cl, library(multiway))
        registerDoParallel(cl)
    }
    for (w in 1:lnfac) {
       optmodel.new <- vector("list", 6)                                        
       if (model == "parafac") {
         if (verbose == T) {
           cat("nfac =", nfac[w], "method = parafac", fill = T)
         }
         tic <- proc.time()
         pfac <- parafac(X = x, nfac = nfac[w], parallel = parallel, cl = cl,   
                         verbose = verbose, ...)
         toc <- proc.time()
         if (w == 1) {const <- pfac$control$const}
         time.pfac <- toc[3] - tic[3]
         Aweights[[w]] <- pfac$A
         Bweights[[w]] <- pfac$B                                               
         train <- train.weights[[w]] <- as.matrix(pfac$C)
         if (lxdim == 4L) {
           Cweights[[w]] <- pfac$C                                               
           train <- train.weights[[w]] <- as.matrix(pfac$D)
         }
       } else {
         if (verbose == T) {
           cat("nfac =", nfac[w], "method = parafac2", fill = T)
         }
         tic <- proc.time()
         pfac <- parafac2(X = x, nfac = nfac[w], parallel = parallel, cl = cl,  
                          verbose = verbose, ...)
         toc <- proc.time()
         if (w == 1) {const <- pfac$control$const}
         time.pfac <- toc[3] - tic[3]
         Aweights[[w]] <- pfac$A
         Bweights[[w]] <- pfac$B                                                  
         train <- train.weights[[w]] <- as.matrix(pfac$C)
         if (lxdim == 4L) {
           Cweights[[w]] <- pfac$C                                              
           train <- train.weights[[w]] <- as.matrix(pfac$D)
         }
         Phi[[w]] <- pfac$Phi
       }
       if ('1' %in% method) {
         if (verbose == T) {
           cat("nfac =", nfac[w], "method = plr", fill = T)
         }
         tic <- proc.time()
         if (nfac[w] == 1 || nfac[w] == 1L) {
           train.plr <- cbind(train, 0)
           if (family == "multinomial") {
             stop("Input combination nfac = 1, family = 'multinomial', and \n
                  method = 'PLR' gives an error due to an issue interfacing \n
                  with glmnet::cv.glmnet. Until resolved, this combination of \n
                  arguments is not permitted.")
           }
         } else {
           train.plr <- train
         }
         plr.results <- kcv.plr(x = train.plr,  y = y, nfolds = nfolds,  
                                foldid = foldid, alpha = alpha,
                                family = family, lambda = lambda,
                                weights = plr.weights, parallel = parallel)
         toc <- proc.time()
         time.plr <- toc[3] - tic[3]
         error.plr <- plr.results$error
         plr.id <- plr.results$alpha.id
         plr.opt <- alpha[plr.id]
         optmodel.new[[1]] <- plr.fit <- plr.results$plr.fit
         if (nfolds == 2) {lambda.min <- plr.results$lambda.min}
         if (nfolds >= 3) {lambda.min <- plr.fit$lambda.min}
       } else {
         lambda.min <- plr.opt <- error.plr <- time.plr <- NA
         optmodel.new[[1]] <- NULL
       }
       if ('2' %in% method) {
         if (verbose == T) {
           cat("nfac =", nfac[w], "method = svm", fill = T)
         }
         tic <- proc.time()
         svm.results <- kcv.svm(x = train, y = y, nfolds = nfolds,
                                foldid = foldid, svm.grid = svm.grid, 
                                class.weights = frac, parallel = parallel)
         toc <- proc.time()
         time.svm <- toc[3] - tic[3]
         error.svm <- svm.results$error
         svm.id <- svm.results$svm.grid.id
         svm.opt <- svm.grid[svm.id, ]
         optmodel.new[[2]] <- svm.fit <- svm.results$svm.fit
       } else {
         error.svm <- time.svm <- NA
         svm.opt <- data.frame(NA, NA)
         optmodel.new[[2]] <- NULL
       }
       if ('3' %in% method) {
         if (verbose == T) {
           cat("nfac =", nfac[w], "method = rf", fill = T)
         }
         tic <- proc.time()
         rf.results <- kcv.rf(x = train, y = y, nfolds = nfolds,
                              foldid = foldid, rf.grid = rf.grid, 
                              classwt = prior, parallel = parallel)
         toc <- proc.time()
         time.rf <- toc[3] - tic[3]
         error.rf <- rf.results$error
         rf.id <- rf.results$rf.grid.id
         rf.opt <- rf.grid[rf.id, ]
         optmodel.new[[3]] <- rf.results$rf.fit
       } else {
         error.rf <- time.rf <- NA
         rf.opt <- data.frame(NA, NA)
         optmodel.new[[3]] <- NULL
       }
       if ('4' %in% method) {
         if (verbose == T) {
           cat("nfac =", nfac[w], "method = nn", fill = T)
         }
         tic <- proc.time()
         nn.results <- kcv.nn(x = train, y = y, nfolds = nfolds,
                              foldid = foldid, nn.grid = nn.grid, 
                              weights = weight, parallel = parallel)
         toc <- proc.time()
         time.nn <- toc[3] - tic[3]
         error.nn <- nn.results$error
         nn.id <- nn.results$nn.grid.id
         nn.opt <- nn.grid[nn.id, ]
         optmodel.new[[4]] <- nn.fit <- nn.results$nn.fit
       } else {
         nn.opt <- error.nn <- time.nn <- NA
         nn.opt <- data.frame(NA, NA)
         optmodel.new[[4]] <- NULL
       }
       if ('5' %in% method) {
         if (verbose == T) {
           cat("nfac =", nfac[w], "method = rda", fill = T)
         }
         tic <- proc.time()
         rda.results <- kcv.rda(x = train, y = as.numeric(y) - 1, 
                                nfolds = nfolds, foldid = foldid, 
                                rda.grid = rda.grid, prior = frac, 
                                parallel = parallel)
         toc <- proc.time()
         time.rda <- toc[3] - tic[3]
         error.rda <- rda.results$error
         rda.id <- rda.results$rda.grid.id
         rda.opt <- rda.grid[rda.id, ]
         optmodel.new[[5]] <- rda.fit <- rda.results$rda.fit                  
       } else {
         error.rda <- time.rda <- NA
         rda.opt <- data.frame(NA, NA)
         optmodel.new[[5]] <- NULL
       }
       if ('6' %in% method) {
         if (verbose == T) {
           cat("nfac =", nfac[w], "method = gbm", fill = T)
         }
         tic <- proc.time()
         gbm.results <- kcv.gbm(x = train, y = y, nfolds = nfolds, 
                                foldid = foldid, family = family, 
                                gbm.grid = gbm.grid, prior = frac, 
                                parallel = parallel)
         toc <- proc.time()
         time.gbm <- toc[3] - tic[3]
         error.gbm <- gbm.results$error
         gbm.id <- gbm.results$gbm.grid.id
         gbm.opt <- gbm.grid[gbm.id, ]
         optmodel.new[[6]] <- gbm.fit <- gbm.results$gbm.fit                  
       } else {
         error.gbm <- time.gbm <- NA
         gbm.opt <- data.frame(NA, NA, NA, NA)
         optmodel.new[[6]] <- NULL
       } 
       opt.model[[w]] <- optmodel.new
       optparam.new <- data.frame(nfac = nfac[w], alpha = plr.opt,
                                  lambda = lambda.min, gamma = svm.opt[1, 1], 
                                  cost = svm.opt[1, 2], ntree = rf.opt[1, 1], 
                                  nodesize = rf.opt[1, 2], size = nn.opt[1, 1], 
                                  decay = nn.opt[1, 2], 
                                  rda.alpha = rda.opt[1, 1], 
                                  delta = rda.opt[1, 2], eta = gbm.opt[1, 1],
                                  max.depth = gbm.opt[1, 2], 
                                  subsample = gbm.opt[1, 3], 
                                  nrounds = gbm.opt[1, 4])
       esttime.new <- data.frame(nfac = nfac[w], time.pfac = time.pfac,
                                 time.plr = time.plr, time.svm = time.svm,
                                 time.rf = time.rf, time.nn = time.nn,
                                 time.rda = time.rda, time.gbm)
       kcv.error.new <- data.frame(nfac = nfac[w], error.plr = error.plr,
                                   error.svm = error.svm, error.rf = error.rf,
                                   error.nn = error.nn, error.rda = error.rda,
                                   error.gbm)
       opt.param <- rbind(opt.param, optparam.new)
       est.time <- rbind(est.time, esttime.new)
       kcv.error <- rbind(kcv.error, kcv.error.new)
    }                                                                      
    if (parallel == T) {stopCluster(cl)}
    levels(y) <- names(frac)
    if (lxdim == 3L) {
      tcpfalist <- list(opt.model = opt.model, opt.param = opt.param, 
                        kcv.error = kcv.error, est.time = est.time, 
                        model = model, method = method, x = x, y = y, 
                        Aweights = Aweights, Bweights = Bweights, 
                        Cweights = NULL, Phi = Phi, const = const,
                        cmode = cmode, family = family, xdim = xdim, 
                        lxdim = lxdim, train.weights = train.weights)
    }
    if (lxdim == 4L) {
      tcpfalist <- list(opt.model = opt.model, opt.param = opt.param, 
                        kcv.error = kcv.error, est.time = est.time, 
                        model = model, method = method, x = x, y = y, 
                        Aweights = Aweights, Bweights = Bweights, 
                        Cweights = Cweights, Phi = Phi, const = const, 
                        cmode = cmode, family = family, xdim = xdim, 
                        lxdim = lxdim, train.weights = train.weights)
    }
    class(tcpfalist) <- "tunecpfa"
    return(tcpfalist)
}