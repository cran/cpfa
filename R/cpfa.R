cpfa <-
  function(x, y, z = NULL, model = c("parafac", "parafac2", "pca"), nfac = 1,
           nrep = 5, ratio = 0.8, nfolds = 10,
           method = c("PLR", "SVM", "RF", "NN", "RDA", "GBM"),
           family = c("binomial", "multinomial"), parameters = list(),
           type.out = c("measures", "descriptives"), foldid = NULL,
           prior = NULL, cmode = NULL, seeds = NULL, plot.out = FALSE,
           plot.measures = NULL, parallel = FALSE, cl = NULL,
           verbose = TRUE, compscale = TRUE, pcarot = c("unrotated", "varimax"),
           ...)
{
    permflag <- FALSE
    models <- c("parafac", "parafac2", "pca")
    model0 <- sum(tolower(model) %in% models)
    if ((model0 == 0L) || (model0 > 1L)) {
      stop("Input 'model' not specified correctly. Input must be one of \n
           either 'parafac', 'parafac2', or 'pca'.")
    }
    model <- tolower(model)
    flattened <- NULL
    if (is.array(x) && (!(is.matrix(x))) && (model %in% c("parafac", "pca"))) {
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
          stop("Input 'cmode' must be 1, 2, or 3 (or 4 if 'x' is 4-way).")
        }
        if (model == "parafac") {
          modeval <- 1:lxdim
          mode.re <- c(modeval[-cmode], cmode)
          x <- aperm(x, mode.re)
          xdim <- dim(x)
          permflag <- TRUE
        } else {
          modeval <- 1:lxdim
          mode.re <- c(cmode, modeval[-cmode])
          x <- matrix(aperm(x, mode.re), nrow = xdim[cmode])
          xdim <- dim(x)
          lxdim <- length(xdim)
          flattened <- TRUE
        }
      } else {
        cmode <- lxdim
        if (model == "pca") {
          modeval <- 1:lxdim
          mode.re <- c(cmode, modeval[-cmode])
          x <- matrix(aperm(x, mode.re), nrow = xdim[cmode])
          xdim <- dim(x)
          lxdim <- length(xdim)
          flattened <- TRUE
          if (length(y) != dim(x)[1]) {
            stop("Input 'x' was a 3-way or 4-way array, and model was 'pca'. \n 
                 Input 'y' had a length different from the last mode of 'x'. \n 
                 When 'x' is a 3-way or 4-way array, and when model is 'pca', \n 
                 the last mode must be the classification mode. Consider \n 
                 permuting the array.")
          }
        }
      }
    } else if (is.array(x) && (!(is.matrix(x))) && (model == "parafac2")) {
      xdim <- dim(x)                             
      lxdim <- length(xdim)
      if (!((lxdim == 3L) || (lxdim == 4L))) {
        stop("Input 'x' must be a 3-way or 4-way array.")
      }
      if (any(is.nan(x)) || any(is.infinite(x))) {
        stop("Input 'x' cannot contain NaN or Inf values.")
      }
      if (any(is.na(x))) {stop("Input 'x' cannot contain missing values.")}
      if (!(is.null(cmode))) {
        cmode <- lxdim
        warning("Input 'cmode' is ignored when 'model = parafac2'. Last mode \n
                is classification mode by default.")
      } else {
        cmode <- lxdim
      }
      storlist <- vector("list", xdim[cmode])
      if (lxdim == 3L) {
        for (k in 1:xdim[cmode]) {storlist[[k]] <- x[, , k]}
      } else {
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
          stop("Input 'x' must be list of matrices with same 
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
    } else if (is.list(x) && (model != "parafac2")) {
      stop("Input 'x' cannot be of class 'list' if model is 'parafac' or \n
           'pca'.")
    } else if (is.matrix(x) && (model != "pca")) {
      stop("Input 'x' cannot be of class 'matrix' if model is 'parafac' or \n
           'parafac2'.")
    } else if (is.matrix(x) && (model == "pca")) {
      xdim <- origdim <- dim(x)                             
      lxdim <- length(xdim)
      if (!((lxdim == 2L))) {
        stop("If 'x' has class 'matrix' and model is 'pca', 'x' must be a \n
             2-way matrix.")
      }
      if (any(is.nan(x)) || any(is.infinite(x))) {
        stop("Input 'x' cannot contain NaN or Inf values.")
      }
      if (any(is.na(x))) {stop("Input 'x' cannot contain missing values.")}
      if (!(is.null(cmode))) {
        if (!(cmode %in% (1:lxdim))) {
          stop("Input 'cmode' must be 1 or 2 when model is 'pca' and when \n
               'x' is of class 'matrix'.")
        }
        if (cmode == 1) {
          x <- x
        } else {
          modeval <- 1:lxdim
          mode.re <- c(cmode, modeval[-cmode])
          x <- aperm(x, mode.re)
          xdim <- dim(x)
          permflag <- TRUE
        }
      } else {
        cmode <- 1
      }
      flattened <- FALSE
    } else {
      stop("Input 'x' can be of class 'array' for any model, can be of class \n
           'list' for 'parafac2', and can be of class 'matrix' for 'pca'. \n
           Else, input 'x' cannot be of a different class.")
    }
    if (!(is.factor(y))) {stop("Input 'y' must be of class 'factor'.")}
    ylev <- length(levels(y))
    if (!(sum(sort(as.numeric(levels(y))) == 0:(ylev - 1)) == ylev)) {
      stop("Input 'y' must contain labels of 0 and 1 for binary problems, or \n
           0, 1, 2, ..., for multiclass problems.")
    }
    if (model == "parafac") {
      if (!(length(y) == origdim[cmode])) {
        stop("Length of 'y' must match number of levels in classification \n
             mode of 'x'.")
      }
    } else if (model == "pca") {
      if (!(length(y) == xdim[1])) {
        stop("Length of 'y' must match number of levels in classification \n
             mode of 'x'.")
      }
    } else {
      if (!(length(y) == xdim[cmode])) {
        stop("Length of 'y' must match number of levels in classification \n
             mode of 'x'.")
      }
    }
    numcheck(nrep)
    if ((!(ceiling(nrep) == nrep)) || (nrep < 1)) {
      stop("Input 'nrep' must be an integer greater than 0.")
    }
    numcheck(ratio)
    if ((ratio >= 1) || (ratio <= 0)) {
      stop("Input 'ratio' must be a number between 0 and 1, exclusive.")
    }
    if (is.null(seeds)) {
      seeds <- 1:nrep
    } else {
      if (!((is.numeric(seeds)) || (is.integer(seeds)))) {
        stop("Input 'seeds' must be of class 'numeric' or 'integer'.")
      }
      if (length(seeds) != nrep) {
        stop("Input 'seeds' must have length equal to input 'nrep'.")
      }
      if (length(unique(seeds)) != length(seeds)) {
        warning("Not all seeds are unique.")
      }
    }
    types <- c("measures", "descriptives")
    numtype <- sum(tolower(type.out) %in% types)
    if (numtype == 0) {
      stop("Input 'type.out' does not contain a valid value. Must specify \n
           either 'measures' or 'descriptives' for 'type.out'.")
    } else if (numtype > 2L) {
      stop("Input 'type.out' contains three or more values. Must specify \n
           either 'measures' or 'descriptives' for 'type.out'.")
    } else if (numtype == 2L) {
      type.out <- "descriptives"
    } else {
      type.out <- tolower(type.out)
    }
    if (model == "parafac") {
      nobs <- dim(x)[lxdim]
    } else if (model == "pca") {
      nobs <- nrow(x)
    } else {
      nobs <- length(x)
    }
    if (!(is.null(z))) {
      zdim <- dim(z)                                                            
      lzdim <- length(zdim)
      if (!((lzdim == 2L))) {stop("Input 'z' must have two modes.")}
      if (!(is.matrix(z))) {stop("Input 'z' must be of class 'matrix'.")}
      if (any(is.nan(z)) || any(is.infinite(z))) {
        stop("Input 'z' cannot contain NaN or Inf values.")
      }
      if (any(is.na(z))) {stop("Input 'z' cannot contain missing values.")}
      if (zdim[1] != length(y)) {
        stop("Input 'z', when provided, must have number of rows equal to \n
             the length of input 'y' (i.e., the number of class labels).")
      }
    }
    ntrain <- ceiling(nobs * ratio)
    if (!(is.null(foldid))) {
      if (!(is.list(foldid))) {
        stop("Input 'foldid' must be of class 'list' when provided.")
      }
      if (length(foldid) != nrep) {
        stop("Input 'foldid' must be a list of length equal to input 'nrep'.")
      }
      if (any(unlist(lapply(foldid, length)) != ntrain)) {
        stop("Each vector in 'foldid' must have length equal to \n
             ceiling(nobs * ratio).")
      }
    }
    logicheck(plot.out)
    if (plot.out) {
      if (is.null(plot.measures)) {
        plottype <- 5
      } else {
        cmeasures <- c("err", "acc", "tpr", "fpr", "tnr", "fnr", "ppv", "npv",
                       "fdr", "fom", "fs")
        plottype.num <- sum(cmeasures %in% plot.measures)
        if (plottype.num == 0) {
          stop("Input 'plot.out' is true, but input 'plot.measures' does not \n
               contain any accepted values. See help file and argument \n
               'plot.measures' for a list of accepted values.")
        }
        plottype <- which(cmeasures %in% plot.measures == TRUE) + 3
      }
    }
    nfac <- sort(nfac)
    stor <- array(0, dim = c(length(nfac) * length(method), 11, nrep))
    predstor <- Aw <- Bw <- Cw <- Pw <- vector(mode = "list", length = nrep)
    trainIDs <- testIDs <- opara <- predstor
    cmode0 <- cmode
    if (cmode == lxdim) {cmode <- NULL}
    logicheck(verbose)
    ccreated <- FALSE
    if ((parallel == TRUE) && (is.null(cl))) {
      cl <- makeCluster(detectCores())
      ccreated <- TRUE
      registerDoParallel(cl)
      clusterEvalQ(cl, library(multiway))
    }
    for (i in 1:nrep) {
      if (verbose == TRUE) {cat("nrep =", i, " \n")}
      set.seed(seed = seeds[i])
      train.id <- trainIDs[[i]] <- sample.int(nobs, size = ntrain)
      alllevels <- 1:nobs
      testIDs[[i]] <- alllevels[-train.id]
      y.train <- y[train.id]
      y.test <- as.numeric(y[-train.id]) - 1
      if (!(is.null(z))) {
        z.train <- z[train.id, , drop = FALSE]
        z.test <- z[-train.id, , drop = FALSE]
      } else {
        z.train <- z.test <- NULL
      }
      if (model == "parafac") {
        if (lxdim == 3L) {
          X.train <- x[, , train.id]                       
          X.test <- x[, , -train.id]
        } else {
          X.train <- x[, , , train.id]
          X.test <- x[, , , -train.id]
        }
      } else if (model == "pca") {
        X.train <- x[train.id, , drop = FALSE]
        X.test <- x[-train.id, , drop = FALSE]
      } else {
        X.train <- x[train.id]
        X.test <- x[-train.id]
      }
      if (is.null(foldid)) {
        cfoldid <- NULL
      } else {
        cfoldid <- foldid[[i]]
      }
      tcpfalist <- tunecpfa(x = X.train, y = y.train, z = z.train, 
                            model = model, nfac = nfac, nfolds = nfolds, 
                            method = method, family = family, 
                            parameters = parameters, foldid = cfoldid, 
                            prior = prior, cmode = NULL, parallel = parallel, 
                            cl = cl, verbose = verbose, compscale = compscale, 
                            pcarot = pcarot, ...)
      Aw[[i]] <- tcpfalist$Aweights
      Bw[[i]] <- tcpfalist$Bweights
      Cw[[i]] <- tcpfalist$Cweights
      Pw[[i]] <- tcpfalist$Phi
      opara[[i]] <- tcpfalist$opt.param
      yhat <- predict(object = tcpfalist, newdata = X.test, newdata.z = z.test, 
                      type = "response")           
      out <- cpm.all(x = yhat, y = y.test, level = levels(y))
      stor[ , , i] <- as.matrix(out$cpms)
      predstor[[i]] <- predict(object = tcpfalist, newdata = X.test,
                               newdata.z = z.test, type = "classify.weights")
    }
    mconst <- tcpfalist$const
    rnam <- rownames(out$cpms)
    cnam <- colnames(out$cpms)
    dimnames(stor)[[1]] <- rnam
    dimnames(stor)[[2]] <- cnam
    train.weights <- list(Atrain.weights = Aw, Btrain.weights = Bw,
                          Ctrain.weights = Cw, Phitrain = Pw)
    mean.tune.param <- Reduce("+", opara) / length(opara)
    if (plot.out == TRUE) {
      ncomps <- length(nfac)
      nmethods <- length(method)
      plotstor <- data.frame(matrix(0, nrow = (ncomps * nmethods * nrep),
                                    ncol = 14))
      plotcname <- c("method", "nfac", "rep", colnames(stor))
      colnames(plotstor) <- plotcname
      matnum <- ncomps * nmethods
      methnames0 <- sapply(strsplit(rownames(stor), split = '.', fixed = TRUE),
                           function(x) (x[2]))
      methnames <- gsub('[[:digit:]]+', '', methnames0)
      nfacnames <- as.numeric(gsub(".*?([0-9]+).*", "\\1", rownames(stor)))
      for (i in 1:nrep) {
        indl <- matnum * (i - 1) + 1
        indu <- matnum * i
        plotstor[indl:indu, 1] <- methnames
        plotstor[indl:indu, 2] <- nfacnames
        plotstor[indl:indu, 3] <- i
        plotstor[indl:indu, 4:14] <- stor[, , i]
      }
      toplot <- colnames(plotstor)[plottype]
      for (j in 1:length(plottype)) {
        pformula <- formula(paste0(toplot[j], " ~ ", "method * nfac"))
        boxplot(pformula, data = plotstor, ylim = c(0, 1),
                xlab = "Method and Number of Components", na.rm = FALSE,
                ylab = toupper(toplot[j]), main = "Performance Measure")
      }
    }
    if (ccreated == TRUE) {stopCluster(cl)}
    if (permflag == TRUE) {x <- aperm(x, perm = order(mode.re))}
    if (type.out == "measures") {
      cpfalist <- list(measure = stor, predweights = predstor,
                       train.weights = train.weights, opt.tune = opara,
                       mean.opt.tune = mean.tune.param, X = x, y = y, z = z,
                       nfac = nfac, model = model, method = method,
                       const = mconst, cmode = cmode0, family = family, 
                       lxdim = lxdim, trainIDs = trainIDs, testIDs = testIDs,
                       flattened = flattened)
      class(cpfalist) <- "wrapcpfa"
      return(cpfalist)                              
    } else {
      dfun <- c("mean", "median", "sd")
      output <- vector(mode = "list", length = length(dfun))
      for (j in seq_along(dfun)) {
         output[[j]] <- apply(stor, 1:2,
                              FUN = function(x){return(get(dfun[j])(x,
                                                                na.rm = TRUE))})
        rownames(output[[j]]) <- rnam
        colnames(output[[j]]) <- cnam
      }
      names(output) <- dfun  
      cpfalist <- list(descriptive = output, predweights = predstor,
                       train.weights = train.weights, opt.tune = opara,
                       mean.opt.tune = mean.tune.param, X = x, y = y, z = z,
                       nfac = nfac, model = model, method = method,
                       const = mconst, cmode = cmode0, family = family, 
                       lxdim = lxdim, trainIDs = trainIDs, testIDs = testIDs,
                       flattened = flattened)
      class(cpfalist) <- "wrapcpfa"
      return(cpfalist)
    }
}