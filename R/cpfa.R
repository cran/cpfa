cpfa <- 
  function(x, y, model = c("parafac", "parafac2"), nfac = 1, 
           nrep = 5, ratio = 0.8, nfolds = 10,
           method = c("PLR", "SVM", "RF", "NN", "RDA", "GBM"),
           family = c("binomial", "multinomial"), parameters = list(), 
           type.out = c("measures", "descriptives"), foldid = NULL, 
           prior = NULL, cmode = NULL, seeds = NULL, plot.out = FALSE, 
           plot.measures = NULL, parallel = FALSE, cl = NULL, 
           verbose = TRUE, ...) 
{   
    models <- c("parafac", "parafac2")
    model0 <- sum(tolower(model) %in% models)
    if ((model0 == 0) || (model0 > 1L)) {
      stop("Input 'model' not specified correctly. Input must be only one of \n
           either 'parafac' or 'parafac2'.")
    }
    model <- tolower(model)
    if (is.array(x) && (model == "parafac")) {
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
        if (!(cmode %in% (1:lxdim))) {
          stop("Input 'cmode' must be 1, 2, or 3 \n 
               (or 4 if 'x' is a 4-way array).")
        }
        modeval <- 1:lxdim
        mode.re <- c(modeval[-cmode], cmode)
        x <- aperm(x, mode.re)
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
      if (any(as.logical(lapply(x, function(a){return(any(is.infinite(a)))})))){
        stop("Input 'x' cannot contain Inf values")
      }
      if (any(as.logical(lapply(x, function(a){return(any(is.na(a)))})))) {
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
    if (!is.factor(y)) {stop("Input 'y' must be of class 'factor'.")}
    if (!(length(y) == xdim[cmode])) {
      stop("Length of 'y' must match number of levels in classification \n 
           mode of 'x'.")
    }
    if ((!(ceiling(nrep) == nrep)) || (nrep < 1) || (length(ratio) != 1L)) {
      stop("Input 'nrep' must be a single integer greater than 0.")
    }
    if (!(is.numeric(ratio))) {stop("Input 'ratio' must be of class numeric.")}
    if ((ratio > 1) || (ratio < 0) || (length(ratio) != 1L)) {
      stop("Input 'ratio' must be a number between 0 and 1, inclusive.")
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
    } else {
      nobs <- length(x)
    }
    ntrain <- ceiling(nobs * ratio)
    if (!(is.null(foldid))) {
      if (length(foldid) != ntrain) {
        stop("Input 'foldid' must have length equal to ceiling(nobs * ratio) \n
             where 'nobs' is the number of observations in the classification \n
             mode.")
      }
    }
    if (!(is.logical(plot.out))) {
      stop("Input 'plot.out' must be a logical value.")
    }
    if (length(plot.out) != 1L) {
      stop("Input 'plot.out' must be a single value.")
    }
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
        plottype <- which(cmeasures %in% plot.measures == T) + 3
      } 
    }
    stor <- array(0, dim = c(length(nfac) * length(method), 11, nrep))
    predstor <- Aw <- Bw <- Cw <- Pw <- vector(mode = "list", length = nrep)
    opara <- predstor
    cmode0 <- cmode
    if (cmode == lxdim) {cmode <- NULL}
    for (i in 1:nrep) {
       if (verbose == T) {cat("nrep =", i, " \n")}
       set.seed(seed = seeds[i])
       train.id <- sample.int(nobs, size = ntrain)
       y.train <- y[train.id]
       y.test <- as.numeric(y[-train.id]) - 1 
       if (model == "parafac") {
         if (lxdim == 3L) { 
           X.train <- x[, , train.id]                                               
           X.test <- x[, , -train.id]
         } else {
           X.train <- x[, , , train.id]
           X.test <- x[, , , -train.id]
         }
       } else {
         X.train <- x[train.id]
         X.test <- x[-train.id]
       }
       tcpfalist <- tunecpfa(x = X.train, y = y.train, nfac = nfac,            
                             nfolds = nfolds, method = method, foldid = foldid, 
                             prior = prior, model = model, family = family,
                             parameters = parameters, parallel = parallel, 
                             cl = cl, verbose = verbose, cmode = cmode, ...)
       Aw[[i]] <- tcpfalist$Aweights
       Bw[[i]] <- tcpfalist$Bweights
       Cw[[i]] <- tcpfalist$Cweights
       Pw[[i]] <- tcpfalist$Phi
       opara[[i]] <- tcpfalist$opt.param
       yhat <- predict(object = tcpfalist, newdata = X.test, type = "response")                        
       out <- cpm.all(x = yhat, y = y.test, level = levels(y))
       stor[ , , i] <- as.matrix(out$cpms)
       predstor[[i]] <- predict(object = tcpfalist, newdata = X.test, 
                                type = "classify.weights")
       if ((plot.out) && (i == 1)) {plot.mind <- tcpfalist$method}
    }
    mconst <- tcpfalist$const
    rnam <- rownames(out$cpms)
    cnam <- colnames(out$cpms)
    dimnames(stor)[[1]] <- rnam
    dimnames(stor)[[2]] <- cnam
    train.weights <- list(Atrain.weights = Aw, Btrain.weights = Bw, 
                          Ctrain.weights = Cw, Phitrain = Pw)
    mean.tune.param <- Reduce("+", opara) / length(opara)
    if (plot.out == T) {
      ncomps <- length(nfac)
      nmethods <- length(method)
      plotstor <- data.frame(matrix(0, nrow = (ncomps * nmethods * nrep), 
                                    ncol = 14))
      plotcname <- c("method", "nfac", "rep", colnames(stor))
      colnames(plotstor) <- plotcname
      matnum <- ncomps * nmethods
      methnames0 <- sapply(strsplit(rownames(stor), split = '.', fixed = T), 
                           function(x) (x[2]))
      methnames <- gsub('[[:digit:]]+', '', methnames0)
      nfacnames <- as.numeric(gsub(".*?([0-9]+).*", "\\1", rownames(stor)))
      for (i in 1:nrep) {
         indl <- matnum * (i - 1) + 1
         indu <- matnum * i
         plotstor[indl:indu, 1] <- methnames
         plotstor[indl:indu, 2] <- nfacnames
         plotstor[indl:indu, 3] <- i
         plotstor[indl:indu, 4:14] <- stor[,,i]
      }
      toplot <- colnames(plotstor)[plottype]
      for (j in 1:length(plottype)) {
         pformula <- formula(paste0(toplot[j], " ~ ", "method * nfac"))
         boxplot(pformula, data = plotstor, ylim = c(0, 1), 
                 xlab = "Method and Number of Components", na.rm = F,
                 ylab = toupper(toplot[j]), main = "Performance Measure")
      }
    }
    if (type.out == "measures") {
      cpfalist <- list(measure = stor, predweights = predstor,
                       train.weights = train.weights, opt.tune = opara,
                       mean.opt.tune = mean.tune.param, X = x, nfac = nfac,
                       model = model, method = method, const = mconst, 
                       cmode = cmode0)
      class(cpfalist) <- "wrapcpfa"
      return(cpfalist)                                                              
    } else {
      dfun <- c("mean", "median", "sd")
      output <- vector(mode = "list", length = length(dfun))
      for (j in seq_along(dfun)) {
         output[[j]] <- apply(stor, 1:2, 
                             FUN = function(x){return(get(dfun[j])(x, 
                                                                na.rm = T))})
         rownames(output[[j]]) <- rnam
         colnames(output[[j]]) <- cnam
      }
      names(output) <- dfun      
      cpfalist <- list(descriptive = output, predweights = predstor,
                       train.weights = train.weights, opt.tune = opara,
                       mean.opt.tune = mean.tune.param, X = x, nfac = nfac,
                       model = model, method = method, const = mconst, 
                       cmode = cmode0)
      class(cpfalist) <- "wrapcpfa"
      return(cpfalist)
    }
}