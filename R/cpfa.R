cpfa <-
  function(x, y, nfac = 1, nfolds = 10, foldid = NULL, prior = NULL,
           method = c("PLR", "SVM", "RF", "NN"), 
           family = c("binomial", "multinomial"),
           alpha = NULL, lambda = NULL, cost = NULL, gamma = NULL, 
           ntree = NULL, nodesize = NULL, size = NULL, decay = NULL,
           parallel = FALSE, cl = NULL, verbose = TRUE, cmode = NULL, ...) 
{
    xdim <- dim(x)
    lxdim <- length(xdim)
    if (!((lxdim == 3L) | (lxdim == 4L)))
      stop("Input 'x' must be a 3-way or 4-way array.")
    if (any(is.nan(x)) | any(is.infinite(x))) 
      stop("Input 'x' cannot contain NaN or Inf values.")
    if (any(is.na(x)))
      stop("Input 'x' cannot contain missing values.")
    xold <- x
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
    if (class(y) != "factor")
      stop("Input 'y' must be of class 'factor'.")
    if (!(length(y) == xdim[cmode]))
      stop("Length of 'y' must match number of levels in 
           classification mode/dimension of 'x'.")
    nfac <- sort(nfac)
    lnfac <- length(nfac)
    if (!is.numeric(nfac))
      stop("Input 'nfac' must contain integer(s).")
    if (length(nfolds) != 1 | nfolds < 2 | nfolds != floor(nfolds))
      stop("Input 'nfolds' must be a single whole 
           number equal to or greater than 2.")
    if (nfolds > length(y))
      stop("Input 'nfolds' must be a single whole
           number equal to or less than the number of labels in 'y'.")
    if (is.null(foldid)) {
      foldid <- sample(rep(1:nfolds, length.out = length(y)))
    }
    if (length(foldid) != xdim[cmode]) {
      stop("Input 'foldid' must match number of levels in classification mode")
    }
    if (!is.integer(foldid))
      stop("Input 'foldid' must be of class integer.")
    if (length(unique(foldid)) < 2L)
      stop("Input 'foldid' must contain IDs for two or more folds.")
    if (!all(foldid == floor(foldid)))
      stop("Input 'foldid' must contain integers.")
    methods <- c("PLR", "SVM", "RF", "NN")
    method <- pmatch(toupper(method), methods)
    if (length(method) == 0) {
      method <- c(1)
    }
    if ('1' %in% method) {
      if (is.null(alpha)) {
        alpha <- seq(0, 1, length = 6)
      }
      if (any(alpha < 0L) | any(alpha > 1L))
        stop("Input 'alpha' must contain real numbers 
                  between zero and one (inclusive).")
      if (!is.numeric(alpha))
        stop("Input 'alpha' must be numeric.")
      alpha <- sort(alpha)
    }
    if ('2' %in% method) {
      if (is.null(gamma)) {
        gamma <- c(0, 0.01, 0.1, 1, 10, 100, 1000)
      }
      if (any(gamma < 0L)) 
        stop("Input 'gamma' must contain real numbers equal 
                  to or greater than zero.")
      if (!is.numeric(gamma))
        stop("Input 'gamma' must be numeric.")
      gamma <- sort(gamma)
      if (is.null(cost)) {
        cost <- c(1, 2, 4, 8, 16, 32, 64)
      }
      if (any(cost <= 0L))
        stop("Input 'cost' must be a real number greater than zero.")
      if (!is.numeric(cost))
        stop("Input 'cost' must be numeric.")
      cost <- sort(cost)
      svm.grid <- expand.grid(gamma, cost)
    }
    if ('3' %in% method) {
      if (is.null(ntree)) {
        ntree <- c(100, 200, 400, 600, 800, 1600, 3200)
      } 
      if (any(ntree < 1L))
        stop("Input 'ntree' must be greater than or equal to one.")
      if (!all(ntree == floor(ntree)))
        stop("Input 'ntree' must contain only integers.")
      if (!is.numeric(ntree))
        stop("Input 'ntree' must be numeric.")
      ntree <- sort(ntree)
      if (is.null(nodesize)) {
        nodesize <- c(1, 2, 4, 8, 16, 32, 64)
      }
      if (any(nodesize < 1L))
        stop("Input 'nodesize' must be greater than or equal to one.") 
      if (!all(nodesize == floor(nodesize)))
        stop("Input 'nodesize' must be integer.")
      if (!is.numeric(nodesize))
        stop("Input 'nodesize' must be numeric.")
      nodesize <- sort(nodesize)
      rf.grid <- expand.grid(ntree, nodesize)
    }
    if ('4' %in% method) {
      if (is.null(size)) {
        size = c(1, 2, 4, 8, 16, 32, 64)
      } 
      if (any(size < 0L))
        stop("Input 'size' must be greater than or equal to zero.")
      if (!all(size == floor(size)))
        stop("Input 'size' must contain only integers.")
      if (!is.numeric(size))
        stop("Input 'size' must be numeric.")
      size <- sort(size)
      if (is.null(decay)) {
        decay = c(0.001, 0.01, 0.1, 1, 2, 4, 8, 16)
      }
      if (!is.numeric(decay))
        stop("Input 'decay' must be numeric.")
      decay <- sort(decay)
      nn.grid <- expand.grid(size, decay)
    }
    families <- c("binomial", "multinomial")
    family <- pmatch(tolower(family), families)
    lfam <- length(family)
    if (lfam == 0L) {
      stop("Input 'family' must not be NULL.")
    }
    if (any(is.na(family))) {
      stop("Input 'family' must be a character value, either 'binomial'
           or 'multinomial'.")
    }
    if (lfam == 2L & (length(levels(y)) == 2L)) {
      family <- "binomial"
    } 
    if (lfam == 2L & (length(levels(y)) > 2L)) {
      family <- "multinomial"
    }
    if (lfam == 1L & family == 1L) {
      family <- "binomial"
    }
    if (lfam == 1L & family == 2L) {
      family <- "multinomial"
    }
    if (!is.null(prior)) {
      if (sum(prior) != 1)
        stop("Values within input 'prior' must sum to one.")
      if (!is.numeric(prior)) 
        stop("Input 'prior' must be numeric.")
      if (family == "binomial") {
        if (!(length(prior) == 2L))
          stop("Input 'prior' must contain two values for 
               family of 'binomial'.")
      }
      if (family == "multinomial") {
        if (!(length(prior) >= 3L))
          stop("Input 'prior' must contain three or more values for 
             family of 'multinomial'.")
      }
      prior <- prior
      frac <- as.table(prior)
      if (family == "binomial") {
        names(frac) <- c(0, 1)
        weight <- as.numeric(frac[y])
      }
      if (family == "multinomial") {
        names(frac) <- 1:length(levels(y)) - 1
        weight <- as.numeric(frac[y])
      }
    }
    if (is.null(prior)) {
      frac <- table(y)/length(y)
      prior <- as.numeric(frac)
      if (family == "binomial") {
        names(frac) <- c(0, 1)
        weight <- as.numeric(frac[y])
      }
      if (family == "multinomial") {
        names(frac) <- 1:length(levels(y)) - 1
        weight <- as.numeric(frac[y])
      }
    }
    opt.model <- vector("list", lnfac)
    Aweights <- vector("list", lnfac)
    Bweights <- vector("list", lnfac)
    Bweights <- vector("list", lnfac)
    Cweights <- vector("list", lnfac)
    opt.param <- NULL
    est.time <- NULL
    kcv.error <- NULL
    clstop <- FALSE
    if (parallel == TRUE) {
      if (is.null(cl)) {
        clstop <- TRUE
        cl <- makeCluster(detectCores())
      }
        ce <- clusterEvalQ(cl, library(multiway))
        registerDoParallel(cl)
    }
    if (('1' %in% method) && (nfolds == 2))
      warning("For 'nfolds = 2', method 'PLR' might give warnings.")
    for (w in 1:lnfac) {
       optmodel.new <- vector("list", 4)
       if (verbose == TRUE) {
         cat("nfac =", nfac[w], "method = parafac", fill = TRUE)
       }
       tic <- proc.time()
       pfac <- parafac(X = x, nfac = nfac[w], parallel = parallel, cl = cl, 
                       verbose = verbose, ...)
       toc <- proc.time()
       if (w == 1) {const <- pfac$control$const}
       time.pfac <- toc[3] - tic[3]
       Aweights[[w]] <- pfac$A
       Bweights[[w]] <- pfac$B
       train <- as.matrix(pfac$C)
       if (lxdim == 4L) {
         Cweights[[w]] <- pfac$C
         train <- as.matrix(pfac$D)
       }
       if ('1' %in% method) {
         if (verbose == TRUE) {
           cat("nfac =", nfac[w], "method = plr", fill = TRUE)
         }
           tic <- proc.time()
           plr.results <- kcv.plr(x = train,  y = y,  nfolds = nfolds,  
                                  foldid = foldid, alpha = alpha,
                                  family = family, lambda = lambda, 
                                  type.measure = "class", 
                                  weights = weight, standardize = FALSE, 
                                  parallel = parallel)
           toc <- proc.time()
           time.plr <- toc[3] - tic[3]
           error.plr <- plr.results$error
           plr.id <- plr.results$alpha.id
           plr.opt <- alpha[plr.id]
           optmodel.new[[1]] <- plr.fit <- plr.results$plr.fit
           if (nfolds == 2) {
             lambda.min <- plr.results$lambda.min
           }
           if (nfolds >= 3) {
             lambda.min <- plr.fit$lambda.min
           }
       } 
       else {
           time.plr <- NA
           error.plr <- NA
           plr.opt <- NA
           lambda.min <- NA
           optmodel.new[[1]] <- NULL
       }
       if ('2' %in% method) {
         if (verbose == TRUE) {
           cat("nfac =", nfac[w], "method = svm", fill = TRUE)
         }
           tic <- proc.time()
           svm.results <- kcv.svm(x = train, y = y, nfolds = nfolds,
                                  foldid = foldid, svm.grid = svm.grid, 
                                  class.weights = frac, parallel = parallel, 
                                  probability = TRUE)
           toc <- proc.time()
           time.svm <- toc[3] - tic[3]
           error.svm <- svm.results$error
           svm.id <- svm.results$svm.grid.id
           svm.opt <- svm.grid[svm.id,]
           optmodel.new[[2]] <- svm.fit <- svm.results$svm.fit
       } 
       else {
           time.svm <- NA
           error.svm <- NA
           svm.opt <- data.frame(NA, NA)
           optmodel.new[[2]] <- NULL
       }
       if ('3' %in% method) {
         if (verbose == TRUE) {
           cat("nfac =", nfac[w], "method = rf", fill = TRUE)
         }
           tic <- proc.time()
           rf.results <- kcv.rf(x = train, y = y, nfolds = nfolds,
                                foldid = foldid, rf.grid = rf.grid, 
                                classwt = prior, parallel = parallel)
           toc <- proc.time()
           time.rf <- toc[3] - tic[3]
           error.rf <- rf.results$error
           rf.id <- rf.results$rf.grid.id
           rf.opt <- rf.grid[rf.id,]
           optmodel.new[[3]] <- rf.fit <- rf.results$rf.fit
       } 
       else {
           time.rf <- NA
           error.rf <- NA
           rf.opt <- data.frame(NA, NA)
           optmodel.new[[3]] <- NULL
       }
       if ('4' %in% method) {
         if (verbose == TRUE) {
           cat("nfac =", nfac[w], "method = nn", fill = TRUE)
         }
         tic <- proc.time()
         nn.results <- kcv.nn(x = train, y = y, nfolds = nfolds,
                              foldid = foldid, nn.grid = nn.grid, 
                              weights = weight, parallel = parallel)
         toc <- proc.time()
         time.nn <- toc[3] - tic[3]
         error.nn <- nn.results$error
         nn.id <- nn.results$nn.grid.id
         nn.opt <- nn.grid[nn.id,]
         optmodel.new[[4]] <- nn.fit <- nn.results$nn.fit
       }
       else {
         time.nn <- NA
         error.nn <- NA
         nn.opt <- NA
         nn.opt <- data.frame(NA, NA)
         optmodel.new[[4]] <- NULL
       }
       opt.model[[w]] <- optmodel.new
       optparam.new <- data.frame(nfac = nfac[w], alpha = plr.opt,
                                  lambda = lambda.min, 
                                  gamma = svm.opt[1,1], cost = svm.opt[1,2],
                                  ntree = rf.opt[1,1], nodesize = rf.opt[1,2],
                                  size = nn.opt[1,1], decay = nn.opt[1,2])
       esttime.new <- data.frame(nfac = nfac[w], time.pfac = time.pfac,
                                 time.plr = time.plr, time.svm = time.svm,
                                 time.rf = time.rf, time.nn = time.nn)
       kcv.error.new <- data.frame(nfac = nfac[w], error.plr = error.plr,
                               error.svm = error.svm, error.rf = error.rf,
                               error.nn = error.nn)
       opt.param <- rbind(opt.param, optparam.new)
       est.time <- rbind(est.time, esttime.new)
       kcv.error <- rbind(kcv.error, kcv.error.new)
    } 
    if (clstop == TRUE) {
      stopCluster(cl)
    }
    levels(y) <- names(frac)
    if (lxdim == 3L) {
      cpfalist <- list(opt.model = opt.model, opt.param = opt.param, 
                       kcv.error = kcv.error, est.time = est.time, 
                       method = method, x = xold, y = y, Aweights = Aweights, 
                       Bweights = Bweights, Cweights = NULL, const = const,
                       cmode = cmode, family = family)
    }
    if (lxdim == 4L) {
      cpfalist <- list(opt.model = opt.model, opt.param = opt.param, 
                       kcv.error = kcv.error, est.time = est.time, 
                       method = method, x = xold, y = y, Aweights = Aweights, 
                       Bweights = Bweights, Cweights = Cweights, const = const,
                       cmode = cmode, family = family)
    }
    class(cpfalist) <- "cpfa"
    return(cpfalist)
}