kcv.svm <-
  function(x, y, foldid, svm.grid, nfolds = NULL, scale = TRUE, 
           kernel = "radial", degree = 3, coef0 = 0, nu = 0.5,
           class.weights = NULL, cachesize = 40, tolerance = 0.001,
           shrinking = TRUE, cross = 0, probability = TRUE, fitted = TRUE,
           na.action = na.omit, parallel = FALSE) 
{
    if (!is.factor(y)) {
      y <- factor(y)
    }
    if (length(unique(y)) == 1L) {
      stop("Input 'y' must contain some variation (i.e., cannot contain \n
           only a single type of label).")
    }
    if (is.null(nfolds)) {
      nfolds <- 10
    } else {
      if (!(is.integer(nfolds) || is.numeric(nfolds))) {
        stop("Input 'nfolds' must be of class 'integer' or 'numeric'.")
      }
      if (length(unique(nfolds)) != 1L) {
        stop("Input 'nfolds' must contain only a single value.")
      }
      if ((nfolds %% 1) != 0) {
        stop("Input 'nfolds' must be an integer.")
      }
      if ((nfolds < 2) | (nfolds > length(y))) {
        stop("Input 'nfolds' must be an integer between 2 and the \n
             number of observations, inclusive. Alternatively, the \n
             number of observations cannot be less than 'nfolds'.")
      }
    }
    if (is.null(parallel)) {
      parallel <- FALSE
    } else {
      if (!(is.logical(parallel))) {
        stop("Input 'parallel' must be of class 'logical'.")
      }
      if (length(unique(parallel)) != 1L) {
        stop("Input 'parallel' must contain only a single value.")
      }  
    }
    if (is.null(foldid)) {
      foldid <- sample(rep(1:nfolds, length.out = length(y)))
    }
    if (length(unique(foldid)) != as.integer(nfolds)) {
      stop("Input 'foldid' must contain the number of unique values equal to \n
           input 'nfolds'.")
    }
    grid.row <- nrow(svm.grid)
    cv.svm <- matrix(rep(0, grid.row * nfolds), ncol = nfolds)
    if (parallel == TRUE) {
      cv.svm <- foreach (gg = 1:nfolds, .combine = cbind, 
                         .packages = 'e1071') %dopar% {
                         x.train <- as.matrix(x[which(foldid != gg), ])
                         y.train <- y[which(foldid != gg)]
                         x.test <- as.matrix(x[which(foldid == gg), ])
                         y.test <- y[which(foldid == gg)]
                         stortune <- matrix(rep(0, grid.row), ncol = 1)
                         for (yy in 1:grid.row) {
                            svm.fit <- svm(x.train, y.train, 
                                           gamma = svm.grid[yy, 1], 
                                           coef0 = coef0, 
                                           cost = svm.grid[yy, 2], 
                                           nu = nu, 
                                           class.weights = class.weights,
                                           cachesize = cachesize, 
                                           tolerance = tolerance, 
                                           shrinking = shrinking, 
                                           cross = cross,
                                           probability = probability, 
                                           fitted = fitted, 
                                           na.action = na.action)
                           svm.pred <- predict(svm.fit, x.test, 
                                               type = 'response')
                           stortune[yy, 1] <- 1 - mean(svm.pred == y.test) 
                         }
                         cv.svm[, gg] <- stortune
                }
    } else {
      for (gg in 1:nfolds) {
         x.train <- as.matrix(x[which(foldid != gg), ])
         y.train <- y[which(foldid != gg)]
         x.test <- as.matrix(x[which(foldid == gg), ])
         y.test <- y[which(foldid == gg)]
         stortune <- matrix(rep(0, grid.row), ncol = 1)
         for (yy in 1:grid.row) {
            svm.fit <- svm(x.train, y.train, gamma = svm.grid[yy, 1],
                           coef0 = coef0, cost = svm.grid[yy, 2], nu = nu,
                           class.weights = class.weights, cachesize = cachesize,
                           tolerance = tolerance, shrinking = shrinking,
                           cross = cross, probability = probability,
                           fitted = fitted, na.action = na.action)
            svm.pred <- predict(svm.fit, x.test, type = 'response')
            stortune[yy, 1] <- 1 - mean(svm.pred == y.test) 
         }
         cv.svm[, gg] <- stortune
      }
    }
    svm.mean <- apply(cv.svm, 1, mean)
    minid <- which.min(svm.mean)
    svm.fit.best <- svm(x, y, gamma = svm.grid[minid, 1], coef0 = coef0,
                        cost = svm.grid[minid, 2], nu = nu, 
                        class.weights = class.weights, cachesize = cachesize,
                        tolerance = tolerance, shrinking = shrinking,
                        cross = cross, probability = probability, 
                        fitted = fitted, na.action = na.action)
    return(list(svm.grid.id = minid, svm.fit = svm.fit.best,
                error = svm.mean[minid]))
}  