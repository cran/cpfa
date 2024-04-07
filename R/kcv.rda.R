kcv.rda <-
  function(x, y, foldid = NULL, rda.grid, nfolds = NULL, prior = NULL,
           regularization = "S", genelist = FALSE, trace = FALSE, 
           parallel = FALSE) 
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
    if (is.null(prior)) {
      prior <- table(y) / length(y)
    }
    grid.row <- nrow(rda.grid)
    cv.rda <- matrix(rep(0, grid.row * nfolds), ncol = nfolds)
    if (parallel == TRUE) {
      cv.rda <- foreach(gg = 1:nfolds, .combine = cbind, 
                        .packages = 'rda') %dopar% {
                        x.train <- as.matrix(x[which(foldid != gg), ])
                        y.train <- y[which(foldid != gg)]
                        x.test <- as.matrix(x[which(foldid == gg), ])
                        y.test <- y[which(foldid == gg)]
                        stortune <- matrix(rep(0, grid.row), ncol = 1)
                        for (yy in 1:grid.row) {
                           rda.fit <- rda(x = t(x.train), y = y.train, 
                                          prior = prior, 
                                          alpha = rda.grid[yy, 1],
                                          delta = rda.grid[yy, 2],
                                          regularization = regularization,
                                          genelist = genelist, 
                                          trace = trace)
                           rda.pred <- predict(rda.fit, x = t(x.train), 
                                               y = y.train, xnew = t(x.test), 
                                               alpha = rda.grid[yy, 1], 
                                               delta = rda.grid[yy, 2], 
                                               type = "class") - 1
                           stortune[yy, 1] <- 1 - mean(rda.pred == y.test) 
                        }
                        cv.rda[, gg] <- stortune
                }
    } else {
      for (gg in 1:nfolds) {
         x.train <- as.matrix(x[which(foldid != gg), ])
         y.train <- y[which(foldid != gg)]
         x.test <- as.matrix(x[which(foldid == gg), ])
         y.test <- y[which(foldid == gg)]
         stortune <- matrix(rep(0, grid.row), ncol = 1)
         for (yy in 1:grid.row) {
            rda.fit <- rda(x = t(unlist(x.train)), y = unlist(y.train), 
                           prior = prior, alpha = rda.grid[yy, 1], 
                           delta = rda.grid[yy, 2], 
                           regularization = regularization, 
                           genelist = genelist, trace = trace)
            rda.pred <- predict(rda.fit, x = t(x.train), y = y.train, 
                                xnew = t(x.test), alpha = rda.grid[yy, 1], 
                                delta = rda.grid[yy, 2], type = "class") - 1
            stortune[yy, 1] <- 1 - mean(rda.pred == y.test) 
         }
         cv.rda[, gg] <- stortune
      }
    }
    rda.mean <- apply(cv.rda, 1, mean)
    minid <- which.min(rda.mean)
    rda.fit.best <- rda(x = t(x), y = y, prior = prior, 
                        alpha = rda.grid[minid, 1], delta = rda.grid[minid, 2],
                        regularization = regularization, genelist = genelist, 
                        trace = trace)
    return(list(rda.grid.id = minid, rda.fit = rda.fit.best,
                error = rda.mean[minid]))
}