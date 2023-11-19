kcv.nn <-
  function(x, y, foldid = NULL, nn.grid, nfolds = NULL, weights = NULL,
           linout = FALSE, censored = FALSE, skip = FALSE, rang = 0.5, 
           parallel = FALSE) 
{
    if (!is.factor(y)) {
      y <- factor(y)
    }
    if (length(unique(y)) == 1L)
      stop("Input 'y' must contain some variation (i.e., cannot contain \n
           only a single type of label).")
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
    if (is.null(weights)) {
      weights <- rep(1, length(y))
    }
    x <- as.matrix(x)
    weights <- as.matrix(weights)
    grid.row <- nrow(nn.grid)
    cv.nn <- matrix(rep(0, grid.row * nfolds), ncol = nfolds)
    if (parallel == TRUE) {
      cv.nn <- foreach (gg = 1:nfolds, .combine = cbind, 
                        .packages = 'nnet') %dopar% {
                        x.train <- as.matrix(x[which(foldid != gg),])
                        y.train <- y[which(foldid != gg)]
                        x.test <- as.matrix(x[which(foldid == gg),])
                        y.test <- as.numeric(y[which(foldid == gg)]) - 1
                        trweights <- as.matrix(weights[which(foldid != gg),])
                        data.ytrain <- data.frame(y.train)
                        con.ytrain <- model.matrix(~y.train - 1, data.ytrain)
                        stortune <- matrix(rep(0, grid.row), ncol = 1)
                        for (yy in 1:grid.row) {
                           nn.fit <- nnet(x = x.train, y = con.ytrain, 
                                          trace = FALSE, size = nn.grid[yy, 1], 
                                          decay = nn.grid[yy, 2], 
                                          linout = linout, weights = trweights,
                                          censored = censored, skip = skip, 
                                          rang = rang, MaxNWts = 10000)
                           nn.pred <- predict(nn.fit, newdata = x.test, 
                                              type = 'raw')
                           y.pred <- as.numeric((apply(nn.pred, 1, 
                                                 which.max))) - 1
                           stortune[yy, 1] <- 1 - mean(y.pred == y.test)
                        }
                        cv.nn[, gg] <- stortune
               }
    } else {
      for (gg in 1:nfolds) {
         x.train <- as.matrix(x[which(foldid != gg),])
         y.train <- y[which(foldid != gg)]
         x.test <- as.matrix(x[which(foldid == gg),])
         y.test <- as.numeric(y[which(foldid == gg)]) - 1
         trweights <- as.matrix(weights[which(foldid != gg),])
         data.ytrain <- data.frame(y.train)
         con.ytrain <- model.matrix(~y.train - 1, data.ytrain)
         stortune <- matrix(rep(0, grid.row), ncol = 1)
         for (yy in 1:grid.row) {
            nn.fit <- nnet(x = x.train, y = con.ytrain, trace = FALSE,
                           size = nn.grid[yy, 1], decay = nn.grid[yy, 2],
                           linout = linout, censored = censored, skip = skip, 
                           weights = trweights, rang = rang, MaxNWts = 10000)
            nn.pred <- predict(nn.fit, newdata = x.test, type = 'raw')
            y.pred <- as.numeric((apply(nn.pred, 1, which.max))) - 1
            stortune[yy, 1] <- 1 - mean(y.pred == y.test) 
         }
         cv.nn[, gg] <- stortune
      }
    }
    nn.mean <- apply(cv.nn, 1, mean)
    minid <- which.min(nn.mean)
    data.y <- data.frame(y)
    con.y <- model.matrix(~y - 1, data.y)
    nn.fit.best <- nnet(x = x, y = con.y, trace = FALSE, 
                        size = nn.grid[minid, 1], decay = nn.grid[minid, 2],
                        linout = linout, censored = censored, skip = skip, 
                        weights = weights, rang = rang, MaxNWts = 10000)
    return(list(nn.grid.id = minid, nn.fit = nn.fit.best, 
                error = nn.mean[minid]))
}