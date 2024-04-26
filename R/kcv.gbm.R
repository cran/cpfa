kcv.gbm <-
  function(x, y, foldid = NULL, gbm.grid, nfolds = NULL, prior = NULL,          
           family = c("binomial", "multinomial"), min_child_weight = 1,
           colsample_bytree = 1, booster = "gbtree", parallel = FALSE) 
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
    if (family == "binomial") {
      objective <- "binary:logistic"
    }
    if (family == "multinomial") {
      objective <- "multi:softprob"
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
      weights <- as.numeric((table(y)[y] / length(y)))
      if (family == "binomial") {
        threshold <- (table(y) / length(y))[1]
      }
    } else {
      weights <- as.numeric((prior * length(y))[y] / length(y))
      if (family == "binomial") {
        threshold <- prior[1]
      }
    }
    y <- as.numeric(y) - 1
    num_classes <- length(unique(y))
    grid.row <- nrow(gbm.grid)
    cv.gbm <- matrix(rep(0, grid.row * nfolds), ncol = nfolds)
    if (parallel == TRUE) {
      cv.gbm <- foreach(gg = 1:nfolds, .combine = cbind, 
                        .packages = 'xgboost') %dopar% {
                        x.train <- as.matrix(x[which(foldid != gg), ])
                        y.train <- y[which(foldid != gg)]
                        trweights <- as.matrix(weights[which(foldid != gg)])
                        xgtrain <- xgb.DMatrix(data = x.train, 
                                               label = y.train, 
                                               weight = trweights)
                        x.test <- as.matrix(x[which(foldid == gg), ])
                        y.test <- y[which(foldid == gg)]
                        xgtest <- xgb.DMatrix(data = x.test)
                        stortune <- matrix(rep(0, grid.row), ncol = 1)
                        for (yy in 1:grid.row) {
                           params <- list(booster = booster,
                                          objective = objective,
                                          eta = gbm.grid[yy, 1],
                                          max_depth = gbm.grid[yy, 2],
                                          min_child_weight = min_child_weight,
                                          subsample = gbm.grid[yy, 3],
                                          colsample_bytree = colsample_bytree)
                           if (family == "multinomial") {
                             params$num_class <- num_classes
                           }                                                    
                           gbm.fit <- xgboost(params = params, data = xgtrain,
                                              nrounds = gbm.grid[yy, 4], 
                                              verbose = 0)
                           if (family == "binomial") {
                             gbm.pred <- predict(gbm.fit, xgtest)
                             gbm.class <- as.numeric(gbm.pred > threshold)
                           } else {
                             gbm.pred <- t(matrix(predict(gbm.fit, xgtest), 
                                                  nrow = num_classes))
                             gbm.class <- as.numeric((apply(gbm.pred, 1, 
                                                            which.max))) - 1
                           }
                           stortune[yy, 1] <- 1 - mean(gbm.class == y.test) 
                        }
                        cv.gbm[, gg] <- stortune
                }
    } else {
      for (gg in 1:nfolds) {
         x.train <- as.matrix(x[which(foldid != gg), ])
         y.train <- y[which(foldid != gg)]
         trweights <- as.matrix(weights[which(foldid != gg)])
         xgtrain <- xgb.DMatrix(data = x.train, 
                                label = y.train, 
                                weight = trweights)
         x.test <- as.matrix(x[which(foldid == gg), ])
         y.test <- y[which(foldid == gg)]
         xgtest <- xgb.DMatrix(data = x.test)
         stortune <- matrix(rep(0, grid.row), ncol = 1)
         for (yy in 1:grid.row) {
            params <- list(booster = booster,
                           objective = objective,
                           eta = gbm.grid[yy, 1],
                           max_depth = gbm.grid[yy, 2],
                           min_child_weight = min_child_weight,
                           subsample = gbm.grid[yy, 3],
                           colsample_bytree = colsample_bytree)
            if (family == "multinomial") {
              params$num_class <- num_classes
            }
            gbm.fit <- xgboost(params = params, data = xgtrain,
                               nrounds = gbm.grid[yy, 4], verbose = 0)
            if (family == "binomial") {
              gbm.pred <- predict(gbm.fit, xgtest)
              gbm.class <- as.numeric(gbm.pred > threshold)
            } else {
              gbm.pred <- t(matrix(predict(gbm.fit, xgtest), 
                                   nrow = num_classes))
              gbm.class <- as.numeric((apply(gbm.pred, 1, which.max))) - 1
            }
            stortune[yy, 1] <- 1 - mean(gbm.class == y.test) 
         }
         cv.gbm[, gg] <- stortune
      }
    }
    gbm.mean <- apply(cv.gbm, 1, mean)
    minid <- which.min(gbm.mean)
    params <- list(booster = booster,
                   objective = objective,
                   eta = gbm.grid[minid, 1],
                   max_depth = gbm.grid[minid, 2],
                   min_child_weight = min_child_weight,
                   subsample = gbm.grid[minid, 3],
                   colsample_bytree = colsample_bytree)
    if (family == "multinomial") {
      params$num_class <- num_classes
    }
    xgdata <- xgb.DMatrix(data = x, label = y, weight = weights)
    gbm.fit.best <- xgboost(params = params, data = xgdata,
                            nrounds = gbm.grid[minid, 4], verbose = 0)
    return(list(gbm.grid.id = minid, gbm.fit = gbm.fit.best, 
                error = gbm.mean[minid]))
}