kcv.rf <-
  function(x, y, foldid = NULL, rf.grid, nfolds = NULL, xtest = NULL, 
           ytest = NULL, classwt = NULL, maxnodes = NULL, importance = FALSE, 
           localImp = FALSE, nPerm = 1, norm.votes = TRUE, do.trace = FALSE, 
           corr.bias = FALSE, keep.inbag = FALSE, parallel = FALSE) 
{
    kcvcheck(y = y, nfolds = nfolds, parallel = parallel, foldid = foldid)
    if (is.logical(unlist(x)[1])) {
      stop("Input 'x' in random forest contains at least one logical value. \n
           Convert logical values to numeric or integer values to proceed.")
    }
    grid.row <- nrow(rf.grid)
    cv.rf <- matrix(rep(0, grid.row * nfolds), ncol = nfolds)
    if (parallel == T) {
      cv.rf <- foreach (gg = 1:nfolds, .combine = cbind, 
                        .packages = 'randomForest') %dopar% {
                        x.train <- as.matrix(x[which(foldid != gg), ])
                        y.train <- y[which(foldid != gg)]
                        x.test <- as.matrix(x[which(foldid == gg), ])
                        y.test <- y[which(foldid == gg)]
                        stortune <- matrix(rep(0, grid.row), ncol = 1)
                        for (yy in 1:grid.row) {
                           rf.fit <- randomForest(x.train, y.train,
                                                  ntree = rf.grid[yy, 1],
                                                  nodesize = rf.grid[yy, 2], 
                                                  xtest = xtest, ytest = ytest, 
                                                  classwt = classwt, 
                                                  maxnodes = maxnodes, 
                                                  importance = importance, 
                                                  localImp = localImp, 
                                                  nPerm = nPerm, 
                                                  norm.votes = norm.votes, 
                                                  do.trace = do.trace, 
                                                  corr.bias = corr.bias, 
                                                  keep.inbag = keep.inbag)
                           rf.pred <- predict(rf.fit, x.test, type = 'class')
                           stortune[yy, 1] <- 1 - mean(rf.pred == y.test) 
                        }
                        cv.rf[, gg] <- stortune
               }
    } else {
      for (gg in 1:nfolds) {
         x.train <- as.matrix(x[which(foldid != gg), ])
         y.train <- y[which(foldid != gg)]
         x.test <- as.matrix(x[which(foldid == gg), ])
         y.test <- y[which(foldid == gg)]
         stortune <- matrix(rep(0, grid.row), ncol = 1)
         for (yy in 1:grid.row) {
            rf.fit <- randomForest(x.train, y.train, ntree = rf.grid[yy, 1],
                                   nodesize = rf.grid[yy, 2], xtest = xtest, 
                                   ytest = ytest, classwt = classwt, 
                                   maxnodes = maxnodes, importance = importance, 
                                   localImp = localImp, nPerm = nPerm, 
                                   norm.votes = norm.votes, do.trace = do.trace, 
                                   corr.bias = corr.bias, 
                                   keep.inbag = keep.inbag)
            rf.pred <- predict(rf.fit, x.test, type = 'class')
            stortune[yy, 1] <- 1 - mean(rf.pred == y.test) 
         }
         cv.rf[, gg] <- stortune
      }
    }
    rf.mean <- apply(cv.rf, 1, mean)
    minid <- which.min(rf.mean)
    rf.fit.best <- randomForest(x, y, ntree = rf.grid[minid, 1],
                                nodesize = rf.grid[minid, 2], xtest = xtest, 
                                ytest = ytest, classwt = classwt, 
                                maxnodes = maxnodes, importance = importance, 
                                localImp = localImp, nPerm = nPerm, 
                                norm.votes = norm.votes, do.trace = do.trace, 
                                corr.bias = corr.bias, keep.inbag = keep.inbag)
    return(list(rf.grid.id = minid, rf.fit = rf.fit.best, 
                error = rf.mean[minid]))
}