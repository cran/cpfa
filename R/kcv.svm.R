kcv.svm <-
  function(x, y, foldid, svm.grid, nfolds = 10, scale = TRUE,
           type = NULL, kernel = "radial", degree = 3, coef0 = 0, nu = 0.5,
           class.weights = NULL, cachesize = 40, tolerance = 0.001,
           shrinking = TRUE, cross = 0, probability = TRUE, fitted = TRUE,
           na.action = na.omit, parallel = FALSE) 
{
    if (class(y) != "factor") {
      y <- factor(y)
    }
    if (is.null(nfolds)) {
      nfolds <- 10
    }
    if (is.null(foldid)) {
      foldid <- sample(rep(1:nfolds, length.out = length(y)))
    }
    npred <- ncol(x)
    grid.row <- nrow(svm.grid)
    cv.svm <- matrix(rep(0, grid.row * nfolds), ncol = nfolds)
    if (parallel == TRUE) {
      cv.svm <- foreach (gg = 1:nfolds, .combine = cbind, 
                       .packages = 'e1071') %dopar% {
                        x.train <- as.matrix(x[which(foldid != gg),])
                        y.train <- y[which(foldid != gg)]
                        x.test <- as.matrix(x[which(foldid == gg),])
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
    } 
    else {
        for (gg in 1:nfolds) {
           x.train <- as.matrix(x[which(foldid != gg),])
           y.train <- y[which(foldid != gg)]
           x.test <- as.matrix(x[which(foldid == gg),])
           y.test <- y[which(foldid == gg)]
           stortune <- matrix(rep(0, grid.row), ncol = 1)
           for (yy in 1:grid.row) {
              svm.fit <- svm(x.train, y.train, gamma = svm.grid[yy, 1],
                              coef0 = coef0, cost = svm.grid[yy, 2], nu = nu,
                              class.weights = class.weights, 
                              cachesize = cachesize,
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