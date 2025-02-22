kcv.nn <-
  function(x, y, foldid = NULL, nn.grid, nfolds = NULL, weights = NULL,
           linout = FALSE, censored = FALSE, skip = FALSE, rang = 0.5, 
           parallel = FALSE) 
{
    kcvcheck(y = y, nfolds = nfolds, parallel = parallel, foldid = foldid)
    if (is.null(weights)) {weights <- as.numeric((table(y)[y] / length(y)))}
    x <- as.matrix(x)
    weights <- as.matrix(weights)                                               
    grid.row <- nrow(nn.grid)
    cv.nn <- matrix(rep(0, grid.row * nfolds), ncol = nfolds)
    if (parallel == T) {
      cv.nn <- foreach (gg = 1:nfolds, .combine = cbind, 
                        .packages = 'nnet') %dopar% {
                        x.train <- as.matrix(x[which(foldid != gg), ])
                        y.train <- y[which(foldid != gg)]
                        x.test <- as.matrix(x[which(foldid == gg), ])
                        y.test <- as.numeric(y[which(foldid == gg)]) - 1
                        trweights <- as.matrix(weights[which(foldid != gg), ])
                        data.ytrain <- data.frame(y.train)
                        con.ytrain <- model.matrix(~y.train - 1, data.ytrain)
                        stortune <- matrix(rep(0, grid.row), ncol = 1)
                        for (yy in 1:grid.row) {
                           nn.fit <- nnet(x = x.train, y = con.ytrain, 
                                          trace = F, size = nn.grid[yy, 1], 
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
         x.train <- as.matrix(x[which(foldid != gg), ])
         y.train <- y[which(foldid != gg)]
         x.test <- as.matrix(x[which(foldid == gg), ])
         y.test <- as.numeric(y[which(foldid == gg)]) - 1
         trweights <- as.matrix(weights[which(foldid != gg), ])
         data.ytrain <- data.frame(y.train)
         con.ytrain <- model.matrix(~y.train - 1, data.ytrain)
         stortune <- matrix(rep(0, grid.row), ncol = 1)
         for (yy in 1:grid.row) {
            nn.fit <- nnet(x = x.train, y = con.ytrain, trace = F,
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
    nn.fit.best <- nnet(x = x, y = con.y, trace = F, size = nn.grid[minid, 1], 
                        decay = nn.grid[minid, 2], linout = linout, 
                        censored = censored, skip = skip, weights = weights, 
                        rang = rang, MaxNWts = 1e4)
    return(list(nn.grid.id = minid, nn.fit = nn.fit.best, 
                error = nn.mean[minid]))
}