kcv.plr <-
  function(x, y, foldid = NULL, alpha, nfolds = NULL, 
           family = c("binomial", "multinomial"), offset = NULL, lambda = NULL, 
           weights = NULL, standardize = FALSE, grouped = TRUE, keep = FALSE, 
           parallel = FALSE, maxit = 1e+06) 
{
    if (!is.factor(y)) {y <- factor(y)} 
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
    lalpha <- length(alpha)
    cvlist <- vector("list", lalpha)
    if (ncol(x) == 1) {x <- cbind(0, x)}
    if (nfolds == 2) {
      if (is.null(weights)) {
        weights <- as.numeric((table(y)[y] / length(y)))
        if (family == "binomial") {
          threshold <- (table(y) / length(y))[1]
        }
      } else {
        if (family == "binomial") {
          threshold <- (table(y) / length(y))[1]
        }
      }
      stor.cvm <- rep(NA, lalpha)
      stor.minlam <- rep(NA, lalpha)
      for (h in 1:lalpha) {
         v.alpha <- alpha[h]
         x1 <- x[foldid == 1, ]
         x2 <- x[foldid == 2, ]
         y1 <- y[foldid == 1]
         y2 <- y[foldid == 2]
         fweight <- cbind(foldid, weights)
         weight1 <- fweight[which(fweight[, 1] == 1), 2]
         weight2 <- fweight[which(fweight[, 1] == 2), 2]
         fit0 <- glmnet(x, y, family = family, weights = weights, 
                        alpha = v.alpha)
         if (is.null(lambda)) {
           lam <- fit0$lambda
           nlam <- length(fit0$lambda)
         }
         if (!is.null(lambda)) {
           lam <- lambda
           nlam <- length(lambda)
         }
         fit1 <- glmnet(x1, y1, family = family, lambda = lam, alpha = v.alpha, 
                        weights = weight1)
         fit2 <- glmnet(x2, y2, family = family, lambda = lam, alpha = v.alpha, 
                        weights = weight2)
         pred1 <- predict(fit2, newx = x1, type = "response")
         pred2 <- predict(fit1, newx = x2, type = "response")
         cvm <- rep(NA, nlam)
         for (i in 1:nlam) {
            if (family == "binomial") {
              class1 <- as.numeric(pred1[, i] > threshold)
              class2 <- as.numeric(pred2[, i] > threshold)
              precvm <- ((1 - mean(class1 == y1)) + (1 - mean(class2 == y2)))
              cvm[i] <- precvm / 2
            }
            if (family == "multinomial") {
              class1 <- as.numeric((apply(pred1[, , i], 1, which.max))) - 1
              class2 <- as.numeric((apply(pred2[, , i], 1, which.max))) - 1
              precvm <- ((1 - mean(class1 == y1)) + (1 - mean(class2 == y2)))
              cvm[i] <- precvm / 2
            }
         }
         minid <- which.min(cvm)
         mincv <- cvm[minid]
         stor.cvm[h] <- mincv
         stor.minlam[h] <- minlam <- fit0$lambda[minid]
      }
      minid <- which.min(stor.cvm)
      minlam <- stor.minlam[minid]
      mincv <- stor.cvm
    } else {
      for (h in 1:lalpha) {
         cvlist[[h]] <- cv.glmnet(x, y, family = family, offset = offset,
                                  weights = weights, type.measure = "class", 
                                  nfolds = nfolds, foldid = foldid, 
                                  alpha = alpha[h], lambda = lambda, 
                                  grouped = grouped, keep = keep, 
                                  parallel = parallel, maxit = maxit, 
                                  standardize = standardize)
      }
      mincv <- sapply(cvlist, function(x) min(x$cvm))
      minid <- which.min(mincv)
    } 
    if (nfolds == 2) {
      return(list(alpha.id = minid, plr.fit = fit0,
                  error = mincv[minid], lambda.min = minlam))
    } else {
      return(list(alpha.id = minid, plr.fit = cvlist[[minid]],
                  error = mincv[minid]))
    }
}