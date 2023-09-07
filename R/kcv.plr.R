kcv.plr <-
  function(x, y, foldid = NULL, alpha, nfolds = 10,
           family = c("binomial", "multinomial"), 
           offset = NULL, lambda = NULL, weights = NULL, standardize = FALSE,
           grouped = TRUE, keep = FALSE, parallel = FALSE, maxit = 1e+06) 
{
    if (!is.factor(y)) {
      y <- factor(y)
    } 
    if (is.null(nfolds)) {
      nfolds <- 10
    }
    if (is.null(foldid)) {
      foldid <- sample(rep(1:nfolds, length.out = length(y)))
    }
    lalpha <- length(alpha)
    yweight <- cbind(y, weights)
    cvlist <- vector("list", lalpha)
    if (ncol(x) == 1) {
      x <- cbind(0, x)
    }
    if (nfolds == 2) {
      if (family == "binomial") {
        threshold <- yweight[which(yweight[,1] == 2), 2]
      }
      stor.cvm <- rep(NA, lalpha)
      stor.minlam <- rep(NA, lalpha)
      for (h in 1:lalpha) {
         v.alpha <- alpha[h]
         x1 <- x[foldid == 1,]
         x2 <- x[foldid == 2,]
         y1 <- y[foldid == 1]
         y2 <- y[foldid == 2]
         fweight <- cbind(foldid, weights)
         weight1 <- fweight[which(fweight[,1] == 1), 2]
         weight2 <- fweight[which(fweight[,1] == 2), 2]
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
         fit1 <- glmnet(x1, y1, family = family, lambda = lam,
                        alpha = v.alpha, weights = weight1)
         fit2 <- glmnet(x2, y2, family = family, lambda = lam,
                        alpha = v.alpha, weights = weight2)
         pred1 <- predict(fit2, newx = x1, type = "response")
         pred2 <- predict(fit1, newx = x2, type = "response")
         cvm <- rep(NA, nlam)
         for (i in 1:nlam) {
           if (family == "binomial") {
             class1 <- as.numeric(pred1[,i] > threshold)
             class2 <- as.numeric(pred2[,i] > threshold)
             cvm[i] <- ((1 - mean(class1 == y1)) + (1 - mean(class2 == y2))) / 2
           }
           if (family == "multinomial") {
             class1 <- as.numeric((apply(pred1[,,i], 1, which.max))) - 1
             class2 <- as.numeric((apply(pred2[,,i], 1, which.max))) - 1
             cvm[i] <- ((1 - mean(class1 == y1)) + (1 - mean(class2 == y2))) / 2
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
    }
    if (nfolds >= 3) {
      for (h in 1:lalpha) {
         cvlist[[h]] <- cv.glmnet(x, y, family = family, offset = offset,
                             weights = weights, type.measure = "class", 
                             nfolds = nfolds, foldid = foldid, alpha = alpha[h], 
                             lambda = lambda, grouped = grouped, keep = keep, 
                             parallel = parallel, maxit = maxit, 
                             standardize = standardize)
      }
      mincv <- sapply(cvlist, function(x) min(x$cvm))
      minid <- which.min(mincv)
    }
    if (nfolds == 2) {
      return(list(alpha.id = minid, plr.fit = fit0,
                  error = mincv[minid], lambda.min = minlam))
    }
    if (nfolds >= 3) {
      return(list(alpha.id = minid, plr.fit = cvlist[[minid]],
             error = mincv[minid]))
    }
}