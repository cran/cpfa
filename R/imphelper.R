imphelper <- 
  function(w, method, opt.model, C.pred, nfac, threshold, storrows, family, 
           y, ytrain, ylab, opt.param, train.weights) 
{
    lmethod <- length(method)
    type <- "response" 
    storfac <- matrix(NA, nrow = storrows, ncol = lmethod)
    colcount <- 1
    if ('1' %in% method) {
      plr.fit <- opt.model[[w]][[1]]
      if (!(is.null(plr.fit$glmnet.fit))) {
        plr.fit.class <- plr.fit$glmnet.fit$classnames
        lambda.min <- plr.fit$lambda.min
      }
      if (is.null(plr.fit$glmnet.fit)) {
        plr.fit.class <- plr.fit$classnames
        lambda.min <- opt.param[which(opt.param$nfac == nfac[w]), ]$lambda
      }
      if ((nfac[w] == 1) || (nfac[w] == 1L)) {
        C.pred.plr <- cbind(C.pred, 0)
      } else {
        C.pred.plr <- C.pred
      }
      if (family == "binomial") {
        vals <- as.numeric(predict(plr.fit, newx = C.pred.plr, type = type, 
                                   s = lambda.min, levels = plr.fit.class))
        storfac[, colcount] <- as.numeric(vals > threshold)
        colcount <- colcount + 1
      }
      if (family == "multinomial") {
        vals <- predict(plr.fit, newx = C.pred.plr, type = type, 
                        s = lambda.min, levels = plr.fit.class)
        storfac[, colcount] <- as.numeric((apply(vals, 1, which.max))) - 1
        colcount <- colcount + 1
      }
    }
    if ('2' %in% method) {
      svm.fit <- opt.model[[w]][[2]]
      if (family == "binomial") {
        svm.prob <- attr(predict(svm.fit, C.pred, type = type,
                                 probability = TRUE), "probabilities")
        svm.prob <- svm.prob[, which(colnames(svm.prob) == "1")]
        storfac[, colcount] <- as.numeric(svm.prob > threshold)
        colcount <- colcount + 1
      }
      if (family == "multinomial") {
        svm.prob <- attr(predict(svm.fit, C.pred, type = type, 
                                 probability = TRUE), "probabilities")
        sord <- cbind(1:ncol(svm.prob), as.numeric(colnames(svm.prob)) + 1)
        svmord <- sord[order(sord[,2]), ]
        svm.prob <- svm.prob[, svmord[, 1]]
        storfac[, colcount] <- as.numeric((apply(svm.prob, 1, which.max))) - 1
        colcount <- colcount + 1
      }
    }
    if ('3' %in% method) {
      rf.fit <- opt.model[[w]][[3]]
      if (family == "binomial") {
        rf.prob <- predict(rf.fit, C.pred, type = "prob")
        rf.prob <- rf.prob[, which(colnames(rf.prob) == "1")]
        storfac[, colcount] <- as.numeric(rf.prob > threshold)
        colcount <- colcount + 1
      }
      if (family == "multinomial") {
        rf.prob <- predict(rf.fit, C.pred, type = "prob")
        storfac[, colcount] <- as.numeric((apply(rf.prob, 1, which.max))) - 1
        colcount <- colcount + 1
      }
    }
    if ('4' %in% method) {
      nn.fit <- opt.model[[w]][[4]]
      if (family == "binomial") {
        nn.prob <- predict(nn.fit, newdata = C.pred, type = "raw")
        nn.prob <- nn.prob[, which(colnames(nn.prob) == "y1")]
        storfac[, colcount] <- as.numeric(nn.prob > threshold)
        colcount <- colcount + 1
      }
      if (family == "multinomial") {
        nn.prob <- predict(nn.fit, newdata = C.pred, type = "raw")
        storfac[, colcount] <- as.numeric((apply(nn.prob, 1, which.max))) - 1
        colcount <- colcount + 1
      }
    }
    if ('5' %in% method) {
      rda.fit <- opt.model[[w]][[5]]
      if (family == "binomial") {
        rda.prob <- predict(rda.fit, x = t(train.weights), 
                            y = as.numeric(ytrain) - 1, 
                            xnew = t(C.pred), type = "posterior")
        colnames(rda.prob) <- c(0, 1)
        rda.prob <- rda.prob[, which(colnames(rda.prob) == "1")]
        storfac[, colcount] <- as.numeric(rda.prob > threshold)
        colcount <- colcount + 1
      }
      if (family == "multinomial") {
        rda.prob <- predict(rda.fit, x = t(train.weights), 
                            y = as.numeric(ytrain) - 1, 
                            xnew = t(C.pred), type = "posterior")
        storfac[, colcount] <- as.numeric((apply(rda.prob, 1, which.max))) - 1
        colcount <- colcount + 1
      }
    }
    if ('6' %in% method) {
      gbm.fit <- opt.model[[w]][[6]]
      xgdata <- xgb.DMatrix(data = C.pred)
      num_classes <- length(unique(y))
      if (family == "binomial") {
        gbm.pred <- predict(gbm.fit, xgdata)
        storfac[, colcount] <- as.numeric(gbm.pred > threshold)       
      }
      if (family == "multinomial") {
        gbm.pred <- matrix(predict(gbm.fit, xgdata), ncol = num_classes)
        storfac[, colcount] <- as.numeric((apply(gbm.pred, 1, which.max))) - 1
      }
    } 
    return(storfac)
}