cpm <- 
  function(x, y, level = NULL, fbeta = NULL, prior = NULL) 
{   
    if (class(x) != class(y)) {
      stop("Class of input 'x' must match class of input 'y'.")
    }
    if (!(class(x) %in% c("numeric", "factor", "integer"))) {
      stop("Input 'x' must be of class 'numeric', 'factor', or 'integer'.")
    }
    if (!(class(y) %in% c("numeric", "factor", "integer"))) {
      stop("Input 'y' must be of class 'numeric', 'factor', or 'integer'.")
    }
    if (any(is.na(x))) {stop("Input 'x' cannot contain missing values.")}
    if (any(is.na(y))) {stop("Input 'y' cannot contain missing values.")}
    if (length(x) != length(y)) {
      stop("Input 'x' must be same length as input 'y'.") 
    }
    if (!(is.null(level))) {
      llev <- length(level)
      if (!(class(level) %in% c("numeric", "integer", "character"))) {
        stop ("Input 'level' must be of class 'numeric', 'integer' \n 
              or 'character'.")
      }
      if (llev == 1) {
        stop ("Input 'level' must contain two or more values when supplied.")
      }
    } else {
      luni <- length(unique(c(x, y)))
      if (luni == 1L) {
        stop("Inputs 'x' and 'y' contain only one unique value. Unclear if \n
              binary or multiclass classification is occurring. Must input \n
              argument 'level' to specify.")
      }
      level <- sort(unique(c(x, y)))
    }
    x <- as.integer(factor(x, levels = level)) - 1
    y <- as.integer(factor(y, levels = level)) - 1
    llev <- length(level)
    cm <- matrix(0, nrow = llev, ncol = llev)
    colnames(cm) <- rownames(cm) <- 1:llev - 1
    for (k in 1:length(x)) {
       roww <- which(rownames(cm) == x[k])
       colum <- which(colnames(cm) == y[k])
       cm[roww, colum] <- cm[roww, colum] + 1
    }
    colnames(cm) <- rownames(cm) <- level
    if (is.null(fbeta)) {fbeta <- 1}
    if (!(is.null(fbeta))) {
      if ((length(fbeta) != 1) || (!(is.numeric(fbeta)))) {
        stop("Input 'fbeta' must be a single real number when inputted and of \n
             class numeric.")
      }
    }
    if (!(is.null(prior))) {
      if (!is.numeric(prior)) {
        stop("Input 'prior' must be of class 'numeric' when provided.")
      }
      if (length(prior) != llev) {
        stop("When inputted, 'prior' must contain as many values as the \n
              number of classes.")
      }
      if ((any(prior < 0)) || (any(prior > 1))) {
        stop("When inputted, 'prior' must only contain values between 0 \n
             and 1, inclusive.")
      }
      if (sum(prior) != 1) {
        stop("When inputted, 'prior' must contain values that sum to 1.")
      }
    }
    if (is.null(prior)) {prior <- rep((1 / llev), llev)}
    acc <- (sum(diag(cm))) / sum(cm)
    err <- 1 - acc
    if (llev == 2) {
      tpr <- cm[2,2] / (cm[2,2] + cm[2,1])
      fnr <- 1 - tpr
      fpr <- cm[1,2] / (cm[1,2] + cm[1,1])
      tnr <- 1 - fpr
      ppv <- cm[2,2] / (cm[2,2] + cm[1,2])
      fdr <- 1 - ppv
      npv <- cm[1,1] / (cm[1,1] + cm[2,1])
      fom <- 1 - npv
    } else {
      tpr.sum <- ppv.sum <- fpr.sum <- npv.sum <- rep(0, llev)
      fom.sum <- fdr.sum <- fnr.sum <- tnr.sum <- tpr.sum
      for (h in 1:llev) {
         tpr.sum[h] <- (cm[h,h] / sum(cm[h, ]))
         ppv.sum[h] <- (cm[h,h] / sum(cm[, h]))
         fpr.sum[h] <- ((sum(cm[,h]) - cm[h,h]) / 
                        (sum(cm[,h]) + sum(diag(cm)) - 2 * (cm[h,h])))
         npv.sum[h] <- ((sum(diag(cm)) - cm[h,h]) / 
                        (sum(diag(cm)) + sum(cm[h,]) - 2 * (cm[h,h])))
         tnr.sum[h] <- 1 - fpr.sum[h]
         fnr.sum[h] <- 1 - tpr.sum[h]
         fdr.sum[h] <- 1 - ppv.sum[h]
         fom.sum[h] <- 1 - npv.sum[h]
      }
      tpr <- sum(prior * tpr.sum)
      ppv <- sum(prior * ppv.sum)
      fpr <- sum(prior * fpr.sum)
      npv <- sum(prior * npv.sum)
      tnr <- sum(prior * tnr.sum)
      fnr <- sum(prior * fnr.sum)
      fdr <- sum(prior * fdr.sum)
      fom <- sum(prior * fom.sum)
    }
    fs <- (1 + (fbeta^2)) * ((ppv * tpr) / (((fbeta^2) * ppv) + tpr))
    class.eval <- data.frame(err = err, acc = acc, tpr = tpr, fpr = fpr, 
                             tnr = tnr, fnr = fnr, ppv = ppv, npv = npv,
                             fdr = fdr, fom = fom, fs = fs)
    output <- list(cm = cm, class.eval = class.eval)
    return(output)
}