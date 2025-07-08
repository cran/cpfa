kcvcheck <-
  function(y, nfolds, parallel, foldid) 
{
    if (!(is.factor(y))) {y <- factor(y)}
    if (length(unique(y)) == 1L) {
      stop("Input 'y' must contain some variation (i.e., cannot contain only \n
           a single type of label).")
    }
    if (is.null(nfolds)) {
      nfolds <- 10
    } else {
      numcheck(nfolds)
      if ((nfolds %% 1) != 0) {stop("Input 'nfolds' must be an integer.")}
      if ((nfolds < 2) || (nfolds > length(y))) {
        stop("Input 'nfolds' must be an integer between 2 and the \n
             number of observations, inclusive. In addition, the \n
             number of observations cannot be less than 'nfolds'.")
      }
    }
    if (is.null(parallel)) {
      parallel <- FALSE
    } else {
      logicheck(parallel)
    }
    if (is.null(foldid)) {
      foldid <- sample(rep(1:nfolds, length.out = length(y)))
    }
    if (length(unique(foldid)) != as.integer(nfolds)) {
      stop("Input 'foldid' must contain the number of unique values equal to \n
           input 'nfolds'.")
    }
}