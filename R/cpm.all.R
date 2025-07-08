cpm.all <- 
  function(x, y, ...) 
{   
    if (!(is.data.frame(x))) {stop("Input 'x' must be of class 'data.frame'. ")}
    if (!(is.numeric(y))) {stop("Input 'y' must be of class 'numeric'.")}
    if (any(is.na(x))) {stop("Input 'x' cannot contain missing values.")}
    if (any(is.na(y))) {stop("Input 'y' cannot contain missing values.")}
    if (nrow(x) != length(y)) {
      stop("Input 'x' must be same length as input 'y'.") 
    }
    nvar <- ncol(x)
    cmat <- vector('list', nvar)
    cpms <- NULL
    values <- apply(x, 2, function(a){cpm(a, y, ...)})
    for (i in 1:nvar) {
       cmat[[i]] <- values[[i]][[1]]
       cpms <- rbind(cpms, values[[i]][[2]]) 
    }
    names(cmat) <- colnames(x)
    rownames(cpms) <- colnames(x)
    output <- list(cm.list = cmat, cpms = cpms)
    return(output)
}