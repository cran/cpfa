plotcpfa <- 
  function(object, cmeasure = "acc", meanvalue = TRUE, supNum = FALSE, 
           parallel = FALSE, cl = NULL, scale.remode = NULL, newscales = 1, 
           scale.abmode = NULL, sign.remode = NULL, newsigns = 1, 
           sign.abmode = NULL, ...) 
{
    if (!inherits(object, "wrapcpfa")) {
      stop("Input 'object' must be of class 'wrapcpfa'.")
    }
    cmeasures <- c("err", "acc", "tpr", "fpr", "tnr", "fnr", "ppv", "npv", 
                   "fdr", "fom", "fs")
    cmeasure <- tolower(cmeasure)
    ctype <- sum(cmeasures %in% cmeasure)
    if (ctype != 1) {
      stop("Input 'cmeasure' must contain a single acceptable value. See help \n
           file for acceptable values.")
    }
    if (!(is.logical(meanvalue))) {
      stop("Input 'meanvalue' must be a logical value.")
    }
    if (length(meanvalue) != 1L) {
      stop("Input 'meanvalue' must be a single value.")
    }
    if (!(is.logical(supNum))) {
      stop("Input 'supNum' must be a logical value.")
    }
    if (length(supNum) != 1L) {
      stop("Input 'supNum' must be a single value.")
    }
    if (!(is.logical(parallel))) {
      stop("Input 'parallel' must be of class logical.")
    }
    if (length(parallel) != 1L) {
      stop("Input 'parallel' must be a single value.")
    }
    model <- object$model
    method <- object$method
    const <- object$const
    X <- object$X
    cmode <- object$cmode
    lxdim <- length(const)
    if (model == "parafac") {
      modeval <- 1:lxdim
      mode.re <- c(modeval[-cmode], cmode)
      X <- aperm(X, mode.re)
    }
    values <- object$descriptive
    if (meanvalue == TRUE) {
      finalvalues <- values$mean
    } else {
      finalvalues <- values$median
    }
    vals <- finalvalues[, which(colnames(finalvalues) == cmeasure)]
    if (length(vals[!is.na(vals)]) == 0) {
      stop("For selected 'cmeasure', all values are 'NA'. Plots cannot be \n
           created. Consider checking output from function 'cpfa'.")
    }
    lowisbest <- c("err", "fpr", "fnr", "fdr", "fom")
    if (cmeasure %in% lowisbest) {
      methodnfac <- names(which(vals == min(vals, na.rm = TRUE)))
    } else {
      methodnfac <- names(which(vals == max(vals, na.rm = TRUE)))
    }
    nfac.opt0 <- as.numeric(gsub("[^0-9]", "", methodnfac))
    nfac.opt <- min(nfac.opt0)
    mapstor <- vector("list", 3)
    if (parallel == TRUE) {
      if (is.null(cl)) {cl <- makeCluster(detectCores())}
      ce <- clusterEvalQ(cl, library(multiway))
      registerDoParallel(cl)
    }
    if (model == "parafac") {
      pfac <- parafac(X = X, nfac = nfac.opt, const = const, cl = cl, 
                      parallel = parallel)#, ...)
      if ((!(is.null(scale.remode))) && (!(is.null(scale.abmode)))) {
        pfac <- rescale(pfac, mode = scale.remode, newscale = newscales, 
                        absorb = scale.abmode)
      }
      if ((!(is.null(sign.remode))) && (!(is.null(sign.abmode)))) {
        pfac <- resign(pfac, mode = sign.remode, newsign = newsigns, 
                       absorb = sign.abmode)
      }
      mapstor[[1]] <- pfac$A
      mapstor[[2]] <- pfac$B                                               
      if (lxdim == 4L) {mapstor[[3]] <- pfac$C}
    } else {
      pfac <- parafac2(X = X, nfac = nfac.opt, const = const, cl = cl, 
                       parallel = parallel)#, ...)
      if ((!(is.null(scale.remode))) && (!(is.null(scale.abmode)))) {
        pfac <- rescale(pfac, mode = scale.remode, newscale = newscales, 
                        absorb = scale.abmode)
      }
      if ((!(is.null(sign.remode))) && (!(is.null(sign.abmode)))) {
        pfac <- resign(pfac, mode = sign.remode, newsign = newsigns, 
                       absorb = sign.abmode)
      }
      mapstor[[2]] <- pfac$B                                             
      if (lxdim == 4L) {mapstor[[3]] <- pfac$C}
    }
    if (parallel == TRUE) {stopCluster(cl)}
    plotlabels <- c("A Weights", "B Weights", "C Weights")
    palette_colors <- colorRampPalette(c("red", "white", "green"))(50)
    par(mfrow = c(1, 1))
    layout(1)
    for (i in 1:3) {
       if (!(is.null(mapstor[[i]]))) {
         mat <- mapstor[[i]]
         mainlab <- plotlabels[i]
         rownames(mat) <- paste0("L", 1:nrow(mat))
         colnames(mat) <- paste0("C", 1:ncol(mat))
         maxabs <- max(abs(mat))
         image(x = 1:ncol(mat), y = 1:nrow(mat), z = t(mat)[, nrow(mat):1],
               col = palette_colors, xlab = "Components", ylab = "Levels", 
               axes = FALSE, zlim = c(-maxabs, maxabs), main = mainlab)
         axis(1, at = 1:ncol(mat), labels = colnames(mat))
         axis(2, at = 1:nrow(mat), labels = rev(rownames(mat)))
         box()
         if (!(supNum == TRUE)) {
           for (j in 1:nrow(mat)) {
              for (k in 1:ncol(mat)) {
                 text(x = k, y = (nrow(mat) - j + 1), 
                      labels = round(mat[j, k], 2), cex = 0.8)
              }
           }
         }
       }
    }
    return(mapstor)
}