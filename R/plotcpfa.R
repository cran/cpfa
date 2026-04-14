plotcpfa <- 
  function(object, cmeasure = "acc", meanvalue = TRUE, supNum = FALSE, 
           cmode = NULL, parallel = FALSE, cl = NULL, scale.remode = NULL, 
           newscales = 1, scale.abmode = NULL, sign.remode = NULL, newsigns = 1, 
           sign.abmode = NULL, ...) 
{
    if (!(inherits(object, "wrapcpfa"))) {
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
    logicheck(meanvalue); logicheck(supNum); logicheck(parallel)
    model <- object$model; const <- object$const
    X <- object$X
    lxdim <- object$lxdim
    flattened <- object$flattened
    if (!(is.null(cmode))) {
      if (model == "parafac") {
        if (!(cmode %in% (1:lxdim))) {
          stop("Input 'cmode' must be 1, 2, or 3 (or 4 if 'x' is four-way).")
        }
        modeval <- 1:lxdim
        mode.re <- c(modeval[-cmode], cmode)
        X <- aperm(X, mode.re)
      } 
      if (model == "pca") {
        if ((flattened == TRUE) && (cmode > 2L)) {
          warning("New input 'cmode' is greater than 2, but original input \n
                  'x' in function 'cpfa' was a 3-way or 4-way array. Because \n 
                  'x' was an array, it was flattened to a 2-way matrix. Thus, \n 
                  new 'cmode' is ignored.")
        } else {
          if (!(cmode %in% (1:lxdim))) {
            stop("Input 'cmode' must be 1 or 2 when model was 'pca'.")
          } else if (cmode == 1L) {
            X <- X
          } else {
            X <- t(X)
          }
        }
      }
    } else {
      if (model == "parafac") {
        cmode <- object$cmode
        modeval <- 1:lxdim
        mode.re <- c(modeval[-cmode], cmode)
        X <- aperm(X, mode.re)
      }
      if ((model == "pca") && (flattened == FALSE)) {
        cmode <- object$cmode
        if (cmode == 2L) {X <- t(X)}
      }
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
      if (nrow(finalvalues) == 1L) {
        methodnfac <- rownames(finalvalues)
      } else {
        methodnfac <- names(which(vals == min(vals, na.rm = TRUE)))
      }
    } else {
      if (nrow(finalvalues) == 1L) {
        methodnfac <- rownames(finalvalues)
      } else {
        methodnfac <- names(which(vals == max(vals, na.rm = TRUE)))
      }
    }
    nfac.opt0 <- as.numeric(gsub("[^0-9]", "", methodnfac))
    nfac.opt <- min(nfac.opt0)
    outstor <- vector("list", 4)
    mapstor <- vector("list", 3)
    if (parallel == TRUE) {
      if (is.null(cl)) {cl <- makeCluster(detectCores())}
      ce <- clusterEvalQ(cl, library(multiway))
      registerDoParallel(cl)
    }
    if (model == "parafac") {
      pfac <- parafac(X = X, nfac = nfac.opt, const = const, cl = cl, 
                      parallel = parallel, ...)
      if ((!(is.null(scale.remode))) && (!(is.null(scale.abmode)))) {
        pfac <- rescale(pfac, mode = scale.remode, newscale = newscales, 
                        absorb = scale.abmode)
      }
      if ((!(is.null(sign.remode))) && (!(is.null(sign.abmode)))) {
        pfac <- resign(pfac, mode = sign.remode, newsign = newsigns, 
                       absorb = sign.abmode)
      }
      mapstor[[1]] <- outstor[[1]] <- pfac$A
      mapstor[[2]] <- outstor[[2]] <- pfac$B   
      outstor[[3]] <- pfac$C
      if (lxdim == 4L) {
        mapstor[[3]] <- pfac$C
        outstor[[4]] <- pfac$D
      }
    } else if (model == "parafac2") {
      pfac <- parafac2(X = X, nfac = nfac.opt, const = const, cl = cl, 
                       parallel = parallel, ...)
      if ((!(is.null(scale.remode))) && (!(is.null(scale.abmode)))) {
        pfac <- rescale(pfac, mode = scale.remode, newscale = newscales, 
                        absorb = scale.abmode)
      }
      if ((!(is.null(sign.remode))) && (!(is.null(sign.abmode)))) {
        pfac <- resign(pfac, mode = sign.remode, newsign = newsigns, 
                       absorb = sign.abmode)
      }
      outstor[[1]] <- pfac$A
      mapstor[[2]] <- outstor[[2]] <- pfac$B   
      outstor[[3]] <- pfac$C
      if (lxdim == 4L) {
        mapstor[[3]] <- pfac$C
        outstor[[4]] <- pfac$D
      }
    } else {
      xcent <- scale(X, center = TRUE, scale = FALSE)
      pcacenter <- attr(xcent, "scaled:center")
      USV <- svd(xcent)
      weights <- USV$v[, 1:nfac.opt, drop = FALSE]
      leftsing <- USV$u[, 1:nfac.opt, drop = FALSE]
      scores <- sweep(leftsing, 2, USV$d[1:nfac.opt], "*")
      if (const == "varimax") {
        if (nfac.opt != 1L) {
          varires <- varimax(weights)
          rotmat <- varires$rotmat
          weights <- unclass(as.matrix(varires$loadings))
          scores <- scores %*% rotmat
        }
      }
      mapstor[[1]] <- outstor[[1]] <- weights
      outstor[[2]] <- scores
    }
    if (parallel == TRUE) {stopCluster(cl)}
    plotlabels <- c("A Weights", "B Weights", "C Weights")
    plotlabels2 <- c(plotlabels, "D Weights")
    palette_colors <- colorRampPalette(c("red", "white", "green"))(50)
    par(mfrow = c(1, 1))
    for (i in 1:3) {
       if (!(is.null(mapstor[[i]]))) {
         mat <- as.matrix(mapstor[[i]])
         mainlab <- plotlabels[i]
         rownames(mat) <- paste0("L", 1:nrow(mat))
         colnames(mat) <- paste0("C", 1:ncol(mat))
         maxabs <- max(abs(mat))
         if (ncol(mat) == 1) {
           zarg <- t(as.matrix(t(mat)[, nrow(mat):1]))
         } else {
           zarg <- as.matrix(t(mat)[, nrow(mat):1])
         }
         image(x = 1:ncol(mat), y = 1:nrow(mat), z = zarg, col = palette_colors, 
               xlab = "Components", ylab = "Levels", axes = FALSE, 
               zlim = c(-maxabs, maxabs), main = mainlab)
         axis(1, at = 1:ncol(mat), labels = colnames(mat))
         axis(2, at = 1:nrow(mat), labels = rev(rownames(mat)))
         box()
         if (!(supNum == TRUE)) {
           for (j in 1:nrow(mat)) {
              for (k in 1:ncol(mat)) {
                 text(x = k, y = (nrow(mat) - j + 1), cex = 0.8,
                      labels = round(mat[j, k], 2))
              }
           }
         }
       }
    }
    names(outstor) <- c("Aweights", "Bweights", "Cweights", "Dweights")
    return(outstor)
}