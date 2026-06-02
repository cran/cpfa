postalign <- 
  function(object) 
{
    if (!(inherits(object, "wrapcpfa"))) {
      stop("Input 'object' must be of class 'wrapcpfa'.")
    }
    if (object$lxdim == 4L) {fourway <- TRUE} else {fourway <- FALSE}
    if (object$model == "parafac2") {isp2 <- TRUE} else {isp2 <- FALSE}
    if (object$model == "pca") {
      stop("Permutation alignment can only be implemented for Parafac or \n 
         Parafac2 models.")
    }
    nreps <- length(object$testIDs)
    splitseq <- 1:nreps
    lnfac <- length(object$nfac)
    pairindex <- combn(nreps, 2)
    rpmat <- NULL
    ftccb <- changestor <- vector(mode = "list", length = lnfac)
    for (ii in 1:lnfac) {
       pmatrix <- matrix(0, nrow = nreps, ncol = nreps)
       cnfac <- object$nfac[ii]
       if (cnfac == 1) {
         rpmat <- c(rpmat, NA)
         next
       }
       Alist <- lapply(object$train.weights$Atrain.weights, `[[`, ii)
       Blist <- lapply(object$train.weights$Btrain.weights, `[[`, ii)
       predw0 <- lapply(object$predweights, `[[`, ii)
       predw <- lapply(predw0, function(x) x[, 1:cnfac, drop = FALSE])
       if ((isp2 == TRUE) && (fourway == TRUE)) {
         Clist <- lapply(object$train.weights$Ctrain.weights, `[[`, ii)
         Philist <- lapply(object$train.weights$Phitrain, `[[`, ii)
         Dlist0 <- lapply(object$train.weights$Classtrain.weights, `[[`, ii)
         Dlist <- lapply(Dlist0, function(x) x[, 1:cnfac, drop = FALSE])
       } else if ((isp2 == TRUE) && (fourway == FALSE)) {
         Philist <- lapply(object$train.weights$Phitrain, `[[`, ii)
         Clist0 <- lapply(object$train.weights$Classtrain.weights, `[[`, ii)
         Clist <- lapply(Clist0, function(x) x[, 1:cnfac, drop = FALSE])
       } else if ((isp2 == FALSE) && (fourway == TRUE)) {
         Clist <- lapply(object$train.weights$Ctrain.weights, `[[`, ii)
         Dlist0 <- lapply(object$train.weights$Classtrain.weights, `[[`, ii)
         Dlist <- lapply(Dlist0, function(x) x[, 1:cnfac, drop = FALSE])
       } else {
         Clist0 <- lapply(object$train.weights$Classtrain.weights, `[[`, ii)
         Clist <- lapply(Clist0, function(x) x[, 1:cnfac, drop = FALSE])
       }
       for (jj in 1:ncol(pairindex)) {
          m1 <- pairindex[1, jj]
          m2 <- pairindex[2, jj]
          if ((isp2 == TRUE) && (fourway == TRUE)) {
            Hmat <- abs(congru(Blist[[m1]], Blist[[m2]])) 
            Hmat <- Hmat * abs(congru(Clist[[m1]], Clist[[m2]]))
            neword <- as.integer(solve_LSAP(Hmat, maximum = TRUE))
          } else if ((isp2 == TRUE) && (fourway == FALSE)) {
            Hmat <- abs(congru(Blist[[m1]], Blist[[m2]])) 
            neword <- as.integer(solve_LSAP(Hmat, maximum = TRUE))
          } else if ((isp2 == FALSE) && (fourway == TRUE)) {
            Hmat <- abs(congru(Alist[[m1]], Alist[[m2]]))
            Hmat <- Hmat * abs(congru(Blist[[m1]], Blist[[m2]]))
            Hmat <- Hmat * abs(congru(Clist[[m1]], Clist[[m2]]))
            neword <- as.integer(solve_LSAP(Hmat, maximum = TRUE))
          } else {
            Hmat <- abs(congru(Alist[[m1]], Alist[[m2]]))
            Hmat <- Hmat * abs(congru(Blist[[m1]], Blist[[m2]]))
            neword <- as.integer(solve_LSAP(Hmat, maximum = TRUE))
          }
          Pmat <- 1 - Hmat
          perform <- sum(Pmat[cbind(seq_along(neword), neword)])
          pmatrix[m1, m2] <- perform
       }
       fmodel <- which.min(apply(pmatrix + t(pmatrix), 1, sum))
       rpmat <- c(rpmat, fmodel)
       if ((isp2 == TRUE) && (fourway == TRUE)) {
         baselist <- list(B = Blist[[fmodel]], C = Clist[[fmodel]])
       } else if ((isp2 == TRUE) && (fourway == FALSE)) {
         baselist <- list(B = Blist[[fmodel]])
       } else if ((isp2 == FALSE) && (fourway == TRUE)) {
         baselist <- list(A = Alist[[fmodel]], B = Blist[[fmodel]], 
                          C = Clist[[fmodel]])
       } else {
         baselist <- list(A = Alist[[fmodel]], B = Blist[[fmodel]])
       }
       toalign <- splitseq[-fmodel]
       TCCB <- changeorder <- matrix(0, nrow = (nreps - 1), ncol = (cnfac + 1))
       for (kk in 1:length(toalign)) {
          cmod <- toalign[kk]
          TCCB[kk, 1] <- changeorder[kk, 1] <- cmod
          if ((isp2 == TRUE) && (fourway == TRUE)) {
            newlist <- list(A = Alist[[cmod]], B = Blist[[cmod]], 
                            C = Clist[[cmod]], Phi = Philist[[cmod]],
                            D = Dlist[[cmod]], predw = predw[[cmod]])
            newlist <- alignmod(baselist = baselist, newlist = newlist, 
                                cnfac = cnfac, isp2 = isp2, fourway = fourway)
            object$train.weights$Atrain.weights[[cmod]][[ii]] <- newlist$A
            object$train.weights$Btrain.weights[[cmod]][[ii]] <- newlist$B
            object$train.weights$Ctrain.weights[[cmod]][[ii]] <- newlist$C
            object$train.weights$Phitrain[[cmod]][[ii]] <- newlist$Phi
            Dsend <- cutback(Dlist0[[cmod]], newlist$D, cnfac) 
            object$train.weights$Classtrain.weights[[cmod]][[ii]] <- Dsend
            predwsend <- cutback(predw0[[cmod]], newlist$predw, cnfac)
            object$predweights[[cmod]][[ii]] <- predwsend
          } else if ((isp2 == TRUE) && (fourway == FALSE)) {
            newlist <- list(A = Alist[[cmod]], B = Blist[[cmod]], 
                            Phi = Philist[[cmod]], C = Clist[[cmod]],
                            predw = predw[[cmod]])
            newlist <- alignmod(baselist = baselist, newlist = newlist, 
                                cnfac = cnfac, isp2 = isp2, fourway = fourway)
            object$train.weights$Atrain.weights[[cmod]][[ii]] <- newlist$A
            object$train.weights$Btrain.weights[[cmod]][[ii]] <- newlist$B
            object$train.weights$Phitrain[[cmod]][[ii]] <- newlist$Phi
            Csend <- cutback(Clist0[[cmod]], newlist$C, cnfac) 
            object$train.weights$Classtrain.weights[[cmod]][[ii]] <- Csend
            predwsend <- cutback(predw0[[cmod]], newlist$predw, cnfac)
            object$predweights[[cmod]][[ii]] <- predwsend
          } else if ((isp2 == FALSE) && (fourway == TRUE)) {
            newlist <- list(A = Alist[[cmod]], B = Blist[[cmod]], 
                            C = Clist[[cmod]], D = Dlist[[cmod]], 
                            predw = predw[[cmod]])
            newlist <- alignmod(baselist = baselist, newlist = newlist, 
                                cnfac = cnfac, isp2 = isp2, fourway = fourway)
            object$train.weights$Atrain.weights[[cmod]][[ii]] <- newlist$A
            object$train.weights$Btrain.weights[[cmod]][[ii]] <- newlist$B
            object$train.weights$Ctrain.weights[[cmod]][[ii]] <- newlist$C
            Dsend <- cutback(Dlist0[[cmod]], newlist$D, cnfac) 
            object$train.weights$Classtrain.weights[[cmod]][[ii]] <- Dsend
            predwsend <- cutback(predw0[[cmod]], newlist$predw, cnfac)
            object$predweights[[cmod]][[ii]] <- predwsend
          } else {
            newlist <- list(A = Alist[[cmod]], B = Blist[[cmod]], 
                            C = Clist[[cmod]], predw = predw[[cmod]])
            newlist <- alignmod(baselist = baselist, newlist = newlist, 
                                cnfac = cnfac, isp2 = isp2, fourway = fourway)
            object$train.weights$Atrain.weights[[cmod]][[ii]] <- newlist$A
            object$train.weights$Btrain.weights[[cmod]][[ii]] <- newlist$B
            Csend <- cutback(Clist0[[cmod]], newlist$C, cnfac) 
            object$train.weights$Classtrain.weights[[cmod]][[ii]] <- Csend
            predwsend <- cutback(predw0[[cmod]], newlist$predw, cnfac)
            object$predweights[[cmod]][[ii]] <- predwsend
          }
          changeorder[kk, 2:(cnfac + 1)] <- newlist$changeorder
          TCCB[kk, 2:(cnfac + 1)] <- newlist$finalB
       }
       changestor[[ii]] <- changeorder
       ftccb[[ii]] <- TCCB
    }
    names(changestor) <- paste0("nfac", object$nfac)
    names(ftccb) <- paste0("nfac", object$nfac)
    rout <- list(object = object, targetmod = rpmat, changeorders = changestor,
                 tccb = ftccb)
    return(rout)
}