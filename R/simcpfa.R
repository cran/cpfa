simcpfa <- function(arraydim = NULL, model = "parafac", 
                    nfac = 2, nclass = 2, nreps = 100, onreps = 10, 
                    corresp = c(0.3, -0.3), meanpred = c(0, 0), modes = 3,
                    corrpred = matrix(c(1, 0.2, 0.2, 1), nrow = 2), 
                    pf2num = NULL, Amat = NULL, Bmat = NULL, Cmat = NULL, 
                    Gmat = NULL, technical = list()) {
if (is.null(nfac)) {
  warning("Input 'nfac' was NULL. 'nfac' was set to 2.")
  nfac <- 2
}
if (any(!is.numeric(nfac)) || (length(nfac) != 1)) {
  stop("Input 'nfac' must be a single numeric value.")
} 
if ((is.infinite(nfac)) || (is.na(nfac)) || (is.nan(nfac))) { 
  stop("Input 'nfac' must be a finite number and cannot be NA or NaN.") 
}
if ((nfac < 1) || (nfac != floor(nfac))) { 
  stop("Input 'nfac' must be an integer value of 1 or greater.") 
}
if (is.null(modes)) {
  warning("Input 'modes' was NULL. 'modes' was set to 3.")
  modes <- 3
}
if (any(!is.numeric(modes)) || (length(modes) != 1)) {
  stop("Input 'modes' must be a single numeric value.")
}
if ((is.infinite(modes)) || (is.na(modes)) || (is.nan(modes))) { 
  stop("Input 'modes' must be a finite number and cannot be NA or NaN.") 
}
if ((!(modes %in% c(3, 4))) || (modes != floor(modes))) { 
  stop("Input 'modes' must be an integer value of either 3 or 4.") 
}
if (is.null(arraydim)) {
  if (modes == 3) { 
    warning("Input 'arraydim' was NULL. 'arraydim' was set to c(10, 10, 100).")
    arraydim <- c(10, 10, 100)
  } else {
    warning("Input 'arraydim' was NULL. 'arraydim' was \n
            set to c(10, 10, 10, 100).")
    arraydim <- c(10, 10, 10, 100)
  }
}
if (length(arraydim) != modes) {
  stop("Input 'arraydim' must have length equal to 'modes' (e.g., 3 or 4).")
}
if (!is.numeric(arraydim)) {stop("Input 'arraydim' must be numeric.")}
if (any(is.na(arraydim)) || any(is.nan(arraydim))) { 
  stop("Input 'arraydim' must not contain NA or NaN values.")
}
if (any(is.infinite(arraydim))) {
  stop("Input 'arraydim' must contain only finite numbers.")
}
if (any(arraydim < 2) || any(arraydim != floor(arraydim))) { 
  stop("Input 'arraydim' must contain integer values of 2 or greater.") 
}
if (is.null(nclass)) {
  warning("Input 'nclass' was NULL. 'nclass' was set to 2.")
  nclass <- 2
}
if (any(!is.numeric(nclass)) || (length(nclass) != 1)) {
  stop("Input 'nclass' must be a single numeric value.")
} 
if ((is.infinite(nclass)) || (is.na(nclass)) || (is.nan(nclass))) { 
  stop("Input 'nclass' must be a finite number and cannot be NA or NaN.") 
}
if ((nclass < 2) || (nclass != floor(nclass))) { 
  stop("Input 'nclass' must be an integer value of 2 or greater.") 
}
if (is.null(corresp)) {
  warning("Input 'corresp' was NULL. 'corresp' was set to rep(0.5, nfac).")
  corresp <- rep(0.5, nfac)
}
if (length(corresp) != nfac) {
  stop("Input 'corresp' must have a length equal to 'nfac'.")
}
if (!is.numeric(corresp)) {stop("Input 'corresp' must be numeric.")}
if (any(is.na(corresp) | is.nan(corresp))) {
  stop("Input 'corresp' must not contain NA or NaN values.")
}
if (any((corresp < -1) | (corresp > 1))) {
  stop("Input 'corresp' must contain values between -1 and 1, inclusive.")
}
n <- arraydim[modes]
if (is.null(nreps)) {
  warning("Input 'nreps' was NULL. 'nreps' was set to 100.")
  nreps <- 100
}
if (any(!is.numeric(nreps)) || (length(nreps) != 1)) {
  stop("Input 'nreps' must be a single numeric value.")
} 
if ((is.infinite(nreps)) || (is.na(nreps)) || (is.nan(nreps))) { 
  stop("Input 'nreps' must be a finite number and cannot be NA or NaN.") 
}
if ((nreps < 1) || (nreps != floor(nreps))) { 
  stop("Input 'nreps' must be an integer value of 1 or greater.") 
}
if (is.null(onreps)) {
  warning("Input 'onreps' was NULL. 'onreps' was set to 10.")
  onreps <- 10
}
if (any(!is.numeric(onreps)) || (length(onreps) != 1)) {
  stop("Input 'onreps' must be a single numeric value.")
}   
if ((is.infinite(onreps)) || (is.na(onreps)) || (is.nan(onreps))) { 
  stop("Input 'onreps' must be a finite number and cannot be NA or NaN.") 
}
if ((onreps < 1) || (onreps != floor(onreps))) { 
  stop("Input 'onreps' must be an integer value of 1 or greater.") 
}
if (is.null(meanpred)) {
  meanpred <- rep(0, nfac)
} else { 
  if (length(meanpred) != nfac) {
    stop("Input 'meanpred' must have a length equal to 'nfac' when provided.")
  }
  if (!is.numeric(meanpred)) {stop("Input 'meanpred' must be numeric.")}
  if (any(is.na(meanpred)) || any(is.nan(meanpred))) { 
    stop("Input 'meanpred' must not contain NA or NaN values.")
  }
  if (any(is.infinite(meanpred))) {
    stop("Input 'meanpred' must contain only finite numbers.")
  }
}
if (is.null(model)) {
  warning("Input 'model' was NULL. 'model' was set to 'parafac'.")
  model <- "parafac"
}
models <- c("parafac", "parafac2")
model0 <- sum(tolower(model) %in% models)
if ((model0 == 0) || (model0 > 1L)) {
  stop("Input 'model' is not specified correctly. Input must be only one of \n
       either 'parafac' or 'parafac2'.")
} 
model <- tolower(model)
if (model == "parafac2") {
  if (is.null(pf2num)) {
    pf2num <- rep(c(nfac, (nfac + 1), (nfac + 2)), length.out = n)
  } else {
    if (length(pf2num) != n) {
      stop("When provided, input 'pf2num' must have length equal to the last \n 
           value of input 'arraydim'.")
    }
    if (!is.numeric(pf2num)) {stop("Input 'pf2num' must be numeric.")}
    if (any(is.na(pf2num) | is.nan(pf2num))) {
      stop("Input 'pf2num' must not contain NA or NaN values.")
    }
    if (any((pf2num < 2))) {
      stop("Input 'pf2num' must contain values of 2 or greater.")
    }
    if (any(pf2num != floor(pf2num))) {
      stop("Input 'pf2num' must contain integer values of 2 or greater.")
    }
    if (any(pf2num < nfac)) { 
      stop("Input 'pf2num' must contain integer values greater \n 
           than or equal to 'nfac'.")
    }
  }
} else {
  if (!(is.null(pf2num))) {
    warning("Input 'pf2num' was provided, but model is 'parafac'. 'pf2num' \n
            was ignored.")
  }
}
if (!is.null(corrpred)) {
  if (!is.matrix(corrpred)) {stop("Input 'corrpred' must be a matrix.")}
  if (any(is.na(corrpred)) || any(is.nan(corrpred))) {
    stop("Input 'corrpred' must not contain NA values.") 
  }
  if (!all(dim(corrpred) == c(nfac, nfac))) {
    stop("Input 'corrpred' must be a square matrix with dimensions equal \n 
         to 'nfac'.") 
  }
  if (!all(diag(corrpred) == 1)) {
    stop("Input 'corrpred' must contain diagonal values of 1.") 
  }
  offdiag <- corrpred[!diag(nrow(corrpred))]
  if (any((offdiag < -1) | (offdiag > 1))) { 
    stop("Input 'corrpred' must contain off-diagonal elements between \n 
         -1 and 1, inclusive.") 
  }
  cpredlow <- corrpred[lower.tri(corrpred)]
  cpredhigh <-  t(corrpred)[lower.tri(t(corrpred))]
  tcheck <- !(all(cpredlow == cpredhigh))
  if (tcheck == TRUE) {
    stop("Input 'corrpred' must be a correlation matrix with lower triangular \n
         elements that correspond to its upper triangular elements.")
  }
  evc <- eigen(corrpred, symmetric = TRUE)
  if (any(evc$values <= 1e-12)) { 
    stop("Input 'corrpred' is not positive definite. Provide a positive \n
         definite correlation matrix for 'corrpred'.") 
  }
} else {
  corrpred <- matrix(c(1, 0.2, 0.2, 1), nrow = 2)
  evc <- eigen(corrpred, symmetric = TRUE)
}
if (!(is.list(technical))) {
  stop("Input 'technical' must be of class 'list' when provided.")
}
if (length(technical) != 0) {
  allowedkeys <- c("distA", "distB", "distC", "distG")
  if (!(all(names(technical) %in% allowedkeys))) {
    stop("When provided, input technical can only contain one or more of \n
         the following: 'distA', 'distB', 'distC', or 'distG'.")
  }
  for (nm in names(technical)) {
     innerlist <- technical[[nm]]
     if (!(is.list(innerlist))) {
       stop(sprintf("For '%s', the inner list for the input 'technical' \n
                    is not of class 'list'.", nm))
     }
     if (!(is.character(names(innerlist)[1]))) {
       stop(sprintf("For '%s', the inner list for the input 'technical' \n
                    does not have names that are of class 'character'.", nm))
     }
     if ((length(innerlist) == 0) || (names(innerlist)[1] != "dname")) {
       stop(sprintf("For '%s', the first element of the inner list for the \n
                    input 'technical' must be named 'dname'.", nm))
     }
  }
  finnam <- lapply(technical, function(x) x[1])
  fintech <- lapply(technical, function(x) x[-1])
  finlets <- sapply(names(technical), function(x) substr(x, nchar(x), nchar(x)))
} else {
  finlets <- finnam <- NULL
  fintech <- list()
}
if (!(is.null(Amat))) {
  if (model == "parafac") {
    if (!(is.matrix(Amat))) {
      stop("Input 'Amat' must be a matrix when provided and when \n
           model = 'parafac'.")
    }
    if (nrow(Amat) != arraydim[1]) {
      stop("Input 'Amat', when provided, must have number of rows equal to \n
           the first value in input 'arraydim' when model = 'parafac'.")
    }
    if (ncol(Amat) != nfac) {
      stop("Input 'Amat', when provided, must have number of columns equal to \n
           input 'nfac' when model = 'parafac'.")
    }
  } else if (model == "parafac2") {
    if (!(is.list(Amat))) {
      stop("Input 'Amat', when provided and when model = 'parafac2', must \n
           be a list.")
    }
    if (length(Amat) < 2) {
      stop("Input 'Amat', when provided and when model = 'parafac2', \n 
           must contain at least two elements.")
    }
    for (i in seq_along(Amat)) {
       if (!is.matrix(Amat[[i]])) {
         stop(sprintf("Element %d of 'Amat' is not a matrix. Provide a \n
                      matrix for all elements of 'Amat' when \n
                      model = 'parafac2' and when 'Amat' is provided.", i))
       }
       if (nrow(Amat[[i]]) < 2) {
         stop(sprintf("Element %d of 'Amat' has less than two rows. Provide \n 
                      a matrix with at least two rows for all elements of \n
                      'Amat' when model = 'parafac2' and when 'Amat' is \n 
                      provided.", i))
       }
       if (ncol(Amat[[i]]) != nfac) {
         stop(sprintf("Element %d of 'Amat' must have number of columns \n
                      equal to input 'nfac'. Provide a matrix with number \n
                      of columns equal to 'nfac' for all elements of \n
                      'Amat' when model = 'parafac2' and when 'Amat' is \n 
                      provided.", i))
       }
       if (!(is.numeric(Amat[[i]]))) {
         stop(sprintf("Element %d of 'Amat' is not numeric. Provide a matrix \n 
                      with only real numbers for all elements of 'Amat' when \n
                      model = 'parafac2' and when 'Amat' is provided.", i))
       }
       if (!(all(is.finite(Amat[[i]])))) {
         stop(sprintf("Element %d of 'Amat' is not finite. Provide a matrix \n 
                      with only finite numbers for all elements of 'Amat' \n
                      when model = 'parafac2' and when 'Amat' is provided.", i))
       }
    }
  }
}
if (!(is.null(Bmat))) {
  if (!(is.matrix(Bmat))) {
    stop("Input 'Bmat', when provided, must be a matrix.")
  }
  if (nrow(Bmat) != arraydim[2]) {
    stop("Input 'Bmat', when provided, must have number of rows equal to \n 
         the second value of input 'arraydim'.")
  }
  if (ncol(Bmat) != nfac) {
    stop("Input 'Bmat', when provided, must have number of columns \n
         equal to input 'nfac'.")
  }
  if (!(is.numeric(Bmat))) {
    stop("Input 'Bmat', when provided, must contain only real numbers.")
  }
  if (!(all(is.finite(Bmat)))) {
    stop("Input 'Bmat', when provided, must contain only finite numbers.")
  }
}
if (!(is.null(Cmat))) {
  if (modes == 4) {
    if (!(is.matrix(Cmat))) {
      stop("Input 'Cmat', when provided and when modes = 4, must \n
           be a matrix.")
    }
    if (nrow(Cmat) != arraydim[3]) {
      stop("Input 'Cmat', when provided and when modes = 4, must have number \n
           of rows equal to the third value in 'arraydim'.")
    }
    if (ncol(Cmat) != nfac) {
      stop("Input 'Cmat', when provided and when modes = 4, must have number \n
           of columns equal to input 'nfac'.")
    }
    if (!(is.numeric(Cmat))) {
      stop("Input 'Cmat', when provided, must contain only real numbers.")
    }
    if (!(all(is.finite(Cmat)))) {
      stop("Input 'Cmat', when provided, must contain only finite numbers.")
    }
  } else {
    warning("Input 'Cmat' was provided when modes = 3. 'Cmat' was ignored.")
  }
}
if (model == "parafac") {
  if (!(is.null(Gmat))) {
    warning("Gmat was provided but was ignored because model = 'parafac'.")
  }
} else {
  if (!(is.null(Gmat))) {
    if (!(is.matrix(Gmat))) {
      stop("Input 'Gmat' must be a matrix when provided and \n 
           model = 'parafac2'.")
    }
    if (nrow(Gmat) != nfac) {
      stop("Input 'Gmat', when provided and model = 'parafac2', must have \n 
           number of rows equal to input 'nfac'.")
    }
    if (ncol(Gmat) != nfac) {
      stop("Input 'Gmat', when provided and model = 'parafac2', must have \n 
           number of columns equal to input 'nfac'.")
    }
    if (!(is.numeric(Gmat))) {
      stop("Input 'Gmat', when provided, must contain only real numbers.")
    }
    if (!(all(is.finite(Gmat)))) {
      stop("Input 'Gmat', when provided, must contain only finite numbers.")
    }
  }
}
finalletters <- c()
if (!(is.null(Amat))) {
  finalletters <- c(finalletters, substr("Amat", 1, 1))
}
if (!(is.null(Bmat))) {
  finalletters <- c(finalletters, substr("Bmat", 1, 1))
}
if ((modes == 4) && (!(is.null(Cmat)))) {
  finalletters <- c(finalletters, substr("Cmat", 1, 1))
}
if ((model == "parafac2") && (!(is.null(Gmat)))) {
  finalletters <- c(finalletters, substr("Gmat", 1, 1))
}
if (!(is.null(finlets))) {
  conflict <- intersect(finalletters, finlets)
  if (length(conflict) > 0) {
    warning(sprintf("Mode(s) '%s' had both weights and distribution \n 
                    information provided. Provided weights were used and \n
                    provided distribution information was ignored.", 
                    paste(conflict, collapse = ", ")))
    outdex <- which(finlets %in% conflict)
    distnam <- finnam[-outdex]
    disttech <- fintech[-outdex]
    distmodes <- finlets[-outdex]
  } else {
    distnam <- finnam
    disttech <- fintech
    distmodes <- finlets
  }
} else {
  distnam <- finnam
  disttech <- fintech
  distmodes <- finlets
}
storXout <- storYout <- NULL
stordatout <- Inf
Sigma.sqrt <- evc$vectors %*% diag(sqrt(evc$values)) %*% t(evc$vectors)
mum <- as.matrix(meanpred)
warnflag <- FALSE
sdfact <- 1
for (j in 1:onreps) {
   storY <- storX <- betabest <- NULL
   stordat <- Inf
   Z <- matrix(rnorm(n * nfac), nrow = n, ncol = nfac)
   Xq <- rep(1, n) %*% t(mum) + Z %*% Sigma.sqrt
   for (i in 1:nreps) {
      if (nclass == 2) {
        if (is.null(betabest)) {
          beta <- matrix(runif(nfac, -1, 1), nrow = nfac)
        } else {
          if (runif(1) < 0.05) {
            beta <- matrix(runif(nfac, -1, 1), nrow = nfac)
          } else {
            beta <- matrix(rnorm(nfac, mean = betabest, sd = sdfact), 
                           nrow = nfac)
          }
        }
        linpred <- scale(Xq %*% beta)
        pro <- 1 / (1 + exp(-linpred))
        Y <- rbinom(n, 1, pro)
      } else {
        if (is.null(betabest)) {
          beta <- matrix(runif(nfac * nclass, -1, 1), nrow = nfac)
        } else {
          if (runif(1) < 0.05) {
            beta <- matrix(runif(nfac * nclass, -1, 1), nrow = nfac)
          } else {
            beta <- matrix(rnorm(nfac * nclass, mean = as.vector(betabest), 
                                 sd = sdfact), nrow = nfac)
          }
        }
        linpred <- scale(Xq %*% beta)
        pro <- t(apply(linpred, 1, function(x) exp(x) / sum(exp(x))))
        Y <- apply(pro, 1, function(p) sample(1:nclass, size = 1, prob = p)) - 1
      }
      if (length(unique(Y)) == 1) {
        warnflag <- TRUE
        next
      }
      outcome <- sum(abs(corresp - cor(Xq, Y)))
      simdat <- list(outcome = outcome, Y = Y, Xq = Xq)
      outdat <- simdat$outcome
      if (outdat < stordat) {
        stordat <- outdat
        storY <- simdat$Y
        storX <- simdat$Xq
        betabest <- beta
      }
   }
   sdfact <- sdfact * 0.95
   if (stordat < stordatout) {
     stordatout <- stordat
     storYout <- as.matrix(storY)
     storXout <- storX
   }
}
if (is.null(storXout)) {
  stop("Component weights for classification mode were not simulated. \n
       Consider modifying 'corresp' or 'corrpred'.")
}
if (is.null(storYout)) {
  stop("Class labels were not simulated. Consider modifying 'corresp' \n
       or 'corrpred'.")
}
if (warnflag == TRUE) {
  warning("At least one simulation had zero variance in the outcome variable.")
}
y <- storYout
if (is.null(Bmat)) {
  if ("B" %in% distmodes) {
    dsupply <- sapply(names(disttech), 
                      function(x) substr(x, nchar(x), nchar(x)))
    techindex <- which("B" == dsupply)
    bvals <- distdraw(dname = distnam[techindex]$distB$dname, 
                      n = (arraydim[2] * nfac), modes = "B",
                      params = disttech[techindex]$distB)
  } else {
    bvals <- distdraw(dname = "normal", n = (arraydim[2] * nfac), 
                      modes = "B", params = NULL)
  }
  Bmat <- matrix(bvals, nrow = arraydim[2], ncol = nfac)
}
if ((modes == 3) && (model == "parafac")) {
  if (is.null(Amat)) {
    if ("A" %in% distmodes) {
      dsupply <- sapply(names(disttech), 
                        function(x) substr(x, nchar(x), nchar(x)))
      techindex <- which("A" == dsupply)
      avals <- distdraw(dname = distnam[techindex]$distA$dname, 
                        n = (arraydim[1] * nfac), modes = "A",
                        params = disttech[techindex]$distA)
    } else {
      avals <- distdraw(dname = "normal", n = (arraydim[1] * nfac), 
                        modes = "A", params = NULL)
    }
    Amat <- matrix(avals, nrow = arraydim[1], ncol = nfac)
  }
  Cmat <- storXout
  Xmat <- tcrossprod(Amat, krprod(Cmat, Bmat))
  Xmat <- array(Xmat, dim = arraydim)
  Emat <- array(rnorm(prod(arraydim)), dim = arraydim)
  Emat <- nscale(Emat, 0, ssnew = sumsq(Xmat)) 
  X <- Xmat + Emat
  dataout <- list(X = X, y = y, model = model, Amat = Amat, Bmat = Bmat, 
                  Cmat = Cmat, Emat = Emat)
} else if ((modes == 3) && (model == "parafac2")) {
  if (is.null(Gmat)) {
    if ("G" %in% distmodes) {
      dsupply <- sapply(names(disttech), 
                        function(x) substr(x, nchar(x), nchar(x)))
      techindex <- which("G" == dsupply)
      gvals <- distdraw(dname = distnam[techindex]$distG$dname, 
                        n = (nfac * nfac), modes = "G",
                        params = disttech[techindex]$distG)
    } else {
      gvals <- distdraw(dname = "normal", n = (nfac * nfac), 
                        modes = "G", params = NULL)
    }
    Gmat <- matrix(gvals, nrow = nfac, ncol = nfac)
  }
  Cmat <- storXout
  nDd <- pf2num
  if (is.null(Amat)) { 
    Amat <- Amat0 <- vector("list", n)
    for (Dd in 1:n) {
       if ("A" %in% distmodes) {
         dsupply <- sapply(names(disttech), 
                           function(x) substr(x, nchar(x), nchar(x)))
         techindex <- which("A" == dsupply)
         avals <- distdraw(dname = distnam[techindex]$distA$dname, 
                           n = (nDd[Dd] * nfac), modes = "A",
                           params = disttech[techindex]$distA)
       } else {
         avals <- distdraw(dname = "normal", n = (nDd[Dd] * nfac), 
                           modes = "A", params = NULL)
       }
       Amat[[Dd]] <- Amat0[[Dd]] <- matrix(avals, nrow = nDd[Dd], ncol = nfac)
       Amat[[Dd]] <- svd(Amat[[Dd]], nv = 0)$u %*% Gmat
    }
  }
  Xmat <- Emat <- vector("list", n)
  for (Dd in 1:n) {
     Xmat[[Dd]] <- tcrossprod(Amat[[Dd]] %*% diag(Cmat[Dd, ]), Bmat)
     Emat[[Dd]] <- matrix(rnorm(nDd[Dd] * arraydim[2]), nrow = nDd[Dd], 
                          ncol = arraydim[2])
  }
  Emat <- nscale(Emat, 0, ssnew = sumsq(Xmat))
  X <- mapply("+", Xmat, Emat, SIMPLIFY = FALSE)
  dataout <- list(X = X, y = y, model = model, Gmat = Gmat, Amat = Amat0, 
                  Bmat = Bmat, Cmat = Cmat, Emat = Emat)
} else if ((modes == 4) && (model == "parafac")) { 
  if (is.null(Amat)) {
    if ("A" %in% distmodes) {
      dsupply <- sapply(names(disttech), 
                        function(x) substr(x, nchar(x), nchar(x)))
      techindex <- which("A" == dsupply)
      avals <- distdraw(dname = distnam[techindex]$distA$dname, 
                        n = (arraydim[1] * nfac), modes = "A",
                        params = disttech[techindex]$distA)
    } else {
      avals <- distdraw(dname = "normal", n = (arraydim[1] * nfac), 
                        modes = "A", params = NULL)
    }
    Amat <- matrix(avals, nrow = arraydim[1], ncol = nfac)
  }
  if (is.null(Cmat)) {
    if ("C" %in% distmodes) {
      dsupply <- sapply(names(disttech), 
                        function(x) substr(x, nchar(x), nchar(x)))
      techindex <- which("C" == dsupply)
      cvals <- distdraw(dname = distnam[techindex]$distC$dname, 
                        n = (arraydim[3] * nfac), modes = "C",
                        params = disttech[techindex]$distC)
    } else {
      cvals <- distdraw(dname = "normal", n = (arraydim[3] * nfac), 
                        modes = "C", params = NULL)
    }
    Cmat <- matrix(cvals, nrow = arraydim[3], ncol = nfac)
  }
  Dmat <- storXout
  Xmat <- tcrossprod(Amat, krprod(Dmat, krprod(Cmat, Bmat)))
  Xmat <- array(Xmat, dim = arraydim)
  Emat <- array(rnorm(prod(arraydim)), dim = arraydim)
  Emat <- nscale(Emat, 0, ssnew = sumsq(Xmat))
  X <- Xmat + Emat
  dataout <- list(X = X, y = y, model = model, Amat = Amat, Bmat = Bmat, 
                  Cmat = Cmat, Dmat = Dmat, Emat = Emat)
} else {
  if (is.null(Gmat)) {
    if ("G" %in% distmodes) {
      dsupply <- sapply(names(disttech), 
                        function(x) substr(x, nchar(x), nchar(x)))
      techindex <- which("G" == dsupply)
      gvals <- distdraw(dname = distnam[techindex]$distG$dname, 
                        n = (nfac * nfac), modes = "G",
                        params = disttech[techindex]$distG)
    } else {
      gvals <- distdraw(dname = "normal", n = (nfac * nfac), 
                        modes = "G", params = NULL)
    }
    Gmat <- matrix(gvals, nrow = nfac, ncol = nfac)
  }
  if (is.null(Cmat)) {
    if ("C" %in% distmodes) {
      dsupply <- sapply(names(disttech), 
                        function(x) substr(x, nchar(x), nchar(x)))
      techindex <- which("C" == dsupply)
      cvals <- distdraw(dname = distnam[techindex]$distC$dname, 
                        n = (arraydim[3] * nfac), modes = "C",
                        params = disttech[techindex]$distC)
    } else {
      cvals <- distdraw(dname = "normal", n = (arraydim[3] * nfac), 
                        modes = "C", params = NULL)
    }
    Cmat <- matrix(cvals, nrow = arraydim[3], ncol = nfac)
  }
  Dmat <- storXout
  nDd <- pf2num
  if (is.null(Amat)) { 
    Amat <- Amat0 <- vector("list", n)
    for (Dd in 1:n) {
       if ("A" %in% distmodes) {
         dsupply <- sapply(names(disttech), 
                           function(x) substr(x, nchar(x), nchar(x)))
         techindex <- which("A" == dsupply)
         avals <- distdraw(dname = distnam[techindex]$distA$dname, 
                           n = (nDd[Dd] * nfac), modes = "A",
                           params = disttech[techindex]$distA)
       } else {
         avals <- distdraw(dname = "normal", n = (nDd[Dd] * nfac), 
                           modes = "A", params = NULL)
       }
       Amat[[Dd]] <- Amat0[[Dd]] <- matrix(avals, nrow = nDd[Dd], ncol = nfac)
       Amat[[Dd]] <- svd(Amat[[Dd]], nv = 0)$u %*% Gmat
    }
  }
  X <- Xmat <- Emat <- vector("list", n)
  for (Dd in 1:n) {
     leftMat <- Amat[[Dd]] %*% diag(Dmat[Dd, ])
     Xmat[[Dd]] <- array(tcrossprod(leftMat, krprod(Cmat, Bmat)), 
                         dim = c(nDd[Dd], arraydim[2], arraydim[3]))
     Emat[[Dd]] <- array(rnorm(nDd[Dd] * arraydim[2] * arraydim[3]), 
                         dim = c(nDd[Dd], arraydim[2], arraydim[3]))
     X[[Dd]] <- Xmat[[Dd]] + Emat[[Dd]]
  }
  dataout <- list(X = X, y = y, model = model, Gmat = Gmat, Amat = Amat0, 
                  Bmat = Bmat, Cmat = Cmat, Dmat = Dmat, Emat = Emat)
}
return(dataout)
}