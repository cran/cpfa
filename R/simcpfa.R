simcpfa <- 
  function(arraydim = NULL, model = "parafac", nfac = 2, nclass = 2, 
           nreps = 100, onreps = 10, corresp = c(0.3, -0.3), meanpred = c(0, 0), 
           modes = 3, corrpred = matrix(c(1, 0.2, 0.2, 1), nrow = 2), 
           pf2num = NULL, Amat = NULL, Bmat = NULL, Cmat = NULL, Dmat = NULL, 
           Gmat = NULL, Emat = NULL, technical = list())
{
if (is.null(nfac)) {
  warning("Input 'nfac' was NULL. 'nfac' was set to 2."); nfac <- 2
}
numcheck(nfac)
if ((nfac < 1) || (nfac != floor(nfac))) { 
  stop("Input 'nfac' must be an integer value of 1 or greater.") 
}
if (is.null(modes)) {
  warning("Input 'modes' was NULL. 'modes' was set to 3.")
  modes <- 3
}
numcheck(modes)
if ((!(modes %in% c(3, 4))) || (modes != floor(modes))) { 
  stop("Input 'modes' must be an integer value of either 3 or 4.") 
}
if (is.null(arraydim)) {
  if (modes == 3) { 
    warning("Input 'arraydim' was NULL. 'arraydim' was set to c(10, 10, 100).")
    arraydim <- c(10, 10, 100)
  } else {
    warning("Input 'arraydim' was NULL. 'arraydim' was set to \n 
            c(10, 10, 10, 100).")
    arraydim <- c(10, 10, 10, 100)
  }
}
if (length(arraydim) != modes) {
  stop("Input 'arraydim' must have length equal to 'modes' (e.g., 3 or 4).")
}
if (!(is.numeric(arraydim))) {stop("Input 'arraydim' must be numeric.")}
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
  warning("Input 'nclass' was NULL. 'nclass' was set to 2."); nclass <- 2
}
numcheck(nclass)
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
if (!(is.numeric(corresp))) {stop("Input 'corresp' must be numeric.")}
if (any(is.na(corresp) | is.nan(corresp))) {
  stop("Input 'corresp' must not contain NA or NaN values.")
}
if (any((corresp < -1) | (corresp > 1))) {
  stop("Input 'corresp' must contain values between -1 and 1, inclusive.")
}
n <- arraydim[modes]
if (is.null(nreps)) {
  warning("Input 'nreps' was NULL. 'nreps' was set to 100."); nreps <- 100
}
numcheck(nreps)
if ((nreps < 1) || (nreps != floor(nreps))) { 
  stop("Input 'nreps' must be an integer value of 1 or greater.") 
}
if (is.null(onreps)) {
  warning("Input 'onreps' was NULL. 'onreps' was set to 10."); onreps <- 10
}
numcheck(onreps)
if ((onreps < 1) || (onreps != floor(onreps))) { 
  stop("Input 'onreps' must be an integer value of 1 or greater.") 
}
if (is.null(meanpred)) {
  meanpred <- rep(0, nfac)
} else { 
  if (length(meanpred) != nfac) {
    stop("Input 'meanpred' must have a length equal to 'nfac' when provided.")
  }
  if (!(is.numeric(meanpred))) {stop("Input 'meanpred' must be numeric.")}
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
if ((!(is.character(model))) || (length(model) != 1L)) {
  stop("Input 'model' must be a single character value of \n 
       'parafac' or 'parafac2'.")
}
if (!(tolower(model) %in% c("parafac", "parafac2"))) {
  stop("Input 'model' is not specified correctly. Input must be only one of \n
       either 'parafac' or 'parafac2'.")
}
model <- tolower(model)
if (model == "parafac2") {
  if (is.null(pf2num)) {
    if (!(is.null(Amat))) {
      if (length(Amat) != n) {
        stop(sprintf("When model = 'parafac2' and pf2num = NULL, 'Amat' \n 
                     must be a list of length %d (the last element of \n 
                     input 'arraydim').", n))
      }
      pf2num <- vapply(Amat, nrow, integer(1))
    } else {
      warning("Input 'pf2num' was NULL. 'pf2num' was set \n 
              to 'rep(c(nfac + 1, nfac + 2, nfac + 3), \n 
              length.out = arraydim[modes])'.")
      pf2num <- rep(c(nfac + 1, nfac + 2, nfac + 3), length.out = n)
    }
  } else {
    if (length(pf2num) != n) {
      stop("When provided, input 'pf2num' must have length equal to the last \n 
           value of input 'arraydim'.")
    }
    if (!(is.numeric(pf2num))) {stop("Input 'pf2num' must be numeric.")}
    if (any(is.na(pf2num) | is.nan(pf2num))) {
      stop("Input 'pf2num' must not contain NA or NaN values.")
    }
    if (any(!(is.finite(pf2num)))) {
      stop("Input 'pf2num' must contain finite values.")
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
if (!(is.null(corrpred))) {
  if (!(is.matrix(corrpred))) {stop("Input 'corrpred' must be a matrix.")}
  if (any(is.na(corrpred)) || any(is.nan(corrpred))) {
    stop("Input 'corrpred' must not contain NA values.") 
  }
  if (!(all(dim(corrpred) == c(nfac, nfac)))) {
    stop("Input 'corrpred' must be a square matrix with dimensions equal \n 
         to 'nfac'.") 
  }
  if (!(all(diag(corrpred) == 1))) {
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
  allowedkeys <- c("distA", "distB", "distC", "distG", "distE")
  if (is.null(names(technical))) {
    stop("Input 'technical', when provided, must contain valid named lists.")
  }
  if (!(all(names(technical) %in% allowedkeys))) {
    stop("Input 'technical', when provided, can only contain one or more of \n
         the following: 'distA', 'distB', 'distC', 'distG' or 'distE'.")
  }
  if (length(technical) != length(unique(names(technical)))) {
    stop("Input 'technical', when provided, can only contain valid lists once.")
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
  finlets <- finnam <- NULL; fintech <- list()
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
    if (!(is.numeric(Amat))) {
      stop("Input 'Amat', when provided, must be numeric.")
    }
  } else {
    if (!(is.list(Amat))) {
      stop("Input 'Amat', when provided and when model = 'parafac2', must \n
           be a list.")
    }
    if (length(Amat) < 2) {
      stop("Input 'Amat', when provided and when model = 'parafac2', \n 
           must contain at least two elements.")
    }
    if (length(Amat) != arraydim[modes]) {
      stop("Input 'Amat', when provided, must have length equal to \n
           the last value in input 'arraydim' when model = 'parafac2'.")
    }
    for (i in seq_along(Amat)) {
       if (!(is.matrix(Amat[[i]]))) {
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
       if (nrow(Amat[[i]]) != pf2num[i]) {
         stop(sprintf("Element %d of 'Amat' has number of rows not equal \n
                      to the corresponding value in 'pf2num'. Provide \n 
                      a matrix in each element of 'Amat' with number of \n
                      rows that matches the corresponding value in input \n
                      'pf2num' when model = 'parafac2' and when 'Amat' is \n 
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
errprov <- FALSE
if (!(is.null(Emat))) {
  errprov <- TRUE
  if (model == "parafac") {
    if (!(is.array(Emat))) {
      stop("Input 'Emat' must be an array when provided and when \n 
           model = 'parafac'.")
    }
    if (!(all(dim(Emat) == arraydim))) {
      stop("Input 'Emat', when provided, must have dimensions equal to those \n
           in input 'arraydim' when model = 'parafac'.")
    }
    if (!(is.numeric(Emat))) {
      stop("Input 'Emat', when provided, must be numeric.")
    }
  } else {
    if (!(is.list(Emat))) {
      stop("Input 'Emat', when provided and when model = 'parafac2', must \n
           be a list.")
    }
    if (length(Emat) < 2) {
      stop("Input 'Emat', when provided and when model = 'parafac2', \n 
           must contain at least two elements.")
    }
    for (i in seq_along(Emat)) {
       if (modes == 3) {
         if (!(is.matrix(Emat[[i]]))) {
           stop(sprintf("Element %d of 'Emat' is not a matrix. Provide a \n
                        matrix for all elements of 'Emat' when \n
                        model = 'parafac2', when 'Emat' is provided, and \n
                        when modes = 3.", i))
         }
         if (nrow(Emat[[i]]) != pf2num[i]) {
           stop(sprintf("Element %d of 'Emat' has number of rows not equal \n
                        to the corresponding value in 'pf2num'. Provide \n 
                        a matrix in each element of 'Emat' with number of \n
                        rows that matches the corresponding value in input \n
                        'pf2num' when model = 'parafac2', when 'Emat' is \n 
                        provided, and when modes = 3.", i))
         }
         if (ncol(Emat[[i]]) != arraydim[2]) {
           stop(sprintf("Element %d of 'Emat' must have number of columns \n
                        equal to the second element of input 'arraydim'. \n 
                        Provide a matrix with number of columns equal to \n 
                        the second element of input 'arraydim' for all \n 
                        elements of 'Emat' when model = 'parafac2', when \n 
                        'Emat' is provided, and when modes = 3.", i))
         }
         if (!(is.numeric(Emat[[i]]))) {
           stop(sprintf("Element %d of 'Emat' is not numeric. Provide a \n
                        matrix with only real numbers for all elements of \n
                        'Emat' when model = 'parafac2', when 'Emat' is \n 
                        provided, and when modes = 3.", i))
         }
         if (!(all(is.finite(Emat[[i]])))) {
           stop(sprintf("Element %d of 'Emat' is not finite. Provide a matrix \n 
                        with only finite numbers for all elements of 'Emat' \n
                        when model = 'parafac2', when 'Emat' is provided, \n 
                        and when modes = 3.", i))
         }
       } else {
         if (!(is.array(Emat[[i]]))) {
           stop(sprintf("Element %d of 'Emat' is not an array. Provide a \n
                        three-way array for all elements of 'Emat' when \n
                        model = 'parafac2', when 'Emat' is provided, and \n
                        when modes = 4.", i))
         }
         if (length(dim(Emat[[i]])) != 3) {
           stop(sprintf("Element %d of 'Emat' must be a three-way array \n
                        when model = 'parafac2', when 'Emat' is provided, and \n
                        when modes = 4.", i))
         }
         if (any(dim(Emat[[i]]) < 2)) {
           stop(sprintf("Element %d of 'Emat' has an array dimension of \n 
                        less than two levels. Provide an array with \n
                        dimensions of at least 2 for all dimensions of \n 
                        'Emat' when model = 'parafac2', when 'Emat' is \n 
                        provided, and when modes = 4.", i))
         } 
         if (!(is.numeric(Emat[[i]]))) {
           stop(sprintf("Element %d of 'Emat' is not numeric. Provide an \n
                        array with only real numbers for all elements of \n
                        'Emat' when model = 'parafac2', when 'Emat' is \n 
                        provided, and when modes = 4.", i))
         }
         if (!(all(is.finite(Emat[[i]])))) {
           stop(sprintf("Element %d of 'Emat' is not finite. Provide an array \n 
                        with only finite numbers for all elements of 'Emat' \n
                        when model = 'parafac2', when 'Emat' is provided, \n 
                        and when modes = 4.", i))
         }
         if (dim(Emat[[i]])[1] != pf2num[i]) {
           stop(sprintf("Element %d of 'Emat' has first dimension number of \n
                        levels not equal to the corresponding value in \n 
                        'pf2num'. Provide an array in each element of 'Emat' \n
                        whose first dimension has number of levels that \n 
                        matches the corresponding value in input 'pf2num' \n 
                        when model = 'parafac2', when 'Emat' is provided, \n 
                        and when modes = 4.", i))
         }
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
cmodesup <- FALSE
if (!(is.null(Cmat))) {
  if (!(is.matrix(Cmat))) {
    stop("Input 'Cmat', when provided, must be a matrix.")
  }
  if (nrow(Cmat) != arraydim[3]) {
    stop("Input 'Cmat', when provided, must have number of rows equal to \n 
         the third value in 'arraydim'.")
  }
  if (ncol(Cmat) != nfac) {
    stop("Input 'Cmat', when provided, must have number of columns equal \n 
         to input 'nfac'.")
  }
  if (!(is.numeric(Cmat))) {
    stop("Input 'Cmat', when provided, must be numeric.")
  }
  if (!(all(is.finite(Cmat)))) {
    stop("Input 'Cmat', when provided, must contain only finite numbers.")
  }
  if (modes == 3) {
    cmodesup <- TRUE
    if (any(apply(Cmat, 2, var) == 0)) {
      stop("Input 'Cmat', when provided as a classification mode weights \n
           matrix, must have a non-zero variance for each column.")
    }
  }
}
if (!(is.null(Dmat))) {
  if (modes == 4) {
    cmodesup <- TRUE
    if (!(is.matrix(Dmat))) {
      stop("Input 'Dmat', when provided, must be a matrix.")
    }
    if (nrow(Dmat) != arraydim[4]) {
      stop("Input 'Dmat', when provided, must have number of rows equal to \n 
         the fourth value in 'arraydim'.")
    }
    if (ncol(Dmat) != nfac) {
      stop("Input 'Dmat', when provided, must have number of columns equal \n 
         to input 'nfac'.")
    }
    if (!(is.numeric(Dmat))) {
      stop("Input 'Dmat', when provided, must be numeric.")
    }
    if (!(all(is.finite(Dmat)))) {
      stop("Input 'Dmat', when provided, must contain only finite numbers.")
    }
    if (any(apply(Dmat, 2, var) == 0)) {
      stop("Input 'Dmat', when provided as a classification mode weights \n
           matrix, must have a non-zero variance for each column.")
    }
  } else {
    Dmat <- NULL
    warning("Input 'Dmat' was provided but was ignored because modes = 3.")
  }
}
if (model == "parafac") {
  if (!(is.null(Gmat))) {
    Gmat <- NULL
    warning("Gmat was provided but was ignored because model = 'parafac'.")
  }
} else {
  if (!(is.null(Gmat))) {
    if (!(is.matrix(Gmat))) {
      stop("Input 'Gmat' must be a matrix when provided and when \n
           model = 'parafac2'.")
    }
    if (nrow(Gmat) != nfac) {
      stop("Input 'Gmat', when provided and when model = 'parafac2', must \n 
           have number of rows equal to input 'nfac'.")
    }
    if (ncol(Gmat) != nfac) {
      stop("Input 'Gmat', when provided and when model = 'parafac2', must \n 
           have number of columns equal to input 'nfac'.")
    }
    if (!(is.numeric(Gmat))) {
      stop("Input 'Gmat', when provided, must be numeric.")
    }
    if (!(all(is.finite(Gmat)))) {
      stop("Input 'Gmat', when provided, must contain only finite numbers.")
    }
  }
}
finalletters <- c()
if (!(is.null(Amat))) {finalletters <- c(finalletters, substr("Amat", 1, 1))}
if (!(is.null(Bmat))) {finalletters <- c(finalletters, substr("Bmat", 1, 1))}
if (!(is.null(Cmat))) {finalletters <- c(finalletters, substr("Cmat", 1, 1))}
if (!(is.null(Dmat))) {finalletters <- c(finalletters, substr("Dmat", 1, 1))}
if ((model == "parafac2") && (!(is.null(Gmat)))) {
  finalletters <- c(finalletters, substr("Gmat", 1, 1))
}
if (!(is.null(Emat))) {finalletters <- c(finalletters, substr("Emat", 1, 1))}
if (!(is.null(finlets))) {
  conflict <- intersect(finalletters, finlets)
  if (length(conflict) > 0) {
    warning(sprintf("Mode(s) '%s' had both weights and distribution \n 
                    information provided. Provided weights were used and \n
                    provided distribution information was ignored.", 
                    paste(conflict, collapse = ", ")))
    outdex <- which(finlets %in% conflict)
    distnam <- finnam[-outdex]; disttech <- fintech[-outdex]
    distmodes <- finlets[-outdex]
  } else {
    distnam <- finnam; disttech <- fintech; distmodes <- finlets
  }
} else {
  distnam <- finnam; disttech <- fintech; distmodes <- finlets
} 
storXout <- storYout <- NULL; stordatout <- Inf; warnflag <- FALSE; sdfact <- 1
Sigma.sqrt <- evc$vectors %*% diag(sqrt(evc$values)) %*% t(evc$vectors) 
mum <- as.matrix(meanpred)
if (cmodesup == TRUE) {
  if (onreps != 1) {
    if (modes == 3) {
      warning("Input 'onreps' was reduced to 1 because a classification mode \n
              weights matrix 'Cmat' was provided.")
    } else {
      warning("Input 'onreps' was reduced to 1 because a classification mode \n
              weights matrix 'Dmat' was provided.")
    }
    onreps <- 1
  }
}
for (j in 1:onreps) {
   storY <- storX <- bbest <- NULL; stordat <- Inf
   if (cmodesup == FALSE) {
     Z <- matrix(rnorm(n * nfac), nrow = n, ncol = nfac)
     Xq <- matrix(1, n, 1) %*% t(mum) + Z %*% Sigma.sqrt
   } else {if (modes == 3) {Xq <- Cmat} else {Xq <- Dmat}}
   for (i in 1:nreps) {
      if (nclass == 2) {
        if (is.null(bbest)) {
          beta <- matrix(runif(nfac, -1, 1), nrow = nfac)
        } else {
          if (runif(1) < 0.05) {
            beta <- matrix(runif(nfac, -1, 1), nrow = nfac)
          } else {
            beta <- matrix(rnorm(nfac, mean = bbest, sd = sdfact), nrow = nfac)
          }
        }
        linpred <- scale(Xq %*% beta); pro <- 1 / (1 + exp(-linpred))
        Y <- rbinom(n, 1, pro)
      } else {
        if (is.null(bbest)) {
          beta <- matrix(runif((nfac * nclass), -1, 1), nrow = nfac)
        } else {
          if (runif(1) < 0.05) {
            beta <- matrix(runif((nfac * nclass), -1, 1), nrow = nfac)
          } else {
            beta <- matrix(rnorm((nfac * nclass), mean = as.vector(bbest), 
                                 sd = sdfact), nrow = nfac)
          }
        }
        linpred <- scale(Xq %*% beta) 
        pro <- t(apply(linpred, 1, function(x) exp(x) / sum(exp(x))))
        Y <- apply(pro, 1, function(p) sample(1:nclass, size = 1, prob = p)) - 1
      }
      if (length(unique(Y)) == 1) {warnflag <- TRUE; next}
      outcome <- sum(abs(corresp - cor(Xq, Y)))
      simdat <- list(outcome = outcome, Y = Y, Xq = Xq)
      outdat <- simdat$outcome
      if (outdat < stordat) {
        stordat <- outdat; storY <- simdat$Y; storX <- simdat$Xq; bbest <- beta
      }
   }
   sdfact <- sdfact * 0.95
   if (stordat < stordatout) {
     stordatout <- stordat; storYout <- as.matrix(storY); storXout <- storX
   }
}
if (warnflag == TRUE) {
  warning("At least one simulation had zero variance in the outcome variable.")
}
if (is.null(storXout)) {
  stop("Component weights for classification mode were not simulated. \n
       Consider modifying 'corresp' or 'corrpred'. Alternatively, consider \n
       increasing the number of levels of the classification mode in \n
       input 'arraydim'.")
}
if (is.null(storYout)) {
  stop("Class labels were not simulated. Consider modifying 'corresp' \n
       or 'corrpred'. Alternatively, consider increasing the number of \n
       levels of the classification mode in input 'arraydim'.")
}
y <- storYout
if (is.null(Bmat)) {
  if ("B" %in% distmodes) {
    dsupply <- sapply(names(disttech), 
                      function(x) substr(x, nchar(x), nchar(x)))
    techindex <- which("B" == dsupply)
    bvals <- distdraw(dname = distnam[[techindex]]$dname, 
                      n = (arraydim[2] * nfac), modes = "B",
                      params = disttech[[techindex]])
  } else {
    bvals <- distdraw(dname = "normal", n = (arraydim[2] * nfac), modes = "B", 
                      params = NULL)
  }
  Bmat <- matrix(bvals, nrow = arraydim[2], ncol = nfac)
}
if ((modes == 3) && (model == "parafac")) {
  if (is.null(Amat)) {
    if ("A" %in% distmodes) {
      dsupply <- sapply(names(disttech), 
                        function(x) substr(x, nchar(x), nchar(x)))
      techindex <- which("A" == dsupply)
      avals <- distdraw(dname = distnam[[techindex]]$dname,
                        n = (arraydim[1] * nfac), modes = "A",
                        params = disttech[[techindex]])
    } else {
      avals <- distdraw(dname = "normal", n = (arraydim[1] * nfac), modes = "A", 
                        params = NULL)
    }
    Amat <- matrix(avals, nrow = arraydim[1], ncol = nfac)
  }
  Cmat <- storXout; Xmat <- tcrossprod(Amat, krprod(Cmat, Bmat))
  Xmat <- array(Xmat, dim = arraydim)
  if (errprov == FALSE) {
    if ("E" %in% distmodes) {
      dsupply <- sapply(names(disttech), 
                        function(x) substr(x, nchar(x), nchar(x)))
      techindex <- which("E" == dsupply)
      evals <- distdraw(dname = distnam[[techindex]]$dname, 
                        n = (prod(arraydim)), modes = "E",
                        params = disttech[[techindex]])
    } else {
      evals <- distdraw(dname = "normal", n = (prod(arraydim)), modes = "E", 
                        params = NULL)
    }
    Emat <- array(evals, dim = arraydim)
  }
  X <- Xmat + Emat
  dataout <- list(X = X, y = y, model = model, Amat = Amat, Bmat = Bmat, 
                  Cmat = Cmat, Emat = Emat)
} else if ((modes == 3) && (model == "parafac2")) {
  if (is.null(Gmat)) {
    if ("G" %in% distmodes) {
      dsupply <- sapply(names(disttech), 
                        function(x) substr(x, nchar(x), nchar(x)))
      techindex <- which("G" == dsupply)
      gvals <- distdraw(dname = distnam[[techindex]]$dname, modes = "G",
                        n = (nfac * nfac), params = disttech[[techindex]])
    } else {
      gvals <- distdraw(dname = "normal", n = (nfac * nfac), modes = "G", 
                        params = NULL)
    }
    Gmat <- matrix(gvals, nrow = nfac, ncol = nfac)
  }
  Cmat <- storXout; nDd <- pf2num
  if (is.null(Amat)) { 
    Amat <- Amat0 <- vector("list", n)
    for (Dd in 1:n) {
       if ("A" %in% distmodes) {
         dsupply <- sapply(names(disttech), 
                           function(x) substr(x, nchar(x), nchar(x)))
         techindex <- which("A" == dsupply)
         avals <- distdraw(dname = distnam[[techindex]]$dname, 
                           n = (nDd[Dd] * nfac), modes = "A",
                           params = disttech[[techindex]])
       } else {
         avals <- distdraw(dname = "normal", n = (nDd[Dd] * nfac), modes = "A", 
                           params = NULL)
       }
       Amat[[Dd]] <- Amat0[[Dd]] <- matrix(avals, nrow = nDd[Dd], ncol = nfac)
       Amat[[Dd]] <- svd(Amat[[Dd]], nv = 0)$u %*% Gmat
    }
  }
  Xmat <- vector("list", n); if (errprov == FALSE) {Emat <- vector("list", n)}
  for (Dd in 1:n) {
     if (nfac == 1) {
       Xmat[[Dd]] <- tcrossprod(Amat[[Dd]] %*% as.matrix(Cmat[Dd, ]), Bmat)
     } else {
       Xmat[[Dd]] <- tcrossprod(Amat[[Dd]] %*% diag(Cmat[Dd, ]), Bmat)
     }
     if (errprov == FALSE) {
       if ("E" %in% distmodes) {
         dsupply <- sapply(names(disttech), 
                           function(x) substr(x, nchar(x), nchar(x)))
         techindex <- which("E" == dsupply)
         evals <- distdraw(dname = distnam[[techindex]]$dname, 
                           n = (nDd[Dd] * arraydim[2]), modes = "E",
                           params = disttech[[techindex]])
       } else {
         evals <- distdraw(dname = "normal", n = (nDd[Dd] * arraydim[2]), 
                           modes = "E", params = NULL)
       }
       Emat[[Dd]] <- matrix(evals, nrow = nDd[Dd], ncol = arraydim[2])
     }
  }
  X <- mapply("+", Xmat, Emat, SIMPLIFY = FALSE)
  dataout <- list(X = X, y = y, model = model, Gmat = Gmat, Amat = Amat0, 
                  Bmat = Bmat, Cmat = Cmat, Emat = Emat)
} else if ((modes == 4) && (model == "parafac")) { 
  if (is.null(Amat)) {
    if ("A" %in% distmodes) {
      dsupply <- sapply(names(disttech), 
                        function(x) substr(x, nchar(x), nchar(x)))
      techindex <- which("A" == dsupply)
      avals <- distdraw(dname = distnam[[techindex]]$dname, 
                        n = (arraydim[1] * nfac), modes = "A",
                        params = disttech[[techindex]])
    } else {
      avals <- distdraw(dname = "normal", n = (arraydim[1] * nfac), modes = "A",
                        params = NULL)
    }
    Amat <- matrix(avals, nrow = arraydim[1], ncol = nfac)
  }
  if (is.null(Cmat)) {
    if ("C" %in% distmodes) {
      dsupply <- sapply(names(disttech), 
                        function(x) substr(x, nchar(x), nchar(x)))
      techindex <- which("C" == dsupply)
      cvals <- distdraw(dname = distnam[[techindex]]$dname, 
                        n = (arraydim[3] * nfac), modes = "C",
                        params = disttech[[techindex]])
    } else {
      cvals <- distdraw(dname = "normal", n = (arraydim[3] * nfac), 
                        modes = "C", params = NULL)
    }
    Cmat <- matrix(cvals, nrow = arraydim[3], ncol = nfac)
  }
  Dmat <- storXout; Xmat <- tcrossprod(Amat, krprod(Dmat, krprod(Cmat, Bmat)))
  Xmat <- array(Xmat, dim = arraydim)
  if (errprov == FALSE) {
    if ("E" %in% distmodes) {
      dsupply <- sapply(names(disttech), 
                        function(x) substr(x, nchar(x), nchar(x)))
      techindex <- which("E" == dsupply)
      evals <- distdraw(dname = distnam[[techindex]]$dname, 
                        n = (prod(arraydim)), modes = "E",
                        params = disttech[[techindex]])
    } else {
      evals <- distdraw(dname = "normal", n = (prod(arraydim)), modes = "E", 
                        params = NULL)
    }
    Emat <- array(evals, dim = arraydim)
  }
  X <- Xmat + Emat
  dataout <- list(X = X, y = y, model = model, Amat = Amat, Bmat = Bmat, 
                  Cmat = Cmat, Dmat = Dmat, Emat = Emat)
} else {
  if (is.null(Gmat)) {
    if ("G" %in% distmodes) {
      dsupply <- sapply(names(disttech), 
                        function(x) substr(x, nchar(x), nchar(x)))
      techindex <- which("G" == dsupply)
      gvals <- distdraw(dname = distnam[[techindex]]$dname, modes = "G",
                        n = (nfac * nfac), params = disttech[[techindex]]) 
    } else {
      gvals <- distdraw(dname = "normal", n = (nfac * nfac), modes = "G",
                        params = NULL)
    }
    Gmat <- matrix(gvals, nrow = nfac, ncol = nfac)
  }
  if (is.null(Cmat)) {
    if ("C" %in% distmodes) {
      dsupply <- sapply(names(disttech), 
                        function(x) substr(x, nchar(x), nchar(x)))
      techindex <- which("C" == dsupply)
      cvals <- distdraw(dname = distnam[[techindex]]$dname,  
                        n = (arraydim[3] * nfac), modes = "C",
                        params = disttech[[techindex]])
    } else {
      cvals <- distdraw(dname = "normal", n = (arraydim[3] * nfac), 
                        modes = "C", params = NULL)
    }
    Cmat <- matrix(cvals, nrow = arraydim[3], ncol = nfac)
  }
  Dmat <- storXout; nDd <- pf2num
  if (is.null(Amat)) { 
    Amat <- Amat0 <- vector("list", n)
    for (Dd in 1:n) {
       if ("A" %in% distmodes) {
         dsupply <- sapply(names(disttech), 
                           function(x) substr(x, nchar(x), nchar(x)))
         techindex <- which("A" == dsupply)
         avals <- distdraw(dname = distnam[[techindex]]$dname, 
                           n = (nDd[Dd] * nfac), modes = "A",
                           params = disttech[[techindex]])
       } else {
         avals <- distdraw(dname = "normal", n = (nDd[Dd] * nfac), modes = "A",
                           params = NULL)
       }
       Amat[[Dd]] <- Amat0[[Dd]] <- matrix(avals, nrow = nDd[Dd], ncol = nfac)
       Amat[[Dd]] <- svd(Amat[[Dd]], nv = 0)$u %*% Gmat
    }
  }
  X <- Xmat <- vector("list", n)
  if (errprov == FALSE) {Emat <- X}
  for (Dd in 1:n) {
     if (nfac == 1) {
       leftMat <- Amat[[Dd]] %*% as.matrix(Dmat[Dd, ])
     } else {
       leftMat <- Amat[[Dd]] %*% diag(Dmat[Dd, ])
     }
     Xmat[[Dd]] <- array(tcrossprod(leftMat, krprod(Cmat, Bmat)), 
                         dim = c(nDd[Dd], arraydim[2], arraydim[3]))
     if (errprov == FALSE) {
       if ("E" %in% distmodes) {
         dsupply <- sapply(names(disttech), 
                           function(x) substr(x, nchar(x), nchar(x)))
         techindex <- which("E" == dsupply)
         evals <- distdraw(dname = distnam[[techindex]]$dname, 
                           n = (nDd[Dd] * arraydim[2] * arraydim[3]), 
                           modes = "E", params = disttech[[techindex]])
       } else {
         evals <- distdraw(dname = "normal", 
                           n = (nDd[Dd] * arraydim[2] * arraydim[3]), 
                           modes = "E", params = NULL)
       }
       Emat[[Dd]] <- array(evals, dim = c(nDd[Dd], arraydim[2], arraydim[3]))
     }
     X[[Dd]] <- Xmat[[Dd]] + Emat[[Dd]]
  }
  dataout <- list(X = X, y = y, model = model, Gmat = Gmat, Amat = Amat0, 
                  Bmat = Bmat, Cmat = Cmat, Dmat = Dmat, Emat = Emat)
}
return(dataout)
}