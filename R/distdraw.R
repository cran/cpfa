distdraw <- 
  function(dname = "normal", n = 1e2, params = NULL, modes = "A") 
{
  if (is.null(n)) {n <- 1e2}
  if (length(n) != 1) {stop("Input 'n' must be a single value.")}
  if ((is.infinite(n)) || (is.na(n)) || (is.nan(n))) {
    stop("Input 'n' must be a finite number and cannot be NA or NaN.")
  }
  if ((!is.numeric(n)) || (n < 1) || (n != floor(n))) { 
    stop("Input 'n' must be an integer value of 1 or greater.") 
  }
  allowmodes <- c("A", "B", "C", "G", "E")
  if (length(modes) != 1) {
    stop("Input 'modes' must be specified as a single letter from: ", 
         paste(allowmodes, collapse = ", "))
  }
  modeselected <- toupper(modes)
  if (!(modeselected %in% allowmodes)) {
    stop("Input 'modes' must be one of the following: ", 
         paste(allowmodes, collapse = ", "))
  }
  modes <- modeselected; dname <- tolower(dname)
  normalparams <- list(mean = 0, sd = 1)
  uniformparams <- list(min = 0, max = 1)
  gammaparams <- list(shape = 1, scale = 1)
  betaparams <- list(shape1 = 1, shape2 = 1)
  binomialparams <- list(size = 1, prob = 0.5)
  poissonparams <- list(lambda = 1)
  exponentialparams <- list(rate = 1)
  geometricparams <- list(prob = 0.5)
  negbinomialparams <- list(size = 1, prob = 0.5)
  hypergeoparams <- list(m = 1, n = 1, k = 1)
  lognormalparams <- list(meanlog = 0, sdlog = 1)
  cauchyparams <- list(location = 0, scale = 1)
  distfuncs <- list(
     normal = function(n, params) do.call(rnorm, c(list(n = n), params)),
     uniform = function(n, params) do.call(runif, c(list(n = n), params)),
     gamma = function(n, params) do.call(rgamma, c(list(n = n), params)),
     beta = function(n, params) do.call(rbeta, c(list(n = n), params)),
     binomial = function(n, params) do.call(rbinom, c(list(n = n), params)),
     poisson = function(n, params) do.call(rpois, c(list(n = n), params)),
     exponential = function(n, params) do.call(rexp, c(list(n = n), params)),
     geometric = function(n, params) do.call(rgeom, c(list(n = n), params)),
     negbinomial = function(n, params) do.call(rnbinom, c(list(n = n), params)),
     hypergeo = function(n, params) do.call(rhyper, c(list(nn = n), params)),
     lognormal = function(n, params) do.call(rlnorm, c(list(n = n), params)),
     cauchy = function(n, params) do.call(rcauchy, c(list(n = n), params)))
  if (!(dname %in% names(distfuncs))) {
    stop(sprintf("For mode '%s', input 'dname' must be the name of a \n
                 valid distribution.", modes))
  }
  if (is.null(params)) {
    params <- switch(dname, normal = normalparams, uniform = uniformparams, 
                     gamma = gammaparams, beta = betaparams,
                     binomial = binomialparams, poisson = poissonparams, 
                     exponential = exponentialparams, 
                     geometric = geometricparams,
                     negbinomial = negbinomialparams, hypergeo = hypergeoparams, 
                     lognormal = lognormalparams, cauchy = cauchyparams)
  } else {
    if ((length(params) > 0) && (!(is.numeric(unlist(params))))) {
      stop(sprintf("For mode '%s', input 'params' must contain numeric \n
                    values.", modes))
    }
    if ((length(params) > 0) && (is.null(names(params)))) {
      stop(sprintf("For mode '%s', input 'params' must be a named list.", 
                   modes))
    }
  }
  allowedparams <- list(normal = c("mean", "sd"), uniform = c("min", "max"), 
                        gamma = c("shape", "scale"), 
                        beta = c("shape1", "shape2"), 
                        binomial = c("size", "prob"), poisson = c("lambda"), 
                        exponential = c("rate"), geometric = c("prob"), 
                        negbinomial = c("size", "prob"), 
                        hypergeo = c("m", "n", "k"), 
                        lognormal = c("meanlog", "sdlog"), 
                        cauchy = c("location", "scale"))
  if (!(all(tolower(names(params)) %in% allowedparams[[dname]]))) {
    invalidnames <- names(params)[!(tolower(names(params)) %in% 
                                      allowedparams[[dname]])]
    stop(sprintf("For mode '%s', invalid parameter name(s) for %s \n
                 distribution: %s", modes, dname, paste(invalidnames, 
                                                           collapse = ", ")))
  }
  names(params) <- tolower(names(params))
  defaultparams <- list(normal = normalparams, uniform = uniformparams, 
                        gamma = gammaparams, beta = betaparams, 
                        binomial = binomialparams, poisson = poissonparams, 
                        exponential = exponentialparams, 
                        geometric = geometricparams, 
                        negbinomial = negbinomialparams, 
                        hypergeo = hypergeoparams, lognormal = lognormalparams, 
                        cauchy = cauchyparams)
  for (param in allowedparams[[dname]]) {
    if ((!(param %in% names(params)))) {
      params[[param]] <- defaultparams[[dname]][[param]]
    }
  }
  switch(dname,
         normal = {
           if (params$sd <= 0) {
             stop(sprintf("For mode '%s', for normal distribution, \n
                          'sd' must be greater than 0.", modes))
           }
         },
         uniform = {
           if (params$min >= params$max) {
             stop(sprintf("For mode '%s', for uniform distribution, 'min' \n
                          must be less than 'max'.", modes))
           }
         },
         gamma = {
           if (params$shape <= 0) {
             stop(sprintf("For mode '%s', for gamma distribution, 'shape' \n 
                          must be greater than 0.", modes))
           }
           if (params$scale <= 0) {
             stop(sprintf("For mode '%s', for gamma distribution, 'scale' \n 
                          must be greater than 0.", modes))
           }
         },
         beta = {
           if (params$shape1 <= 0) {
             stop(sprintf("For mode '%s', for beta distribution, 'shape1' \n 
                          must be greater than 0.", modes))
           }
           if (params$shape2 <= 0) {
             stop(sprintf("For mode '%s', for beta distribution, 'shape2' \n 
                          must be greater than 0.", modes))
           }
         },
         binomial = {
           if (params$size < 0) {
             stop(sprintf("For mode '%s', for binomial distribution, 'size' \n 
                          must be zero or positive.", modes))
           }
           mnumb <- .Machine$double.eps^0.5
           if ((!(abs(params$size - round(params$size)) < mnumb))) {
             stop(sprintf("For mode '%s', for binomial distribution, 'size' \n 
                          must be an integer.", modes))
           }
           if ((params$prob < 0) || (params$prob > 1)) {
             stop(sprintf("For mode '%s', for binomial distribution, 'prob' \n 
                          must be between 0 and 1, inclusive.", modes))
           }
         },
         poisson = {
           if (params$lambda < 0) {
             stop(sprintf("For mode '%s', for poisson distribution, 'lambda' \n 
                          must be greater than or equal to 0.", modes))
           }
         },
         exponential = {
           if (params$rate <= 0) {
             stop(sprintf("For mode '%s', for exponential distribution, \n 
                          'rate' must be greater than 0.", modes))
           }
         },
         geometric = {
           if ((params$prob <= 0) || (params$prob > 1)) {
             stop(sprintf("For mode '%s', for geometric distribution, 'prob' \n 
                          must be greater than 0 and less than or equal to 1", 
                          modes))
           }
         },
         negbinomial = {
           if (params$size <= 0) {
             stop(sprintf("For mode '%s', for negative binomial distribution, \n
                          'size' must be a positive integer.", modes))
           }
           mnumb <- .Machine$double.eps^0.5
           if ((!(abs(params$size - round(params$size)) < mnumb))) {
             stop(sprintf("For mode '%s', for negative binomial distribution, \n 
                          'size' must be an integer.", modes))
           }
           if ((params$prob <= 0) || (params$prob > 1)) {
             stop(sprintf("For mode '%s', for negative binomial distribution, \n 
                          'prob' must be greater than 0 and less than or \n
                          equal to 1.", modes))
           }
         },
         hypergeo = {
           if (params$m < 0) {
             stop(sprintf("For mode '%s', for hypergeometric distribution, \n 
                          'm' must be greater than or equal to 0.", modes))
           }
           if (params$n < 0) {
             stop(sprintf("For mode '%s', for hypergeometric distribution, \n 
                          'n' must be greater than or equal to 0.", modes))
           }
           if (params$k < 0) {
             stop(sprintf("For mode '%s', for hypergeometric distribution, \n 
                          'k' must be greater than or equal to 0.", modes))
           }
           mnumb <- .Machine$double.eps^0.5
           if ((!(abs(params$m - round(params$m)) < mnumb))) {
             stop(sprintf("For mode '%s', for hypergeometric distribution, \n 
                          'm' must be an integer.", modes))
           }
           if ((!(abs(params$n - round(params$n)) < mnumb))) {
             stop(sprintf("For mode '%s', for hypergeometric distribution, \n 
                          'n' must be an integer.", modes))
           }
           if ((!(abs(params$k - round(params$k)) < mnumb))) {
             stop(sprintf("For mode '%s', for hypergeometric distribution, \n 
                          'k' must be an integer.", modes))
           }
           if (params$k > (params$m + params$n)) {
             stop(sprintf("For mode '%s', for hypergeometric distribution, \n 
                          'k' cannot be greater than the sum of 'm' and \n
                          'n'.", modes))
           }
         },
         lognormal = {
           if (params$sdlog <= 0) {
             stop(sprintf("For mode '%s', for lognormal distribution, \n
                          'sdlog' must be greater than 0.", modes))
           }
         },
         cauchy = {
           if (params$scale <= 0) {
             stop(sprintf("For mode '%s', for Cauchy distribution, 'scale' \n 
                          must be greater than 0.", modes))
           }
         })
  result <- distfuncs[[dname]](n, params)
  return(result)
}