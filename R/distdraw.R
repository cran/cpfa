distdraw <- function(dname = "normal", n = 1e2, params = NULL, 
                     modes = "A") {
  if (is.null(n)) {
    n <- 100
  }
  if (length(n) != 1) {stop("Input 'n' must be a single value.")}
  if ((is.infinite(n)) || (is.na(n)) || (is.nan(n))) {
    stop("Input 'n' must be a finite number and cannot be NA or NaN.")
  }
  if ((!is.numeric(n)) || (n < 2) || (n != floor(n))) { 
    stop("Input 'n' must be an integer value of 2 or greater.") 
  }
  allowmodes <- c("A", "B", "C", "G")
  if (length(modes) != 1) {
    stop("Input 'modes' must be specified as a single letter from: ", 
         paste(allowmodes, collapse = ", "))
  }
  modeselected <- toupper(modes)
  if (!modeselected %in% allowmodes) {
    stop("Input 'modes' must be one of the following: ", 
         paste(allowmodes, collapse = ", "))
  }
  modes <- modeselected
  normalparams <- list(mean = 0, sd = 1)
  uniformparams <- list(min = 0, max = 1)
  gammaparams <- list(shape = 1, scale = 1)
  betaparams <- list(shape1 = 1, shape2 = 1)
  binomialparams <- list(size = 1, prob = 0.5)
  poissonparams <- list(lambda = 1)
  exponentialparams <- list(rate = 1)
  dname <- tolower(dname)
  distfuncs <- list(
        normal = function(n, params) do.call(rnorm, c(list(n = n), params)),
        uniform = function(n, params) do.call(runif, c(list(n = n), params)),
        gamma = function(n, params) do.call(rgamma, c(list(n = n), params)),
        beta = function(n, params) do.call(rbeta, c(list(n = n), params)),
        binomial = function(n, params) do.call(rbinom, c(list(n = n), params)),
        poisson = function(n, params) do.call(rpois, c(list(n = n), params)),
        exponential = function(n, params) do.call(rexp, c(list(n = n), params)))
  if (!(dname %in% names(distfuncs))) {
    stop(sprintf("For mode '%s', input 'dname' must be the name of a \n
                 valid distribution.", modes))
  }
  if (is.null(params)) {
    params <- switch(dname, normal = normalparams, uniform = uniformparams, 
                     gamma = gammaparams, beta = betaparams,
                     binomial = binomialparams, poisson = poissonparams, 
                     exponential = exponentialparams)
  } else {
    if (!(is.numeric(unlist(params)))) {
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
                        exponential = c("rate"))
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
                        exponential = exponentialparams)
  for (param in allowedparams[[dname]]) {
    if ((!(param %in% names(params)))) {
      params[[param]] <- defaultparams[[dname]][[param]]
    }
  }
  switch(dname,
         normal = {
           if ((!(is.null(params$sd))) && (params$sd <= 0)) {
             stop(sprintf("For mode '%s', for normal distribution, \n
                          'sd' must be greater than 0.", modes))
           }
         },
         uniform = {
           if ((!(is.null(params$min))) && (!(is.null(params$max))) && 
               (params$min >= params$max)) {
             stop(sprintf("For mode '%s', for uniform distribution, 'min' \n
                          must be less than 'max'.", modes))
           }
         },
         gamma = {
           if ((!(is.null(params$shape))) && (params$shape <= 0)) {
             stop(sprintf("For mode '%s', for gamma distribution, 'shape' \n 
                          must be greater than 0.", modes))
           }
           if ((!(is.null(params$scale))) && (params$scale <= 0)) {
             stop(sprintf("For mode '%s', for gamma distribution, 'scale' \n 
                          must be greater than 0.", modes))
           }
         },
         beta = {
           if ((!(is.null(params$shape1))) && (params$shape1 < 0)) {
             stop(sprintf("For mode '%s', for beta distribution, 'shape1' \n 
                          must be greater than or equal to 0.", modes))
           }
           if ((!(is.null(params$shape2))) && (params$shape2 < 0)) {
             stop(sprintf("For mode '%s', for beta distribution, 'shape2' \n 
                          must be greater than or equal to 0.", modes))
           }
         },
         binomial = {
           if ((!(is.null(params$size))) && (params$size < 0)) {
             stop(sprintf("For mode '%s', for binomial distribution, 'size' \n 
                          must be zero or positive.", modes))
           }
           if ((!(is.null(params$prob))) && 
               ((params$prob < 0) || (params$prob > 1))) {
             stop(sprintf("For mode '%s', for binomial distribution, 'prob' \n 
                          must be between 0 and 1.", modes))
           }
         },
         poisson = {
           if ((!(is.null(params$lambda))) && (params$lambda <= 0)) {
             stop(sprintf("For mode '%s', for poisson distribution, 'lambda' \n 
                          must be greater than 0.", modes))
           }
         },
         exponential = {
           if ((!(is.null(params$rate))) && (params$rate <= 0)) {
             stop(sprintf("For mode '%s', for exponential distribution, \n 
                          'rate' must be greater than 0.", modes))
           }
         })
  result <- distfuncs[[dname]](n, params)
  return(result)
}