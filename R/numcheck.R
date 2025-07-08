numcheck <- 
  function(argum) 
{
    argnam <- deparse(substitute(argum))
    if (length(argum) != 1L) {
      stop(sprintf("Input '%s' must contain only a single value.", argnam), 
           call. = FALSE)
    }
    if (!(is.numeric(argum))) {
      stop(sprintf("Input '%s' must be of class 'numeric'.", argnam), 
           call. = FALSE)
    }
    if ((is.infinite(argum)) || (is.na(argum)) || (is.nan(argum))) { 
      stop(sprintf("Input '%s' must be a finite number and cannot be NA \n 
                   or NaN.", argnam), call. = FALSE)
    }
}