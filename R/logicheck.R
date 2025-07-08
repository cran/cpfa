logicheck <- 
  function(argu) 
{
    argnam <- deparse(substitute(argu))
    if (length(argu) != 1L) {
      stop(sprintf("Input '%s' must contain only a single value.", argnam), 
           call. = FALSE)
    }
    if (!(is.logical(argu))) {
      stop(sprintf("Input '%s' must be of class 'logical'.", argnam), 
           call. = FALSE)
    }
}