repsums <- 
  function(mats) 
{
  one <- mats[[1L]]
  d <- dim(one)
  arr <- array(unlist(mats, use.names = FALSE), 
               dim = c(if (is.null(d)) length(one) else d, length(mats)))
  mar <- seq_len(length(dim(arr)) - 1L)
  red <- function(f) {
    z <- apply(arr, mar, f) 
    if (is.null(d)) names(z) <- names(one) else dimnames(z) <- dimnames(one)
    return(z)
  }
  flist <- list(mean = red(function(z) mean(z, na.rm = TRUE)), 
                median = red(function(z) median(z, na.rm = TRUE)), 
                sd = red(function(z) sd(z, na.rm = TRUE)), 
                n = length(mats), nvalid = red(function(z) sum(!(is.na(z)))))
  return(flist)
}