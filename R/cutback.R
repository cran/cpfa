cutback <- 
  function(fullmatrix, newcomps, cnfac) 
{
    fullmatrix[, 1:cnfac] <- newcomps
    return(fullmatrix)
}