alignmod <- 
  function(baselist, newlist, cnfac, isp2 = NULL, fourway = NULL) 
{
    if ((isp2 == TRUE) && (fourway == TRUE)) {
      Hmat <- abs(congru(baselist$B, newlist$B)) 
      Hmat <- Hmat * abs(congru(baselist$C, newlist$C))
      neword <- as.integer(solve_LSAP(Hmat, maximum = TRUE))
      newlist$A <- lapply(newlist$A, function(mat) mat[, neword, drop = FALSE])
      newlist$B <- newlist$B[, neword]; newlist$C <- newlist$C[, neword]
      newlist$D <- newlist$D[, neword]; newlist$predw <- newlist$predw[, neword]
      newlist$Phi <- newlist$Phi[neword, neword, drop = FALSE]
    } else if ((isp2 == TRUE) && (fourway == FALSE)) {
      Hmat <- abs(congru(baselist$B, newlist$B))
      neword <- as.integer(solve_LSAP(Hmat, maximum = TRUE))
      newlist$A <- lapply(newlist$A, function(mat) mat[, neword, drop = FALSE])
      newlist$B <- newlist$B[, neword]; newlist$C <- newlist$C[, neword]
      newlist$predw <- newlist$predw[, neword]
      newlist$Phi <- newlist$Phi[neword, neword, drop = FALSE]
    } else if ((isp2 == FALSE) && (fourway == TRUE)) {
      Hmat <- abs(congru(baselist$A, newlist$A))
      Hmat <- Hmat * abs(congru(baselist$B, newlist$B)) 
      Hmat <- Hmat * abs(congru(baselist$C, newlist$C))
      neword <- as.integer(solve_LSAP(Hmat, maximum = TRUE))
      newlist$A <- newlist$A[, neword]; newlist$B <- newlist$B[, neword]
      newlist$C <- newlist$C[, neword]; newlist$D <- newlist$D[, neword]
      newlist$predw <- newlist$predw[, neword]
    } else {
      Hmat <- abs(congru(baselist$A, newlist$A))
      Hmat <- Hmat * abs(congru(baselist$B, newlist$B))
      neword <- as.integer(solve_LSAP(Hmat, maximum = TRUE))
      newlist$A <- newlist$A[, neword]; newlist$B <- newlist$B[, neword]
      newlist$C <- newlist$C[, neword]; newlist$predw <- newlist$predw[, neword]
    }
    finalB <- diag(abs(congru(baselist$B, newlist$B)))
    newlist$changeorder <- neword
    newlist$finalB <- finalB
    return(newlist)
}