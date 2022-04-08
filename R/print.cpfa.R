print.cpfa <-
  function(x, ...)
{
   kcv.error <- x$kcv.error
   est.time <- x$est.time
   method <- x$method
   nfac <- x$opt.param$nfac
   nway <- length(dim(x$x))
   method.char <- NULL
   if ('1' %in% method) {method.char <- c(method.char, "PLR")}
   if ('2' %in% method) {method.char <- c(method.char, "SVM")}
   if ('3' %in% method) {method.char <- c(method.char, "RF")}
   if ('4' %in% method) {method.char <- c(method.char, "NN")}
   cat(paste0("Parafac Models Estimated:"))
   cat(paste0("\n", nway, "-way Parafac with ", nfac, " factors"), "\n")
   cat(paste0("\n", "Classification Methods Tuned:"))
   cat(paste0("\n",method.char), "\n")
   cat(paste0("\n", "KCV Misclassification Error (estimation time in seconds) by 
Model and Method:"), "\n")
   for (w in 1:length(nfac)) {
      cnfac <- nfac[w]
      cerror <- kcv.error[which(nfac == cnfac),]
      ctime <- est.time[which(nfac == cnfac),]
      cat(paste0("\n","Parafac with ", cnfac, " factors:"), "\n")
      if ('1' %in% method) {
        error <- round(cerror$error.plr, 4)
        time <- round(ctime$time.plr, 4)
        cat(paste0("  PLR:"))
        cat(paste0("  ", "Error = ", error, " (", time, ")"))
      }
      if ('2' %in% method) {
        error <- round(cerror$error.svm, 4)
        time <- round(ctime$time.svm, 4)
        cat(paste0("\n", "  SVM:"))
        cat(paste0("  ", "Error = ", error, " (", time, ")"))
      }
      if ('3' %in% method) {
        error <- round(cerror$error.rf, 4)
        time <- round(ctime$time.rf, 4)
        cat(paste0("\n", "  RF:"))
        cat(paste0("  ", " Error = ", error, " (", time, ")"))
      }
      if ('4' %in% method) {
        error <- round(cerror$error.nn, 4)
        time <- round(ctime$time.nn, 4)
        cat(paste0("\n", "  NN:"))
        cat(paste0("  ", " Error = ", error, " (", time, ")"))
      }
   }
}