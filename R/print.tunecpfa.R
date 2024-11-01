print.tunecpfa <-
  function(x, ...)
{
    if (!inherits(x, "tunecpfa")) {
      stop("Input 'x' must be of class 'tunecpfa'.")
    }
    kcv.error <- x$kcv.error
    est.time <- x$est.time
    method <- x$method
    nfac <- x$opt.param$nfac
    nway <- x$lxdim
    model <- x$model
    method.char <- NULL
    mname <- c("PLR", "SVM", "RF", "NN", "RDA", "GBM")
    nus <- as.character(1:6)
    for (hh in 1:length(nus)) {
       if (nus[hh] %in% method) {method.char <- c(method.char, mname[hh])}
    }
    cat(paste0("Parafac Models Estimated:"))
    cat(paste0("\n", nway, "-way ", toupper(model), " with ", nfac, " factors"),
        "\n")
    cat(paste0("\n", "Classification Methods Tuned:"))
    cat(paste0("\n", method.char), "\n")
    cat(paste0("\n", "KCV Misclassification Error (est. time in seconds) \n 
               by Model and Method:"), "\n")
    for (w in seq_along(nfac)) {
       cnfac <- nfac[w]
       cerror <- kcv.error[which(nfac == cnfac), ]
       ctime <- est.time[which(nfac == cnfac), ]
       cat(paste0("\n", toupper(model), " with ", cnfac, " factors:"), "\n")
       for (hh in 1:length(nus)) {
          if (nus[hh] %in% method) {
            error <- round(eval(parse(text = paste0("cerror$error.", 
                                                    tolower(mname[hh])))), 4)
            time <- round(eval(parse(text = paste0("ctime$time.", 
                                                   tolower(mname[hh])))), 4)
            cat(paste0("  ", mname[hh], ":"))
            cat(paste0("  ", "Error = ", error, " (", time, ")\n"))
          }
       }
    }
}