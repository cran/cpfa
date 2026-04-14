print.tunecpfa <-
  function(x, ...)
{
    if (!(inherits(x, "tunecpfa"))) {
      stop("Input 'x' must be of class 'tunecpfa'.")
    }
    kcv.error <- x$kcv.error; est.time <- x$est.time
    method <- x$method; nfac <- x$opt.param$nfac
    nway <- x$lxdim; model <- x$model
    consts0 <- x$const; zout <- x$z
    if (!(is.null(zout))) {ncovari <- dim(zout)[2]}
    consts <- paste(consts0, collapse = "-")
    method.char <- NULL
    mname <- c("PLR", "SVM", "RF", "NN", "RDA", "GBM")
    nus <- as.character(1:6)
    for (hh in 1:length(nus)) {
       if (nus[hh] %in% method) {method.char <- c(method.char, mname[hh])}
    }
    cat(paste0("Component Models Fit:"))
    cat(paste0("\n", nway, "-way ", toupper(model), " with ", nfac, 
               " component(s)", " and constraint(s): ", consts), "\n")
    cat(paste0("\n", "Classification Methods Tuned:"))
    cat(paste0("\n", method.char), "\n")
    if (!(is.null(zout))) {
      cat(paste0("\n", "Number of Additional Features Included:"))
      cat(paste0("\n", ncovari), "\n")
    }
    cat(paste0("\n", "KCV Misclassification Error (est. time in seconds) \n 
               by Model and Method:"), "\n")
    for (w in seq_along(nfac)) {
       cnfac <- nfac[w]
       cerror <- kcv.error[which(nfac == cnfac), ]
       ctime <- est.time[which(nfac == cnfac), ]
       cat(paste0("\n", toupper(model), " with ", cnfac, " component(s):"), 
           "\n")
       for (hh in 1:length(nus)) {
          if (nus[hh] %in% method) {
            colzname <- paste0("error.", tolower(mname[hh]))
            error <- round(cerror[[colzname]], 4)
            colznametime <- paste0("time.", tolower(mname[hh]))
            time <- round(ctime[[colznametime]], 4)
            cat(paste0("  ", mname[hh], ":"))
            cat(paste0("  ", "Error = ", error, " (", time, ")\n"))
          }
       }
    }
    invisible(x)
}