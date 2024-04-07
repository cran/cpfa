\name{tune.cpfa}
\alias{tune.cpfa}
\title{
Deprecated - Tuning for Classification with Parallel Factor Analysis
}
\description{
Deprecated. Use function \code{tunecpfa} instead. See help file for \code{tunecpfa}.
}
\usage{
tune.cpfa(x, y, model = c("parafac", "parafac2"), nfac = 1, nfolds = 10, 
          method = c("PLR", "SVM", "RF", "NN", "RDA"), 
          family = c("binomial", "multinomial"), parameters = list(), 
          foldid = NULL, prior = NULL, cmode = NULL, parallel = FALSE, 
          cl = NULL, verbose = TRUE, ...)
}
\arguments{
\item{x}{
  For Parafac or Parafac2, a three-way or four-way data array. For Parafac2, can be a list of length \code{K} where the \code{k}-th element is a matrix or three-way array associated with the \code{k}-th element. Array or list must contain real numbers.
} 
\item{y}{
  A vector containing at least two unique class labels. Should be a factor that contains two or more levels . For binary case, ensure the order of factor levels (left to right) is such that negative class is first and positive class is second.
}
\item{model}{
  Character designating the Parafac model to use, either \code{model = "parafac"} to fit the Parafac model or \code{model = "parafac2"} to fit the Parafac2 model.
}
\item{nfac}{
  Number of components for each Parafac or Parafac2 model to fit. Default is \code{nfac = 1}.
}
\item{nfolds}{
  Numeric setting number of folds for k-fold cross-validation. Must be 2 or greater. Default is \code{nfolds = 10}.
}
\item{method}{
  Character vector indicating classification methods to use. Possible methods include penalized logistic regression (PLR); support vector machine (SVM); random forest (RF); feed-forward neural network (NN); and regularized discriminant analysis (RDA). If none selected, default is to use all methods.
}
\item{family}{
  Character value specifying binary classification (\code{family = "binomial"}) or multiclass classification (\code{family = "multinomial"}). If not provided, number of levels of input \code{y} is used, where two levels is binary, and where three or more levels is multiclass.
}
\item{parameters}{
  List containing arguments related to classification methods. When specified, must contain one or more of the following:
  \describe{
  \item{alpha}{
  Values for penalized logistic regression alpha parameter; default is \code{alpha = seq(0, 1, length = 6)}. Must be numeric and contain only real numbers between 0 and 1, inclusive.
}
  \item{lambda}{
  Optional user-supplied lambda sequence for \code{cv.glmnet} for penalized logistic regression. Default is NULL.
}
  \item{cost}{
  Values for support vector machine cost parameter; default is \code{cost = c(1, 2, 4, 8, 16, 32, 64)}. Must be numeric and contain only real numbers greater than or equal to zero.
}
  \item{gamma}{
  Values for support vector machine gamma parameter; default is \code{gamma = c(0, 0.01, 0.1, 1, 10, 100, 1000)}. Must be numeric and greater than or equal to 0.
}
  \item{ntree}{
  Values for random forest number of trees parameter; default is \code{ntree = c(100, 200, 400, 600, 800, 1600, 3200)}. Must be numeric and contain only integers greater than or equal to 1.
}
  \item{nodesize}{
  Values for random forest node size parameter; default is \code{nodesize = c(1, 2, 4, 8, 16, 32, 64)}. Must be numeric and contain only integers greater than or equal to 1.
}
  \item{size}{
  Values for neural network size parameter; default is \code{size = c(1, 2, 4, 8, 16, 32, 64)}. Must be numeric and contain only integers greater than or equal to 0.
}
  \item{decay}{
  Values for neural network decay parameter; default is \code{decay = c(0.001, 0.01, 0.1, 1, 2, 4, 8, 16)}. Must be numeric and contain only real numbers.
}
  \item{rda.alpha}{
  Values for regularized discriminant analysis alpha parameter; default is \code{rda.alpha = seq(0, 0.999, length = 6)}. Must be numeric and contain only real numbers between 0 (inclusive) and 1 (exclusive).
}
  \item{delta}{
  Values for regularized discriminant analysis delta parameter; default is \code{delta = c(0, 0.1, 1, 2, 3, 4)}. Must be numeric and contain only real numbers greater than or equal to 0.
}
}
}
\item{foldid}{
  Vector containing fold IDs for k-fold cross-validation. Can be of class integer, numeric, or data frame. Should contain integers from 1 through the number of folds. If not provided, fold IDs are generated randomly for observations using 1 through the number of folds \code{nfolds}.
}
\item{prior}{
  Prior probabilities of class membership. If unspecified, the class proportions for input \code{y} are used. If specified, the probabilities should be in the order of the factor levels of input \code{y}.
}
\item{cmode}{
  Integer value of 1, 2, or 3 (or 4 if \code{x} is a four-way array) specifying the mode whose component weights will be predictors for classification. Defaults to the last mode of the inputted array (i.e., defaults to 3 for three-way array, and to 4 for four-way array). If \code{model = "parafac2"}, last mode will be used.
}
\item{parallel}{
  Logical indicating if parallel computing should be implemented. If TRUE, the package \bold{parallel} is used for parallel computing. For all classification methods except penalized logistic regression, the \bold{doParallel} package is used as a wrapper. Defaults to FALSE, which implements sequential computing.
}
\item{cl}{
  Cluster for parallel computing, which is used when \code{parallel = TRUE}. Note that if \code{parallel = TRUE} and \code{cl = NULL}, then the cluster is defined as \code{makeCluster(detectCores())}.
}
\item{verbose}{
  If TRUE, progress is printed.
}
\item{...}{
  Additional arguments to be passed to function \code{parafac} for fitting a Parafac model or function \code{parafac2} for fitting a Parafac2 model. Example: can impose different constraints on different modes of the input array using the argument \code{const}. See help file for function \code{parafac} or for function \code{parafac2} for additional details.
}
}
\details{
See help file for \code{tunecpfa}.
}
\value{
Returns an object of class \code{cpfa} with the following elements:
  
\item{opt.model}{
  List containing optimal model for tuned classification methods for each Parafac or Parafac2 model that was fit.
}
\item{opt.param}{
  Data frame containing optimal parameters for tuned classification methods.
}
\item{kcv.error}{
  Data frame containing KCV misclassification error for optimal parameters for tuned classification methods.
}
\item{est.time}{
  Data frame containing times for fitting Parafac or Parafac2 model and for tuning classification methods.
}
\item{method}{
  Numeric indicating classification methods used. Value of '1' indicates 'PLR'; value of '2' indicates 'SVM'; value of '3' indicates 'RF'; value of '4' indicates 'NN'; and value of '5' indicates 'RDA'.
}
\item{x}{
  Three-way or four-way array used. If a list was used with \code{model = "parafac2"}, returns list of matrices or three-way arrays used.
}
\item{y}{
  Factor containing class labels used. Note that output \code{y} is recoded such that the input labels of \code{y} are converted to numeric integers from 0 through the number of levels, which are then applied as labels for output \code{y}. 
}
\item{Aweights}{
  List containing estimated A weights for each Parafac or Parafac2 model that was fit.
}
\item{Bweights}{
  List containing estimated B weights for each Parafac or Parafac2 model that was fit.
}
\item{Cweights}{
  List containing estimated C weights for each Parafac or Parafac2 model that was fit. Null if inputted argument \code{x} was a three-way array.
}
\item{Phi}{
  If \code{model = "parafac2"}, a list containing estimated \code{Phi} from the Parafac2 model. \code{Phi} is the common cross product matrix shared by all levels of the last mode (see help file for function \code{parafac2} in package \bold{multiway} for additional details). NULL if \code{model = "parafac"}.
}
\item{const}{
  Constraints used in fitting Parafac or Parafac2 models. If argument \code{const} was not inputted, no constraints will be used.
}
\item{cmode}{
  Integer value of 1, 2, or 3 (or 4 if \code{x} is a four-way array) specifying mode whose component weights were predictors for classification.
}
\item{family}{
  Character value specifying whether classification was binary (\code{family = "binomial"}) or multiclass (\code{family = "multinomial"}).
}
\item{xdim}{
  Numeric value specifying number of levels for each mode of input \code{x}. If \code{model = "parafac2"}, number of levels for first mode is designated as \code{NA} because the number of levels can differ across levels of the last mode.
}
\item{lxdim}{
  Numeric value specifying number of modes of input \code{x}.
}
\item{train.weights}{
  List containing classification component weights for each fit Parafac or Parafac2 model, for possibly different numbers of components. The weights used to train classifiers.
}
}
\author{
Matthew A. Snodgress <snodg031@umn.edu>
}