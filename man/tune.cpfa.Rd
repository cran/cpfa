\name{tune.cpfa}
\alias{tune.cpfa}
\title{
Tuning for Classification with Parallel Factor Analysis
}
\description{
Fits Richard A. Harshman's Parallel Factor Analysis (Parafac) model-1 to a three-way or four-way data tensor/array. Uses Parafac factor weights from a single mode of this model as predictors to tune parameters for one or more classification methods via a k-fold cross-validation procedure. Supports binary and multiclass classification.
}
\usage{
tune.cpfa(x, y, nfac = 1, nfolds = 10, foldid = NULL, prior = NULL,
           method = c("PLR", "SVM", "RF", "NN"), 
           family = c("binomial", "multinomial"),
           alpha = NULL, lambda = NULL, cost = NULL, gamma = NULL, 
           ntree = NULL, nodesize = NULL, size = NULL, decay = NULL, 
           parallel = FALSE, cl = NULL, 
           verbose = TRUE, cmode = NULL, ...)
}
\arguments{
\item{x}{
  Three-way or four-way data array. See note below.
} 
\item{y}{
  A factor with two or more levels containing class labels. For binary case, ensure the order of factor levels (left to right) is such that negative class is first and positive class is second.
}
\item{nfac}{
  Number of factors for each Parafac model to estimate. Default is nfac = 1.
}
\item{nfolds}{
  Numeric setting number of folds for k-fold cross-validation. Must be 2 or greater. Default is 10 folds.
}
\item{foldid}{
  Integer vector containing fold IDs for k-fold cross-validation. If not provided, fold IDs are generated randomly for number of folds 'nfolds'.
}
\item{prior}{
  Prior probabilities of class membership. If unspecified, the class proportions for input 'y' are used. If present, the probabilities should be specified in the order of the factor levels of input 'y'.
}
\item{method}{
  Character vector indicating classification methods to use. Possible methods include penalized logistic regression (PLR), support vector machine (SVM), and random forest (RF). If none selected, default is to use all methods. See example.
}
\item{family}{
  Character value specifying binary classification (\code{family = "binomial"}) or multiclass classification (\code{family = "multinomial"}). If not provided, number of levels of input 'y' is used, where two levels is binary, and where three or more levels is multiclass.
}
\item{alpha}{
   Values for penalized logistic regression alpha parameter; default is alpha = seq(0, 1, length = 6). Must be numeric and contain only real numbers between 0 and 1, inclusive.
}
\item{lambda}{
   Optional user-supplied lambda sequence for \code{cv.glmnet}. Default is NULL.
}
\item{cost}{
  Values for support vector machine cost parameter; default is cost = c(1, 2, 4, 8, 16, 32, 64). Must be numeric and contain only real numbers greater than or equal to zero.
}
\item{gamma}{
  Values for support vector machine gamma parameter; default is
gamma = c(0, 0.01, 0.1, 1, 10, 100, 1000). Must be numeric and greater than or equal to 0.
}
\item{ntree}{
  Values for random forest number of trees parameter; default is ntree = c(100, 200, 400, 600, 800, 1600, 3200). Must be numeric and contain only integers greater than or equal to 1.
}
\item{nodesize}{
  Values for random forest node size parameter; default is nodesize = c(1, 2, 4, 8, 16, 32, 64). Must be numeric and contain only integers greater than or equal to 1.
}
\item{size}{
  Values for single-layer neural network size parameter; default is size = c(1, 2, 4, 8, 16, 32, 64). Must be numeric and contain only integers greater than or equal to 0.
}
\item{decay}{
  Values for neural network decay parameter; default is decay = c(0.001, 0.01, 0.1, 1, 2, 4, 8, 16). Must be numeric and contain only real numbers.
}
\item{parallel}{
  Logical indicating if package \code{parallel} should be used for parallel computing. For support vector machine and random forest, \code{doParallel} package is used as a wrapper. Defaults to FALSE, which implements sequential computing.
}
\item{cl}{
  Cluster for parallel computing, which is used when \code{parallel = TRUE}. Note that if \code{parallel = TRUE} and \code{cl = NULL}, then the cluster is defined as makeCluster(detectCores()).
}
\item{verbose}{
  If TRUE, progress is printed.
}
\item{cmode}{
  Integer value of 1, 2, or 3 (or 4 if 'x' is a four-way array) specifying mode whose factor weights will be predictors for classification. Defaults to last mode of inputted array (i.e. 3 for three-way array, and 4 for four-way array). 
}
\item{...}{
  Additional arguments to be passed to function \code{parafac} for Parafac model estimation.
}
}
\details{
After fitting a Parafac model with package \code{multiway} (see \code{parafac} in \code{multiway} for details), estimated classification mode weight matrix is passed to one or several of four classification methods--including penalized logistic regression (PLR); support vector machine (SVM); random forest (RF); and feed-forward neural network (NN). 

Package \code{glmnet} fits models for PLR. PLR tunes penalty parameter lambda; elastic net parameter alpha is set by user (see \code{cv.glmnet} in package \code{glmnet}). For SVM, package \code{e1071} is used with a radial basis kernel. Penalty parameter cost and radial basis parameter gamma are used (see \code{svm} in package \code{e1071}). For RF, package \code{randomForest} is used and implements Breiman's random forest algorithm. Number of predictors sampled at each node split is set at default of sqrt(R), where R is the number of Parafac factors. Two tuning parameters allowed are ntree, the number of trees to be grown, and nodesize, the minimum size of terminal nodes (see \code{randomForest} in package \code{randomForest}). For NN, package \code{nnet} fits a single-hidden-layer, feed-forward neural network model. Penalty parameters size (number of hidden layer units) and decay (weight decay) are used (see \code{nnet}).

For all four methods, k-fold cross-validation is implemented to tune classification parameters where the number of folds is set by argument \code{nfolds}.
}
\value{
  Returns an object of class 'cpfa' with the following elements:

  \item{opt.model}{
  List containing optimal model for tuned classification methods for each Parafac model estimated.
}
  \item{opt.param}{
  Data frame containing optimal parameters for tuned classification methods.
}
  \item{kcv.error}{
  Data frame containing KCV misclassification error for optimal parameters for tuned classification methods.
}
  \item{est.time}{
  Data frame containing estimation times for fitting Parafac model and for tuning classification methods.
}
  \item{method}{
  Numeric indicating classification methods used. Value of '1' indicates PLR; value of '2' indicates SVM; and value of '3' indicates 'RF'.
}
  \item{x}{
  three-way or four-way array used. 
}
  \item{y}{
  Factor containing class labels used. Note that output \code{y} is recoded such that the input labels of \code{y} are converted to numeric integers from 0 through the number of levels, which are then applied as labels for output \code{y}. 
}
  \item{Aweights}{
  List containing estimated A weights for each Parafac model estimated.
}
  \item{Bweights}{
  List containing estimated B weights for each Parafac model estimated.
}
  \item{Cweights}{
  List containing estimated C weights for each Parafac model estimated. Null if inputted argument \code{x} was a 3 array.
}
  \item{const}{
  Constraints used in Parafac model estimation. If argument \code{const} was not inputted, output value will be unconstrained for all modes.
}
  \item{cmode}{
  Integer value of 1, 2, or 3 (or 4 if \code{x} is a four-way array) specifying mode whose factor weights were predictors for classification.
}
  \item{family}{
  Character value specifying whether classification was binary  (\code{family = "binomial"}) or multiclass (\code{family = "multinomial"}).
}

}
\note{
If argument \code{cmode} is not null, input array \code{x} is reshaped with function \code{aperm} such that the \code{cmode} dimension of \code{x} is ordered last. Estimated mode A and B (and mode C for a four-way array \code{x}) weights that are outputted as \code{Aweights} and \code{Bweights} (and \code{Cweights}) reflect this permutation. For example, if \code{x} is a four-way array and \code{cmode = 2}, original modes/dimensions 1, 2, 3, and 4 will correspond to output order 1, 3, 4, 2. Here, output A = input 1; B = 3, and C = 4 (i.e. second mode specified by \code{cmode} has been moved to the last/D mode).  
}
\author{
Matthew A. Snodgress <snodg031@umn.edu>
}
\references{
Breiman, L. (2001). Random forests. Machine Learning, 45(1), 5-32.

Cortes, C. and Vapnik, V. (1995). Support-vector networks. Machine Learning, 20(3), 273-297.

Friedman, J. Hastie, T., and Tibshirani, R. (2010). Regularization Paths for Generalized Linear Models via Coordinate Descent. Journal of Statistical Software, 33(1), 1-22.

Harshman, R. A. (1970). Foundations of the PARAFAC procedure: Models and conditions for an "explanatory" multimodal factor analysis. UCLA Working Papers in Phonetics, 16, 1-84.

Harshman, R. A. and Lundy, M. E. (1994). PARAFAC: Parallel factor analysis. Computational Statistics and Data Analysis, 18, 39-72.

Helwig, N. E. (2017). Estimating latent trends in multivariate longitudinal data via Parafac2 with functional and structural constraints. Biometrical Journal, 59(4), 783-803.

Helwig, N. E. (2019). multiway: Component Models for Multi-Way Data. R package version 1.0-6.

Liaw, A. and Wiener, M. (2002). Classification and Regression by randomForest. R News 2(3), 18--22.

Meyer, D., Dimitriadou, E., Hornik, K., Weingessel, A., and Leisch, F. (2021). e1071: Misc Functions of the Department of Statistics, Probability Theory Group (Formerly: E1071), TU Wien. R package version 1.7-6.

Ripley, B. (1994). Neural networks and related methods for classification. Journal of the Royal Statistical Society: Series B (Methodological), 56(3), 409-437.

Venables, W. and Ripley, B. (2002) Modern Applied Statistics with S. Fourth Edition. Springer, New York. ISBN 0-387-95457-0.

Zou, H. and Hastie, T. (2005). Regularization and variable selection via the elastic net. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 67(2), 301-320.
}
\examples{
########## Random binary example with 3-way array ##########

# create random data array with Parafac structure
set.seed(3)
mydim <- c(60, 16, 80)
nf <- 3
Amat <- matrix(rnorm(mydim[1]*nf), nrow = mydim[1], ncol = nf)
Bmat <- matrix(runif(mydim[2]*nf), nrow = mydim[2], ncol = nf)
Cmat <- matrix(runif(mydim[3]*nf), nrow = mydim[3], ncol = nf)
Xmat <- tcrossprod(Amat, krprod(Cmat, Bmat))
Xmat <- array(Xmat, dim = mydim)
Emat <- array(rnorm(prod(mydim)), dim = mydim)
X <- Xmat + Emat
y <- factor(rbinom(mydim[3], 1, 0.4))

# initialize
alpha <- seq(0, 1, length = 2)
gamma <- c(0, 0.01)
cost <- c(1, 2)
ntree <- c(100, 200)
nodesize <- c(1, 2)
size <- c(1, 2)
decay <- c(0, 1)
method <- c("PLR", "SVM", "RF", "NN")
nfolds <- 3
const <- c("orthog", "uncons", "uncons")
foldid <- sample(rep(1:nfolds, length.out = length(y)))

# estimate Parafac models and use third mode to tune classification methods
tune.object <- tune.cpfa(x = X, y = y, nfac = 3, nfolds = nfolds,
                         foldid = foldid, method = method, alpha = alpha,
                         gamma = gamma, cost = cost, ntree = ntree,
                         nodesize = nodesize, size = size, decay = decay,
                         parallel = FALSE, const = const)
print.cpfa(tune.object)
}