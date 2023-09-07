\name{tune.cpfa}
\alias{tune.cpfa}
\title{
Tuning for Classification with Parallel Factor Analysis
}
\description{
Fits Richard A. Harshman's Parallel Factor Analysis-1 (Parafac) model or Parallel Factor Analysis-2 (Parafac2) model to a three-way or four-way data tensor/array. Uses Parafac component weights from a single mode of this model as predictors to tune parameters for one or more classification methods via a k-fold cross-validation procedure. Supports binary and multiclass classification.
}
\usage{
tune.cpfa(x, y, nfac = 1, nfolds = 10, foldid = NULL, prior = NULL,
           model = c("parafac", "parafac2"),
           method = c("PLR", "SVM", "RF", "NN"), 
           family = c("binomial", "multinomial"),
           alpha = NULL, lambda = NULL, cost = NULL, gamma = NULL, 
           ntree = NULL, nodesize = NULL, size = NULL, decay = NULL, 
           parallel = FALSE, cl = NULL, 
           verbose = TRUE, cmode = NULL, ...)
}
\arguments{
\item{x}{
  Three-way or four-way data array. When \code{model = "parafac2"}, can also be a list containing either matrices or three-way arrays. See note below.
} 
\item{y}{
  A factor with two or more levels containing class labels. For binary case, ensure the order of factor levels (left to right) is such that negative class is first and positive class is second.
}
\item{nfac}{
  Number of components for each Parafac or Parafac2 model to estimate. Default is \code{nfac = 1}.
}
\item{nfolds}{
  Numeric setting number of folds for k-fold cross-validation. Must be 2 or greater. Default is \code{nfolds = 10}.
}
\item{foldid}{
  Integer vector containing fold IDs for k-fold cross-validation. If not provided, fold IDs are generated randomly for number of folds \code{nfolds}.
}
\item{prior}{
  Prior probabilities of class membership. If unspecified, the class proportions for input \code{y} are used. If present, the probabilities should be specified in the order of the factor levels of input \code{y}.
}
\item{model}{
  Character designating the Parafac model to use, either \code{model = "parafac"} to fit the Parafac model or \code{model = "parafac2"} to fit the Parafac2 model.
}
\item{method}{
  Character vector indicating classification methods to use. Possible methods include penalized logistic regression (PLR), support vector machine (SVM), random forest (RF), and feed-forward neural networks (NN). If none selected, default is to use all methods. See example.
}
\item{family}{
  Character value specifying binary classification (\code{family = "binomial"}) or multiclass classification (\code{family = "multinomial"}). If not provided, number of levels of input \code{y} is used, where two levels is binary, and where three or more levels is multiclass.
}
\item{alpha}{
   Values for penalized logistic regression alpha parameter; default is \code{alpha = seq(0, 1, length = 6)}. Must be numeric and contain only real numbers between 0 and 1, inclusive.
}
\item{lambda}{
   Optional user-supplied lambda sequence for \code{cv.glmnet}. Default is NULL.
}
\item{cost}{
  Values for support vector machine cost parameter; default is \code{cost = c(1, 2, 4, 8, 16, 32, 64)}. Must be numeric and contain only real numbers greater than or equal to zero.
}
\item{gamma}{
  Values for support vector machine gamma parameter; default is
\code{gamma = c(0, 0.01, 0.1, 1, 10, 100, 1000)}. Must be numeric and greater than or equal to 0.
}
\item{ntree}{
  Values for random forest number of trees parameter; default is \code{ntree = c(100, 200, 400, 600, 800, 1600, 3200)}. Must be numeric and contain only integers greater than or equal to 1.
}
\item{nodesize}{
  Values for random forest node size parameter; default is \code{nodesize = c(1, 2, 4, 8, 16, 32, 64)}. Must be numeric and contain only integers greater than or equal to 1.
}
\item{size}{
  Values for feed-forward neural network size parameter; default is \code{size = c(1, 2, 4, 8, 16, 32, 64)}. Must be numeric and contain only integers greater than or equal to 0.
}
\item{decay}{
  Values for feed-forward neural network decay parameter; default is \code{decay = c(0.001, 0.01, 0.1, 1, 2, 4, 8, 16)}. Must be numeric and contain only real numbers.
}
\item{parallel}{
  Logical indicating if package \bold{parallel} should be used for parallel computing. For support vector machine and random forest, \bold{doParallel} package is used as a wrapper. Defaults to FALSE, which implements sequential computing.
}
\item{cl}{
  Cluster for parallel computing, which is used when \code{parallel = TRUE}. Note that if \code{parallel = TRUE} and \code{cl = NULL}, then the cluster is defined as \code{makeCluster(detectCores())}.
}
\item{verbose}{
  If TRUE, progress is printed.
}
\item{cmode}{
  Integer value of 1, 2, or 3 (or 4 if 'x' is a four-way array) specifying mode whose component weights will be predictors for classification. Defaults to last mode of inputted array (i.e., 3 for three-way array, and 4 for four-way array). If \code{model = "parafac2"}, last mode will be used.
}
\item{...}{
  Additional arguments to be passed to function \code{parafac} for fitting a Parafac model or function \code{parafac2} for fitting a Parafac2 model.
}
}
\details{
After fitting a Parafac or Parafac2 model with package \bold{multiway} (see \code{parafac} or \code{parafac2} in \bold{multiway} for details), estimated classification mode weight matrix is passed to one or several of four classification methods--including penalized logistic regression (PLR); support vector machine (SVM); random forest (RF); and feed-forward neural network (NN). 

Package \bold{glmnet} fits models for PLR. PLR tunes penalty parameter lambda; elastic net parameter alpha is set by user (see \code{cv.glmnet} in package \bold{glmnet}). For SVM, package \bold{e1071} is used with a radial basis kernel. Penalty parameter cost and radial basis parameter gamma are used (see \code{svm} in package \bold{e1071}). For RF, package \bold{randomForest} is used and implements Breiman's random forest algorithm. Number of predictors sampled at each node split is set at default of sqrt(R), where R is the number of Parafac or Parafac2 components. Two tuning parameters allowed are ntree, the number of trees to be grown, and nodesize, the minimum size of terminal nodes (see \code{randomForest} in package \bold{randomForest}). For NN, package \bold{nnet} fits a single-hidden-layer, feed-forward neural network model. Penalty parameters size (i.e., number of hidden layer units) and decay (i.e., weight decay) are used (see \bold{nnet}).

For all four methods, k-fold cross-validation is implemented to tune classification parameters where the number of folds is set by argument \code{nfolds}.
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
  Numeric indicating classification methods used. Value of '1' indicates PLR; value of '2' indicates SVM; value of '3' indicates 'RF'; and value of '4' indicates 'NN'.
}
  \item{x}{
  Three-way or four-way array used. If a list was used with \code{model = "parafac2"}, list of matrices or three-way arrays used.
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
  List containing estimated C weights for each Parafac or Parafac2 model that was fit. Null if inputted argument \code{x} was a 3-way array.
}
  \item{Phi}{
  If \code{model = "parafac2"}, a list containing estimated \code{Phi} from the Parafac2 model. \code{Phi} is the common crossproduct matrix shared by all levels of the last mode (see help file for function \code{parafac2} in package \bold{multiway} for additional details). NULL if \code{model = "parafac"}.
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

}
\note{
For fitting the Parafac model, if argument \code{cmode} is not null, input array \code{x} is reshaped with function \code{aperm} such that the \code{cmode} dimension of \code{x} is ordered last. Estimated mode A and B (and mode C for a four-way array \code{x}) weights that are outputted as \code{Aweights} and \code{Bweights} (and \code{Cweights}) reflect this permutation. For example, if \code{x} is a four-way array and \code{cmode = 2}, original modes/dimensions 1, 2, 3, and 4 will correspond to output order 1, 3, 4, 2. Here, output A = input 1; B = 3, and C = 4 (i.e. second mode specified by \code{cmode} has been moved to the last/D mode). For \code{model = "parafac2"}, classification mode is assumed to be the last mode (i.e., mode C for three-way array and mode D for four-way array).

In addition, note that there are cases where the following combination of arguments will give an error: \code{nfac = 1, family = "multinomial", method = "PLR"}.
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
model <- "parafac"
method <- c("PLR", "SVM", "RF", "NN")
nfolds <- 3
const <- c("orthog", "uncons", "uncons")
foldid <- sample(rep(1:nfolds, length.out = length(y)))

# estimate Parafac models and use third mode to tune classification methods
tune.object <- tune.cpfa(x = X, y = y, nfac = 3, nfolds = nfolds,
                         foldid = foldid, model = model, method = method, 
                         alpha = alpha, gamma = gamma, cost = cost, 
                         ntree = ntree, nodesize = nodesize, size = size, 
                         decay = decay, parallel = FALSE, const = const)
print(tune.object)
}