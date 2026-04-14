\name{tunecpfa}
\alias{tunecpfa}
\title{
Tuning for Classification with Parallel Factor Analysis
}
\description{
Fits Richard A. Harshman's Parallel Factor Analysis-1 (Parafac) model or 
Parallel Factor Analysis-2 (Parafac2) model to a three-way or four-way data 
array. Also fits a principal component analysis model (PCA) to a two-way matrix.
Allows for different constraint options on multiple tensor modes. For PCA, 
allows for either an unrotated or varimax rotated solution. Uses component 
weights from a single mode of the model as predictors to tune parameters for 
one or more classification methods via a k-fold cross-validation procedure. 
Allows for additional features to be included. Supports binary and multiclass 
classification.
}
\usage{
tunecpfa(x, y, z = NULL, model = c("parafac", "parafac2", "pca"), nfac = 1, 
         nfolds = 10, method = c("PLR", "SVM", "RF", "NN", "RDA", "GBM"), 
         family = c("binomial", "multinomial"), parameters = list(), 
         foldid = NULL, prior = NULL, cmode = NULL, parallel = FALSE, 
         cl = NULL, verbose = TRUE, compscale = TRUE, 
         pcarot = c("unrotated", "varimax"), ...)
}
\arguments{
\item{x}{
  A three-way or four-way data array. For Parafac2, can be a list where each 
  element is a matrix or three-way array. Array or list must contain only real 
  numbers. See note below. For PCA, can be a two-way matrix, a three-way array,
  or a four-way array. If a three-way or four-way array for PCA, array is 
  unfolded across all modes except the classification mode.
} 
\item{y}{
  A vector containing at least two unique class labels. Should be a factor that 
  contains two or more levels. For binary case, ensure the order of factor 
  levels (left to right) is such that negative class is first and positive 
  class is second.
}
\item{z}{
  A matrix containing additional features for classification. Must have number
  of rows equal to the number of class labels (i.e., equal to the length of 
  \code{y}). 
}
\item{model}{
  Character designating the component model to use, including: 
  \code{model = "parafac"} to fit the Parafac model, \code{model = "parafac2"} 
  to fit the Parafac2 model, or \code{model = "pca"} to fit the PCA model. One
  must be provided.
}
\item{nfac}{
  Vector containing integers that specify the number of components for each 
  component model to fit. Default is \code{nfac = 1}.
}
\item{nfolds}{
  Numeric value specifying the number of folds for k-fold cross-validation. 
  Must be 2 or greater. Default is \code{nfolds = 10}.
}
\item{method}{
  Character vector indicating classification methods to use. Possible methods 
  include penalized logistic regression (PLR); support vector machine (SVM); 
  random forest (RF); feed-forward neural network (NN); regularized 
  discriminant analysis (RDA); and gradient boosting machine (GBM). If none 
  selected, default is to use all methods.
}
\item{family}{
  Character value specifying binary classification (\code{family = "binomial"}) 
  or multiclass classification (\code{family = "multinomial"}). If not 
  provided, number of levels of input \code{y} is used, where two levels is 
  binary, and where three or more levels is multiclass.
}
\item{parameters}{
  List containing arguments related to classification methods. When specified, 
  must contain one or more of the following:
  \describe{
  \item{alpha}{
  Values for penalized logistic regression alpha parameter; default is 
  \code{alpha = seq(0, 1, length = 6)}. Must be numeric and contain only real 
  numbers between 0 and 1, inclusive.
}
  \item{lambda}{
  Optional user-supplied lambda sequence for \code{cv.glmnet} for penalized 
  logistic regression. Default is NULL.
}
  \item{cost}{
  Values for support vector machine cost parameter; default is 
  \code{cost = c(1, 2, 4, 8, 16, 32, 64)}. Must be numeric and contain only 
  real numbers greater than 0.
}
  \item{gamma}{
  Values for support vector machine gamma parameter; default is 
  \code{gamma = c(0, 0.01, 0.1, 1, 10, 100, 1000)}. Must be numeric and greater 
  than or equal to 0.
}
  \item{ntree}{
  Values for random forest number of trees parameter; default is 
  \code{ntree = c(100, 200, 400, 600, 800, 1600, 3200)}. Must be numeric and 
  contain only integers greater than or equal to 1.
}
  \item{nodesize}{
  Values for random forest node size parameter; default is 
  \code{nodesize = c(1, 2, 4, 8, 16, 32, 64)}. Must be numeric and contain 
  only integers greater than or equal to 1.
}
  \item{size}{
  Values for neural network size parameter; default is 
  \code{size = c(1, 2, 4, 8, 16, 32, 64)}. Must be numeric and contain only 
  integers greater than or equal to 0.
}
  \item{decay}{
  Values for neural network decay parameter; default is 
  \code{decay = c(0.001, 0.01, 0.1, 1, 2, 4, 8, 16)}. Must be numeric and 
  contain only real numbers.
}
  \item{rda.alpha}{
  Values for regularized discriminant analysis alpha parameter; default is 
  \code{rda.alpha = seq(0, 0.999, length = 6)}. Must be numeric and contain 
  only real numbers between 0 (inclusive) and 1 (exclusive).
}
  \item{delta}{
  Values for regularized discriminant analysis delta parameter; default is 
  \code{delta = c(0, 0.1, 1, 2, 3, 4)}. Must be numeric and contain only real 
  numbers greater than or equal to 0.
}
\item{eta}{
  Values for gradient boosting machine eta parameter; default is 
  \code{eta = c(0.1, 0.3, 0.5, 0.7, 0.9)}. Must be numeric and contain only 
  real numbers greater than 0 and less than 1.
}
\item{max.depth}{
  Values for gradient boosting machine max.depth parameter; default is 
  \code{max.depth = c(1, 2, 3, 4)}. Must be numeric and contain only integers 
  greater than or equal to 1.
}
\item{subsample}{
  Values for gradient boosting machine subsample parameter; default is 
  \code{subsample = c(0.6, 0.7, 0.8, 0.9)}. Must be numeric and contain only 
  real numbers greater than 0 and less than or equal to 1.
}
\item{nrounds}{
  Values for gradient boosting machine nrounds parameter; default is 
  \code{nrounds = c(100, 200, 300, 500)}. Must be numeric and contain only 
  integers greater than or equal to 1.
}
}
}
\item{foldid}{
  Vector containing fold IDs for k-fold cross-validation. Must be of class 
  integer. Should contain integers from 1 through the number of folds. If not 
  provided, fold IDs are generated randomly for observations using integers of 
  1 through the number of folds \code{nfolds}.
}
\item{prior}{
  Prior probabilities of class membership. If specified, the probabilities 
  should be in the order of the factor levels of input \code{y}. If unspecified, 
  the observed class proportions for input \code{y} are used. Based on 
  \code{prior}, inverse probability weights are calculated to account for class 
  imbalance. Note that RF and RDA ignore \code{prior} and use uniform priors 
  to handle imbalance.
}
\item{cmode}{
  Integer value of 1, 2, or 3 (or 4 if \code{x} is a four-way array) specifying 
  the mode whose component weights will be predictors for classification. 
  Defaults to the last mode of the input array (i.e., defaults to 3 for 
  a three-way array, and to 4 for a four-way array). 
  If \code{model = "parafac2"}, last mode will be used. If \code{model = "pca"},
  cmode is set to the first mode. 
}
\item{parallel}{
  Logical indicating whether to use parallel computing. If TRUE, the package 
  \bold{parallel} is used for parallel computing. For all classification 
  methods except penalized logistic regression, the \bold{doParallel} package 
  is used as a wrapper. Defaults to FALSE, which implements sequential 
  computing.
}
\item{cl}{
  Cluster for parallel computing, which is used when \code{parallel = TRUE}. 
  Note that if \code{parallel = TRUE} and \code{cl = NULL}, then the cluster 
  is defined as \code{makeCluster(detectCores())}.
}
\item{verbose}{
  If TRUE, progress is printed.
}
\item{compscale}{
  Logical indicating whether to scale each column of the estimated 
  classification component weights matrix (i.e., the features/predictors). If
  \code{TRUE}, each column is scaled to have mean zero and unit variance. If 
  \code{FALSE}, no scaling is performed.
}
\item{pcarot}{
  Character indicating whether to apply a varimax rotation or leave PCA loadings
  unrotated when \code{model} is set to \code{pca}. Ignored when \code{model} is
  not \code{pca}. Defaults to \code{"unrotated"} if not specified.
}
\item{...}{
  Additional arguments to be passed to function \code{parafac} for fitting a 
  Parafac model or function \code{parafac2} for fitting a Parafac2 model. 
  Example: can impose different constraints on different modes of the input 
  array using the argument \code{const}. See help file for function 
  \code{parafac} or for function \code{parafac2} for additional details.
}
}
\details{
After fitting a Parafac or Parafac2 model with the training set using package 
\bold{multiway} (see \code{parafac} or \code{parafac2} in \bold{multiway} for 
details), or after fitting a PCA model using the singular value decomposition, 
the estimated classification mode weight matrix is passed to one or more 
classification methods. The methods include: penalized logistic regression 
(PLR); support vector machine (SVM); random forest (RF); feed-forward neural 
network (NN); regularized discriminant analysis (RDA); and gradient boosting 
machine (GBM).

Package \bold{glmnet} fits models for PLR. PLR tunes penalty parameter lambda 
while the elastic net parameter alpha is set by the user (see the help file for 
function \code{cv.glmnet} in package \bold{glmnet}). For SVM, package 
\bold{e1071} is used with a radial basis kernel. Penalty parameter cost and 
radial basis parameter gamma are used (see \code{svm} in package \bold{e1071}). 
For RF, package \bold{randomForest} is used and implements Breiman's random 
forest algorithm. The number of predictors sampled at each node split is set 
at the default of sqrt(R), where R is the number of Parafac or Parafac2 
components. Two tuning parameters allowed are ntree, the number of trees to 
be grown, and nodesize, the minimum size of terminal nodes 
(see \code{randomForest} in package \bold{randomForest}). For NN, package 
\bold{nnet} fits a single-hidden-layer, feed-forward neural network model. 
Penalty parameters size (i.e., number of hidden layer units) and decay 
(i.e., weight decay) are used (see \bold{nnet}). For RDA, package \bold{rda} 
fits a shrunken centroids regularized discriminant analysis model. Tuning 
parameters include rda.alpha, the shrinkage penalty for the within-class 
covariance matrix, and delta, the shrinkage penalty of class centroids towards 
the overall dataset centroid. For GBM, package \bold{xgboost} fits a gradient 
boosting machine model. Four tuning parameters are allowed: (1) eta, the 
learning rate; (2) max.depth, the maximum tree depth; (3) subsample, the 
fraction of samples per tree; and (4) nrounds, the number of boosting trees 
to build.

For all six methods, k-fold cross-validation is implemented to tune 
classification parameters where the number of folds is set by argument 
\code{nfolds}.
}
\value{
Returns an object of class \code{tunecpfa} with the following elements:
  
\item{opt.model}{
  List containing optimal model for tuned classification methods for each 
  component model that was fit.
}
\item{opt.param}{
  Data frame containing optimal parameters for tuned classification methods.
}
\item{kcv.error}{
  Data frame containing KCV misclassification error for optimal parameters for 
  tuned classification methods.
}
\item{est.time}{
  Data frame containing times for fitting the component model and for 
  tuning classification methods.
}
\item{model}{
  Character designating the component model that was used, including:
  \code{model = "parafac"} for the Parafac model, \code{model = "parafac2"} 
  for the Parafac2 model, or \code{model = "pca"} for the PCA model.
}
\item{method}{
  Numeric indicating classification methods used. Value of '1' indicates 
  'PLR'; value of '2' indicates 'SVM'; value of '3' indicates 'RF'; value of 
  '4' indicates 'NN'; value of '5' indicates 'RDA'; and value of '6' indicates 
  'GBM'.
}
\item{x}{
  Two-way matrix, or three-way or four-way data array or list used in argument 
  \code{x}. If \code{x} was a three-way or four-way array, and if \code{model}
  was \code{pca}, returns the flattened two-way matrix.
}
\item{y}{
  Factor containing class labels used. Note that output \code{y} is recoded 
  such that the input labels of \code{y} are converted to numeric integers 
  from 0 through the number of levels, which are then applied as labels for 
  output \code{y}. 
}
\item{z}{
  Matrix containing additional features for classification.
}
\item{Aweights}{
  List containing estimated A weights for each component model that was fit. For
  PCA, contains the loadings.
}
\item{Bweights}{
  List containing estimated B weights for each component model that was fit. 
  Null if \code{model} was \code{pca}.
}
\item{Cweights}{
  List containing estimated C weights for each Parafac or Parafac2 model that 
  was fit. Null if input argument \code{x} was a three-way array or if 
  \code{model} was \code{pca}.
}
\item{Phi}{
  If \code{model = "parafac2"}, a list containing estimated \code{Phi} from 
  the Parafac2 model. \code{Phi} is the common cross product matrix shared by 
  all levels of the last mode (see help file for function \code{parafac2} in 
  package \bold{multiway} for additional details). NULL if 
  \code{model = "parafac"} or if \code{model = "pca"}.
}
\item{const}{
  Constraints used in fitting the Parafac, Parafac2, or PCA model. When 
  \code{model} was \code{pca}, contains the value used in input \code{pcarot}
  to indicate whether a rotation was applied to the solution.
}
\item{cmode}{
  Integer value of 1, 2, or 3 (or 4 if \code{x} is a four-way array) specifying 
  mode whose component weights were predictors for classification.
}
\item{family}{
  Character value specifying whether classification was binary 
  (\code{family = "binomial"}) or multiclass (\code{family = "multinomial"}).
}
\item{xdim}{
  Numeric value specifying number of levels for each mode of input \code{x}. 
  If \code{model = "parafac2"}, the number of levels for the first mode is 
  designated as \code{NA} because the number of levels can differ across levels 
  of the last mode. If \code{model = "pca"} and if input \code{x} was a 3-way or 
  4-way array, \code{xdim} contains the number of levels of the matrix that is 
  the unfolded version of the input array.
}
\item{lxdim}{
  Integer specifying the number of modes of the output \code{x}. Identical to 
  the number of modes of input \code{x}, except when input \code{x} was a 
  3-way or 4-way array and when \code{model} was \code{pca}. In this case,
  \code{lxdim} is 2.
}
\item{train.weights}{
  List containing classification component weights for each fit Parafac or 
  Parafac2 model, for possibly different numbers of components. The weights 
  used to train classifiers.
}
\item{priorweights}{
  List containing three elements based on input \code{prior}: (1) \code{weight}, 
  inverse probability weights per observation, used to balance classes for PLR, 
  NN, and GBM; (2) \code{frac}, inverse probability weights per class, used to 
  balance classes for SVM; and (3) \code{pricorrect}, uniform priors used to 
  balance classes for RF and RDA. 
}
\item{scenters}{
  List containing center of the scale for each Parafac or Parafac2 model that 
  was fit. Note that, for each component, estimated classification weights are 
  mean-centered before being passed to classification methods. Returns 
  \code{FALSE} if argument \code{compscale = FALSE}.
}
\item{sscales}{
  List containing standard devition of the scale for each Parafac or Parafac2 
  model that was fit. Note that, for each component, estimated classification 
  weights are scaled to have a variance of 1 before being passed to 
  classification methods. Returns \code{FALSE} if argument 
  \code{compscale = FALSE}.
}
\item{pcacenter}{
  Numeric vector containing centers for each feature used in PCA. NULL if 
  \code{model} is not \code{pca}.
}
}
\note{
For fitting the Parafac model, if argument \code{cmode} is not NULL, input 
array \code{x} is reshaped with function \code{aperm} such that the 
\code{cmode} dimension of \code{x} is ordered last. Estimated mode A and B 
(and mode C for a four-way array) weights that are outputted as \code{Aweights} 
and \code{Bweights} (and \code{Cweights}) reflect this permutation. For example, 
if \code{x} is a four-way array and \code{cmode = 2}, the original input 
modes 1, 2, 3, and 4 will correspond to output modes 1, 3, 4, 2. Here, 
output A = input 1; B = 3, and C = 4 (i.e., the second mode specified by 
\code{cmode} has been moved to the D mode/last mode). For 
\code{model = "parafac2"}, classification mode is assumed to be the last mode 
(i.e., mode C for three-way array and mode D for four-way array). For PCA, if 
argument \code{cmode} is not NULL, and if \code{x} is a 3-way or 4-way array, 
the array is reshaped with \code{aperm} such that the \code{cmode} dimension
of \code{x} is ordered first. Then, the array is unfolded into a two-way matrix.
If \code{x} is already input as a two-way matrix, the matrix is reshaped if 
\code{cmode} is 2, such that the matrix is transposed. Note that for PCA, 
\code{Aweights} contains the PCA model loadings.

In addition, note that the following combination of arguments will give an 
error: \code{nfac = 1, family = "multinomial", method = "PLR"}. The issue 
arises from providing \code{glmnet::cv.glmnet} input \code{x} with a matrix 
that has a single column. The issue is resolved for \code{family = "binomial"} 
because a column of 0s is appended to the single column, but this solution 
does not appear to work for the multiclass case. As such, this combination of 
arguments is not currently allowed. This issue will be resolved in a future 
update.

Applications of this function to real datasets can be explored at the following
repository: \url{https://github.com/matthewasisgress/multiway-classification/}.
}
\author{
Matthew A. Asisgress <mattgress@protonmail.ch>
}
\references{
Breiman, L. (2001). Random forests. Machine Learning, 45(1), 5-32.

Chen, T., He, T., Benesty, M., Khotilovich, V., Tang, Y., Cho, H., Chen, K., 
Mitchell, R., Cano, I., Zhou, T., Li, M., Xie, J., Lin, M., Geng, Y., Li, Y., 
Yuan, J. (2025). xgboost: Extreme gradient boosting. R Package Version 1.7.9.1.

Cortes, C. and Vapnik, V. (1995). Support-vector networks. Machine Learning, 
20(3), 273-297.

Friedman, J. H. (2001). Greedy function approximation: a gradient boosting 
machine. Annals of Statistics, 29(5), 1189-1232.

Friedman, J. H. (1989). Regularized discriminant analysis. Journal of the 
American Statistical Association, 84(405), 165-175.

Friedman, J., Hastie, T., and Tibshirani, R. (2010). Regularization paths for 
generalized linear models via coordinate descent. Journal of Statistical 
Software, 33(1), 1-22.

Gaujoux, R. (2025). doRNG: Generic reproducible parallel backend for 'foreach' 
loops. R Package Version 1.8.6.2.

Guo, Y., Hastie, T., and Tibshirani, R. (2007). Regularized linear discriminant 
analysis and its application in microarrays. Biostatistics, 8(1), 86-100.

Guo, Y., Hastie, T., and Tibshirani, R. (2023). rda: Shrunken centroids 
regularized discriminant analysis. R Package Version 1.2-1.

Harshman, R. (1970). Foundations of the PARAFAC procedure: Models and conditions 
for an "explanatory" multimodal factor analysis. UCLA Working Papers in 
Phonetics, 16, 1-84.

Harshman, R. (1972). PARAFAC2: Mathematical and technical notes. UCLA Working 
Papers in Phonetics, 22, 30-44.

Harshman, R. and Lundy, M. (1994). PARAFAC: Parallel factor analysis. 
Computational Statistics and Data Analysis, 18, 39-72.

Helwig, N. (2017). Estimating latent trends in multivariate longitudinal data 
via Parafac2 with functional and structural constraints. Biometrical Journal, 
59(4), 783-803.

Helwig, N. (2025). multiway: Component models for multi-way data. R Package 
Version 1.0-7.

Liaw, A. and Wiener, M. (2002). Classification and regression by randomForest. 
R News 2(3), 18--22.

Meyer, D., Dimitriadou, E., Hornik, K., Weingessel, A., and Leisch, F. (2024). 
e1071: Misc functions of the Department of Statistics, Probability Theory Group 
(Formerly: E1071), TU Wien. R Package Version 1.7-16.

Ripley, B. (1994). Neural networks and related methods for classification. 
Journal of the Royal Statistical Society: Series B (Methodological), 56(3), 
409-437.

Venables, W. and Ripley, B. (2002). Modern applied statistics with S. Fourth 
Edition. Springer, New York. ISBN 0-387-95457-0.

Zou, H. and Hastie, T. (2005). Regularization and variable selection via the 
elastic net. Journal of the Royal Statistical Society: Series B (Statistical 
Methodology), 67(2), 301-320.
}

\examples{
########## Parafac example with 3-way array and binary response ##########
\dontrun{
# set seed and simulate a three-way array connected to a binary response
set.seed(5)

# define list of arguments specifying distributions for A and G weights
techlist <- list(distA = list(dname = "poisson", 
                              lambda = 3),                 # for A weights
                 distG = list(dname = "gamma", shape = 2, 
                              scale = 4))                  # for G weights

# define target correlation matrix for columns of C mode weights matrix
cormat <- matrix(c(1, .6, .6, .6, 1, .6, .6, .6, 1), nrow = 3, ncol = 3)

# simulate a three-way array connected to a response
data <- simcpfa(arraydim = c(11, 12, 100), model = "parafac", nfac = 3, 
                nclass = 2, nreps = 1e2, onreps = 10, corresp = rep(.6, 3), 
                meanpred = rep(2, 3), modes = 3, corrpred = cormat,
                technical = techlist, smethod = "eigende")
                
# initialize
alpha <- seq(0, 1, length = 2)
gamma <- c(0, 0.01)
cost <- c(1, 2)
ntree <- c(100, 200)
nodesize <- c(1, 2)
size <- c(1, 2)
decay <- c(0, 1)
rda.alpha <- c(0.1, 0.6)
delta <- c(0.1, 2)
eta <- c(0.3, 0.7)
max.depth <- c(1, 2)
subsample <- c(0.75)
nrounds <- c(100)
method <- c("PLR", "SVM", "RF", "NN", "RDA", "GBM")
family <- "binomial"
parameters <- list(alpha = alpha, gamma = gamma, cost = cost, ntree = ntree,
                   nodesize = nodesize, size = size, decay = decay, 
                   rda.alpha = rda.alpha, delta = delta, eta = eta,
                   max.depth = max.depth, subsample = subsample,
                   nrounds = nrounds)
model <- "parafac"
nfolds <- 3
nstart <- 3

# constrain first mode weights to be orthogonal
const <- c("orthog", "uncons", "uncons")

# fit Parafac models and use third mode to tune classification methods
tune.object <- tunecpfa(x = data$X, y = data$y, z = data$Cmat, model = model, 
                        nfac = 3, nfolds = nfolds, method = method, 
                        family = family, parameters = parameters, 
                        parallel = FALSE, const = const, nstart = nstart)

# print tuning object
tune.object
}
}