\name{cpfa}
\alias{cpfa}
\title{
Classification with Parallel Factor Analysis
}
\description{
Fits Richard A. Harshman's Parallel Factor Analysis-1 (Parafac) model or Parallel Factor Analysis-2 (Parafac2) model to a three-way or four-way data array. Allows for different constraint options on multiple tensor modes. Uses Parafac component weights from a single mode of this model as predictors to tune parameters for one or more classification methods via a k-fold cross-validation procedure. Predicts class labels and calculates multiple performance measures for binary or multiclass classification over some number of replications with different train-test splits. Provides descriptive statistics to pool output across replications.
}
\usage{
cpfa(x, y, model = c("parafac", "parafac2"), nfac = 1, nrep = 5, ratio = 0.8,
     nfolds = 10, method = c("PLR", "SVM", "RF", "NN", "RDA", "GBM"), 
     family = c("binomial", "multinomial"), parameters = list(), 
     type.out = c("measures", "descriptives"), foldid = NULL, 
     prior = NULL, cmode = NULL, seeds = NULL, plot.out = FALSE, 
     plot.measures = NULL, parallel = FALSE, cl = NULL, verbose = TRUE, ...)
}
\arguments{
\item{x}{
  A three-way or four-way data array. For Parafac2, can be a list of length \code{K} where the \code{k}-th element is a matrix or three-way array associated with the \code{k}-th element. Array or list must contain only real numbers. See note below.
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
\item{nrep}{
  Number of replications to repeat the procedure. Default is \code{nrep = 5}.
}
\item{ratio}{
  Split ratio for dividing data into train and test sets. Default is \code{ratio = 0.8}.
}
\item{nfolds}{
  Numeric setting number of folds for k-fold cross-validation. Must be 2 or greater. Default is \code{nfolds = 10}.
}
\item{method}{
  Character vector indicating classification methods to use. Possible methods include penalized logistic regression (PLR); support vector machine (SVM); random forest (RF); feed-forward neural network (NN); regularized discriminant analysis (RDA); and gradient boosting machine (GBM). If none are selected, default is to use all methods with \code{method = c("PLR", "SVM", "RF", "NN", "RDA", "GBM")}.
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
\item{eta}{
  Values for gradient boosting machine eta parameter; default is \code{eta = c(0.1, 0.3, 0.5, 0.7, 0.9)}. Must be numeric and contain only real numbers greater than 0 and less than 1.
}
\item{max.depth}{
  Values for gradient boosting machine max.depth parameter; default is \code{max.depth = c(1, 2, 3, 4)}. Must be numeric and contain only integers greater than or equal to 1.
}
\item{subsample}{
  Values for gradient boosting machine subsample parameter; default is \code{subsample = c(0.6, 0.7, 0.8, 0.9)}. Must be numeric and contain only real numbers greater than 0 and less than or equal to 1.
}
\item{nrounds}{
  Values for gradient boosting machine nrounds parameter; default is \code{nrounds = c(100, 200, 300, 500)}. Must be numeric and contain only integers greater than or equal to 1.
}
}
}
\item{type.out}{
  Type of output desired: \code{type.out = "measures"} gives array containing classification performance measures for all replications while \code{type.out = "descriptives"} gives list of descriptive statistics calculated across all replications for each performance measure. Both options also provide the estimated training weights and classification weights. Defaults to \code{type.out = "descriptives"}.
}
\item{foldid}{
  Integer vector containing fold IDs for k-fold cross-validation. If not provided, fold IDs are generated randomly for number of folds \code{nfolds}.
}
\item{prior}{
  Prior probabilities of class membership. If unspecified, the class proportions for input \code{y} are used. If specified, the probabilities should be in the order of the factor levels of input \code{y}.
}
\item{cmode}{
  Integer value of 1, 2, or 3 (or 4 if \code{x} is a four-way array) specifying the mode whose component weights will be predictors for classification. Defaults to the last mode of the inputted array (i.e., defaults to 3 for three-way array, and to 4 for four-way array). If \code{model = "parafac2"}, last mode will be used.
}
\item{seeds}{
  Random seeds to be associated with each replication. Default is \code{seeds = 1:nrep}.
}
\item{plot.out}{
  Logical indicating whether to output one or more box plots of classification performance measures that are plotted across classification methods and number of components.
}
\item{plot.measures}{
  Character vector containing values that specify for plotting one or more of 11 possible classification performance measures. Only relevant when \code{plot.out = TRUE}. Should contain one or more of the following labels: \code{c("err", "acc", "tpr", "fpr", "tnr", "fnr", "ppv", "npv", "fdr", "fom", "fs")}. A box plot will be created for each measure that is specified, summarizing output across replications. Note that additional information about each label is available in the Details section of the help file for function \code{cpm}. Note also that there are a few cases where the x-axis tick labels for a plot might not appear. This issue will be resolved in a future update.
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
Data are split into a training set and a testing set. After fitting a Parafac or Parafac2 model with the training set using package \bold{multiway} (see \code{parafac} or \code{parafac2} in \bold{multiway} for details), the estimated classification mode weight matrix is passed to one or several of six classification methods. The methods include: penalized logistic regression (PLR); support vector machine (SVM); random forest (RF); feed-forward neural network (NN); regularized discriminant analysis (RDA); and gradient boosting machine (GBM).

Package \bold{glmnet} fits models for PLR. PLR tunes penalty parameter lambda while the elastic net parameter alpha is set by the user (see the help file for function \code{cv.glmnet} in package \bold{glmnet}). For SVM, package \bold{e1071} is used with a radial basis kernel. Penalty parameter cost and radial basis parameter gamma are used (see \code{svm} in package \bold{e1071}). For RF, package \bold{randomForest} is used and implements Breiman's random forest algorithm. The number of predictors sampled at each node split is set at the default of sqrt(R), where R is the number of Parafac or Parafac2 components. Two tuning parameters allowed are ntree, the number of trees to be grown, and nodesize, the minimum size of terminal nodes (see \code{randomForest} in package \bold{randomForest}). For NN, package \bold{nnet} fits a single-hidden-layer, feed-forward neural network model. Penalty parameters size (i.e., number of hidden layer units) and decay (i.e., weight decay) are used (see \bold{nnet}). For RDA, package \bold{rda} fits a shrunken centroids regularized discriminant analysis model. Tuning parameters include rda.alpha, the shrinkage penalty for the within-class covariance matrix, and delta, the shrinkage penalty of class centroids towards the overall dataset centroid. For GBM, package \bold{xgboost} fits a gradient boosting machine model. Four tuning parameters are allowed: (1) eta, the learning rate; (2) max.depth, the maximum tree depth; (3) subsample, the fraction of samples per tree; and (4) nrounds, the number of boosting trees to build.

For all six methods, k-fold cross-validation is implemented to tune classification parameters where the number of folds is set by argument \code{nfolds}. Separately, the trained Parafac or Parafac2 model is used to predict the classification mode's component weights using the testing set data. The predicted component weights and the optimized classification method are then used to predict class labels. Finally, classification performance measures are calculated. The process is repeated over a number of replications with different random splits of the input array and of the class labels at each replication.
}
\value{
  Returns an object of class \code{wrapcpfa} either with a three-way array with classification performance measures for each model and for each replication, or with a list containing matrices with descriptive statistics for performance measures calculated across all replications. Specify \code{type.out = "measures"} to output the array of performance measures. Specify \code{type.out = "descriptives"} to output descriptive statistics across replications. In addition, for both options, the following are also provided:
  \item{predweights}{
    List of predicted classification weights for each Parafac or Parafac2 model and for each replication.
  }
  \item{train.weights}{
    List of lists of training weights for each Parafac or Parafac2 model and for each replication.
  }
  \item{opt.tune}{
    List of optimal tuning parameters for classification methods for each Parafac or Parafac2 model and for each replication.
  }
  \item{mean.opt.tune}{
    Mean across all replications of optimal tuning parameters for classification methods for each Parafac or Parafac2 model. 
  }
  \item{X}{
    Three-way or four-way data array or list used in argument \code{x}.
  }
  \item{nfac}{
    Number of components used to fit each Parafac or Parafac2 model.
  }
  \item{model}{
    Character designating the Parafac model that was used, either \code{model = "parafac"} for the Parafac model or \code{model = "parafac2"} for the Parafac2 model.
  }
  \item{method}{
    Classification methods used.
  }
  \item{const}{
    Constraints used in fitting Parafac or Parafac2 models.
  }
  \item{cmode}{
    Integer value used to specify the mode whose component weights were predictors for classification. 
}
}
\note{
If argument \code{cmode} is not null, input array \code{x} is reshaped with function \code{aperm} such that the \code{cmode} dimension of \code{x} is ordered last. Estimated mode A and B (and mode C for a four-way array) weights that are outputted as \code{Aweights} and \code{Bweights} (and \code{Cweights}) reflect this permutation. For example, if \code{x} is a four-way array and \code{cmode = 2}, the original input modes 1, 2, 3, and 4 will correspond to output modes 1, 3, 4, 2. Here, output A = input 1; B = 3, and C = 4 (i.e., the second mode specified by \code{cmode} has been moved to the D mode/last mode). For \code{model = "parafac2"}, classification mode is assumed to be the last mode (i.e., mode C for three-way array and mode D for four-way array). 

In addition, note that the following combination of arguments will give an error: \code{nfac = 1, family = "multinomial", method = "PLR"}. The issue arises from providing \code{glmnet::cv.glmnet} input \code{x} a matrix with a single column. The issue is resolved for \code{family = "binomial"} because a column of 0s is appended to the single column, but this solution does not appear to work for the multiclass case. As such, this combination of arguments is not currently allowed. This issue will be resolved in a future update.
}
\author{
Matthew A. Snodgress <snodg031@umn.edu>
}
\references{
Breiman, L. (2001). Random forests. Machine Learning, 45(1), 5-32.

Chen, T., He, T., Benesty, M., Khotilovich, V., Tang, Y., Cho, H., Chen, K., Mitchell, R., Cano, I., Zhou, T., Li, M., Xie, J., Lin, M., Geng, Y., Li, Y., Yuan, J. (2024). xgboost: Extreme gradient boosting. R Package Version 1.7.7.1.

Cortes, C. and Vapnik, V. (1995). Support-vector networks. Machine Learning, 20(3), 273-297.

Friedman, J. H. (2001). Greedy function approximation: a gradient boosting machine. Annals of Statistics, 29(5), 1189-1232.

Friedman, J. H. (1989). Regularized discriminant analysis. Journal of the American Statistical Association, 84(405), 165-175.

Friedman, J. Hastie, T., and Tibshirani, R. (2010). Regularization paths for generalized linear models via coordinate descent. Journal of Statistical Software, 33(1), 1-22.

Guo, Y., Hastie, T., and Tibshirani, R. (2007). Regularized linear discriminant analysis and its application in microarrays. Biostatistics, 8(1), 86-100.

Guo Y., Hastie T., and Tibshirani, R. (2023). rda: Shrunken centroids regularized discriminant analysis. R Package Version 1.2-1.

Harshman, R. (1970). Foundations of the PARAFAC procedure: Models and conditions for an "explanatory" multimodal factor analysis. UCLA Working Papers in Phonetics, 16, 1-84.

Harshman, R. (1972). PARAFAC2: Mathematical and technical notes. UCLA Working Papers in Phonetics, 22, 30-44.

Harshman, R. and Lundy, M. (1994). PARAFAC: Parallel factor analysis. Computational Statistics and Data Analysis, 18, 39-72.

Helwig, N. (2017). Estimating latent trends in multivariate longitudinal data via Parafac2 with functional and structural constraints. Biometrical Journal, 59(4), 783-803.

Helwig, N. (2019). multiway: Component models for multi-way data. R Package Version 1.0-6.

Liaw, A. and Wiener, M. (2002). Classification and regression by randomForest. R News 2(3), 18--22.

Meyer, D., Dimitriadou, E., Hornik, K., Weingessel, A., and Leisch, F. (2023). e1071: Misc functions of the Department of Statistics, Probability Theory Group (Formerly: E1071), TU Wien. R Package Version 1.7-13.

Ripley, B. (1994). Neural networks and related methods for classification. Journal of the Royal Statistical Society: Series B (Methodological), 56(3), 409-437.

Venables, W. and Ripley, B. (2002). Modern applied statistics with S. Fourth Edition. Springer, New York. ISBN 0-387-95457-0.

Zou, H. and Hastie, T. (2005). Regularization and variable selection via the elastic net. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 67(2), 301-320.
}
\examples{
########## Parafac2 example with 4-way array and multiclass response ##########

# set seed and specify dimensions of a four-way tensor
set.seed(5)
mydim <- c(10, 11, 12, 100)
nf <- 3

# create correlation matrix between response and fourth mode's weights 
rho.dd <- .35 
rho.dy <- .75 
cormat.values <- c(1, rho.dd, rho.dd, rho.dy, rho.dd, 1, rho.dd, rho.dy, 
                   rho.dd, rho.dd, 1, rho.dy, rho.dy, rho.dy, rho.dy, 1)
cormat <- matrix(cormat.values, nrow = (nf + 1), ncol = (nf + 1))

# sample from a multivariate normal with specified correlation structure
ymean <- Dmean <- 2
mu <- as.matrix(c(Dmean, Dmean, Dmean, ymean))
eidecomp <- eigen(cormat, symmetric = TRUE)
L.sqrt <- diag(eidecomp$values^0.5)
cormat.sqrt <- eidecomp$vectors \%*\% L.sqrt \%*\% t(eidecomp$vectors)
Z <- matrix(rnorm(mydim[4] * (nf + 1)), nrow = mydim[4], ncol = (nf + 1))
Xw <- rep(1, mydim[4]) \%*\% t(mu) + Z \%*\% cormat.sqrt
Dmat <- Xw[, 1:nf]

# create a random four-way data tensor with D weights related to a response
Bmat <- matrix(runif(mydim[2] * nf), nrow = mydim[2], ncol = nf)
Cmat <- matrix(runif(mydim[3] * nf), nrow = mydim[3], ncol = nf)
nDd <- rep(c(10, 12, 14), length.out = mydim[4])
Gmat <- matrix(rnorm(nf * nf), nrow = nf)
Amat <- vector("list", mydim[4])
X <- Xmat <- Emat <- Amat
for (Dd in 1:mydim[4]) {
   Amat[[Dd]] <- matrix(nf * rnorm(nDd[Dd]), nrow = nDd[Dd], ncol = nf)
   Amat[[Dd]] <- svd(Amat[[Dd]], nv = 0)$u \%*\% Gmat
   leftMat <- Amat[[Dd]] \%*\% diag(Dmat[Dd,])
   Xmat[[Dd]] <- array(tcrossprod(leftMat, krprod(Cmat, Bmat)), 
                       dim = c(nDd[Dd], mydim[2], mydim[3]))
   Emat[[Dd]] <- array(rnorm(nDd[Dd] * mydim[2] * mydim[3]), 
                       dim = c(nDd[Dd], mydim[2], mydim[3]))
   X[[Dd]] <- Xmat[[Dd]] + Emat[[Dd]]
}

# create a multiclass response
stor <- matrix(rep(1, nrow(Xw)), nrow = nrow(Xw))
stor[which(Xw[, (nf + 1)] < (ymean - 0.4 * sd(Xw[, (nf + 1)])))] <- 2
stor[which(Xw[, (nf + 1)] > (ymean + 0.4 * sd(Xw[, (nf + 1)])))] <- 0
y <- factor(stor)

# initialize
alpha <- seq(0, 1, length = 2)
gamma <- c(0, 1)
cost <- c(0.1, 5)
ntree <- c(200, 300)
nodesize <- c(1, 2)
size <- c(1, 2)
decay <- c(0, 1)
rda.alpha <- seq(0.1, 0.9, length = 2)
delta <- c(0.1, 2)
eta <- c(0.3, 0.7)
max.depth <- c(1, 2)
subsample <- c(0.75)
nrounds <- c(100)
method <- c("PLR", "SVM", "RF", "NN", "RDA", "GBM")
family <- "multinomial"
parameters <- list(alpha = alpha, gamma = gamma, cost = cost, ntree = ntree,
                   nodesize = nodesize, size = size, decay = decay, 
                   rda.alpha = rda.alpha, delta = delta, eta = eta,
                   max.depth = max.depth, subsample = subsample,
                   nrounds = nrounds)
model <- "parafac2"
nfolds <- 3
nstart <- 3

# constrain first mode weights to be orthogonal, fourth mode to be nonnegative
const <- c("orthog", "uncons", "uncons", "nonneg")

# fit Parafac2 model and use fourth mode weights to tune classification
# methods, to predict class labels, and to return classificaiton 
# performance measures pooled across multiple train-test splits
output <- cpfa(x = X, y = y, model = model, nfac = nf, nrep = 2, ratio = 0.8, 
               nfolds = nfolds, method = method, family = family, 
               parameters = parameters, type.out = "descriptives", 
               seeds = NULL, plot.out = TRUE, parallel = FALSE, const = const,
               nstart = nstart)

# print performance measure means across train-test splits
output$descriptive$mean
}