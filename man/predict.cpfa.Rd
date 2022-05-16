\name{predict.cpfa}
\alias{predict.cpfa}
\title{
Predict Method for Classification with Parallel Factor Analysis
}
\description{
\code{predict} method for class "cpfa".
}
\usage{
\method{predict}{cpfa}(object, newdata = NULL, nfac = NULL, method = NULL,
        type = c("response", "prob", "classify.weights"),
        threshold = NULL, ...)
}
\arguments{
  \item{object}{
    A fit object of class 'cpfa' from function 'cpfa'.
}
  \item{newdata}{
    An optional three-way or four-way data array used to predict Parafac factor weights using estimated Parafac model factor weights from inputted object. Dimensions must match dimensions of original data for all modes except classification mode. If omitted, the original data are used.
}
  \item{nfac}{
    Number of factors in models for prediction. Defaults to number of factors for each Parafac model from input 'object'.
}
  \item{method}{
   Classification methods used for prediction. Defaults to methods used in input 'object'.
}
  \item{type}{
  Character vector indicating type of prediction to return. Possible values include: (1) \code{"response"}, returning predicted class labels; (2) \code{"prob"}, returning predicted class probabilities; or (3) \code{"classify.weights"}, returning predicted factor weights used in classification from Parafac models specified. Defaults to \code{"response"}.
}
  \item{threshold}{
  For binary classification, value indicating prediction threshold over which observations are classified as the positive class. If not provided, calculates threshold using class proportions in original data. For multiclass classification, \code{threshold} is not currently implemented. 
}
  \item{...}{
additional arguments affecting the prediction produced (currently ignored).  
}
}

\details{
Predicts class labels for a binary or a multiclass outcome. Specifically, predicts factor weights for one mode of a Parallel Factor Analysis (Parafac) model-1 using newdata and previously estimated mode weights from original data. Passes predicted factor weights to one or several classification methods as new data for predicting class labels.

Tuning parameters optimized by k-fold cross-validation are used for each classification method (see help for \code{cpfa}). If not supplied in argument 'threshold', prediction threshold for all classification methods is calculated using proportions of class labels for original data in the binary case (and positive proportion is set as threshold). For multiclass case, class with highest probability is chosen. Calculates and returns only predicted probabilities without class assignment by specifying type = "prob". Returns only predicted factor weights if type = "classify.weights".
}
\value{

Returns one of the following: (1) a data frame 'storfac'; (2) a list 'storprob',
or (3) a list 'classify.weights':

  \item{storfac}{A data frame containing predicted class labels or probabilities (binary case) for each Parafac model and classification method selected (see argument 'type'). Number of columns is equal to number of methods times number of Parafac models. Number of rows is equal to number of predicted observations.
}

  \item{storprob}{A list containing predicted probabilities for each Parafac model and classification method selected (see argument 'type'). Only returned if original response was multiclass (3 or more class labels). Number of elements is equal to number of methods times number of Parafac models.
}
  \item{cweights}{List containing predicted factor weights for each Parafac model. Length is equal to number of Parafac models estimated.
}

}
\author{
 Matthew Snodgress <snodg031@umn.edu>
}
\references{

Breiman, L. (2001). Random forests. Machine learning, 45(1), 5-32.

Cortes, C., and Vapnik, V. (1995). Support-vector networks. Machine learning, 20(3), 273-297.

Friedman, J. Hastie, T., and Tibshirani, R. (2010). Regularization Paths for Generalized Linear Models via Coordinate Descent. Journal of Statistical Software, 33(1), 1-22.

Harshman, R. A. (1970). Foundations of the PARAFAC procedure: Models and conditions for an "explanatory" multimodal factor analysis. UCLA Working Papers in Phonetics, 16, 1-84.

Harshman, R. A. and Lundy, M. E. (1994). PARAFAC: Parallel factor analysis. Computational Statistics and Data Analysis, 18, 39-72.

Helwig, N. E. (2017). Estimating latent trends in multivariate longitudinal data via Parafac2 with functional and structural constraints. Biometrical Journal, 59(4), 783-803.

Helwig, N. E. (2019). multiway: Component Models for Multi-Way Data. R Package version 1.0-6.

Liaw, A., and Wiener, M. (2002). Classification and Regression by randomForest. R News 2(3), 18-22.

Meyer, D., Dimitriadou, E., Hornik, K., Weingessel, A., and Leisch, F. (2021). e1071: Misc Functions of the Department of Statistics, Probability Theory Group (Formerly: E1071), TU Wien. R package version 1.7-6.

Ripley, B. (1994). Neural networks and related methods for classification. Journal of the Royal Statistical Society: Series B (Methodological), 56(3), 409-437.

Venables, W. and Ripley, B. (2002) Modern Applied Statistics with S. Fourth Edition. Springer, New York. ISBN 0-387-95457-0.

Zou, H. and Hastie, T. (2005). Regularization and variable selection via the elastic net. Journal of the Royal Statistical Society: Series B (Statistical Methodology). 67(2), 301-320.

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
tune.object <- cpfa(x = X, y = y, nfac = 3, nfolds = nfolds,
                         foldid = foldid, method = method, alpha = alpha,
                         gamma = gamma, cost = cost, ntree = ntree,
                         nodesize = nodesize, size = size, decay = decay,
                         parallel = FALSE, const = const)
                         
# create random data array with Parafac structure and same A and B weights
mydim.new <- c(60, 16, 20)
Cmat <- matrix(runif(mydim.new[3]*nf), nrow = mydim.new[3], ncol = nf)
Xmat <- tcrossprod(Amat, krprod(Cmat, Bmat))
Xmat <- array(Xmat, dim = mydim.new)
Emat <- array(rnorm(prod(mydim.new)), dim = mydim.new)
Xnew <- Xmat + Emat

# predict class labels
predict.labels <- predict(object = tune.object, newdata = Xnew,
                                    type = "response")
head(predict.labels)

# Note: see package vignette for more examples
}