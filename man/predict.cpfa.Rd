\name{predict.cpfa}
\alias{predict.cpfa}
\title{
Predict Method for Classification with Parallel Factor Analysis
}
\description{
Obtains predictions for class labels from a 'cpfa' model object obtained using function \code{tune.cpfa}.
}
\usage{
\method{predict}{cpfa}(object, newdata = NULL, method = NULL,
        type = c("response", "prob", "classify.weights"),
        threshold = NULL, ...)
}
\arguments{
  \item{object}{
    A fit object of class 'cpfa' from function \code{tune.cpfa}.
}
  \item{newdata}{
    An optional three-way or four-way data tensor/array used to predict Parafac or Parafac2 component weights using estimated Parafac or Parafac2 model component weights from inputted object. Dimensions must match dimensions of original data for all modes except the classification mode. If omitted, the original data are used.
}
  \item{method}{
  Character vector indicating classification methods to use. Possible methods include penalized logistic regression (PLR); support vector machine (SVM); random forest (RF); feed-forward neural network (NN); and regularized discriminant analysis (RDA). If none selected, default is to use all methods.
}
  \item{type}{
  Character vector indicating type of prediction to return. Possible values include: (1) \code{"response"}, returning predicted class labels; (2) \code{"prob"}, returning predicted class probabilities; or (3) \code{"classify.weights"}, returning predicted component weights used in classification from Parafac models specified. Defaults to \code{"response"}.
}
  \item{threshold}{
  For binary classification, value indicating prediction threshold over which observations are classified as the positive class. If not provided, calculates threshold using class proportions in original data. For multiclass classification, \code{threshold} is not currently implemented. 
}
  \item{...}{
  Currently ignored. Additional predict arguments.
}
}

\details{
Predicts class labels for a binary or a multiclass outcome. Specifically, predicts component weights for one mode of a Parallel Factor Analysis-1 (Parafac) model or a Parallel Factor Analysis-2 (Parafac2) model using new data and previously estimated mode weights from original data. Passes predicted component weights to one or several classification methods as new data for predicting class labels.

Tuning parameters optimized by k-fold cross-validation are used for each classification method (see help for \code{tune.cpfa}). If not supplied in argument \code{threshold}, prediction threshold for all classification methods is calculated using proportions of class labels for original data in the binary case (and the positive class proportion is set as the threshold). For multiclass case, class with highest probability is chosen.
}
\value{
Returns one of the following, depending on the choice for argument \code{type}:

  \item{type = "response"}{A data frame containing predicted class labels or probabilities (binary case) for each Parafac model and classification method selected (see argument \code{type}). Number of columns is equal to number of methods times number of Parafac models. Number of rows is equal to number of predicted observations.
}
  \item{type = "prob"}{A list containing predicted probabilities for each Parafac model and classification method selected (see argument \code{type}). Only returned if original response was multiclass (i.e., contained three or more class labels). The number of list elements is equal to number of methods times the number of Parafac models.
}
  \item{type = "classify.weights"}{List containing predicted component weights for each Parafac or Parafac2 model. Length is equal to number of Parafac models that were fit.
}
}
\author{
 Matthew Snodgress <snodg031@umn.edu>
}
\references{
Breiman, L. (2001). Random forests. Machine Learning, 45(1), 5-32.

Cortes, C. and Vapnik, V. (1995). Support-vector networks. Machine Learning, 20(3), 273-297.

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
rda.alpha <- c(0.1, 0.6)
delta <- c(0.1, 2)
method <- c("PLR", "SVM", "RF", "NN", "RDA")
family <- "binomial"
parameters <- list(alpha = alpha, gamma = gamma, cost = cost, ntree = ntree,
                   nodesize = nodesize, size = size, decay = decay, 
                   rda.alpha = rda.alpha, delta = delta)
model <- "parafac"
nfolds <- 3
const <- c("orthog", "uncons", "uncons")

# fit Parafac models and use third mode to tune classification methods
tune.object <- tune.cpfa(x = X, y = y, model = model, nfac = 3, nfolds = nfolds, 
                         method = method, family = family, 
                         parameters = parameters, parallel = FALSE, 
                         const = const)
                         
# create random data array with Parafac structure and with the A 
# and B weights used above
mydim.new <- c(60, 16, 20)
Cmat <- matrix(runif(mydim.new[3]*nf), nrow = mydim.new[3], ncol = nf)
Xmat <- tcrossprod(Amat, krprod(Cmat, Bmat))
Xmat <- array(Xmat, dim = mydim.new)
Emat <- array(rnorm(prod(mydim.new)), dim = mydim.new)
Xnew <- Xmat + Emat

# predict class labels
predict.labels <- predict(object = tune.object, newdata = Xnew, 
                          type = "response")
                          
# view predicted labels               
head(predict.labels)
}