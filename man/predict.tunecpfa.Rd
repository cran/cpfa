\name{predict.tunecpfa}
\alias{predict.tunecpfa}
\title{
Predict Method for Tuning for Classification with Parallel Factor Analysis
}
\description{
Obtains predictions for class labels from a 'tunecpfa' model object obtained using function \code{tunecpfa}.
}
\usage{
\method{predict}{tunecpfa}(object, newdata = NULL, method = NULL, 
        type = c("response", "prob", "classify.weights"), 
        threshold = NULL, ...)
}
\arguments{
\item{object}{
  A fit object of class 'tunecpfa' produced by function \code{tunecpfa}.
}
\item{newdata}{
  An optional three-way or four-way data array used to predict Parafac or Parafac2 component weights using estimated Parafac or Parafac2 model component weights from inputted object. For Parafac2, can be a list of length \code{K} where the \code{k}-th element is a matrix or three-way array associated with the \code{k}-th element. Array or list must contain only real numbers. Dimensions must match dimensions of original data for all modes except the classification mode. If omitted, the original data are used.
}
\item{method}{
  Character vector indicating classification methods to use. Possible methods include penalized logistic regression (PLR); support vector machine (SVM); random forest (RF); feed-forward neural network (NN); regularized discriminant analysis (RDA); and gradient boosting machine (GBM). If none selected, default is to use all methods.
}
\item{type}{
  Character vector indicating type of prediction to return. Possible values include: (1) \code{"response"}, returning predicted class labels; (2) \code{"prob"}, returning predicted class probabilities; or (3) \code{"classify.weights"}, returning predicted component weights used in classification from Parafac models specified. Defaults to \code{"response"}.
}
\item{threshold}{
  For binary classification, value indicating prediction threshold over which observations are classified as the positive class. If not provided, calculates threshold using class proportions in original data. For multiclass classification, \code{threshold} is not currently implemented.
}
\item{...}{
  Additional predict arguments. Currently ignored.
}
}

\details{
Predicts class labels for a binary or a multiclass outcome. Specifically, predicts component weights for one mode of a Parallel Factor Analysis-1 (Parafac) model or a Parallel Factor Analysis-2 (Parafac2) model using new data and previously estimated mode weights from original data. Passes predicted component weights to one or several classification methods as new data for predicting class labels.

Tuning parameters optimized by k-fold cross-validation are used for each classification method (see help for \code{tunecpfa}). If not supplied in argument \code{threshold}, prediction threshold for all classification methods is calculated using proportions of class labels for original data in the binary case (and the positive class proportion is set as the threshold). For multiclass case, class with highest probability is chosen.
}
\value{
Returns one of the following, depending on the choice for argument \code{type}:

\item{type = "response"}{A data frame containing predicted class labels or probabilities (binary case) for each Parafac model and classification method selected (see argument \code{type}). Number of columns is equal to number of methods times number of Parafac models. Number of rows is equal to number of predicted observations.
}
\item{type = "prob"}{A list containing predicted probabilities for each Parafac model and classification method selected (see argument \code{type}). Only returned if original response was multiclass (i.e., contained three or more class labels). The number of list elements is equal to the number of methods times the number of Parafac models.
}
\item{type = "classify.weights"}{List containing predicted component weights for each Parafac or Parafac2 model. Length is equal to number of Parafac models that were fit.
}
}
\author{
 Matthew Snodgress <snodg031@umn.edu>
}
\references{
See help file for function \code{tunecpfa} for a list of references.
}
\examples{
########## Parafac2 example with 4-way array and multiclass response ##########
\dontrun{
# set seed and specify dimensions of a four-way tensor
set.seed(5)
mydim <- c(10, 11, 12, 90)
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
rda.alpha <- seq(0.1, 0.9, length = 2)
delta <- c(0.1, 2)
eta <- c(0.3, 0.7)
max.depth <- c(1, 2)
subsample <- c(0.75)
nrounds <- c(100)
method <- c("RDA", "GBM")
family <- "multinomial"
parameters <- list(rda.alpha = rda.alpha, delta = delta, eta = eta,
                   max.depth = max.depth, subsample = subsample, 
                   nrounds = nrounds)
model <- "parafac2"
nfolds <- 3
nstart <- 3

# constrain first mode weights to be orthogonal, fourth mode to be nonnegative
const <- c("orthog", "uncons", "uncons", "nonneg")

# fit Parafac2 model and use fourth mode to tune classification methods
tune.object <- tunecpfa(x = X, y = y, model = model, nfac = nf, 
                        nfolds = nfolds, method = method, family = family, 
                        parameters = parameters, parallel = FALSE, 
                        const = const, nstart = nstart)

# create new data with Parafac2 structure and D weights related to response
mydim.new <- c(10, 11, 12, 10)
Znew <- matrix(rnorm(mydim.new[4] * (nf + 1)), nrow = mydim.new[4], 
               ncol = (nf + 1))
Xwnew <- rep(1, mydim.new[4]) \%*\% t(mu) + Znew \%*\% cormat.sqrt
Dmatnew <- Xwnew[, 1:nf]
Amat <- vector("list", mydim.new[4])
Xnew <- Xmat <- Emat <- Amat
for (Dd in 1:mydim.new[4]) {
   Amat[[Dd]] <- matrix(nf * rnorm(nDd[Dd]), nrow = nDd[Dd], ncol = nf)
   Amat[[Dd]] <- svd(Amat[[Dd]], nv = 0)$u \%*\% Gmat
   leftMat <- Amat[[Dd]] \%*\% diag(Dmatnew[Dd, ])
   Xmat[[Dd]] <- array(tcrossprod(leftMat, krprod(Cmat, Bmat)), 
                       dim = c(nDd[Dd], mydim.new[2], mydim.new[3]))
   Emat[[Dd]] <- array(rnorm(nDd[Dd] * mydim.new[2] * mydim.new[3]), 
                       dim = c(nDd[Dd], mydim.new[2], mydim.new[3]))
   Xnew[[Dd]] <- Xmat[[Dd]] + Emat[[Dd]]
}

# create new random class labels for two levels
stor <- matrix(rep(1, nrow(Xwnew)), nrow = nrow(Xwnew))
stor[which(Xwnew[, (nf + 1)] < (ymean - 0.4 * sd(Xwnew[, (nf + 1)])))] <- 2
stor[which(Xwnew[, (nf + 1)] > (ymean + 0.4 * sd(Xwnew[, (nf + 1)])))] <- 0
newlabels <- as.numeric(stor)

# predict class labels
predict.labels <- predict(object = tune.object, newdata = Xnew, 
                          type = "response")

# print predicted labels
predict.labels
}
}