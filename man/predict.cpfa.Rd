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
    An optional three-way or four-way data array used to predict Parafac or Parafac2 component weights using estimated Parafac or Parafac2 model component weights from inputted object. Dimensions must match dimensions of original data for all modes except the classification mode. If omitted, the original data are used.
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
########## Parafac example with 3-way array and binary response ##########

# set seed and specify dimensions of a three-way tensor
set.seed(3)
mydim <- c(10, 11, 80)
nf <- 3

# create correlation matrix between response and third mode's weights 
rho.cc <- .35 
rho.cy <- .75 
cormat.values <- c(1, rho.cc, rho.cc, rho.cy, rho.cc, 1, rho.cc, rho.cy, 
                   rho.cc, rho.cc, 1, rho.cy, rho.cy, rho.cy, rho.cy, 1)
cormat <- matrix(cormat.values, nrow = (nf + 1), ncol = (nf + 1))

# sample from a multivariate normal with specified correlation structure
ymean <- Cmean <- 2
mu <- as.matrix(c(Cmean, Cmean, Cmean, ymean))
eidecomp <- eigen(cormat, symmetric = TRUE)
L.sqrt <- diag(eidecomp$values^0.5)
cormat.sqrt <- eidecomp$vectors \%*\% L.sqrt \%*\% t(eidecomp$vectors)
Z <- matrix(rnorm(mydim[3]*(nf + 1)), nrow = mydim[3], ncol = (nf + 1))
Xw <- rep(1, mydim[3]) \%*\% t(mu) + Z \%*\% cormat.sqrt
Cmat <- Xw[, 1:nf]

# create a random three-way data tensor with C weights related to a response
Amat <- matrix(rnorm(mydim[1]*nf), nrow = mydim[1], ncol = nf)
Bmat <- matrix(runif(mydim[2]*nf), nrow = mydim[2], ncol = nf)
Xmat <- tcrossprod(Amat, krprod(Cmat, Bmat))
Xmat <- array(Xmat, dim = mydim)
Emat <- array(rnorm(prod(mydim)), dim = mydim)
Emat <- nscale(Emat, 0, ssnew = sumsq(Xmat))  
X <- Xmat + Emat

# create a binary response by dichotomizing at the specified response mean
y <- factor(as.numeric(Xw[ , (nf + 1)] > ymean))

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

# constrain first mode weights to be orthogonal
const <- c("orthog", "uncons", "uncons")

# fit Parafac models and use third mode to tune classification methods
tune.object <- tune.cpfa(x = X, y = y, model = model, nfac = nf, 
                         nfolds = nfolds, method = method, family = family, 
                         parameters = parameters, parallel = FALSE, 
                         const = const)
                         
# create new data with Parafac structure and C weights related to response
mydim.new <- c(10, 11, 20)
Znew <- matrix(rnorm(mydim.new[3]*(nf + 1)), 
               nrow = mydim.new[3], ncol = (nf + 1))
Xwnew <- rep(1, mydim.new[3]) \%*\% t(mu) + Znew \%*\% cormat.sqrt
Cmatnew <- Xwnew[, 1:nf]
Xnew0 <- tcrossprod(Amat, krprod(Cmatnew, Bmat))
Xnew0 <- array(Xnew0, dim = mydim.new)
Ematnew <- array(rnorm(prod(mydim.new)), dim = mydim.new)
Ematnew <- nscale(Ematnew, 0, ssnew = sumsq(Xnew0))  
Xnew <- Xnew0 + Ematnew

# predict class labels
predict.labels <- predict(object = tune.object, newdata = Xnew, 
                          type = "response")

# print predicted labels
predict.labels

########## Parafac2 example with 4-way array and multiclass response ##########

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
alpha <- seq(0, 1, length = 2)
gamma <- c(0, 1)
cost <- c(0.1, 5)
ntree <- c(200, 300)
nodesize <- c(1, 2)
size <- c(1, 2)
decay <- c(0, 1)
rda.alpha <- seq(0.1, 0.9, length = 2)
delta <- c(0.1, 1)
method <- c("PLR", "SVM", "RF", "NN", "RDA")
family <- "multinomial"
parameters <- list(alpha = alpha, gamma = gamma, cost = cost, ntree = ntree,
                   nodesize = nodesize, size = size, decay = decay, 
                   rda.alpha = rda.alpha, delta = delta)
model <- "parafac2"
nfolds <- 3
nstart <- 5

# constrain first mode weights to be orthogonal, fourth mode to be nonnegative
const <- c("orthog", "uncons", "uncons", "nonneg")

# fit Parafac2 model and use fourth mode to tune classification methods
tune.object <- tune.cpfa(x = X, y = y, model = model, nfac = nf, 
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