\name{cpm.all}
\alias{cpm.all}
\title{
Wrapper for Calculating Classification Performance Measures
}
\description{
Applies function \code{cpm} to multiple sets of class labels. Each set of class labels is evaluated against the same set of predicted labels. Works with output from function \code{predict.cpfa} and calculates classification performance measures for multiple classifiers or numbers of components.
}
\usage{
cpm.all(x, y,  ...)
}
\arguments{
  \item{x}{
  A data frame where each column contains a set of known class labels of class numeric, factor, or integer. If a set is of class factor, that set is converted to class integer in the order of factor levels with integers beginning at 0 (i.e., for binary classification, factor levels become 0 and 1; for multiclass, levels become 0, 1, 2, etc.). 
}
  \item{y}{
  Predicted class labels of class numeric, factor, or integer. If factor, converted to class integer in order of factor levels with integers beginning at 0 (i.e., for binary classification, factor levels become 0 and 1; for multiclass, 0, 1, 2, etc.).
}
 \item{...}{
  Additional arguments to be passed to function \code{cpm} for calculating classification performance measures.
}
}
\details{
Wrapper function that applies function \code{cpm} to multiple sets of class labels and one set of predicted labels. See help file for function \code{cpm} for additional details.
}
\value{
Returns a list with the following two elements: 
  \item{cm.list}{
  A list of confusion matrices, denoted \code{cm}, where each confusion matrix is associated with one comparison.
}
  \item{cpms}{
  A data frame containing classification performance measures where each row contains measures for one comparison.
}
}
\references{
Sokolova, M. and Lapalme, G. (2009). A systematic analysis of performance measures for classification tasks. Information Processing and Management, 45(4), 427-437.
}
\author{
Matthew Snodgress <snodg031@umn.edu>
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
Z <- matrix(rnorm(mydim[3] * (nf + 1)), nrow = mydim[3], ncol = (nf + 1))
Xw <- rep(1, mydim[3]) \%*\% t(mu) + Z \%*\% cormat.sqrt
Cmat <- Xw[, 1:nf]

# create a random three-way data tensor with C weights related to a response
Amat <- matrix(rnorm(mydim[1] * nf), nrow = mydim[1], ncol = nf)
Bmat <- matrix(runif(mydim[2] * nf), nrow = mydim[2], ncol = nf)
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
Znew <- matrix(rnorm(mydim.new[3] * (nf + 1)), 
               nrow = mydim.new[3], ncol = (nf + 1))
Xwnew <- rep(1, mydim.new[3]) \%*\% t(mu) + Znew \%*\% cormat.sqrt
Cmatnew <- Xwnew[, 1:nf]
Xnew0 <- tcrossprod(Amat, krprod(Cmatnew, Bmat))
Xnew0 <- array(Xnew0, dim = mydim.new)
Ematnew <- array(rnorm(prod(mydim.new)), dim = mydim.new)
Ematnew <- nscale(Ematnew, 0, ssnew = sumsq(Xnew0))  
Xnew <- Xnew0 + Ematnew

# create new random class labels for two levels
newlabel <- as.numeric(Xwnew[, (nf + 1)] > ymean)

# predict class labels
predict.labels <- predict(object = tune.object, newdata = Xnew, 
                          type = "response")
                        
# calculate performance measures for predicted class labels
evalmeasure <- cpm.all(x = predict.labels, y = newlabel)

# print performance measures
evalmeasure
}