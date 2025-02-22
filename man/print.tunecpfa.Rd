\name{print.tunecpfa}
\alias{print.tunecpfa}
\title{
Print Method for Tuning for Classification with Parallel Factor Analysis
}
\description{
Prints summary of results from a 'tunecpfa' model object obtained using function \code{tunecpfa}. 
}
\usage{
\method{print}{tunecpfa}(x, ...)
}
\arguments{
\item{x}{
  A fit object of class 'tunecpfa' from function \code{tunecpfa}.
}
\item{...}{
  Additional print arguments.
}
}
\details{
Prints names of the models and methods used to create the input 'tunecpfa' model object. Prints misclassification error rates and estimation times in seconds.
}
\value{
Returns a summary of the 'tunecpfa' model object.
}
\author{
 Matthew Snodgress <snodg031@umn.edu>
}
\references{
See help file for function \code{tunecpfa} for a list of references.
}
\examples{
########## Parafac example with 3-way array and binary response ##########
\dontrun{
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
method <- c("PLR", "SVM")
family <- "binomial"
parameters <- list(alpha = alpha, gamma = gamma, cost = cost)
model <- "parafac"
nfolds <- 3
nstart <- 3

# constrain first mode weights to be orthogonal
const <- c("orthog", "uncons", "uncons")

# fit Parafac models and use third mode to tune classification methods
tune.object <- tunecpfa(x = X, y = y, model = model, nfac = nf, 
                        nfolds = nfolds, method = method, family = family, 
                        parameters = parameters, parallel = FALSE, 
                        const = const, nstart = nstart)
                         
# print summary of output
print(tune.object)
}
}