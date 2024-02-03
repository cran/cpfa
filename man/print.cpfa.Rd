\name{print.cpfa}
\alias{print.cpfa}
\title{
Print Method for Classification with Parallel Factor Analysis
}
\description{
Prints summary of results from a 'cpfa' model object obtained using function \code{tune.cpfa}. 
}
\usage{
\method{print}{cpfa}(x, ...)
}
\arguments{
  \item{x}{
    A fit object of class 'cpfa' from function \code{tune.cpfa}.
}
  \item{...}{
  Additional print arguments.
}
}

\details{
Prints names of the models and methods used to create the input 'cpfa' model object. Prints misclassification error rates and estimation times in seconds.
}
\value{
Returns a summary of the 'cpfa' model object.
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
mydim <- c(60, 16, 80)
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
const <- c("orthog", "uncons", "uncons")

# fit Parafac models and use third mode to tune classification methods
tune.object <- tune.cpfa(x = X, y = y, model = model, nfac = nf, 
                         nfolds = nfolds, method = method, family = family, 
                         parameters = parameters, parallel = FALSE, 
                         const = const)
                         
# print summary of output
print(tune.object)
}