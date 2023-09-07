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

Returns a summary of 'cpfa' model object.

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