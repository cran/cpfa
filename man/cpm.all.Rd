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

# create new random class labels for two levels
newlabel <- as.numeric(factor(rbinom(mydim.new[3], 1, 0.4))) - 1 
 
# calculate performance measures for predicted class labels
evalmeasure <- cpm.all(x = predict.labels, y = newlabel)
evalmeasure
}