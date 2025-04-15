\name{cpm.all}
\alias{cpm.all}
\title{
Wrapper for Calculating Classification Performance Measures
}
\description{
Applies function \code{cpm} to multiple sets of class labels. Each set of class 
labels is evaluated against the same set of predicted labels. Works with output 
from function \code{predict.tunecpfa} and calculates classification performance 
measures for multiple classifiers or numbers of components.
}
\usage{
cpm.all(x, y, ...)
}
\arguments{
\item{x}{
  A data frame where each column contains a set of class labels of class 
  numeric, factor, or integer. If a set is of class factor, that set is 
  converted to class integer in the order of factor levels with integers 
  beginning at 0 (i.e., for binary classification, factor levels become 0 and 
  1; for multiclass, levels become 0, 1, 2, etc.). 
}
\item{y}{
  Class labels of class numeric, factor, or integer. If factor, converted to 
  class integer in order of factor levels with integers beginning at 0 (i.e., 
  for binary classification, factor levels become 0 and 1; for multiclass, 0, 
  1, 2, etc.).
}
\item{...}{
  Additional arguments passed to function \code{cpm} for calculating 
  classification performance measures.
}
}
\details{
Wrapper function that applies function \code{cpm} to multiple sets of class 
labels and one other set of labels. See help file for function \code{cpm} for 
additional details.
}
\value{
Returns a list with the following two elements: 
\item{cm.list}{
  A list of confusion matrices, denoted \code{cm}, where each confusion matrix 
  is associated with one comparison.
}
\item{cpms}{
  A data frame containing classification performance measures where each row 
  contains measures for one comparison.
}
}
\references{
Sokolova, M. and Lapalme, G. (2009). A systematic analysis of performance 
measures for classification tasks. Information Processing and Management, 
45(4), 427-437.
}
\author{
Matthew Snodgress <mattgress@protonmail.ch>
}
\examples{
########## Parafac example with 3-way array and binary response ##########
\dontrun{
# set seed and simulate a three-way array related to a binary response
set.seed(5)
cormat <- matrix(c(1, .35, .35, .35, 1, .35, .35, .35, 1), nrow = 3, ncol = 3)
data <- simcpfa(arraydim = c(10, 11, 100), model = "parafac", nfac = 3, 
                nclass = 2, nreps = 1e2, onreps = 10, corresp = rep(.75, 3), 
                meanpred = rep(2, 3), modes = 3, corrpred = cormat)

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
tune.object <- tunecpfa(x = data$X[, , 1:80], y = as.factor(data$y[1:80, ]), 
                        model = model, nfac = 3, nfolds = nfolds, 
                        method = method, family = family, 
                        parameters = parameters, parallel = FALSE, 
                        const = const, nstart = nstart)
                    
# predict class labels
predict.labels <- predict(object = tune.object, newdata = data$X[, , 81:100], 
                          type = "response")
                        
# calculate performance measures for predicted class labels
evalmeasure <- cpm.all(x = predict.labels, y = as.numeric(data$y[81:100, ]))

# print performance measures
evalmeasure
}
}