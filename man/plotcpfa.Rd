\name{plotcpfa}
\alias{plotcpfa}
\title{
Plot Optimal Model from Classification with Parallel Factor Analysis
}
\description{
Plots optimal model based on results from a 'wrapcpfa' object obtained using function \code{cpfa}. 
}
\usage{
plotcpfa(object, cmeasure = "acc", meanvalue = TRUE, supNum = FALSE, 
         parallel = FALSE, cl = NULL, scale.remode = NULL, newscales = 1, 
         scale.abmode = NULL, sign.remode = NULL, newsigns = 1, 
         sign.abmode = NULL, ...)
}
\arguments{
\item{object}{
  An object of class 'wrapcpfa' from function \code{cpfa}.
}
\item{cmeasure}{
  Classification performance measure used to select the optimal number of components. Options include \code{c("err", "acc", "tpr", "fpr", "tnr", "fnr", "ppv", "npv", "fdr", "fom", "fs")}. If \code{cmeasure} is in \code{c("err", "fpr", "fnr", "fdr", "fom")}, the number of components that minimized \code{cmeasure} is selected among all classification methods. Otherwise, the number that maximized \code{cmeasure} is selected.
}
\item{meanvalue}{
  Logical indicating whether to find the optimal number of components based on the mean performance across replications from the results generated by \code{cpfa}. If \code{meanvalue = FALSE}, the median is used. 
}
\item{supNum}{
  Logical indicating whether to suppress text displaying component weight values within plot cells. If TRUE, values are not displayed.
}
\item{parallel}{
  Logical indicating if parallel computing should be implemented. If TRUE, parallel computing is used.
}
\item{cl}{
  Cluster for parallel computing, which is used when \code{parallel = TRUE}. Note that if \code{parallel = TRUE} and \code{cl = NULL}, then the cluster is defined as \code{makeCluster(detectCores())}.
}
\item{scale.remode}{
  Character that indicates a mode to rescale. Must be one of \code{c("A", "B", "C", "D")}. Sent directly to argument \code{mode} in function \code{rescale} from package \bold{multiway}. See help file for \code{rescale} for additional details.
}
\item{newscales}{
  The root mean-square for columns of the mode indicated by \code{scale.remode}. See help file for \code{rescale} for additional details.
}
\item{scale.abmode}{
  Character that indicates the mode that absorbs the inverse of rescalings applied to the mode indicated by \code{scale.remode}. Must be one of \code{c("A", "B", "C", "D")}. Sent directly to argument \code{absorb} in function \code{rescale} from package \bold{multiway}. See help file for \code{rescale} for additional details.
}
\item{sign.remode}{
  Character that indicates a mode to resign. Must be one of \code{c("A", "B", "C", "D")}. Sent directly to argument \code{mode} in function \code{resign} from package \bold{multiway}. See help file for \code{resign} for additional details.
}
\item{newsigns}{
  Scalar or vector indicating resignings for columns of the mode indicated by \code{sign.remode}. See help file for \code{resign} for additional details.
}
\item{sign.abmode}{
  Character that indicates the mode that absorbs the negation of the resignings applied to the mode indicated by \code{sign.remode}. Must be one of \code{c("A", "B", "C", "D")}. Sent directly to argument \code{absorb} in function \code{resign} from package \bold{multiway}. See help file for \code{resign} for additional details.
}
\item{...}{
  Additional arguments to be passed to function \code{parafac} for fitting a Parafac model or function \code{parafac2} for fitting a Parafac2 model. See help file for function \code{parafac} or for function \code{parafac2} for additional details.
}
}
\details{
Selects the number of components that optimized a performance measure across all classification methods used by \code{cpfa}. With this optimal number of components, fits the Parafac or Parafac2 model that was used by \code{cpfa} to create the input 'wrapcpfa' object. Uses same constraints used in \code{cpfa}. Plots component weights for this optimal model using heatmaps. Darker red indicates component weights that are more negative while darker green indicates component weights that are more positive. For three-way Parafac, plots A and B weights. For four-way Parafac, plots A, B, and C weights. For three-way Parafac2, plots B weights. For four-way Parafac2, plots B and C weights.
}
\value{
Returns one or more heatmap plots of component weights for the optimal Parafac or Parafac2 model. Returns list of estimated component weights used in the plots.
}
\author{
 Matthew Snodgress <snodg031@umn.edu>
}
\references{
See help file for function \code{cpfa} for a list of references.
}
\examples{
########## Parafac2 example with 4-way array and multiclass response ##########

# set seed and specify dimensions of a four-way tensor
set.seed(5)
mydim <- c(10, 11, 12, 100)
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
rda.alpha <- seq(0.1, 0.9, length = 2)
delta <- c(0.1, 2)
method <- c("PLR", "SVM", "RDA")
family <- "multinomial"
parameters <- list(alpha = alpha, gamma = gamma, cost = cost, 
                   rda.alpha = rda.alpha, delta = delta)
model <- "parafac2"
nfolds <- 3
nstart <- 1

# constrain first mode weights to be orthogonal, fourth mode to be nonnegative
const <- c("orthog", "uncons", "uncons", "nonneg")

# fit Parafac2 model and use fourth mode weights to tune classification
# methods, to predict class labels, and to return classificaiton 
# performance measures pooled across multiple train-test splits
output <- cpfa(x = X, y = y, model = model, nfac = nf, nrep = 2, ratio = 0.8, 
               nfolds = nfolds, method = method, family = family, 
               parameters = parameters, type.out = "descriptives", 
               seeds = NULL, plot.out = TRUE, parallel = FALSE, const = const,
               nstart = nstart, ctol = 1e-2)

# plot heatmap of component weights for optimal model
plotcpfa(output, nstart = nstart, ctol = 1e-2)
}