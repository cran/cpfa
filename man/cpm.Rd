\name{cpm}
\alias{cpm}
\title{
Classification Performance Measures
}
\description{
Calculates multiple performance measures for binary or multiclass classification. Uses known class labels and evaluates against predicted labels.
}
\usage{
cpm(x, y, level = NULL, fbeta = NULL, prior = NULL)
}
\arguments{
  \item{x}{
  Known class labels of class numeric, factor, or integer. If factor, converted to class integer in order of factor levels with integers beginning at 0 (i.e. for binary classification, factor levels become 0 and 1; for multiclass, 0, 1, 2, etc.). 
}
  \item{y}{
  Predicted class labels of class numeric, factor, or integer. If factor, converted to class integer in order of factor levels with integers beginning at 0 (i.e. for binary classification, factor levels become 0 and 1; for multiclass, 0, 1, 2, etc.).
}
  \item{level}{
  Optional argument specifying possible class labels. For cases when \code{x} or \code{y} do not contain all possible classes. Can be of class numeric, integer, or character. Must contain two elements for binary classification, and contain three or more elements for multiclass classification. If integer, integers should be ordered (e.g. \code{c(0, 1)}; or \code{c(0, 1, 2)}). Note: if both \code{x} and \code{y} jointly contain only a single value (e.g. 1), must specify argument \code{level} in order to identify classification as binary or multiclass. 
}

  \item{fbeta}{
  Optional numeric argument specifying beta value for F-score. Defaults to \code{fbeta = 1}, providing an F1-score (balanced harmonic mean between precision and recall).
}
  \item{prior}{
  Optional numeric argument specifying weights for classes. Defaults to \code{prior = c(rep(1/llev, llev))}, where 'llev' is the number of classes, providing equal importance across classes.
  }
}
\details{
Selecting one class as a negative class and one class as positive, binary classification generates four possible outcomes: (1) negative cases classified as positives, called false positives (FP); (2) negative cases classified as negatives, called true negatives (TN); (3) positive cases classified as negatives, called false negatives (FN); and (4) positive cases classified as positives, called true positives (TP). 

Multiple evaluation measures are calculated using these four outcomes. Measures include: overall error (ERR), also called fraction incorrect;  overall accuracy (ACC), also called fraction correct; true positive rate (TPR), also called recall, hit rate, or sensitivity; false negative rate (FNR), also called miss rate; false positive rate (FPR), also called fall-out; true negative rate (TNR), also called specificity or selectivity; positive predictive value (PPV), also called precision; false discovery rate (FDR); negative predictive value (NPV); false omission rate (FOR); and F-score (FS).

In multiclass classification, the four outcomes are possible for each individual class in macro-averaging, and evaluation measures are averaged over classes. Macro-averaging gives equal importance to all classes. For multiclass classification, calculated measures are currently only macro-averaged.

Note that binary classification assumes a positive class and negative class (i.e. contains a reference group) and is ordered. Multiclass classification is currently unordered. 

Computational details:

ERR = (FP + FN) / (TP + TN + FP + FN).

ACC = (TP + TN) / (TP + TN + FP + FN), and ACC = 1 - ERR.

TPR = TP / (TP + FN).

FNR = FN / (FN + TP), and FNR = 1 - TPR.

FPR = FP / (FP + TN).

TNR = TN / (TN + FP), and TNR = 1 - FPR.

PPV = TP / (TP + FP).

FDR = FP / (FP + TP), and FDR = 1 - PPV.

NPV = TN / (TN + FN).

FOR = FN / (FN + TN), and FOR = 1 - NPV.

FS = (1 + beta^2) * ((PPV * TPR) / (((beta^2)*PPV) + TPR)).

All evaluation measures calculated are between 0 and 1, inclusive. For multiclass classification, macro-averaged values are provided for each performance measure. Note that 'beta' in FS represents the relative weight such that recall (TPR) is beta times more important than precision (PPV). See reference in help document for more details.
}
\value{

  Returns list where first element is a full confusion matrix \code{cm} and second element is a data frame containing performance measures. Note that for multiclass classification, macro-averaged values are provided (each measure calculated for each class, then averaged over all classes; average is weighted by argument \code{prior} if provided). Excluding \code{cm}, second list element contains following performance measures:

  \item{cm}{
  A confusion matrix with counts for each of possible outcomes.
}
  \item{err}{
  Overall error (ERR). Also called fraction incorrect.
}
  \item{acc}{
  Overall accuracy (ACC). Also called fraction correct.
}
  \item{tpr}{
  True positive rate (TPR). Also called recall, hit rate, or sensitivity.
}
  \item{fpr}{
  False positive rate (FPR). Also called fall-out.
}
  \item{tnr}{
  True negative rate (TNR). Also called specificity or selectivity.
}
  \item{fnr}{
  False negative rate (FNR). Also called miss rate.
}
  \item{ppv}{
  Positive predictive value (PPV). Also called precision.
}
  \item{npv}{
  Negative predicted value (NPV).
}
  \item{fdr}{
  False discovery rate (FDR).
}
  \item{fom}{
  False omission rate (FOR).
}
  \item{fs}{
  F-score. Mean between TPR (recall) and PPV (precision) varying by importance 
given to recall over precision (see Details section and argument \code{fbeta}). 
}

}
\references{
Sokolova, M., and Lapalme, G. (2009). A systematic analysis of performance measures for classification tasks. Information processing and management, 45(4), 427-437.
}
\author{
Matthew Snodgress <snodg031@umn.edu>
}
\examples{
########## Random binary example with 3-way array ##########

# create random data array with Parafac structure and random class labels
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
method <- c("PLR", "SVM", "RF")
nfolds <- 3
const <- c("orthog", "uncons", "uncons")
foldid <- sample(rep(1:nfolds, length.out = length(y)))

# estimate Parafac models and use third mode to tune classification methods
tune.object <- cpfa(x = X, y = y, nfac = 3, nfolds = nfolds,
                    foldid = foldid, method = method, alpha = alpha,
                    gamma = gamma, cost = cost, ntree = ntree,
                    nodesize = nodesize, parallel = FALSE, const = const)

# create random data array with Parafac structure and same A and B weights
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
 
# calculate evaluation measure for predicted
y.pred <- predict.labels[,1]
evalmeasure <- cpm(x = newlabel, y = y.pred)
evalmeasure
}