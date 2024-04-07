\name{predict.cpfa}
\alias{predict.cpfa}
\title{
Deprecated - Predict Method for Classification with Parallel Factor Analysis
}
\description{
Deprecated. Use function \code{predict.tunecpfa} instead. See help file for \code{predict.tunecpfa}.
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
See help file for \code{predict.tunecpfa}.
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