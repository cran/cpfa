\name{print.cpfa}
\alias{print.cpfa}
\title{
Deprecated - Print Method for Classification with Parallel Factor Analysis
}
\description{
Deprecated. Use function \code{print.tunecpfa} instead. See help file for \code{print.tunecpfa}.
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