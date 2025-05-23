\name{simcpfa}
\alias{simcpfa}
\title{
Simulate Data for Classification with Parallel Factor Analysis
}
\description{
Simulates a three-way or four-way data array and a set of class labels that are 
related to the simulated array through one mode of the array. Data array is 
simulated using either a Parafac or Parafac2 model with no constraints. Weights 
for mode weight matrices can be drawn from 12 common probability 
distributions. Alternatively, custom weights can be provided for any mode.
}
\usage{
simcpfa(arraydim = NULL, model = "parafac", nfac = 2, nclass = 2, nreps = 100, 
        onreps = 10, corresp = c(0.3, -0.3), meanpred = c(0, 0), modes = 3, 
        corrpred = matrix(c(1, 0.2, 0.2, 1), nrow = 2), pf2num = NULL, 
        Amat = NULL, Bmat = NULL, Cmat = NULL, Dmat = NULL, Gmat = NULL, 
        Emat = NULL, technical = list())
}
\arguments{
\item{arraydim}{
  Numeric vector containing the number of dimensions for each mode of the 
  simulated data array. Must contain integers greater than or equal to 2.
}
\item{model}{
  Character specifying the model to use for simulating the data array. Must be 
  either 'parafac' or 'parafac2'.
}
\item{nfac}{
  Number of components in the Parafac or Parafac2 model. Must be an integer 
  greater than or equal to 1. 
}
\item{nclass}{
  Number of classes in simulated class labels. Must be an integer greater than 
  or equal to 2.
}
\item{nreps}{
  Number of replications for simulating class labels for a given set of 
  classification mode component weights. 
}
\item{onreps}{
  Number of replications for simulating a set of classification mode component 
  weights. 
}
\item{corresp}{
  Numeric vector of target correlations between simulated class labels and 
  columns of the classification mode component weight matrix. Must have length 
  equal to 'nfac'.
}
\item{meanpred}{
  Numeric vector of means used to generate the classification mode component 
  weights. Must be real numbers. Operates as the mean vector parameterizing a 
  multivariate normal distribution from which classification mode component 
  weights are generated. Length must be equal to input \code{nfac}.
}
\item{modes}{
  Single integer of either 3 or 4, indicating whether to simulate a three-way 
  or four-way data array, respectively.
}
\item{corrpred}{
  A positive definite correlation matrix containing the target correlations for 
  the classification mode component weights. Must have number of rows and 
  columns equal to input 'nfac'. Operates as the covariance matrix 
  parameterizing a multivariate normal distribution from which classification 
  mode component weights are generated.
}
\item{pf2num}{
  When \code{model = 'parafac2'}, number of rows for each simulated matrix in 
  the list of matrices \code{Amat}. Replaces the first element of input 
  \code{arraydim} because, for the Parafac2 model, the number of rows in each 
  simulated matrix can vary. If not specified when \code{model = 'parafac2'}, 
  defaults to \code{rep(c((nfac + 1), (nfac + 2), (nfac + 3)), 
  length.out = arraydim[modes])}.
}
\item{Amat}{
  When \code{model = 'parafac'}, a matrix of A mode weights with number of rows 
  equal to the first element of input 'arraydim' and with number of columns 
  equal to input 'nfac'. When \code{model = 'parafac2'}, a list with length 
  equal to the last element of input 'arraydim', where each list element 
  contains a matrix with number of rows of at least 2 and with number 
  of columns equal to input 'nfac'. When provided, replaces a 
  simulated \code{Amat}.
}
\item{Bmat}{
  A matrix of B mode weights with number of rows equal to the second element of 
  input 'arraydim' and with number of columns equal to the input 'nfac'. When 
  provided, replaces a simulated \code{Bmat}.
}
\item{Cmat}{
  A matrix of C mode weights with number of rows equal to the third element of 
  input 'arraydim' and with number of columns equal to the input 'nfac'. When 
  provided, replaces a simulated \code{Cmat} when \code{modes = 4}. When
  \code{modes = 3}, replaces the simulated classification mode weight matrix. If
  provided when \code{modes = 3}, \code{onreps} is reduced to one.
}
\item{Dmat}{
  A matrix of D mode weights with number of rows equal to the fourth element of 
  input 'arraydim' and with number of columns equal to the input 'nfac'. When
  \code{modes = 4}, replaces the simulated classification mode weight matrix.
  When \code{modes = 3}, this argument is ignored.  If provided when 
  \code{modes = 4}, \code{onreps} is reduced to one.
}
\item{Gmat}{
  When \code{model = 'parafac2'}, a matrix of G mode weights with number of rows 
  equal to input 'nfac' and with number of columns equal to input 'nfac'. When 
  provided, replaces a simulated \code{Gmat}.
}
\item{Emat}{
  When \code{model = 'parafac'}, an array containing noise to be added to the
  corresponding elements in the simulated data array. Error array dimensions
  must be equal to the values contained in \code{arraydim}. When
  \code{model = 'parafac2'}, a list containing either matrices 
  (i.e., when \code{modes = 3}) or three-way arrays (i.e., when 
  \code{modes = 4}) whose elements contain noise to be added to corresponding 
  elements in the simulated data array. When provided, replaces a simulated
  \code{Emat}.
}
\item{technical}{
  List containing arguments related to distributions from which to simulate 
  data. When specified, must contain one or more of the following:
  \describe{
\item{distA}{
  List containing arguments specifying the distribution from which deviates 
  are drawn for A mode weights contained in \code{Amat}. Defaults to standard 
  normal distribution when not specified. See Details section for additional 
  information on acceptable arguments.
}
\item{distB}{
  List containing arguments specifying the distribution from which deviates 
  are drawn for B mode weights contained in \code{Bmat}. Defaults to standard 
  normal distribution when not specified. See Details section for additional 
  information on acceptable arguments.
}
\item{distC}{
  For when \code{modes = '4'}, list containing arguments specifying the distribution from which deviates 
  are drawn for C mode weights contained in \code{Cmat}. Defaults to standard 
  normal distribution when not specified. See Details section for additional 
  information on acceptable arguments.
}
\item{distG}{
  For when \code{model = 'parafac2'}, list containing arguments specifying the 
  distribution from which deviates are drawn for G weights contained 
  in \code{Gmat}. Defaults to standard normal distribution when not specified. 
  See Details section for additional information on acceptable arguments.
}
\item{distE}{
  List containing arguments specifying the distribution from which deviates 
  are drawn for error contained in \code{Emat}. Defaults to standard 
  normal distribution when not specified. See Details section for additional 
  information on acceptable arguments.
}
}
}
}
\details{
  Data array simulation consists of two steps. First, a Monte Carlo simulation
  is conducted to simulate class labels using a binomial logistic (i.e., 
  in the binary case) or multinomial logistic (i.e., in the multiclass case) 
  regression model. Specifically, columns of the classification mode weights 
  matrix (e.g., \code{Cmat} when \code{modes = 3}) are generated from a 
  multivariate normal distribution with mean vector \code{meanpred} and with 
  covariance matrix \code{corrpred}. Values are then drawn randomly from a 
  uniform or a normal distribution and serve as beta coefficients. A linear 
  combination of these beta coefficients and the generated classification 
  weights produces a linear systematic part, which is passed through the 
  logistic function (i.e., the sigmoid) in the binary case or through the 
  softmax function in the multiclass case. Resulting probabilities are used 
  to assign class labels. The simulation repeats classification weights 
  generation \code{onreps} times and repeats class label generation, within 
  each \code{onreps} iteration, a total of \code{nreps} times. The generated 
  class labels that correlate best with the generated classification weights 
  (i.e., with correlations closest to \code{corresp}) are retained as the final 
  class labels with corresponding final classification weights. An adaptive 
  sampling technique is used during the simulation such that optimal beta 
  coefficients from previous iterations are used to parameterize a normal 
  distribution, from which new coefficients are drawn in subsequent iterations. 
  Note that, if any simulation replicate produces a set of class labels where 
  all labels are the same (i.e., have no variance), that replicate is discarded.
  Note also that \code{onreps} is ignored when the classification mode weight 
  matrix (i.e., \code{Cmat} when \code{modes = 3} or \code{Dmat} when 
  \code{modes = 4}) is provided; in this case, class labels are simulated with 
  respect to the provided classification mode weight matrix.
  
  Second, depending on the chosen model (i.e., Parafac or Parafac2) specified 
  via \code{model}, and depending on the number of modes specified 
  via \code{modes}, component matrices are randomly generated for each mode 
  of the data array. A data array is then constructed using a Parafac or 
  Parafac2 structure from these weight matrices, including the generated 
  classification mode weight matrix (i.e., \code{Cmat} or \code{Dmat}) from the 
  first step. Alternatively, weight matrices can be provided to override random 
  generation for any weight matrix with the exception of the classification 
  mode. When provided, weight matrices are used to form the final data array. 
  Finally, random noise is added to each value in the array. 
  The resulting output is a synthetic multiway data array 
  paired, through one mode of the array, with a simulated binary or multiclass 
  response.
  
  The \code{technical} argument controls the probability distributions used to 
  simulate weights for different modes. Currently, \code{technical} is highly 
  structured. In particular, \code{technical} must be provided as a named list 
  whose elements must be one of 'distA', 'distB', 'distC', 'distG', or 'distE', 
  with the last letter of each name designating a mode or, in the case of 
  'distE', designating error. Each element provided must itself be a list where 
  the first inner list element is named 'dname', specifying the distribution to 
  be used to generate weights for a given mode or for error. There are 12 
  'dname' options: 'normal', 'uniform', 'gamma', 'beta', 'binomial', 'poisson', 
  'exponential', 'geometric', 'negbinomial', 'hypergeo', 'lognormal', and 
  'cauchy'. Additional arguments can be added to each inner list to parameterize 
  the probability distribution being used. These arguments can be one of the 
  following, for each distribution allowed: 
  
  For \code{dname = 'normal'}, allowed arguments are \code{mean} or 
  \code{sd} (i.e., function \code{rnorm} is called).
  
  For \code{dname = 'uniform'}, allowed arguments are \code{min} or 
  \code{max} (i.e., function \code{runif} is called).

  For \code{dname = 'gamma'}, allowed arguments are \code{shape} or 
  \code{scale} (i.e., function \code{rgamma} is called).
  
  For \code{dname = 'beta'}, allowed arguments are \code{shape1} or 
  \code{shape2} (i.e., function \code{rbeta} is called).

  For \code{dname = 'binomial'}, allowed arguments are \code{size} or 
  \code{prob} (i.e., function \code{rbinom} is called).

  For \code{dname = 'poisson'}, allowed argument is \code{lambda} (i.e., 
  function \code{rpois} is called).

  For \code{dname = 'exponential'}, allowed argument is \code{rate} (i.e., 
  function \code{rexp} is called).
  
  For \code{dname = 'geometric'}, allowed argument is \code{prob} (i.e., 
  function \code{rgeom} is called).
  
  For \code{dname = 'negbinomial'}, allowed arguments are \code{size} or 
  \code{prob} (i.e., function \code{rnbinom} is called).

  For \code{dname = 'hypergeo'}, allowed arguments are \code{m}, \code{n}, or 
  \code{k} (i.e., function \code{rhyper} is called).

  For \code{dname = 'lognormal'}, allowed arguments are \code{meanlog} or 
  \code{sdlog} (i.e., function \code{rlnorm} is called).

  For \code{dname = 'cauchy'}, allowed arguments are \code{location} or 
  \code{scale} (i.e., function \code{rcauchy} is called).

  Note that if a weight matrix and technical information are both provided 
  for a given mode (or for error), the weight matrix is used while technical 
  information is ignored. See Examples below for an example of how to set up 
  \code{technical}.
}
\value{
  \item{X}{
    Simulated data array with dimensions specified by \code{arraydim} and, when 
    \code{model = 'parafac2'}, also by \code{pf2num}. When \code{model = 
    'parafac'}, \code{X} is an object of class 'array'. When 
    \code{model = 'parafac2'}, \code{X} is an object of class 'list'.
  }
  \item{y}{
    Simulated class labels provided as an object of class 'matrix', with number 
    of rows equal to the last element of \code{arraydim} and with number of 
    columns equal to 1.
  }
  \item{model}{
    Character value indicating whether Parafac or Parafac2 model was used to 
    simulate the data array.
  }
  \item{Amat}{
    Simulated A mode weights. When \code{model = 'parafac'}, output is a matrix
    with number of rows equal to the first element of \code{arraydim} and with 
    number of columns equal to the number of components \code{nfac}. When 
    \code{model = 'parafac2'}, output is a list of matrices with number of rows
    for each matrix equal to those specified by \code{pf2num} and with number of
    columns equal to \code{nfac}. If \code{Amat} was supplied, returns original 
    \code{Amat} instead of a simulated \code{Amat}.
  }
  \item{Bmat}{
    Simulated B mode weights provided as a matrix with number of rows equal to 
    the second element of \code{arraydim} and with number of columns equal to 
    the number of components \code{nfac}. If \code{Bmat} was supplied,
    returns original \code{Bmat} instead of a simulated \code{Bmat}.
  }
   \item{Cmat}{
    Simulated C mode weights provided as a matrix with number of rows equal to 
    the third element of \code{arraydim} and with number of columns equal to 
    the number of components \code{nfac}. If \code{Cmat} was supplied when
    \code{modes = 4}, returns original \code{Cmat} instead of a simulated 
    \code{Cmat}.
  }
   \item{Dmat}{
    Simulated D mode weights provided when \code{modes = 4}. Output is a matrix 
    with number of rows equal to the fourth element of \code{arraydim} and with 
    number of columns equal to the number of components \code{nfac}.
  }
  \item{Gmat}{
    Simulated G weights provided when \code{model = 'parafac2'}. Provided as 
    a matrix with number of rows and columns equal to \code{nfac}. If 
    \code{Gmat} was supplied, returns original \code{Gmat} instead of a 
    simulated \code{Gmat}.
  }
  \item{Emat}{
    Error array or list containing noise added to corresponding elements of 
    simulated data array. Output has dimensions specified by 
    \code{arraydim} and, when \code{model = 'parafac2'}, also by \code{pf2num}. 
    When \code{model = 'parafac'}, \code{Emat} is an object of class 'array'. 
    When \code{model = 'parafac2'}, \code{Emat} is an object of class 'list'. 
  }
}
\note{
This simulation implementation contains at least two limitations. First, there 
is currently no argument to control the proportions of generated class labels. 
Second, the covariance matrix parameterizing the multivariate normal 
distribution generating classification mode weights is restricted to a 
correlation matrix. Future updates are planned to address these limitations.

In addition, the simulation could be expanded in at least two ways. First, 
the Monte Carlo simulation, as a brute-force strategy, is simple but 
not optimal and could be replaced with a more efficient approach. Note that 
the correlations between simulated class labels and classification mode 
weights, while ideally close to \code{corresp}, represent a best-case 
scenario---given values for \code{onreps}, \code{nreps}, and \code{corrpred}. 
Second, class labels are currently connected to the data array through only one 
mode but could be simulated such that they are connected through two or more 
modes. Future updates are planned to implement these enhancements.
}
\author{
 Matthew Snodgress <mattgress@protonmail.ch>
}
\references{
See help file for function \code{cpfa} for a list of references.
}
\examples{
########## Parafac2 example with 4-way array and multiclass response ##########
\dontrun{
# set seed for reproducibility
set.seed(5)

# define list of arguments specifying distributions for A and G weights
techlist <- list(distA = list(dname = "poisson", 
                              lambda = 3),                 # for A weights
                 distG = list(dname = "gamma", shape = 2, 
                              scale = 4))                  # for G weights

# define target correlation matrix for columns of D mode weights matrix
cormat <- matrix(c(1, .35, .35, .35, 1, .35, .35, .35, 1), nrow = 3, ncol = 3)

# simulate a four-way ragged array connected to a response
data <- simcpfa(arraydim = c(10, 11, 12, 100), model = "parafac2", nfac = 3, 
                nclass = 3, nreps = 1e2, onreps = 10, corresp = rep(.75, 3), 
                meanpred = rep(2, 3), modes = 4, corrpred = cormat,
                technical = techlist)
                
# examine correlations among columns of classification mode matrix Dmat
cor(data$Dmat)

# examine correlations between columns of classification mode matrix Dmat and
# simulated class labels
cor(data$Dmat, data$y)
}
}