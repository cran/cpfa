## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  cache     = FALSE,
  warning   = FALSE,
  message   = FALSE,
  seed      = 500
)

## ----message = FALSE----------------------------------------------------------
library(cpfa)

## ----message = FALSE, warning = FALSE-----------------------------------------
# set seed for reproducibility
set.seed(500)

# specify correlation
cp <- 0.1

# define target correlation matrix for columns of fourth mode weight matrix
corrpred <- matrix(c(1, cp, cp, cp, 1, cp, cp, cp, 1), nrow = 3, ncol = 3)

# define correlations between fourth mode weight matrix and response vector
corresp <- rep(.85, 3)

# specify number of rows in the three-way array for each level of fourth mode
pf2num <- rep(c(7, 8, 9), length.out = 100)

# simulate a four-way ragged array connected to a response
data <- simcpfa(arraydim = c(10, 11, 12, 100), model = "parafac2", nfac = 3, 
                nclass = 3, nreps = 10, onreps = 10, corresp = corresp,
                pf2num = pf2num, modes = 4, corrpred = corrpred, 
                meanpred = c(10, 20, 30))

# define simulated array 'X' and response vector 'y' from the output
X <- data$X
y <- data$y

## -----------------------------------------------------------------------------
# examine data object X
class(X)
length(X)
dim(X[[1]])
dim(X[[2]])

# examine data object y
class(y)
length(y)
table(y)

## -----------------------------------------------------------------------------
# examine correlations between columns of fourth mode weights 'Dmat' and 
# simulated response vector 'y'
cor(data$Dmat, data$y)

## -----------------------------------------------------------------------------
# set seed
set.seed(500)

# initialize alpha and store within a list called 'parameters'
alpha <- seq(0, 1, length.out = 11)
parameters <- list(alpha = alpha)

# initialize inputs
method <- "PLR"
model <- "parafac2"
nfolds <- 3
nstart <- 3
nfac <- c(2, 3)
family <- "multinomial"
nrep <- 3
ratio <- 0.9
plot.out <- TRUE
const <- c("uncons", "uncons", "uncons", "nonneg")
foldid <- rep(1:nfolds, length.out = ratio * length(y))

# implement train-test splits with inner k-fold CV to optimize classification
output <- cpfa(x = X, y = as.factor(y), model = model, nfac = nfac, 
               nrep = nrep, ratio = ratio, nfolds = nfolds, method = method, 
               family = family, parameters = parameters, plot.out = plot.out, 
               parallel = FALSE, const = const, foldid = foldid, 
               nstart = nstart, verbose = FALSE)

## -----------------------------------------------------------------------------
# examine classification performance measures - median across train-test splits
output$descriptive$median[, 1:2]

## -----------------------------------------------------------------------------
# examine optimal tuning parameters averaged across train-test splits
output$mean.opt.tune

## -----------------------------------------------------------------------------
# set seed
set.seed(500)

# plot heat maps of component weights for optimal model
results <- plotcpfa(output, nstart = 3, ctol = 1e-1, verbose = FALSE)

## ----message = FALSE, warning = FALSE-----------------------------------------
# set seed for reproducibility
set.seed(400)

# specify correlation
cp <- 0.1

# define target correlation matrix for columns of third mode weight matrix
corrpred <- matrix(c(1, cp, cp, 1), nrow = 2, ncol = 2)

# define correlations between third mode weight matrix and response vector
corresp <- rep(.9, 2)

# simulate a three-way array connected to a binary response
data <- simcpfa(arraydim = c(10, 11, 100), model = "parafac", nfac = 2, 
                nclass = 2, nreps = 10, onreps = 10, corresp = corresp,
                modes = 3, corrpred = corrpred, meanpred = c(10, 20))

# define simulated array 'X' and response vector 'y' from the output
X <- data$X
y <- data$y

## -----------------------------------------------------------------------------
# examine data object X
class(X)
dim(X)

# examine data object y
class(y)
length(y)
table(y)

## -----------------------------------------------------------------------------
# examine correlations between columns of third mode weights 'Cmat' and 
# simulated response vector 'y'
cor(data$Cmat, data$y)

## -----------------------------------------------------------------------------
# set seed
set.seed(300)

# initialize tuning parameters and store within a list called 'parameters'
alpha <- seq(0, 1, length.out = 3)
ntree <- c(200, 400)
nodesize <- c(2, 4)
parameters <- list(alpha = alpha, ntree = ntree, nodesize = nodesize)

# initialize inputs
method <- c("PLR", "RF")
model <- "parafac"
nfolds <- 3
nstart <- 3
nfac <- c(2, 3)
family <- "binomial"
nrep <- 3
ratio <- 0.9
plot.out <- TRUE
const <- c("uncons", "orthog", "uncons")

# implement train-test splits with inner k-fold CV to optimize classification
output <- cpfa(x = X, y = as.factor(y), model = model, nfac = nfac, 
               nrep = nrep, ratio = ratio, nfolds = nfolds, method = method, 
               family = family, parameters = parameters, plot.out = plot.out, 
               parallel = FALSE, const = const, nstart = nstart, 
               verbose = FALSE)

## -----------------------------------------------------------------------------
# examine classification performance measures - median across train-test splits
output$descriptive$mean[, 1:2]

## -----------------------------------------------------------------------------
# examine optimal tuning parameters averaged across train-test splits
output$mean.opt.tune

## ----eval = FALSE-------------------------------------------------------------
# # set seed
# set.seed(400)
# 
# # plot heat maps of component weights for optimal model
# results <- plotcpfa(output, nstart = 3, ctol = 1e-1, verbose = FALSE)

