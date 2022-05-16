## ---- eval = FALSE------------------------------------------------------------
#  install.packages("cpfa", repos = "https://cran.us.r-project.org")

## -----------------------------------------------------------------------------
library(cpfa)

## -----------------------------------------------------------------------------

# create random data for three-way tensor with Parafac structure and response
set.seed(123)
mydim <- c(32, 25, 240)
nfac <- nf <- 3
rho.c.c <- .35
rho.c.y <- .75
R <- matrix(c(  1, rho.c.c, rho.c.c, rho.c.y,
                rho.c.c, 1, rho.c.c, rho.c.y, 
                rho.c.c, rho.c.c, 1, rho.c.y,
                rho.c.y, rho.c.y, rho.c.y, 1), nrow = nfac+1, ncol = nfac+1)
Nsubj <- mydim[3]
Nvar <- nfac + 1
C.col <- runif(Nsubj*nfac)
y.col <- rbinom(Nsubj, 1, 0.5)
values <- c(C.col, y.col)
Y <- matrix(values, nrow = Nsubj, ncol = Nvar)
Y <- Y - matrix(1, Nsubj, 1) %*% apply(Y, 2, mean)
S <- t(Y) %*% Y
M <- t(chol(S))
Minv <- solve(M)
L <- t(chol(R))
Xq <- Y %*% t(Minv) %*% t(L)
Cmat <- Xq[, 1:3]*2
Amat <- matrix(rnorm(mydim[1]*nfac), nrow = mydim[1], ncol = nfac)
Bmat <- matrix(runif(mydim[2]*nfac), nrow = mydim[2], ncol = nfac)
Xmat <- tcrossprod(Amat, krprod(Cmat, Bmat))
Xmat <- array(Xmat, dim = mydim)
Emat <- array(rnorm(prod(mydim)), dim = mydim)
Emat <- nscale(Emat, 0, ssnew = sumsq(Xmat))  
X <- Xmat + Emat
y <- factor(as.numeric(Xq[, 4] > 0))


## -----------------------------------------------------------------------------

# examine data object X
dim(X)
class(X)

# examine data object y
class(y)
length(y)
table(y)


## -----------------------------------------------------------------------------

# split data into training and testing sets
split <- 0.8
nlev <- length(y)
index <- round(split*nlev)
X.train <- X[, , 1:index]
X.test <- X[, , (index+1):nlev]
y.train <- y[1:index]
y.test <- as.numeric(y[(index + 1):nlev]) - 1


## -----------------------------------------------------------------------------

# initialize inputs for cpfa
nfac <- 1:3
nfolds <- 5
foldid <- sample(rep(1:nfolds, length.out = length(y.train)))
method <- c("PLR", "SVM", "RF", "NN")
ntree <- c(100, 300, 500)
family <- "binomial"
parallel <- FALSE

# initialize inputs passed directly to 'parafac' within package 'multiway'
const <- c("orthog", "smooth", "uncons")
ctol <- 1e-02
nstart <- 5


## -----------------------------------------------------------------------------

tune.object <- cpfa(x = X.train, y = y.train, nfac = nfac, nfolds = nfolds, 
                    foldid = foldid, method = method, ntree = ntree, 
                    family = family, parallel = parallel, const = const, 
                    ctol = ctol, nstart = nstart)


## -----------------------------------------------------------------------------

tune.object


## -----------------------------------------------------------------------------

yhat <- predict(tune.object, newdata = X.test, type = "response") 


## -----------------------------------------------------------------------------

perform.output <- cpm(yhat$fac.3rf, y.test)
perform.output$cm
perform.output$class.eval


## -----------------------------------------------------------------------------

apply(yhat, 2, function(x){cpm(x, y.test)})


## -----------------------------------------------------------------------------

# create random data for four-way tensor with Parafac structure and response
set.seed(123)
nfac <- nf <- 3
mydim <- c(32, 25, 240, 20)
aseq <- seq(-3, 3, length.out = mydim[1])
Amat <- cbind(dnorm(aseq), dchisq(aseq+3.1, df=3),
              dt(aseq-2, df=4), dgamma(aseq+3.1, shape=3, rate=1))[,1:3]
Bmat <- svd(matrix(runif(mydim[2]*nf), nrow = mydim[2], ncol = nf), nv = 0)$u
rho.c.c <- .725
rho.c.y <- .9
R <- matrix(c(  1, rho.c.c, rho.c.c, rho.c.y,
                rho.c.c, 1, rho.c.c, rho.c.y, 
                rho.c.c, rho.c.c, 1, rho.c.y,
                rho.c.y, rho.c.y, rho.c.y, 1), nrow = nfac+1, ncol = nfac+1)
Nsubj <- mydim[3]
Nvar <- nfac + 1
C.col <- runif(Nsubj*nfac)
y.col <- rbinom(Nsubj, 1, 0.5)
values <- c(C.col, y.col)
Y <- matrix(values, nrow = Nsubj, ncol = Nvar)
Y <- Y - matrix(1, Nsubj, 1) %*% apply(Y, 2, mean)
S <- t(Y) %*% Y
M <- t(chol(S))
Minv <- solve(M)
L <- t(chol(R))
Xq <- Y %*% t(Minv) %*% t(L)
Cmat <- Xq[, 1:3]*2
y <- Xq[,4]
for(i in 1:length(y)){
  if(abs(y[i]) > 0.09){
    y[i] <- 2
  } 
  if((y[i] < 0.09) & (y[i] > 0)){
    y[i] <- 1
  }
  if((y[i] > -0.09) & (y[i] < 0)){
    y[i] <- 0
  }
}
Dmat <- matrix(runif(mydim[4]*nf), nrow = mydim[4], ncol = nf)
Xmat <- tcrossprod(Amat, krprod(Dmat, krprod(Cmat, Bmat)))
Xmat <- array(Xmat, dim = mydim)
Emat <- array(rnorm(prod(mydim)), dim = mydim)
Emat <- nscale(Emat, 0, ssnew = sumsq(Xmat))   # SNR = 1
X <- Xmat + Emat
X2 <- X
y2 <- factor(y)


## -----------------------------------------------------------------------------

dim(X2)
class(X2)


## -----------------------------------------------------------------------------

class(y2)
length(y2)
table(y2)


## -----------------------------------------------------------------------------

# split data into training and testing sets
split <- 0.8
nlev <- length(y2)
index <- round(split*nlev)
X.train <- X2[, , 1:index, ]
X.test <- X2[, , (index+1):nlev, ]
y.train <- y2[1:index]
y.test <- as.numeric(y2[(index + 1):nlev]) - 1


## -----------------------------------------------------------------------------

# initialize inputs for cpfa
nfac <- 3
nfolds <- 3
foldid <- sample(rep(1:nfolds, length.out = length(y.train)))
method <- c("PLR", "SVM", "RF", "NN")
ntree <- c(100, 300, 500)
family <- "multinomial"
parallel <- FALSE

# initialize inputs passed directly to 'parafac' within package 'multiway'
const <- c("orthog", "smooth", "uncons", "uncons")
ctol <- 1e-02
nstart <- 5


## -----------------------------------------------------------------------------

cmode <- 3


## -----------------------------------------------------------------------------

tune.object2 <- cpfa(x = X.train, y = y.train, nfac = nfac, nfolds = nfolds, 
                    foldid = foldid, method = method, ntree = ntree, 
                    family = family, parallel = parallel, const = const, 
                    ctol = ctol, nstart = nstart, cmode = cmode)


## -----------------------------------------------------------------------------

yhat <- predict(tune.object2, newdata = X.test, type = "response") 
apply(yhat, 2, function(x){cpm(x, y.test)})


