#################################################################
### Lec2.r  Script file for lecture notes pp. 100 - 162
#################################################################
##  Illustrate matrix computations using the Toluca dataset
##   Several methods are illustrated for finding Yhat and residuals vectors.
toluca <- read.table("Datasets/CH01TA01.txt",header=FALSE)
names(toluca) <- c("size","hours")
# It is convenient to attach the dataset, but this is not good
#  programming.  Be sure to detach at the end.
attach(toluca)
dim(toluca) # 25 x 2 dataframe
n <- nrow(toluca)
Y <- toluca$hours
Y[1:3]
X <- matrix(c(rep(1,n),toluca$size), nrow=n)
X[1:3,]
XpX <- t(X) %*% X; XpX
XpX.inv <- solve(XpX); XpX.inv
XpX %*% XpX.inv
XpY <- t(X) %*% Y; XpY
betahat <- XpX.inv %*% XpY; betahat
# Get vectors of fits and residuals directly from definition:
Yhat <- X %*% betahat; Yhat[1:3]
resids <- Y - Yhat; resids[1:3]
# Matrix formula for SSE:
SSE <- as.numeric(t(Y) %*% Y - t(betahat) %*% t(X) %*% Y); SSE
#  SSE directly as squared length of the residuals vector.
SSE <- resids %*% resids; SSE
#
# The Hat matrix; slick way to get Yhat, residuals vectors, and SSE
H <- X %*% XpX.inv %*% t(X)
I.mat <- diag(nrow=n)
Yhat <- H %*% Y
resids <- (I.mat - H) %*% Y
SSE <- as.numeric(t(Y) %*% (I.mat - H) %*% Y) # quadratic form for SSE
dfE <- n - qr(X)$rank; dfE # df for Error is n - rank(X)
MSE <- SSE/dfE  # mean sq for resids: divide SSE by df(Error) = n-2 in SLR
s2betahat <- MSE * XpX.inv # estimated var-cov matrix of betahat
## Estimate mean response at Xh= (1, 65)
##   Here we get the estimate and its standard error (SE of the fit)
##   These are the key ingredients for CIs for mean response at a given x
Xh <- c(1,65)
Yhhat <- t(Xh) %*% betahat; Yhhat # predicted value
se.fit <- sqrt(t(Xh) %*% s2betahat %*% Xh) # st. error of the fit
se.fit
#
# T test statistic for test of H0: beta1=0
tobs <- betahat[2]/sqrt(s2betahat[2,2]); tobs
# Clean up.
rm(Y,X, SSE, XpX,XpX.inv,XpY,betahat,Yhat,toluca)
detach(toluca)
## Ex. Multiple regr. with two predictors
##   Dwaine Studios.
##   You might start RStudio over to make sure not to mix up
##      the Toluca and Dwaine datasets.
dwaine <- read.table("Datasets/CH07TA05.txt",
                     header=FALSE)
names(dwaine) <- c("Y","x1","x2")
attach(dwaine)
dim(dwaine)
n <- nrow(dwaine)
Y[1:3]
X <- matrix(c(rep(1,n),x1,x2), nrow=n)
X[1:3,]
XpX <- t(X) %*% X; XpX
XpX.inv <- solve(XpX); round(XpX.inv, digits=3)
round(XpX %*% XpX.inv, digits=3)
XpY <- t(X) %*% Y; XpY
betahat <- XpX.inv %*% XpY; betahat
## Get vectors of fits and residuals directly from definition:
Yhat <- X %*% betahat; Yhat[1:3]
resids <- Y - Yhat; resids[1:3]
## Matrix formula for SSE
SSE <- as.numeric(t(Y) %*% Y - t(betahat) %*% t(X) %*% Y); SSE
##  SSE directly as squared length of the residuals vector.
SSE <- t(resids) %*% resids; SSE
##
## The Hat matrix; slick way to get Yhat, residual vectors, and SSE
H <- X %*% XpX.inv %*% t(X)
I.mat <- diag(nrow=n)
Yhat <- H %*% Y
resids <- (I.mat - H) %*% Y
SSE <- as.numeric(t(Y) %*% (I.mat - H) %*% Y)# quadratic form for SSE
dfE <- n - qr(X)$rank; dfE # df for Error is n - rank(X)
MSE <- SSE/dfE  # mean sq for resids: divide SSE by dfE = n-3
s2betahat <- MSE * XpX.inv # estimated var-cov matrix of betahat
round(s2betahat,digits=3)
## Estimate mean response at Xh= (1, 65.4, 17.6)
Xh <- c(1, 65.4, 17.6)
Yhhat <- t(Xh) %*% betahat; Yhhat # predicted value
se.fit <- sqrt(t(Xh) %*% s2betahat %*% Xh) # st. error of the fit
se.fit
##

## T test statistics for tests of H0: beta1=0, controlling for x2,
##                            and H0: beta2=0, controlling for x1
##   (Test about the intercept is the first one listed.)
##
tobs <- betahat/sqrt(diag(s2betahat)); tobs
# Scatterplot matrix
pairs(dwaine[,c(1,2,3)])
## Clean up
detach(dwaine)

