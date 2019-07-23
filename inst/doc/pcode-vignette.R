## ------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----warning=FALSE,message=FALSE-----------------------------------------
#load dependencies
library(pCODE)
library(deSolve)
library(fda)
library(MASS)
library(pracma)
library(Hmisc)
#set seed for reproducibility
set.seed(123)
ode.model <- function(t,state,parameters){
            with(as.list(c(state,parameters)),{
                 dX <- theta*X*(1-X/10)
                return(list(dX))})}

## ------------------------------------------------------------------------
#define model parameters
model.par   <- c(theta = c(0.1))
#define state initial value
state       <- c(X     = 0.1)

## ------------------------------------------------------------------------
times  <- seq(0,100,length.out=101)
mod    <- ode(y=state,times=times,func=ode.model,parms = model.par)
nobs   <- length(times)
scale  <- 0.5
noise  <- scale * rnorm(n = nobs, mean = 0 , sd = 1)
observ <- mod[,2] + noise

## ----fig.align='center',fig.width=6,fig.height=3.5-----------------------
#plot simulated data against generating model
plot(mod,ylab=names(state))      #curve
points(times, observ,pch='*',col='blue')    #observation

## ------------------------------------------------------------------------
#Generate basis object for interpolation and as argument of pcode
#21 konts equally spaced within [0,100]
knots <- seq(0,100,length.out=21)
#order of basis functions
norder <- 4
#number of basis funtions
nbasis <- length(knots) + norder - 2
#creating Bspline basis
basis  <- create.bspline.basis(c(0,100),nbasis,norder,breaks = knots)

## ----results="hide"------------------------------------------------------
#parameter estimation
pcode.result <- pcode(data = observ, time = times, ode.model = ode.model,
                      par.initial = 0.3, par.names = 'theta',state.names = 'X',
                      basis.list = basis, lambda = 1e2)

## ------------------------------------------------------------------------
pcode.result$structural.par
pcode.result$nuisance.par

## ------------------------------------------------------------------------
deltavar(data = observ, time = times, ode.model = ode.model,
                        par.initial = 0.3, par.names = 'theta',state.names = 'X',
                        basis.list = basis, lambda = 1e2,
                        stepsize = 1e-5,y_stepsize = 1e-5)

## ----eval=FALSE----------------------------------------------------------
#  #Tune lambda based on k-fold cross-validation
#  cv.result <- tunelambda(data = observ, time = times, ode.model = ode.model,
#  	               	   par.initial = 0.3, par.names = 'theta',state.names = 'X',
#  	               	   basis.list = basis, lambda_grid = 10^(-3:10),cv_portion = .01,
#                         rep = 20, kfolds = 5)

## ----loaddata, echo=FALSE------------------------------------------------
load('simpleode.Rdata')

## ---- fig.width=7, fig.height=4,fig.show='hold',fig.align='center'-------
cv.means <- apply(cv.result$cv.score, 1, mean)
cv.sd    <- apply(cv.result$cv.score, 1, sd)/sqrt(20)
plot(-2:10, cv.means[2:14], type = "l", lwd = 2, col = gray(0.4),ylim= c(-0,50),ylab='',xlab='')
errbar(-2:10, cv.means[2:14], cv.means[2:14] + cv.sd[2:14], cv.means[2:14] - cv.sd[2:14], add = TRUE, col = "steelblue2", pch = 19,
       lwd = 0.5)

plot(1:10, cv.means[5:14], type = "l", lwd = 2, col = gray(0.4),ylim= c(-0,5),ylab='',xlab='')
errbar(1:10, cv.means[5:14], cv.means[5:14] + cv.sd[5:14], cv.means[5:14] - cv.sd[5:14], add = TRUE, col = "steelblue2", pch = 19,
       lwd = 0.5)

