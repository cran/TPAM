\name{tpam.cv}
\alias{tpam.cv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cross-validation
}
\description{
This function uses a form of cross-validation to estimate the optimal feature threshold.
}
\usage{
tpam.cv(fit, data, nfold = NULL, folds = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
The result of a call to tpam.train
}
  \item{data}{
Data object of form described in tpam.train documentation
}
  \item{nfold}{
Number of cross-validation folds
}
  \item{folds}{
Lists of indices of cross-validation folds (optional)
}
}

\value{
list(threshold=threshold, error=error, loglik=loglik,size=size, yhat=yhat,y=y,prob=prob,folds=folds, cv.objects=cv.objects, call = this.call)
%%  If it is a LIST, use
\item{thresholds }{Vector of thresholds considered}
\item{error }{The number of cross-validation errors for each threshold value}
\item{loglik}{The cross-validated multinomial log-likelihood value for each threshold value}
\item{size}{A vector of the number of genes that survived the thresholding, for each threshold value tried.}
\item{yhat}{A matrix of size n by nthreshold, containing the cross-validated class predictions for each threshold value, in each column.}
\item{y}{A vector of true class labels.}
\item{prob}{A matrix of size n by nthreshold, containing the cross-validated probabilities of class predictions for each threshold value, in each column.}
\item{folds}{The cross-validation folds used}
\item{cv.objects}{Train objects, from each of the CV folds}
\item{call}{The calling sequence used}
%% ...
}
\author{
Yuping Zhang
}

\examples{
x = list()
for(i in 1:2){
	set.seed(i+123)
	x[[i]] = matrix(rnorm(500*100), ncol=100)
}
y = factor(sample(c(1:2), size=100, replace=TRUE))

data = list(x = x, y=y, genenames = as.character(paste("gene", c(1:500), sep="")), geneid = 	as.character(paste("gene", c(1:500), sep="")))

x = list()
for(i in 1:2){
	set.seed(i+133)
	x[[i]] = matrix(rnorm(500*100), ncol=100)
}
y = factor(sample(c(1:2), size=100, replace=TRUE))

data.test = list(x = x, y=y, genenames = as.character(paste("gene", c(1:500), sep="")), geneid = 	as.character(paste("gene", c(1:500), sep="")))

obj = tpam.train(data, data.test)

cv.obj = tpam.cv(obj$fit, data=data, nfold=2)
}

