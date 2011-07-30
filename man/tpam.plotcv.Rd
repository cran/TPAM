\name{tpam.plotcv}
\alias{tpam.plotcv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to plot the cross-validated error curves.
}
\description{
A function to plot the cross-validated error curves.
}
\usage{
tpam.plotcv(fit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
The result of a call to tpam.cv
}
}
\details{
plots the cross-validated misclassification error curves.
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
tpam.plotcv(cv.obj)
}
