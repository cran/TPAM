\name{tpam.predict}
\alias{tpam.predict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function giving prediction information.
}
\description{
Prediction information is provided.
}
\usage{
tpam.predict(fit, newdata, threshold, type = c("class", "posterior", "centroid", "nonzero"), prior = fit$prior, threshold.scale = fit$threshold.scale)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
The result of a call to tpam.train
}
  \item{newdata}{
Weighted test data, the result of a call to tpam.predict
}
  \item{threshold}{
The desired threshold value
}
  \item{type}{
Type of prediction desired: class predictions, posterior probabilities, (unshrunken) class centroids, vector of genes surviving the threshold.
}
  \item{prior}{
Prior probabilities for each class. Default is that specified in "fit".
}
  \item{threshold.scale}{
Additional scaling factors to be applied to the thresholds. Vector of length equal to the number of classes. Default is that specified in "fit".
}
}

\value{
Give a cross-tabulation of true versus predicted classes for the fit returned by tpam.train or tpam.cv, at the specified threshold.
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

predict.obj = tpam.predict(obj$fit, obj$proj.obj$data.test, threshold=0.8)
}
