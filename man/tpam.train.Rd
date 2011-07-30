\name{tpam.train}
\alias{tpam.train}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Classification using time-course gene expression 
}
\description{
Does prediction of a categorical outcome, using the time-course gene expression. 
}
\usage{
tpam.train(data, data.test)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Data object with components x- a list of p by n matrix of features,
one observation per column, one matrix per time point; y- n-vector of outcome measurements;
genenames - a vector of gene names; geneid - a vector of gene identifiers.
}
  \item{data.test}{
Data object with components x- a list of p by n matrix of features,
one observation per column, one matrix per time point; y- n-vector of outcome measurements;
genenames - a vector of gene names; geneid - a vector of gene identifiers.
}
}

\value{
	\item{proj.obj}{projection of training data and test data}
	\item{fit.obj}{fitted object using training data}
}
\author{
Yuping Zhang
}

\examples{
# generate some data

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

}
