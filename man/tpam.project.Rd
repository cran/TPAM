\name{tpam.project}
\alias{tpam.project}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Project time-course gene expression to weighted gene expression
}
\description{
Project time-course gene expression to weighted gene expression
}
\usage{
tpam.project(data, data.test)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
List of training data, of form described in tpam.train documentation
}
  \item{data.test}{
List of test data, of form described in tpam.test documentation
}
}

\value{
%%  ~Describe the value returned
list(data.train = wdata.train, data.test = wdata.test)
\item{data.train }{Projection of training data}
 \item{data.test }{Projection of test data}

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

proj.obj = tpam.project(data, data.test)
}
