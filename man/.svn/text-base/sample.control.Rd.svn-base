% File src/library/base/man/load.Rd
\name{sample.control}
\alias{sample.control}
\title{Sample while controlling for target}
\description{
  sample with sampled data  controlled with the target data
}
\usage{
sample.control(data, size, target,  prob=seq(0,1,.05), replace=F)
}
\arguments{
  \item{data}{ input data}
  \item{size}{  a non-negative integer giving the number of items to choose.}
  \item{replace}{ replace or not}
  \item{target} {target distribution to be matched}
}
\seealso{
  \code{\link{save}}.
}
\examples{
## save all data
save(list = ls(), file= "all.RData")

## restore the saved values to the current environment
load("all.RData")

## restore the saved values to the workspace
load("all.RData", .GlobalEnv)
}
\keyword{file}
