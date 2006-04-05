\name{disscosangle}

\alias{disscosangle}
\alias{disseuclid}
\alias{disscor}
\alias{dissabscosangle}
\alias{dissabseuclid}
\alias{dissabscor}
\alias{vdisscosangle}
\alias{vdisseuclid}
\alias{vdisscor}
\alias{vdissabscosangle}
\alias{vdissabseuclid}
\alias{vdissabscor}

\title{Functions to compute pair-wise distances}

\description{
Given a matrix \code{X}, these functions compute the \code{nrow(X)} by \code{nrow{X}} matrix of pair-wise distances between all variables (rows) in \code{X}, across all observations (columns) of \code{X}. Each function uses a different distance metric, i.e. definition of what it means for two variables to be similar.
}

\usage{
disscosangle(X, na.rm = TRUE)

disseuclid(X, na.rm = TRUE)

disscor(X, na.rm = TRUE)

dissabscosangle(X, na.rm = TRUE)

dissabseuclid(X, na.rm = TRUE)

dissabscor(X, na.rm = TRUE)

vdisscosangle(X, y, na.rm = TRUE)

vdisseuclid(X, y, na.rm = TRUE)

vdisscor(X, y, na.rm = TRUE)

vdissabscosangle(X, y, na.rm = TRUE)

vdissabseuclid(X, y, na.rm = TRUE)

vdissabscor(X, y, na.rm = TRUE)
}

\arguments{
  \item{X}{A numeric data matrix. Each column corresponds to an observation, and each row corresponds to a variable. In the gene expression context, observations are arrays and variables are genes. All values must be numeric. Missing values are ignored.}
  \item{na.rm}{Indicator of whether to remove missing values (i.e. only compute distance over non-missing observations).}
  \item{y}{A numeric data vector of length \code{ncol(X)}.}
}

\details{
Different choices of distance metric are discussed in the references. Briefly, Euclidean distance (\code{disseuclid}) defines two variables to be close if they are similar in magnitude across observations. Correlation distance (\code{disscor}), in contrast, defines similarity to mean having the same pattern, but not necessarily the same magnitude. Cosine-angle (\code{disscosangle}) distance is a correlation distance that also accounts for magnitude. Cosine-angle distance is also known as uncentered correlation distance. The distance metrics with 'abs' in their names are absolute versions of each metric; the absolute value is applied to the data before computing the distance.
}

\value{
A numeric \code{nrow(X)} by \code{nrow{X}} matrix of pair-wise distances between all variables (rows) in \code{X}. For the vector versions (e.g. \code{vdisscosangle}), a numeric vector of \code{nrow(X)} pair-wise distances between each variable (row) in \code{X} and the vector \code{y}.
}

\references{
van der Laan, M.J. and Pollard, K.S. A new algorithm for hybrid hierarchical clustering with visualization and the bootstrap. Journal of Statistical Planning and Inference, 2003, 117, pp. 275-303.

\url{http://www.stat.berkeley.edu/~laan/Research/Research_subpages/Papers/hopach.pdf}

\url{http://www.bepress.com/ucbbiostat/paper107/}

\url{http://www.stat.berkeley.edu/~laan/Research/Research_subpages/Papers/jsmpaper.pdf}
}

\author{Katherine S. Pollard <kpollard@wald.ucdavis.edu> and Mark J. van der Laan <laan@stat.berkeley.edu>}

\seealso{\code{\link{distancematrix}}}

\examples{
data<-matrix(rnorm(50),nr=5)
disscosangle(data)
}

\keyword{multivariate}
\keyword{internal}
