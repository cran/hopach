\name{correlationordering}

\alias{correlationordering}
\alias{improveordering}

\title{function to compute empirical correlation between distance in a list and distance by a metric}

\description{
Given a matrix of pair wise distances based on a choice of distance metric, \code{correlationordering} computes the empirical correlation (over all pairs of elements) between the distance apart in the rows/columns of the matrix and the distance according to the metric. Correlation ordering will be high if elements close to each other in the matrix have small pair wise distances. If the rows/columns of the distance matrix are ordered according to a clustering of the elements, then correlation ordering should be large compared to a matrix with randomly ordered rows/columns.
}

\usage{
correlationordering(dist)

improveordering(dist,echo=FALSE)
}

\arguments{
  \item{dist}{matrix of all pair wise distances between a set of 'p' elements, 
	as produced, for example, by the \code{distancematrix} function. 
	The value in row 'j' and column 'i' is the distance between element 'i'
 	and element 'j'. The matrix must be symmetric. The ordering of the rows/
	columns is compared to the values in the matrix.}

  \item{echo}{indicator of whether the value of correlation ordering before and after rearranging the ordering should be printed.}
}

\details{
Correlation ordering is defined as the empirical correlation between distance in a list and distance according to some other metric. The value in row 'i' and column 'j' of \code{dist} is compared to 'j-i'. The function \code{correlationordering} computes the correlation ordering for a matrix \code{dist}, whereas the function \code{improveordering} swaps the ordering of elements in \code{dist} until doing so no longer improves correlation ordering. The algorithm for \code{improveordering} is not optimized, so that the function can be quite slow for more than 50 elements. These functions are used by the \code{hopach} clustering function to sensibly order the clusters in the first level of the hierarchical tree, and can also be used to order elements within clusters when the number of elements is not too large.
}

\value{
	For \code{correlationordering}, a number between -1 and 1, as returned by the \code{cor} function, equal to the correlation ordering for the matrix \code{dist}.

	For \code{improveordering}, a vector of length 'p' containing the row indices for the new ordering of the rows/columns of \code{dist}, so that dist[improveordering(dist)] now has higher correlation ordering.
}

\references{

van der Laan, M.J. and Pollard, K.S. A new algorithm for hybrid hierarchical clustering with visualization and the bootstrap. Journal of Statistical Planning and Inference, 2003, 117, pp. 275-303.

\url{http://www.stat.berkeley.edu/~laan/Research/Research_subpages/Papers/hopach.pdf}

}

\author{Katherine S. Pollard <kpollard@wald.ucdavis.edu> and Mark J. van der Laan <laan@stat.berkeley.edu>}

\section{Warning }{The function \code{improveordering} can be very slow for more than about 50 elements. The method employed is a greedy, step-wise algorithm, in which sequentially swaps all pairs of elements and accepts any swap that improves correlation ordering. }

\seealso{\code{\link{distancematrix}}, \code{\link{hopach}}}

\examples{
mydata<-matrix(rnorm(50),nrow=10)
mydist<-distancematrix(mydata,d="euclid")
image(mydist)
correlationordering(mydist)
neword<-improveordering(mydist,echo=TRUE)
correlationordering(mydist[neword,neword])
image(mydist[neword,neword])
}

\keyword{multivariate}
\keyword{cluster}
