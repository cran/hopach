\name{hdist-class}
\docType{class}
\alias{hdist-class}
\alias{hdist}
\alias{is.hdist}
\alias{as.hdist}
\alias{as.hdist,matrix-method}
\alias{as.matrix,hdist-method}
\alias{as.vector,hdist,missing-method}
\alias{[,hdist-method}
\alias{coerce,matrix,hdist-method}
\alias{coerce,hdist,matrix-method}
\alias{dim,hdist-method}
\alias{labels,hdist-method}
\alias{length,hdist-method}
\alias{show,hdist-method}

\title{Class "hdist" - S4 class to hold distance matrices.}
\description{	
	Class hdist was created to take advantage of the structure innate to symmetric 
	matrices.  It stores only the lower triangle of the matrix, thus reducing
	the size (and memory usage) from n x n to [n x (n - 1)] / 2.  	

	Like a matrix, a hdist object is subsettable; thus, hdist[i,j] will 
	return the value at row 'i' column 'j'. Most valid indices for a matrix are also
	valid for a hdist object. (See examples below)
}
\section{Slots}{
	 \describe{
    \item{\code{Data}:}{
	Object of class \code{"numeric"} a vector containing the stacked columns of
	the lower triangle of a symmetric matrix -- often the symmetric matrix is a
	distance matrix.}
    \item{\code{Size}:}{Object of class \code{"numeric"} the dimension of the 
	symmetric matrix, from which \code{Data} was constructed.}
    \item{\code{Labels}:}{Object of class \code{"numeric"} a list of values of 
	length \code{Size} to allow for pretty printing.}
    \item{\code{Call}:}{Object of class \code{"character"} a character string 
	specifying the method used to create the distance matrix
	from which \code{Data} was constructed.}
  }
}
\section{Methods}{
  \describe{
    \item{hdist}{\code{signature{Data = "numeric", Size = "numeric", Labels = "numeric", Call = "character"}}: Create a new hdist object. }
    \item{as.hdist}{\code{signature{from = "matrix"}}: Converts a matrix to a hdist object. }
    \item{as.matrix}{\code{signature(x = "hdist")}: Converts a hdist object to a matrix. }
    \item{as.vector}{\code{signature(x = "hdist", mode = "missing")}: Returns the hdist object as a vector. }
    \item{[}{\code{signature(x = "hdist")}: Subsetting function for hdist objects. See examples and warning. }
    \item{coerce}{\code{signature(from = "matrix", to = "hdist")}: Converts a matrix to a hdist object. }
    \item{coerce}{\code{signature(from = "hdist", to = "matrix")}: Converts a hdist object to a matrix. }
    \item{dim}{\code{signature(x = "hdist")}: Returns the dimension of the hdist object if expanded to a square matrix.}
    \item{labels}{\code{signature(object = "hdist")}: Returns the labels used for printing. }
    \item{length}{\code{signature(x = "hdist")}: Returns the number of rows in hdist object. }
    \item{show}{\code{signature(object = "hdist")}: Prints the hdist object. }
	 }
}
\references{

van der Laan, M.J. and Pollard, K.S. A new algorithm for hybrid hierarchical clustering with visualization and the bootstrap. Journal of Statistical Planning and Inference, 2003, 117, pp. 275-303.

\url{http://www.stat.berkeley.edu/~laan/Research/Research_subpages/Papers/hopach.pdf}

\url{http://www.bepress.com/ucbbiostat/paper107/}

\url{http://www.stat.berkeley.edu/~laan/Research/Research_subpages/Papers/jsmpaper.pdf}

Kaufman, L. and Rousseeuw, P.J. (1990). Finding Groups in Data: An Introduction to Cluster Analysis. Wiley, New York.

}

\author{Katherine S. Pollard <kpollard@ucdavis.edu> and Gregory D. Wall <gwall@wald.ucdavis.edu>}
\note{Thank you to Larry Tai for his assistance creating run-time comparisons.}

\section{Warning }{
	A hdist object is NOT closed under the subsetting operation. For instance, if
	a 100 x 100 symmetric matrix is stored as an hdist object, hdist[c(3,4,5),c(7,8,9)]
	will return a 3 x 3 matrix, since the subsetting will not result in a symmetric matrix.
	However, if index i = j, then subsetting a hdist object will result in a symmetric 
	matrix, and thus a hdist object will be returned. 
	(See examples below) 
}
\seealso{
	\code{\link{hopach}} 
}
\examples{
	showClass("hdist")

	library(hopach) 
	X <- matrix(rnorm(60,mean=10,sd=2),nrow=10,ncol=6,byrow=TRUE)
	dmat <- disscosangle(X) 
	dmat
	str(dmat) 
	
	# Examples where subsetting a hdist object returns a matrix...	
	dmat[c(3,4,5),c(5,6,7,8)]
	dmat[c(TRUE,FALSE),c(FALSE,TRUE)]
	dmat[c(4,5,6), ] 

	# Examples where subsetting a hdist object returns a hdist object...
	dmat[c(3,4,5,6,7),c(3,4,5,6,7)]
	dmat[c(TRUE,FALSE),c(TRUE,FALSE)]
	
	# Expand hdist object to a symmetric matrix...
	as.matrix(dmat)

}
\keyword{classes}
