\name{distancematrix}

\alias{distancematrix}
\alias{distancevector}
\alias{dissmatrix}
\alias{dissvector}
\alias{vectmatrix}

\title{functions to compute pair wise distances between vectors}

\description{
The function \code{distancematrix} is applied to a matrix of data to
compute the pair wise distances between all rows of the matrix. The
function \code{distancevector} is applied to a matrix and a vector
to compute the pair wise distances between each row of the matrix and 
the vector. Both functions allow different choices of distance metric.
The functions \code{dissmatrix} and \code{dissvector} allow one to 
convert between a distance matrix and a vector of the upper triangle. 
The function \code{vectmatrix} is used internally.
}

\usage{
distancematrix(X, d, na.rm=TRUE)

distancevector(X, y, d, na.rm=TRUE)

dissmatrix(v)

dissvector(M)

vectmatrix(index, p)
}

\arguments{
  \item{X}{a numeric matrix. Missing values will be ignored if na.rm=TRUE.}
  \item{y}{a numeric vector, possibly a row of X. Missing values will be ignoredif na.rm=TRUE.}
  \item{na.rm}{an indicator of whether or not to remove missing values. If na.rm=TRUE (default), then distances are computed over all pairwise non-missing values. Else missing values are propagated through the distance computation.}
  \item{d}{character string specifying the metric to be used for calculating 
	dissimilarities between vectors. The currently available options are 
	"cosangle" (cosine angle or uncentered correlation distance), "abscosangle" 
	(absolute cosine angle or absolute uncentered correlation distance), 
	"euclid" (Euclidean distance), "abseuclid" (absolute Euclidean distance),
	"cor" (correlation distance), and "abscor" (absolute correlation distance).
	Advanced users can write their own distance functions and add these.}
  \item{M}{a symmetric matrix of pair wise distances.}
  \item{v}{a vector of pair wise distances corresponding to the upper triangle of
	a distance matrix, stored by rows.}
  \item{index}{index in a distance vector, like that returned by \code{dissvector}.}
  \item{p}{number of elements, e.g. the number of rows in a distance matrix.}
}

\value{
	For \code{distancematrix}, a matrix of all pair wise distances between
	rows of 'X'. The value in row 'j' and column 'i' is the distance
	between rows 'i' and 'j'. The matrix is symmetric, and can be converted to
	a vector containing the upper triangle using the function \code{dissvector}. 

	For \code{distancevector}, a vector of all pair wise distances between
	rows of 'X' and the vector 'y'. Entry 'j' is the distance between row 'j'
	of 'X' and the vector 'y'.

	For \code{dissmatrix}, the corresponding distance vector. For 
	\code{dissvector}, the corresponding distance matrix. If 'M' has
	'p' rows (and columns), then 'v' is length 'p*(p-1)/2'.

	For \code{vectmatrix}, the indices of the row and column of a distance
	matrix corresponding to entry \code{index} in the corresponding 
	distance vector.
}

\references{

van der Laan, M.J. and Pollard, K.S. A new algorithm for hybrid hierarchical clustering with visualization and the bootstrap. Journal of Statistical Planning and Inference, 2003, 117, pp. 275-303.

\url{http://www.stat.berkeley.edu/~laan/Research/Research_subpages/Papers/hopach.pdf}
}

\author{Katherine S. Pollard <kpollard@soe.ucsc.edu> and Mark J. van der Laan <laan@stat.berkeley.edu>}

\section{Warning }{The correlation and absolute correlation distance functions call the \code{cor} function, and will therefore fail if there are missing values in the data and na.rm!=TRUE.
}

\seealso{\code{\link{hopach}}, \code{\link{correlationordering}}, \code{\link{disscosangle}}}

\examples{
mydata<-matrix(rnorm(50),nrow=10)
deuclid<-distancematrix(mydata,d="euclid")
vdeuclid<-dissvector(deuclid)
ddaisy<-daisy(mydata)
vdeuclid
ddaisy/sqrt(length(mydata[1,]))

d1<-distancematrix(mydata,d="abscosangle")
d2<-distancevector(mydata,mydata[1,],d="abscosangle")
d1[1,]
d2 #equal to d1[1,]

d3<-dissvector(d1)
pair<-vectmatrix(5,10)
d1[pair[1],pair[2]]
d3[5]
}

\keyword{multivariate}
\keyword{cluster}
