\name{boothopach}
\alias{boothopach}
\alias{bootmedoids}

\title{functions to perform non-parametric bootstrap resampling of hopach clustering results}

\description{
The function \code{boothopach} takes gene expression data and corresponding \code{hopach} gene clustering output and performs non-parametric bootstrap resampling. The medoid genes (cluster profiles) from the original \code{hopach} clustering result are fixed, and in each bootstrap resampled data set, each gene is assigned to the closest medoid. The proportion of bootstrap samples in which each gene appears in each cluster is an estimate of the gene's membership in each cluster. These membership probabilities can be viewed as a "fuzzy" clustering result. The function \code{bootmedoids} take medoids and a distance function, rather than a hopach object, as input.
}

\usage{
boothopach(data, hopachobj, B = 1000, I, hopachlabels = FALSE)

bootmedoids(data, medoids, d = "cosangle", B = 1000, I)
}

\arguments{
  \item{data}{data matrix, data frame or exprSet of gene expression measurements. Each column corresponds to an array, and each row corresponds to a gene. All values must be numeric. Missing values are ignored.}
  \item{hopachobj}{output of the \code{hopach} function.}
  \item{B}{number of bootstrap resampled data sets.}
  \item{I}{number of bootstrap resampled data sets (deprecated, retaining til v1.2 for back compatibility).}
  \item{hopachlabels}{indicator of whether to use the hopach cluster labels 
	\code{hopachobj$clustering$labels} for the row names (TRUE) versus the 
	numbers 0 to 'k-1', where 'k' is the number of clusters (FALSE).}
  \item{medoids}{row indices of \code{data} for the cluster medoids.}
  \item{d}{character string specifying the metric to be used for calculating 
	dissimilarities between vectors. The currently available options are 
	"cosangle" (cosine angle or uncentered correlation distance), "abscosangle" 
	(absolute cosine angle or absolute uncentered correlation distance), 
	"euclid" (Euclidean distance), "abseuclid" (absolute Euclidean distance),
	"cor" (correlation distance), and "abscor" (absolute correlation distance).
	Advanced users can write their own distance functions and add these.}
}

\details{
The function \code{boothopach} requires only data and the corresponding output from the HOPACH clustering algorithm produced by the \code{hopach} function. The function \code{bootmedoids} is designed to work for any clustering result; the user imputs data, medoid row indices, and the distance metric. The supplied distance metrics are the same as for the \code{distancematrix} function. Each non-parametric bootstrap resampled data set consists of resampling the 'n' columns of \code{data} with replacement 'n' times. The distance between each element and each of the medoid elements is computed using \code{d} for each bootstrap data set, and every element is assigned (for that resampled data set) to the cluster whose medoid is closest. These bootstrap cluster assignments are tabulated over all \code{I} bootstrap data sets.
}

\value{
A matrix of bootstrap estimated cluster membership probabilities, which sum to 1 (over the clusters) for each element being clustered. This matrix has one row for each element being clustered and one column for each of the original clusters (one cluster for each medoid). The value in row 'j' and column 'i' is the proportion of the I bootstrap resampled data sets that element 'j' appeared in cluster 'i' (i.e. was closest to medoid 'i').
}

\references{

van der Laan, M.J. and Pollard, K.S. A new algorithm for hybrid hierarchical clustering with visualization and the bootstrap. Journal of Statistical Planning and Inference, 2003, 117, pp. 275-303.

\url{http://www.stat.berkeley.edu/~laan/Research/Research_subpages/Papers/hopach.pdf}

\url{http://www.bepress.com/ucbbiostat/paper107/}

\url{http://www.stat.berkeley.edu/~laan/Research/Research_subpages/Papers/jsmpaper.pdf}

Kaufman, L. and Rousseeuw, P.J. (1990). Finding Groups in Data: An Introduction to Cluster Analysis. Wiley, New York.

}

\author{Katherine S. Pollard <kpollard@gladstone.ucsf.edu> and Mark J. van der Laan <laan@stat.berkeley.edu>}

\seealso{\code{\link{distancematrix}}, \code{\link{hopach}}}

\examples{

#25 variables from two groups with 3 observations per variable
mydata<-rbind(cbind(rnorm(10,0,0.5),rnorm(10,0,0.5),rnorm(10,0,0.5)),cbind(rnorm(15,5,0.5),rnorm(15,5,0.5),rnorm(15,5,0.5)))
dimnames(mydata)<-list(paste("Var",1:25,sep=""),paste("Exp",1:3,sep=""))
mydist<-distancematrix(mydata,d="cosangle") #compute the distance matrix.

#clusters and final tree
clustresult<-hopach(mydata,dmat=mydist)

#bootstrap resampling
myobj<-boothopach(mydata,clustresult)
table(apply(myobj,1,sum)) # all 1
myobj[clustresult$clust$medoids,] # identity matrix
}

\keyword{cluster}
\keyword{nonparametric}
\keyword{multivariate}
