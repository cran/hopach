\name{labelstomss}

\alias{labelstomss}
\alias{labelstosil}
\alias{medstosil}
\alias{msscheck}
\alias{silcheck}

\title{Functions to compute silhouettes and split silhouettes} 

\description{
Silhouettes measure how well an element belongs to its cluster, and the average silhouette measures the strength of cluster membership overall. 
The Median (or Mean) Split Silhouette (MSS) is a measure of cluster 
heterogeneity. Given a partitioning of elements into groups, the MSS
algorithm considers each group separately and computes the split
silhouette for that group, which evaluates evidence in favor of
further splitting the group. If the median (or mean) split silhouette
over all groups in the partition is low, the groups are homogeneous.
}

\usage{
labelstomss(labels, dist, khigh = 9, within = "med", between = "med", 
hierarchical = TRUE)

labelstosil(labels, dist)

medstosil(medoids, dist)

msscheck(dist, kmax = 9, khigh = 9, within = "med", between = "med", 
    force = FALSE, echo = FALSE, graph = FALSE)

silcheck(data, kmax = 9, diss = FALSE, echo = FALSE, graph = FALSE)
}

\arguments{
  \item{labels}{vector of cluster labels for each element in the set.}
  \item{dist}{numeric distance matrix containing the pair wise distances 
	between all elements. All values must be numeric and missing values
	are not allowed.}
  \item{medoids}{a vector indicating the rows/cols of \code{dist} that are the cluster medoids, i.e. profiles (or centroids) for each cluster.}
  \item{data}{a data matrix. Each column corresponds to an observation, and each row corresponds to a variable. In the gene expression context, observations are arrays and variables are genes. All values must be numeric. Missing values are ignored. In \code{silcheck}, \code{data} may also be a distance matrix or dissimilarity object if the argument \code{diss=TRUE}.}
  \item{khigh}{integer between 1 and 9 specifying the maximum number of 
	children for each cluster when computing MSS.}
  \item{kmax}{integer between 1 and 9 specifying the maximum number of clusters to consider. Can be different from khigh, though typically these are the same value.}
  \item{within}{character string indicating how to compute the split silhouette
	for each cluster. The available options are "med" (median over all
	elements in the cluster) or "mean" (mean over all elements in the 
	cluster).}
  \item{between}{character string indicating how to compute the MSS over all
	clusters. The available options are "med" (median over all
	clusters) or "mean" (mean over all clusters). Recommended to use
	the same value as \code{within}.}
  \item{hierarchical}{logical indicating if 'labels' should be treated as
	encoding a hierarchical tree, e.g. from HOAPCH.}
  \item{force}{indicator of whether to require at least 2 clusters, if FALSE (default), one cluster is considered.}
  \item{echo}{indicator of whether to print the selected number of clusters and corresponding MSS.}
  \item{graph}{indicator of whether to generate a plot of MSS (or average silhouette in \code{silcheck}) versus number of clusters.}
  \item{diss}{idicator of whether \code{data} is a dissimilarity matrix (or dissimilarity object), as in the \code{pam} function of the \code{cluster} package. If TRUE then \code{data} will be considered as a dissimilarity matrix. If FALSE, then \code{data} will be considered as a data matrix (observations by variables).} 
}

\details{
	The Median (and mean) Split Silhouette (MSS) criteria is defined in 
	paper107 listed in the references (below). This criteria is based on the
	criteria function 'silhouette', proposed by Kaufman and Rousseeuw (1990).
	While average silhouette is a good global measure of cluster strength, MSS
	was developed to be more "aggressive" for finding small, homogeneous clusters
	in large data sets. MSS is a measure of average cluster homogeneity. The
	Median version is more robust than the Mean.
}

\value{
For \code{labelstomss}, the median (or mean or combination) split silhouette, depending on the values of \code{within} and \code{between}. \cr

For \code{medstosil} and \code{labelstosil}, a list with first component the cluster label for each element and second compenent the silhouette for that element. The average silhouette is simply the mean of the second component.

For \code{msscheck}, a vector with first component the chosen number of clusters (minimizing MSS) and second component the corresponding MSS.

For \code{silcheck}, a vector with first component the chosen number of clusters (maximizing average silhouette) and second component the corresponding average silhouette.

}

\references{

\url{http://www.bepress.com/ucbbiostat/paper107/}

\url{http://www.stat.berkeley.edu/~laan/Research/Research_subpages/Papers/jsmpaper.pdf}

Kaufman, L. and Rousseeuw, P.J. (1990). Finding Groups in Data: An Introduction to Cluster Analysis. Wiley, New York.

}


\author{Katherine S. Pollard <kpollard@wald.ucdavis.edu> and Mark J. van der Laan <laan@stat.berkeley.edu>}

\seealso{\code{pam}, \code{\link{hopach}}, \code{\link{distancematrix}}}

\examples{

mydata<-rbind(cbind(rnorm(10,0,0.5),rnorm(10,0,0.5),rnorm(10,0,0.5)),cbind(rnorm(15,5,0.5),rnorm(15,5,0.5),rnorm(15,5,0.5)))
mydist<-distancematrix(mydata,d="cosangle") #compute the distance matrix.

#pam
result1<-pam(mydata,k=2)
result2<-pam(mydata,k=5)
labelstomss(result1$clust,mydist,hierarchical=FALSE)
labelstomss(result2$clust,mydist,hierarchical=FALSE)

#hopach
result3<-hopach(mydata,dmat=mydist)
labelstomss(result3$clustering$labels,mydist)
labelstomss(result3$clustering$labels,mydist,within="mean",between="mean")

}

\keyword{cluster}
\keyword{multivariate}
