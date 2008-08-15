\name{hopach-internal}

\alias{collap}
\alias{cutdigits}
\alias{cutzeros}
\alias{digits}
\alias{msscollap}
\alias{msscomplete}
\alias{mssinitlevel}
\alias{mssmulticollap}
\alias{mssnextlevel}
\alias{mssrundown}
\alias{msssplitcluster}
\alias{newnextlevel}
\alias{newsplitcluster}
\alias{nonzeros}
\alias{orderelements}
\alias{paircoll}

\title{Functions used internally by the hopach package}

\description{
The Hierarchical Ordered Partitioning and Collapsing Hybrid (HOPACH)
clustering algorithm builds a hierarchical tree of clusters. These functions are used internally to create and order new levels of the tree.
}

\usage{

collap(data, level, d = "cosangle", dmat = NULL, newmed = "medsil")

cutdigits(labels, dig)

cutzeros(labels) 

digits(label)

msscollap(data, level, khigh, d = "cosangle", dmat = NULL, 
newmed = "medsil", within = "med", between = "med", impr = 0)

msscomplete(level, data, K = 16, khigh = 9, d = "cosangle", 
dmat = NULL, within = "med", between = "med", verbose=FALSE)

mssinitlevel(data, kmax = 9, khigh = 9, d = "cosangle", 
dmat = NULL, within = "med", between = "med", ord = "co", verbose=FALSE)

mssmulticollap(data, level, khigh, d = "cosangle", dmat = NULL, 
newmed = "medsil", within = "med", between = "med", impr = 0)

mssnextlevel(data, prevlevel, dmat, kmax = 9, khigh = 9, 
within = "med", between = "med") 

mssrundown(data, K = 16, kmax = 9, khigh = 9, d = "cosangle", 
dmat = NULL, initord = "co", coll = "seq", newmed = "medsil", 
stop = TRUE, finish = FALSE, within = "med", between = "med", 
impr = 0, verbose = FALSE) 

msssplitcluster(clust1, l1, id1, medoid1, med2dist, right, dist1, 
kmax = 9, khigh = 9, within = "med", between = "med")

newnextlevel(data, prevlevel, dmat, klow = 2, khigh = 6) 

newsplitcluster(clust1, l1, id1, klow = 2, khigh = 2, medoid1, 
med2dist, right, dist1) 

nonzeros(labels)

orderelements(level, data, rel = "own", d = "cosangle", dmat = NULL) 

paircoll(i, j, data, level, d = "cosangle", dmat = NULL, 
newmed = "medsil") 

}

\arguments{

\item{data}{A numeric data matrix.}
\item{level, prevlevel}{A list representing a level of the hopach hierarchical tree with the following components: number of clusters, row indices of medoids,cluster sizes, labels, indicator of whether this is the final level, and a matrix containing the label and medoid for each node in the tree.}
\item{d}{A character string specifying the metric to be used for calculating 
	dissimilarities between variables. The currently available options are 
	"cosangle" (cosine angle or uncentered correlation distance), "abscosangle" 
	(absolute cosine angle or absolute uncentered correlation distance), 
	"euclid" (Euclidean distance), "abseuclid" (absolute Euclidean distance),
	"cor" (correlation distance), and "abscor" (absolute correlation distance).}
\item{dmat}{A numeric matrix of pair wise distances. Missing values are not allowed. If NULL, this matrix is computed using the metric specified by \code{d}. If a matrix is provided, the user is responsible for ensuring that the metric used agrees with \code{d}.}
\item{newmed}{A character string specifying how to choose a medoid for the new
	cluster after collapsing a pair of clusters. The options are "medsil" 
	(maximizer of medoid based silhouette, i.e.: (a-b)/max(a,b), where 
	a is distance to medoid and b is distance to next closest medoid),
	"nn" (nearest neighbor of mean of two collapsed cluster medoids weighted
	by cluster size), "uwnn" (unweighted version of nearest neighbor, i.e. 
	each cluster - rather than each element - gets equal weight), "center" 
	(minimizer of average distance to the medoid).}
\item{labels}{A vector of labels encoding the history of each element through the hopach tree, with one digit for each level of the tree.}
\item{dig}{Number of digits to remove from \code{labels}.}
\item{label}{A single label, e.g. one entry in the vector \code{labels}.}
\item{khigh}{An integer between 1 and 9 specifying the maximum number of children at each node in the tree when computing MSS. Can be different from kmax, though typically these are the same value.}
\item{K}{A positive integer specifying the maximum number of levels in the tree. Must be 16 or less, due to computational limitations (overflow).}
\item{kmax}{An integer between 1 and 9 specifying the maximum number of children at each node in the tree.}
\item{klow}{An integer specifying the minimum number of children at each node in the tree.}
\item{within}{A character specifying what summary measure to use when combining split silhouettes within a cluster in computation of the MSS criteria. The options are "med" (median split silhouette) or "mean" (mean split silhouette).}
\item{between}{A character specifying what summary measure to use when combining mean or median split silhouettes across clusters in computation of the MSS criteria. The options are "med" (median) or "mean" (mean). Typically the same value as \code{within}.}
\item{impr}{A number between 0 and 1 specifying the margin of improvement in MSS needed to accept a collapse step. If (MSS.before - MSS.after)/MSS.before is less than \code{impr}, then the collapse is not performed.}
\item{ord, initord}{A character string specifying how to order the clusters in the initial level of the tree. The options are "co" (maximize correlation ordering, i.e. the empirical correlation between distance apart in the ordering and distance between the cluster medoids) or "clust" (apply hopach with binary splits to the cluster medoids and use the final level of that tree as the ordering). In subsequent levels, the clusters are ordered relative to the previous level, so this initial ordering determines the overall structure of the tree.}
\item{coll}{A character string specifying how collapsing steps are performed at each level. The options are "seq" (begin with the closest pair of clusters and collapse pairs sequentially as long as MSS decreases) and "all" (consider all pairs of clusters and collapse any that decrease MSS).}
\item{stop}{An indicator of whether the hopach algorithm should stop and identify "main clusters" as soon as MSS no longer improves by at least \code{impr}.}
\item{finish}{An indicator of whether the hopach algorithm should compute all K levels of the tree and return the one with minimum MSS. When \code{finish=FALSE} and \code{stop=FALSE}, level \code{K} is returned.}
\item{clust1}{A numeric data matrix, typically the data for the elements in a single cluster that is going to be split.}
\item{l1}{A single label shared by the elements in \code{clust1}.}
\item{id1}{The row indices (in the original data matrix) for the elements in \code{clust1}.}
\item{medoid1}{The row index of the medoid for the elements in \code{clust1}.}
\item{med2dist}{Distance between each element in \code{clust1} and the neighboring cluster's medoid. This is the next medoid to the right, unless \code{clust1} is the last cluster (\code{right=FALSE}), in which case it is the next medoid to the left.}
\item{right}{Indicator of whether the cluster's neighbor is to the right.}
\item{dist1}{The distance matrix for elements in \code{clust1}.}
\item{rel}{A character string specifying how to order the elements within clusters. The options are "own" (order based on distance from cluster medoid with medoid first, i.e. leftmost), "neighbor" (order based on distance to the medoid of the next cluster to the right), or "co" (maximize correlation ordering - can be slow for large clusters!).}
\item{i, j}{Indices of two clusters to be collapsed.}
\item{verbose}{Indicates whether or not to print verbose output.}
}

\value{

The following functions return a list encoding a level of the hopach tree: \code{collap}, \code{msscollap}, \code{paircoll}, \code{mssmulticollap}, \code{mssinitlevel}, \code{msscomplete}, \code{mssnextlevel}, \code{mssrundown}, \code{newnextlevel}.

The following functions return the components for a new level of the tree corresponding to splitting a single cluster: \code{msssplitcluster} and \code{newsplitcluster}. These are called by \code{mssnextlevel} and \code{newnextlevel}, respectively.

The function \code{orderelements} returns a list with two components: (i) the ordered data and (ii) the vector of indices for ordering (as returned by \code{order}).

The function \code{cutdigits} returns truncated labels.

The function \code{cutzeros} returns labels with trailing zeros removed.

The function \code{digits} returns an integer denoting the number of digits in \code{label}.

The function \code{nonzeros} returns a vector denoting the number of digits in \code{labels} after removing trailing zeros.

}

\references{

van der Laan, M.J. and Pollard, K.S. A new algorithm for hybrid hierarchical clustering with visualization and the bootstrap. Journal of Statistical Planning and Inference, 2003, 117, pp. 275-303.

\url{http://www.stat.berkeley.edu/~laan/Research/Research_subpages/Papers/hopach.pdf}

\url{http://www.bepress.com/ucbbiostat/paper107/}

\url{http://www.stat.berkeley.edu/~laan/Research/Research_subpages/Papers/jsmpaper.pdf}

Kaufman, L. and Rousseeuw, P.J. (1990). Finding Groups in Data: An Introduction to Cluster Analysis. Wiley, New York.

}

\author{Katherine S. Pollard <kpollard@gladstone.ucsf.edu> and Mark J. van der Laan <laan@stat.berkeley.edu>}

\seealso{\code{\link{hopach}}}

\examples{
#These are not user level functions
#See: hopach examples
#? hopach
}

\keyword{internal}
\keyword{cluster}
\keyword{multivariate}

