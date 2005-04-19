\name{hopach}

\alias{hopach}

\title{function to perform HOPACH hierarchical clustering}

\description{
The Hierarchical Ordered Partitioning and Collapsing Hybrid (HOPACH)
clustering algorithm builds a hierarchical tree by recursively 
partitioning a data set (e.g., gene expression measurements)
with the PAM algorithm, while ordering
and possibly collapsing clusters at each level. The algorithm uses
the Mean/Median Split Silhouette (MSS) criteria to identify the
level of the tree with maximally homogeneous clusters. It also runs
the tree down to produce a final ordered list of the elements.
}

\usage{
hopach(data, dmat = NULL, d = "cosangle", clusters = "best", K = 15, 
kmax = 9, khigh = 9, coll = "seq", newmed = "medsil", mss = "med", 
impr = 0, initord = "co", ord = "own")
}

\arguments{
  \item{data}{data matrix or data frame. All values must be numeric. Missing values are ignored.}
  \item{dmat}{matrix of pair wise distances between all elements. All values must be
	numeric, and missing values are not allowed. If NULL, this matrix is computed 
	using the metric specified by \code{d}. If a matrix is provided, the user is
	responsible for ensuring that the metric used agrees with \code{d}.}
  \item{d}{character string specifying the metric to be used for calculating 
	dissimilarities between variables. The currently available options are 
	"cosangle" (cosine angle or uncentered correlation distance), "abscosangle" 
	(absolute cosine angle or absolute uncentered correlation distance), 
	"euclid" (Euclidean distance), "abseuclid" (absolute Euclidean distance),
	"cor" (correlation distance), and "abscor" (absolute correlation distance).
	Advanced users can write their own distance functions and add these to the
	functions \code{distancematrix()} and \code{distancevector()}.}
  \item{clusters}{character string specifying if clusters are to be identified
	as the level of the tree with the minimum mean/median split 
	silhouette (MSS) ("best"), the first level of the tree below which MSS 
	increases ("greedy"), or not at all ("none").}
  \item{K}{positive integer specifying the maximum number of levels in the tree. 
	Must be 15 or less, due to computational limitations (overflow).}
  \item{kmax}{integer between 1 and 9 specifying the maximum number of children at
	each node in the tree.}
  \item{khigh}{integer between 1 and 9 specifying the maximum number of children
	at each node in the tree when computing MSS. Can be different from kmax, 
	though typically these are the same value.}
  \item{coll}{character string specifying how collapsing steps are performed at each 
	level. The options are "seq" (begin with the closest pair of clusters and
	collapse pairs sequentially as long as MSS decreases) and "all" (consider all
	pairs of clusters and collapse any that decrease MSS).}
  \item{newmed}{character string specifying how to choose a medoid for the new
	cluster after collapsing a pair of clusters. The options are "medsil" 
	(maximizer of medoid based silhouette, i.e.: (a-b)/max(a,b), where 
	a is distance to medoid and b is distance to next closest medoid),
	"nn" (nearest neighbor of mean of two collapsed cluster medoids weighted
	by cluster size), "uwnn" (unweighted version of nearest neighbor, i.e. 
	each cluster - rather than each element - gets equal weight), "center" 
	(minimizer of average distance to the medoid).}
  \item{mss}{character vector specifying what criteria function to use. The options
	are "med" (median split silhouette) or "mean" (mean split silhouette). See
	details for definition of split silhouettes. The MSS criteria is used to
	determine the number of children at each node, to decide what collapsing 
	should be performed at each level, and to determine the main clusters.}
  \item{impr}{number between 0 and 1 specifying the margin of improvement in MSS 
	needed to accept a collapse step. If (MSS.before - MSS.after)/MSS.before 
	is less than \code{impr}, then the collapse is not performed.}
  \item{initord}{character string specifying how to order the clusters in the initial
	level of the tree. The options are "co" (maximize correlation ordering, i.e. 
	the empirical correlation between distance apart in the ordering and distance
	between the cluster medoids) or "clust" (apply hopach with binary splits to 
	the cluster medoids and use the final level of that tree as the ordering). In
	subsequent levels, the clusters are ordered relative to the previous level, so
	this initial ordering determines the overall structure of the tree.}
  \item{ord}{character string specifying how to order the elements within clusters. This
	method is used to create an ordering of all elements at the level of the tree 
	corresponding to the main clusters. The options are "own" (order based on distance
	from cluster medoid with medoid first, i.e. leftmost), "neighbor" (order based on
	distance to the medoid of the next cluster to the right), or "co" (maximize
	correlation ordering - can be slow for large clusters!).}
}

\details{
The HOPACH hierarchical clustering algorithm is a hybrid between an agglomerative (bottom
up) and a divisive (top down) algorithm. The HOPACH tree is built from the root node (all 
elements) down to the leaf nodes, but at each level collapsing steps are used to unite
similar clusters. In addition, the clusters in each level are ordered with a deterministic
algorithm based on the same distance metric that is used in the clustering. In this way,
the ordering produced in the final level of the tree does not depend on the order of the
data in the original data set (as can be the case with algorithms that have a random 
component in their ordering methods). Unlike other hierarchical clustering methods, HOPACH 
builds a tree of clusters in which the nodes need not be binary, i.e. there can be more
than two children at each split. The divisive steps of the HOPACH algorithm are performed
using the PAM algorithm described in chapter 2 of Kaufman and Rousseeuw (1990) and the R
package 'cluster'.

The Median (or Mean) Split Silhouette (MSS) criteria is used by HOPACH to (i) determine 
the optimal number of children at each node, (ii) decide which pairs of clusters to 
collapse at each level, and (iii) identify the first level of the tree with maximally
homogeneous clusters. In each case, the goal is to minimize MSS, which is a measure of
cluster heterogeneity described in http://www.bepress.com/ucbbiostat/paper107/. 
}

\value{A list with the following components:

  \item{clustering}{the partitioning or 'main clusters' with the following components

	'k' is an integer specifying the number of clusters identified by minimizing MSS.

	'medoids' is a vector indicating the rows of \code{data} that are the 'k'
	cluster medoids, i.e. profiles (or centroids) for each cluster.

	'sizes' is a vector containing the 'k' cluster sizes.

	'labels' is a vector containing the main cluster labels for every variable. Each
	label consists of one digit per level of the tree (up to the level identified as
	the main clusters). The digit (1-9) indicates which child cluster the variable
	was in at that level. For example, '124' means the fist (leftmost in the tree)
	cluster in level 1, the second child of cluster '1' in level 2, and the fourth 
	child of cluster '12' in level 3. These can be mapped to the numbers 1:k for
	simplicity, though the tree structure and relationship amongst the clusters is
	then lost, e.g. 1211 is closer to 1212 than to 1221.

	'order' is a vector containing the ordering of variables within the main clusters.
	The clusters are ordered deterministically as the tree is built. The elements within
	each of the main clusters are ordered with the method determined by the value of
	\code{ord}: "own" (relative to own medoid), "neighbor" (relative to next medoid
	to the right), or "co" (maximize correlation ordering).
  }

  \item{final}{the final level of the hierarchical tree with the following components

	'labels' is a vector containing the final labels for every variable. Each label
	consists of one digit per level of the tree (up to the final level), and the 
	format for the labels is the same as for the clustering labels. The final labels
	contain the entire history of the tree. In fact, internal level 'n' can be 
	reproduced by truncating the final labels to 'n' digits. Ordering the final 
	labels produces the final ordering (final level of the tree), while ordering
	internal level labels produces an ordering of the clusters at that level.

	'order' is a vector containing the ordering of variables at the final level of
	the tree. Essentially, this is the numeric ordering of the final labels. Due to
	the limit on the largest possible integer (overflow), the final labels can have 
	at most 16 digits, i.e. the tree can have at most 16 levels. For large data sets,
	this may not be enough partitioning steps to result in final nodes (leaves) with
	only one variable each. Furthermore, PAM can not partition a node of size 3 or 
	less, so that leaves may contain 2 or 3 variables regardless of the number of
	levels in the tree. Hence, the final ordering of variables is completed by
	ordering the variables in any leaf of size 2 or larger with the method determined
	by the value of \code{ord}: "own" (relative to own medoid), "neighbor" (relative 
	to next medoid to the right), or "co" (maximize correlation ordering).


	'medoids' is a matrix containing the labels and corresponding medoids for each
	internal node and leaf of the tree. The number of digits in the label indicates
	the level for that node. The medoid refers to a row of \code{data}
  }

  \item{call}{the matched 'call' generating the HOPACH output}

  \item{metric}{the distance metric}
}

\references{

van der Laan, M.J. and Pollard, K.S. A new algorithm for hybrid hierarchical clustering with visualization and the bootstrap. Journal of Statistical Planning and Inference, 2003, 117, pp. 275-303.

\url{http://www.stat.berkeley.edu/~laan/Research/Research_subpages/Papers/hopach.pdf}

\url{http://www.bepress.com/ucbbiostat/paper107/}

\url{http://www.stat.berkeley.edu/~laan/Research/Research_subpages/Papers/jsmpaper.pdf}

Kaufman, L. and Rousseeuw, P.J. (1990). Finding Groups in Data: An Introduction to Cluster Analysis. Wiley, New York.

}

\author{Katherine S. Pollard <kpollard@soe.ucsc.edu> and Mark J. van der Laan <laan@stat.berkeley.edu>}

\notes{Thank you to Karen Vranizan <vranizan@uclink.berkeley.edu> for her input
}

\seealso{\code{\link{distancematrix}}, \code{\link{labelstomss}}, \code{\link{boothopach}}, \code{pam}}

\examples{

#25 variables from two groups with 3 observations per variable
mydata<-rbind(cbind(rnorm(10,0,0.5),rnorm(10,0,0.5),rnorm(10,0,0.5)),cbind(rnorm(15,5,0.5),rnorm(15,5,0.5),rnorm(15,5,0.5)))
dimnames(mydata)<-list(paste("Var",1:25,sep=""),paste("Exp",1:3,sep=""))
mydist<-distancematrix(mydata,d="cosangle") #compute the distance matrix.

#clusters and final tree
clustresult<-hopach(mydata,dmat=mydist)
clustresult$clustering$k #number of clusters.
dimnames(mydata)[[1]][clustresult$clustering$medoids] #medoids of clusters.
table(clustresult$clustering$labels) #equal to clustresult$clustering$sizes.

#faster, sometimes fewer clusters
greedyresult<-hopach(mydata,clusters="greedy",dmat=mydist) 

#only get the final ordering (no partitioning into clusters)
orderonly<-hopach(mydata,clusters="none",dmat=mydist)

#cluster the columns (rather than rows)
colresult<-hopach(t(mydata),dmat=distancematrix(t(mydata),d="euclid"))

}

\keyword{cluster}
\keyword{multivariate}
