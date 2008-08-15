\name{prune}

\alias{prune}

\title{function to prune a HOPACH tree.}

\description{
The \code{hopach} clustering function identifies a level of the tree with minimum MSS as the main clusters and also runs the tree down all the way to the final level. The \code{prune} function allows one to access a level of the tree other than the main clusters or the final level.}

\usage{
prune(data,hobj,level,dmat=NULL,ord="own")
}

\arguments{
  \item{data}{data matrix, data frame or exprSet of gene expression measurements. Typically, each column corresponds to an array, and each row corresponds to a gene. For clustering arrays, the arrays appear in the rows and the genes in the columns. All values must be numeric. Missing values are ignored.}
  \item{hobj}{output of the \code{hopach} function.}
  \item{level}{an integer specifying the level to which the tree should be pruned - can be greater than or less than the level with the main clusters.}
  \item{dmat}{matrix of pair wise distances between all genes (arrays). All values must be numeric. If NULL, this matrix is computed using the metric specified by the 'metric' given in hobj. If a matrix is provided, the user is responsible for ensuring that the metric used agrees with that used in computing hobj.}
  \item{ord}{character string indicating which of the two orderings produced by \code{hopach} should be used for the plot. If ord="final", the ordering of elements in the final level of the hierarchical tree is used. If ord="cluster", the ordering from the level of the tree corresponding to the main clusters is used. In both cases, the elements from each cluster will be contiguous. If ord="final", then the medoid element will appear in the middle of each cluster. Else, the ordering depends on the value of \code{ord} passed to the \code{hopach} function. If ord="none", then the elements are plotted in the same order as in dist.}
}

\value{A list with the same components as are returned by the \code{hopach} function. The clustering component now contains the specified level instead of the main clusters.}

\references{

van der Laan, M.J. and Pollard, K.S. A new algorithm for hybrid hierarchical clustering with visualization and the bootstrap. Journal of Statistical Planning and Inference, 2003, 117, pp. 275-303.

\url{http://www.stat.berkeley.edu/~laan/Research/Research_subpages/Papers/hopach.pdf}
}

\author{Katherine S. Pollard <kpollard@gladstone.ucsf.edu> and Mark J. van der Laan <laan@stat.berkeley.edu>}

\seealso{\code{\link{hopach}}, \code{\link{makeoutput}}}

\examples{
mydata<-matrix(rnorm(600),nrow=100)
mydist<-distancematrix(mydata,d="cosangle")
clustresult<-hopach(mydata,dmat=mydist)
level2<-prune(mydata,clustresult,level=2,dmat=mydist,ord="own")
clustresult$clustering$k
level2$clustering$k
}

\keyword{cluster}
\keyword{multivariate}

