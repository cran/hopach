\name{bootplot}

\alias{bootplot}

\title{function to make a barplot of bootstrap estimated cluster membership probabilities}

\description{
After clustering, the \code{boothopach} or \code{bootmedoids} function can be used to estimated the membership of each element being clustered in each of the identified clusters (fuzzy clustering). The proportion of bootstrap resampled data sets in which each element is assigned to each cluster is called the "reappearance proportion" for the element and that cluster. This function plots these proportions in a colored barplot.
}

\usage{
bootplot(bootobj, hopachobj, ord = "bootp", main = NULL, labels = NULL, 
showclusters = TRUE, ...)
}

\arguments{
  \item{bootobj}{output of \code{boothopach} or \code{bootmedoids} applied to the genes - a matrix of bootstrap estimated cluster membership probabilities, with a row for each row in \code{data} and a column for each cluster.}
  \item{hopachobj}{output of the \code{hopach} function. If \code{bootobj} was generated using \code{bootmedoids} (i.e. \code{hopach} was not run), then the \code{bootplot} function can be used by creating a \code{hopachobj} which is a list with at least the following two components: \code{hopachobj$clustering$sizes} (number of elements in each cluster - length should be \code{ncol(bootobj)} and \code{hopachobj$clustering$order} (an ordering of the elements so that elements in the same cluster appear next to each other and elements may also be ordered within cluster). By changing the value of \code{hopachobj$clustering$order}, the order of the elements in the barplot can be altered.}
  \item{ord}{character string indicating how to order the elements (rows) in the barplot. If ord="none", then the elements are plotted in the same order as in \code{bootobj}, i.e. the same order as the original data matrix. If ord="final", the ordering of elements in the final level of the \code{hopach} hierarchical tree is used. If ord="cluster", the ordering from the level of the \code{hopach} tree corresponding to the main clusters is used. If ord="bootp", the elements are ordered first by main cluster and then by bootstrap reappearance proportion within cluster, so that elements with the highest membership in the cluster appear at the bottom. In the last three cases, the elements from each cluster will be contiguous. If ord="final", then the medoid element will appear in the middle of each cluster. If ord="clust", the ordering depends on the value of the \code{ord} argument passed to the \code{hopach} function.  For example, when ord="own" in \code{hopach}, the elements are ordered within cluster based on distance to the medoid, so that the medoid appears first (at the bottom) in the cluster.}
  \item{main}{character string to be used as the main title}
  \item{labels}{a vector of labels for the elements being clustered to be used on the axes. If the number of elements is lager than 50, the labels are not shown.}
  \item{showclusters}{indicator of whether or not to show the cluster boundaries on the plot. If show.clusters=TRUE, solid lines are drawn at the edges of the clusters.}
  \item{\dots}{additional arguments to the \code{barplot} plotting function}
}

\details{
Each cluster (column of \code{bootobj}) is represented by a color. The proportion of bootstrap resampled data sets in which an element appeared in that cluster determines the proportion of the bar for that element which is the corresponding color. As a key, the clusters are labeled on the right margin in text of the same color. 
}

\value{The function \code{bootplot} has no value. It does generate a plot.
}

\references{

van der Laan, M.J. and Pollard, K.S. A new algorithm for hybrid hierarchical clustering with visualization and the bootstrap. Journal of Statistical Planning and Inference, 2003, 117, pp. 275-303.

%\begin{verbatim}
\url{http://www.stat.berkeley.edu/~laan/Research/Research_subpages/Papers/hopach.pdf}
%\end{verbatim}
%\emph{}

}

\author{Katherine S. Pollard <kpollard@soe.ucsc.edu>}

\note{Thank you to Sandrine Dudoit <sandrine@stat.berkeley.edu> for her input and to Jenny Bryan for the original clusplot code.}

\seealso{\code{\link{hopach}}, \code{\link{boothopach}}, \code{\link{bootmedoids}}, \code{\link{barplot}}}

\examples{
mydata<-rbind(cbind(rnorm(10,0,0.5),rnorm(10,0,0.5),rnorm(10,0,0.5)),cbind(rnorm(15,5,0.5),rnorm(15,5,0.5),rnorm(15,5,0.5)))
dimnames(mydata)<-list(paste("Var",1:25,sep=""),paste("Exp",1:3,sep=""))
mydist<-distancematrix(mydata,d="euclid")

#hopach clustering
clustresult<-hopach(mydata,dmat=mydist)

#bootstrap
myobj<-boothopach(mydata,clustresult)

#plots
bootplot(myobj,clustresult,showclusters=FALSE)
bootplot(myobj,clustresult,labels=paste("Sample",LETTERS[1:25],sep=" "))
}

\keyword{cluster}
\keyword{nonparametric}
\keyword{multivariate}
