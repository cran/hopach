\name{makeoutput}

\alias{makeoutput}

\title{function to write a text file with hopach output}

\description{
The function \code{makeoutput} takes a data matrix and corresponding \code{hopach} clustering output, plus possibly bootstrap resampling output, and makes a table summarizing the clustering results. The table is written to a tab delimited text file.
}

\usage{
makeoutput(data, hopachobj, bootobj = NULL, file = "HOPACH.out", 
gene.names = NULL)
}

\arguments{
  \item{data}{data matrix, data frame or exprSet of gene expression measurements. Typically, each column corresponds to an array, and each row corresponds to a gene. For clustering arrays, the arrays appear in the rows and the genes in the columns. All values must be numeric. Missing values are ignored.}
  \item{hopachobj}{output of the \code{hopach} function.}
  \item{bootobj}{optional output of \code{boothopach} or \code{bootmedoids} - a matrix of bootstrap estimated cluster membership probabilities, with a row for each row in \code{data} and a column for each cluster.}
  \item{file}{filename for the table produced.}
  \item{gene.names}{optional names or annotations for the genes (arrays), which can be different from the row names of \code{data}}
}

\details{
The output table contains information about the rows of \code{data}: Index, UID and Name; the main cluster results: Cluster Number, Cluster Label (from \code{hopach}), and Cluster Level Order (the ordering of the elements in the level of the \code{hopach} tree at which the main clusters were identified); and the final level of the tree: Final Label, and Final Order (the ordering of elements in the final level of the tree). Sorting this table on Index results in the rows having the same order as in data. Sorting on Cluster Level Order results in the rows being ordered by cluster, and then within cluster based on the value of the argument \code{ord} to \code{hopach} (default is distance to the medoid). Sorting on Final Level Order results in the rows being ordered as in the leaves of the \code{hopach} tree, where clusters are still ordered and elements near each other in the ordering will have small pairwise distances.
}

\value{
The function \code{makeoutput} has no value. It writes a tab delimited text file to the current working directory.
}

\references{

van der Laan, M.J. and Pollard, K.S. A new algorithm for hybrid hierarchical clustering with visualization and the bootstrap. Journal of Statistical Planning and Inference, 2003, 117, pp. 275-303.

\url{http://www.stat.berkeley.edu/~laan/Research/Research_subpages/Papers/hopach.pdf}
}

\author{Katherine S. Pollard <kpollard@soe.ucsc.edu>}

\note{Thank you to Karen Vranizan <vranizan@uclink.berkeley.edu> for helping to write this function.}

\seealso{\code{\link{hopach}}, \code{\link{boothopach}}, \code{\link{bootmedoids}}}

\examples{

#25 variables from two groups with 3 observations per variable
mydata<-rbind(cbind(rnorm(10,0,0.5),rnorm(10,0,0.5),rnorm(10,0,0.5)),cbind(rnorm(15,5,0.5),rnorm(15,5,0.5),rnorm(15,5,0.5)))
dimnames(mydata)<-list(paste("Var",1:25,sep=""),paste("Exp",1:3,sep=""))
mydist<-distancematrix(mydata,d="cosangle") #compute the distance matrix.

#clusters and final tree
clustresult<-hopach(mydata,dmat=mydist)

#bootstrap resampling
myobj<-boothopach(mydata,clustresult)

#write output file
makeoutput(mydata,clustresult,myobj)
}


\keyword{cluster}
\keyword{nonparametric}
\keyword{multivariate}

