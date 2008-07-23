\name{boot2fuzzy}

\alias{boot2fuzzy}

\title{function to write MapleTree files for viewing bootstrap estimated cluster membership probabilities based on hopach clustering results}

\description{
The MapleTree software (http://mapletree.sourceforge.net/) is an open source, cross-platform, visualization tool to graphically browse results of cluster analyses. The \code{boot2fuzzy} function takes a data matrix, plus corresponding \code{hopach} clustering output and bootstrap resampling output, and writes the (.cdt, .fct, and .mb) files needed to view these "fuzzy clustering" results in MapleTree.
}

\usage{
boot2fuzzy(data, bootobj, hopach.genes, hopach.arrays = NULL, 
file="hopach", clust.wts = NULL, gene.wts = NULL, array.wts = NULL, 
gene.names = NULL)
}

\arguments{
  \item{data}{data matrix, data frame or exprSet of gene expression measurements. Each column corresponds to an array, and each row corresponds to a gene. All values must be numeric. Missing values are ignored.}
  \item{bootobj}{output of \code{boothopach} or \code{bootmedoids} applied to the genes - a matrix of bootstrap estimated cluster membership probabilities, with a row for each row in \code{data} and a column for each cluster.}
  \item{hopach.genes}{output of the \code{hopach} function applied to genes (rows of \code{data}.}
  \item{hopach.arrays}{optional output of the \code{hopach} function applied to arrays (columns of \code{data}.}
  \item{file}{name for the output files (the extensions .cdt, .mb and .fct will be added).}
  \item{clust.wts}{an optional vector of numeric weights for the clusters.}
  \item{gene.wts}{an optional vector of numeric weights for the genes.}
  \item{array.wts}{an optional vector of numeric weights for the arrays.}
  \item{gene.names}{optional vector of names or annotations for the genes, which can be different from the row names of \code{data}}
}

\value{
The function \code{boot2fuzzy} has no value. It writes three text files to the current working directory.
}

\references{

van der Laan, M.J. and Pollard, K.S. A new algorithm for hybrid hierarchical clustering with visualization and the bootstrap. Journal of Statistical Planning and Inference, 2003, 117, pp. 275-303.

%\begin{verbatim}
\url{http://www.stat.berkeley.edu/~laan/Research/Research_subpages/Papers/hopach.pdf}
%\end{verbatim}
%\emph{}

}

\author{Katherine S. Pollard <kpollard@soe.ucsc.edu>}

\note{Thank you to Lisa Simirenko <lsimirenko@lbl.gov> for providing HOPACH views in MapleTree, and to Karen Vranizan <vranizan@uclink.berkeley.edu> for her input. The MapleTree software can be downloaded from: http://sourceforge.net/projects/mapletree/}

\seealso{\code{\link{hopach}}, \code{\link{boothopach}}, \code{\link{bootmedoids}}, \code{\link{hopach2tree}}}

\examples{

#25 variables from two groups with 3 observations per variable
mydata<-rbind(cbind(rnorm(10,0,0.5),rnorm(10,0,0.5),rnorm(10,0,0.5)),cbind(rnorm(15,5,0.5),rnorm(15,5,0.5),rnorm(15,5,0.5)))
dimnames(mydata)<-list(paste("Var",1:25,sep=""),paste("Exp",1:3,sep=""))
mydist<-distancematrix(mydata,d="cosangle") #compute the distance matrix.

#clusters and final tree
clustresult<-hopach(mydata,dmat=mydist)

#bootstrap resampling
myobj<-boothopach(mydata,clustresult)

#write MapleTree files
boot2fuzzy(mydata,myobj,clustresult)

}

\keyword{cluster}
\keyword{nonparametric}
\keyword{multivariate}

