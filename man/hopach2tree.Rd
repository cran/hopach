\name{hopach2tree}

\alias{hopach2tree}
\alias{makeTree}

\title{function to write MapleTree files for viewing hopach hierarchical clustering results}

\description{
The MapleTree software (http://mapletree.sourceforge.net/) is an open source, cross-platform, visualization tool to graphically browse results of cluster analyses. The \code{hopach2tree} function takes a data matrix, plus corresponding \code{hopach} clustering output for genes and/or arrays, and writes the (.cdt, .gtr, and .atr) files needed to view these hierarchical clustering results in MapleTree. The function \code{makeTree} is called internally by \code{hopach2tree}.
}

\usage{
hopach2tree(data, file = "HOPACH", hopach.genes = NULL, hopach.arrays = NULL, 
dist.genes = NULL, dist.arrays = NULL, d.genes = "cosangle", 
d.arrays = "euclid", gene.wts = NULL, array.wts = NULL, gene.names = NULL)

makeTree(labels, ord, medoids, dist, side = "GENE")

}

\arguments{
  \item{data}{data matrix, data frame or exprSet of gene expression measurements. Each column corresponds to an array, and each row corresponds to a gene. All values must be numeric. Missing values are ignored.}
  \item{file}{name for the output files (the extensions .cdt, .gtr and .atr will be added).}
  \item{hopach.genes}{output of the \code{hopach} function applied to genes (rows of \code{data}. If only arrays are clustered, hopach.genes can be NULL. There must be at least K=2 levels in the hopach final tree (ie: hopach.genes$final$labels greater than 1 digit) for a gtr file to be generated.}
  \item{hopach.arrays}{optional output of the \code{hopach} function applied to arrays (columns of \code{data}. There must be at least K=2 levels in the hopach final tree (ie: hopach.arrays$final$labels greater than 1 digit) for an atr file to be generated.}
  \item{dist.genes}{matrix of pair wise distances between all genes. All values must be	numeric, and missing values are not allowed. If NULL, this matrix is computed using the metric specified by \code{d.genes}. Only needed if genes are clustered (hopach.genes!=NULL).}
  \item{dist.arrays}{matrix of pair wise distances between all arrays. All values must be numeric, and missing values are not allowed. If NULL, this matrix is computed using the metric specified by \code{d.arrays}. Only needed if arrays are clustered (hopach.arrays!=NULL).}
  \item{d.genes}{character string specifying the metric to be used for calculating dissimilarities between genes. The currently available options are "cosangle" (cosine angle or uncentered correlation distance), "abscosangle" (absolute cosine angle or absolute uncentered correlation distance), "euclid" (Euclidean distance), "abseuclid" (absolute Euclidean distance), "cor" (correlation distance), and "abscor" (absolute correlation distance). Advanced users can write their own distance functions and add these to the functions \code{distancematrix()} and \code{distancevector()}.}
  \item{d.arrays}{character string specifying the metric to be used for calculating dissimilarities between arrays.}
  \item{gene.wts}{an optional vector of numeric weights for the genes.}
  \item{array.wts}{an optional vector of numeric weights for the arrays.}
  \item{gene.names}{optional vector of names or annotations for the genes, which can be different from the row names of \code{data}.}
  \item{labels}{final cluster labels from a hopach object.}
  \item{ord}{final ordering from a hopach object.}
  \item{medoids}{final medoids matrix from a hopach object.}
  \item{dist}{gene or array distance matrix.}
  \item{side}{character string specifying if the tree is for genes ("GENE", default) or arrays ("ARRY").}
}

\value{
The function \code{hopach2tree} has no value. It writes up to three text files to the current working directory. A .cdt file is always produced. This file can be used to visualize the data matrix as a heat map in MapleTree or other viewers such as TreeView (http://rana.lbl.gov/EisenSoftware.htm), jtreeview (http://sourceforge.net/projects/jtreeview/), and GeneXPress (http://genexpress.stanford.edu/). When \code{hopach.genes}!=NULL, a .gtr is produced, and gene clustering results can be viewed, including ordering the genes in the heat map according to the final level of the \code{hopach} tree and drawing the dendogram for hierarchical gene clustering. Similarly, when \code{hopach.arrays}!=NULL, an .atr file is produced and array clustering results can be viewed.

The function \code{makeTree} is called internally by \code{hopach2tree} to make the objects needed to write the MapleTree files for a gene and/or array HOAPCH clustering result.
}

\references{

van der Laan, M.J. and Pollard, K.S. A new algorithm for hybrid hierarchical clustering with visualization and the bootstrap. Journal of Statistical Planning and Inference, 2003, 117, pp. 275-303.

\url{http://www.stat.berkeley.edu/~laan/Research/Research_subpages/Papers/hopach.pdf}
}

\author{Katherine S. Pollard <kpollard@soe.ucsc.edu>}

\note{Thank you to Lisa Simirenko <lsimirenko@lbl.gov> for providing HOPACH views in MapleTree, and to Karen Vranizan <vranizan@uclink.berkeley.edu> for her input. 

The MapleTree software can be downloaded from: http://sourceforge.net/projects/mapletree/}

\section{Warning }{Operating systems use different end of line characters. These characters can cause errors in MapleTree when files generated on one OS are visualized on another OS. Hence, \code{hopach2tree} should be run on the same OS as MapleTree whenever possible.}

\seealso{\code{\link{hopach}}, \code{\link{boothopach}}, \code{\link{bootmedoids}}, \code{\link{boot2fuzzy}}}

\examples{

#25 variables from two groups with 3 observations per variable
mydata<-rbind(cbind(rnorm(10,0,0.5),rnorm(10,0,0.5),rnorm(10,0,0.5)),cbind(rnorm(15,5,0.5),rnorm(15,5,0.5),rnorm(15,5,0.5)))
dimnames(mydata)<-list(paste("Var",1:25,sep=""),paste("Exp",1:3,sep=""))
mydist<-distancematrix(mydata,d="cosangle") #compute the distance matrix.

#clusters and final tree
clustresult<-hopach(mydata,dmat=mydist)

#write MapleTree files
hopach2tree(mydata,hopach.genes=clustresult,dist.genes=mydist)

}

\keyword{cluster}
\keyword{multivariate}
