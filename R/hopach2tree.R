#MakeTree: a function to make the .gtr (or .atr) file for genes (or arrays)
# labels = the final labels from hopach(): $final$labels
# ord = the order of the elements (corresponding to labels) in the .cdt file
# medoids = the final medoids matrix from hopach(): $final$medoids
# dist = the corresponding distance matrix
# side = dimension clustered: "ARRY" or "GENE" (default).

makeTree<-function(labels,ord,medoids,dist,side="GENE"){
	if((!is.vector(labels) | !is.numeric(labels)))
		stop("First arg to makeTree() must be a numeric vector")
	if(!is.matrix(medoids))
		stop("Second arg to makeTree() must be a matrix")
	if(length(medoids[1,])!=2)
		stop("Second arg to makeTree() must have 2 columns")
	if(is.matrix(dist))
		dist <- as.hdist(dist) 
	if(!is.hdist(dist))
		stop("Third arg to makeTree() must be a matrix, or a hdist object")
	if(length(labels)!=length(dist[1,]))
		stop("Labels and distance matrix dimensions do not agree in makeTree()")
	labels<-labels[ord]
	dist<-dist[ord,ord]
	maxnode<-max(c(max(table(labels)),9))
	treematrix<-cbind(data.frame(node=NA,label=0,avgdist=0,branch=0,medoid=0),matrix(NA,nc=maxnode,dimnames=list(NULL,paste("child",1:maxnode,sep=""))))
	D<-digits(labels)
	if(D<2)
		stop("Can not run makeTree() with single digit labels")
	arrays<-paste(side,1:length(labels),"X",sep="")[ord]
	cutlab<-cutzeros(medoids[,1])
	ucutlab<-unique(cutlab)
	cutmedoids<-data.frame(label=0,medoid=0)
	for(i in 1:length(ucutlab)){
		set<-medoids[cutlab==ucutlab[i],]
		if(sum(cutlab==ucutlab[i])==1)
			cutmedoids<-rbind(cutmedoids,set)
		else{
			digset<-apply(as.matrix(set[,1]),1,digits)
			maxdig<-max(digset)
			cutmedoids<-rbind(cutmedoids,set[digset==maxdig,])
		}
	}
	medoids<-cutmedoids[-1,]
	rm(cutmedoids)
	medoids<-cbind(medoids,apply(as.matrix(medoids[,1]),1,digits))
	nodecount<-1
	cat("Working on level",D,"\n")
	block<-medoids[medoids[,3]==D,]
	block<-block[order(block[,1]),]
	nd<-sum(medoids[,3]==D)
	if(nd==1)
		block<-t(matrix(block))
	if(nd>0){
		for(i in 1:nd){
			ni<-sum(labels==block[i,1])
			if(ni>1){
				#avgdist<-mean(dissvector(as.matrix(dist[labels==block[i,1],labels==block[i,1]])))
				avgdist<-mean(as.vector(dist[labels==block[i,1],labels==block[i,1]]))
				treematrix<-rbind(treematrix,c(paste("NODE",nodecount,"X",sep=""),cutzeros(block[i,1]),avgdist,2*block[i,3]/(D+1)-1,block[i,2],arrays[labels==block[i,1]],rep(NA,maxnode-ni)))
				nodecount<-nodecount+1
			}
		}
	}
	labels.d1<-labels
	labels.d<-cutdigits(labels,D-1)
	for(d in (D-1):1){
		cat("Working on level",d,"\n")
		block<-medoids[medoids[,3]==d,]
		block<-block[order(block[,1]),]
		nd<-sum(medoids[,3]==d)
		if(nd==1)
			block<-t(matrix(unlist(block)))
		if(nd>0){
			for(i in 1:nd){
				uchildren<-sort(unique(labels.d1[labels.d==block[i,1]]))
				if(length(uchildren)>1){
					children<-NULL
					for(j in 1:length(uchildren)){
						if(sum(treematrix$label==uchildren[j]))
							children<-c(children,treematrix$node[treematrix$label==uchildren[j]])
						else
							children<-c(children,arrays[labels.d1==uchildren[j]])
					}
					ni<-length(children)
					if(ni>1){
						#avgdist<-mean(dissvector(as.matrix(dist[labels.d==block[i,1],labels.d==block[i,1]])))
						avgdist<-mean(as.vector(dist[labels.d==block[i,1],labels.d==block[i,1]]))
						treematrix<-rbind(treematrix,c(paste("NODE",nodecount,"X",sep=""),cutzeros(block[i,1]),avgdist,2*block[i,3]/(D+1)-1,block[i,2],children,rep(NA,maxnode-ni)))
						nodecount<-nodecount+1
					}
				}
				else{
					treematrix$label[treematrix$label==uchildren]<-cutzeros(block[i,1])
				}
			}
		}
		labels.d1<-labels.d
		labels.d<-cutdigits(labels.d,d-1)
	}
	uchildren<-sort(unique(cutdigits(labels,1)))
	children<-NULL
	for(j in 1:length(uchildren)){
		if(sum(treematrix$label==uchildren[j]))
			children<-c(children,treematrix$node[treematrix$label==uchildren[j]])
		else
			children<-c(children,arrays[labels.d1==uchildren[j]])
	}
	ni<-length(children)
	treematrix<-rbind(treematrix,c(paste("NODE",nodecount,"X",sep=""),0,mean(as.vector(dist)),-1,order(rowSums(as.matrix(dist)))[1],children,rep(NA,maxnode-ni)))
	treematrix$medoid<-paste(side,treematrix$medoid,"X",sep="")
	return(treematrix[-1,-2])
}

#hopach2tree: make TreeView/MapleTree files from HOPACH output

hopach2tree<-function(data,file="HOPACH",hopach.genes=NULL,hopach.arrays=NULL,dist.genes=NULL,dist.arrays=NULL,d.genes="cosangle",d.arrays="euclid",gene.wts=NULL,array.wts=NULL,gene.names=NULL){
	olddig<-options("digits")$digits
	options(digits=10)
	if(inherits(data,"ExpressionSet")) 
		data<-exprs(data)
	data<-as.matrix(data)
	if(!is.matrix(data))
		stop("First arg to hoapch2tree() must be a matrix")
	p<-length(data[,1])
	n<-length(data[1,])
	if(is.null(hopach.genes)){
		row.labels<-1:p
		row.medoids<-matrix(rep(row.labels,2),nc=2,dimnames=list(NULL,c("label","medoid")))
		row.ord<-1:p
	}
	else{
		row.labels<-hopach.genes$final$labels
		row.medoids<-hopach.genes$final$medoids
		row.ord<-hopach.genes$final$order
	}
	if(is.null(hopach.arrays)){
		col.labels<-1:n
		col.medoids<-matrix(rep(col.labels,2),nc=2,dimnames=list(NULL,c("label","medoid")))
		col.ord<-1:n
	}
	else{
		col.labels<-hopach.arrays$final$labels
		col.medoids<-hopach.arrays$final$medoids
		col.ord<-hopach.arrays$final$order
	}
	if(is.null(gene.wts))
		gene.wts<-rep(1,p)
	if(is.null(array.wts))
		array.wts<-rep(1,n)
	if(is.null(gene.names))
		gene.names<-dimnames(data.frame(data))[[1]]
	gene.names<-as.character(gene.names)
	if(length(array.wts)!=n)
		stop("Array weights and data matrix dimensions do not match in hopach2tree()")
	if(length(gene.wts)!=p)
		stop("Gene weights and data matrix dimensions do not match in hopach2tree()")
	if(length(row.labels)!=p)
		stop("Gene labels and data matrix dimensions do not agree in hopach2tree()")
	if(length(col.labels)!=n)
		stop("Array labels and data matrix dimensions do not agree in hopach2tree()")
	if(length(row.ord)!=p)
		stop("Row ordering and data matrix dimensions do not agree in hopach2tree()")
	if(length(col.ord)!=n)
		stop("Column ordering and data matrix dimensions do not agree in hopach2tree()")
	cat("GENE TREE:\n")
	if(digits(row.labels[1])>1){
		if(is.null(dist.genes))
			dist.genes<-distancematrix(data,d=d.genes)
		if(is.hdist(dist.genes))
			dist.genes <- as.matrix(dist.genes) 
		if(length(dist.genes[1,])!=p)
			stop("Data and gene distance matrix dimensions do not agree in hopach2tree()")
		if(length(dist.genes[,1])!=p)
			stop("Gene distance matrix must be square in hopach2tree()")
		gene.result<-makeTree(row.labels,row.ord,row.medoids,dist.genes,"GENE")
		nr<-nrow(gene.result)
		write.table(paste("DistanceMetric=",d.genes,sep=""),file=paste(file,".gtr",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
		for(i in 1:nr){
			write.table(t(gene.result[i,][!is.na(gene.result[i,])]),paste(file,".gtr",sep=""),append=TRUE,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
		}
		cdt<-c("GID","UID","NAME","GWEIGHT",dimnames(data.frame(data))[[2]][col.ord])
	}
	else{
		cat("Only one level in gene tree: no .gtr file created \n")
		cdt<-c("UID","NAME","GWEIGHT",dimnames(data.frame(data))[[2]][col.ord])
	}
	cdtsize<-1
	cat("ARRAY TREE:\n")
	if(digits(col.labels[1])>1){
		if(is.null(dist.arrays))
			dist.arrays<-distancematrix(t(data),d=d.arrays)
		if(is.hdist(dist.arrays))
			dist.arrays <- as.matrix(dist.arrays) 
		if(length(dist.arrays[1,])!=n)
			stop("Data and array distance matrix dimensions do not agree in hopach2tree()")
		if(length(dist.arrays[,1])!=n)
			stop("Array distance matrix must be square in hopach2tree()")
		sample.result<-makeTree(col.labels,col.ord,col.medoids,dist.arrays,"ARRY")[,-4]
		nr<-nrow(sample.result)
		write.table(paste("DistanceMetric=",d.arrays,sep=""),file=paste(file,".atr",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
		for(i in 1:nr){
			write.table(t(sample.result[i,][!is.na(sample.result[i,])]),paste(file,".atr",sep=""),append=TRUE,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
		}
		if(digits(row.labels[1])>1)
			cdt<-rbind(cdt,c("AID","","","",paste("ARRY",1:n,"X",sep="")[col.ord]))
		else
			cdt<-rbind(cdt,c("AID","","",paste("ARRY",1:n,"X",sep="")[col.ord]))
		cdtsize<-cdtsize+1
	}
	else
		cat("Only one level in array tree: no .atr file created \n")
	if(digits(row.labels[1])>1)
		write.table(matrix(rbind(cdt,c("EWEIGHT","","","",array.wts[col.ord]),cbind(paste("GENE",1:p,"X",sep="")[row.ord],dimnames(data.frame(data))[[1]][row.ord],gene.names[row.ord],gene.wts[row.ord],data[row.ord,col.ord])),nr=(cdtsize+p+1),dimnames=list(NULL,NULL)),paste(file,".cdt",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
	else
		write.table(matrix(rbind(cdt,c("EWEIGHT","","",array.wts[col.ord]),cbind(dimnames(data.frame(data))[[1]][row.ord],gene.names[row.ord],gene.wts[row.ord],data[row.ord,col.ord])),nr=(cdtsize+p+1),dimnames=list(NULL,NULL)),paste(file,".cdt",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
	options(digits=olddig)
}
