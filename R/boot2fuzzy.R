boot2fuzzy<-function(data,bootobj,hopach.genes,hopach.arrays=NULL,file="hopach",clust.wts=NULL,gene.wts=NULL,array.wts=NULL,gene.names=NULL){
	olddig<-options("digits")$digits
	options(digits=10)
	if(inherits(data,"ExpressionSet")) 
		data<-exprs(data)
	data<-as.matrix(data)
	p<-length(data[,1])
	n<-length(data[1,])
	if(dim(bootobj)[1]!=p)
			stop("Bootstrap and data matrix dimensions do not agree in boot2fuzzy()")
	if(is.null(dimnames(data)[[1]]))
		r.names<-paste("row",1:p,sep="")
	else
		r.names<-dimnames(data)[[1]]
	if(is.null(gene.names))
		gene.names<-r.names
	gene.names<-as.character(gene.names)
	if(length(gene.names)!=p)
		stop("Gene names and data dimensions do not match in boot2fuzzy()")
	if(is.null(dimnames(data)[[2]]))
		c.names<-paste("col",1:n,sep="")
	else
		c.names<-dimnames(data)[[2]]
	medoids<-hopach.genes$clust$medoids
	row.ord<-hopach.genes$final$ord
	if(is.null(hopach.arrays))
		col.ord<-1:n
	else
		col.ord<-hopach.arrays$final$ord
	k<-length(medoids)
	if(dim(bootobj)[2]!=k)
		stop("Number of clusters in hopach and bootstrap do not match in boot2fuzzy()")
	if(is.null(clust.wts))
		clust.wts<-rep(1,k)
	if(is.null(gene.wts))
		gene.wts<-rep(1,p)
	if(is.null(array.wts))
		array.wts<-rep(1,n)
	if(length(array.wts)!=n)
		stop("Array weights and data matrix dimensions do not match in boot2fuzzy()")
	if(length(gene.wts)!=p)
		stop("Gene weights and data matrix dimensions do not match in boot2fuzzy()")
	if(length(clust.wts)!=k)
		stop("Cluster weights and number of clusters do not match in boot2fuzzy()")
	if(length(col.ord)!=n)
		stop("Column ordering and data matrix dimensions do not match in boot2fuzzy()")
	if(length(row.ord)!=p)
		stop("Row ordering and data matrix dimensions do not match in boot2fuzzy()")
	write.table(cbind(0:(k-1),paste("Centroid",0:(k-1)),clust.wts,data[medoids,col.ord]),paste(file,"fct",sep="."),quote=FALSE,sep="\t",row.names=FALSE,col.names=c("UID","NAME","GWEIGHT",c.names[col.ord]))
	write.table(cbind(r.names[row.ord],gene.names[row.ord],bootobj[row.ord,]),paste(file,"mb",sep="."),quote=FALSE,sep="\t",row.names=FALSE,col.names=c("UID","NAME",paste(0:(k-1))))
	write.table(rbind(c("UID","NAME","GWEIGHT",c.names[col.ord]),c("EWEIGHT","","",array.wts[col.ord]),cbind(r.names[row.ord],gene.names[row.ord],gene.wts[row.ord],data[row.ord,col.ord])),paste(file,"cdt",sep="."),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
	options(digits=olddig)
}

