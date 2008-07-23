makeoutput<-function(data,hopachobj,bootobj=NULL,file="HOPACH.out",gene.names=NULL){
	olddig<-options("digits")$digits
	options(digits=16)
	if(inherits(data,"ExpressionSet")) 
		data<-exprs(data)
	data<-as.matrix(data)
	p<-nrow(data)
	out<-row.names(data)
	if(is.null(out))
		out<-1:p
	out<-data.frame(UID=out)
	if(!is.null(gene.names)){
		gene.names<-as.vector(gene.names)
		if(length(gene.names)!=p)
			stop("Row names and data matrix dimensions do not agree in makeoutput().")
		out<-data.frame(out,Name=gene.names)
	}
	out<-data.frame(Index=(1:p),out)	
	uclust<-sort(unique(hopachobj$clust$labels))
	clust<-NULL
	for(i in 1:length(uclust)){
		clust[hopachobj$clust$label==uclust[i]]<-(i-1)
	}
	out<-data.frame(out,Cluster.Number=clust)
	out<-data.frame(out,Cluster.Label=hopachobj$clustering$labels)
	out<-out[hopachobj$clustering$order,]
	out<-data.frame(out,Cluster.Level.Order=(1:p))
	out<-out[order(out[,1]),]
	out<-data.frame(out,Final.Label=hopachobj$final$label)
	out<-out[hopachobj$final$order,]
	out<-data.frame(out,Final.Level.Order=1:p)
	out<-out[order(out[,1]),]
	if(!is.null(bootobj)){
		bootobj<-as.matrix(bootobj)
		if(dim(bootobj)[1]!=p)
			stop("Bootstrap and data matrix dimensions do not agree in makeoutput().")
		k<-ncol(bootobj)
		if(length(hopachobj$clust$medoids)!=k)
			stop("Number of clusters in hopach and bootstrap do not match in boot2fuzzy()")
		colnames(bootobj)<-paste("Cluster",0:(k-1),".Membership",sep="")
		out<-data.frame(out,bootobj)
	}
	write.table(out,file,sep="\t",row.names=FALSE)
	options(digits=olddig)
}

