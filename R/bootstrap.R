bootmedoids<-function(data,medoids,d="cosangle",I=1000){
	if(inherits(data,"exprSet")) 
		data<-exprs(data)
	data<-as.matrix(data)
	p<-length(data[,1])
	n<-length(data[1,])
	k<-length(medoids)
	blabs<-matrix(0,nrow=p,ncol=I) #holds the bootstrap cluster labels
	bdist<-matrix(0,nrow=p,ncol=k) #holds distance to medoids (recycled)
	props<-matrix(0,nrow=p,ncol=k) #holds the cluster probabilities
	for(i in 1:I){
		samp<-sample(1:n,replace=TRUE)
		for(j in 1:k){
			bdist[,j]<-distancevector(data[,samp],data[medoids[j],samp],d)
		}
       	blabs[,i]<-apply(bdist,1,order)[1,]
	}
	for(i in 1:k)
		props[,i]<-apply(blabs==i,1,mean,na.rm=TRUE)
	props[medoids,]<-diag(k)
	dimnames(props)<-list(dimnames(data)[[1]],paste("Cluster",0:(k-1),sep=""))
	return(props)
}

boothopach<-function(data,hopachobj,I=1000,hopachlabels=FALSE){
	if(inherits(data,"exprSet")) 
		data<-exprs(data)
	data<-as.matrix(data)
	p<-length(data[,1])
	n<-length(data[1,])
	medoids<-hopachobj$clust$medoids
	d<-hopachobj$metric
	k<-length(medoids)
	blabs<-matrix(0,nrow=p,ncol=I) #holds the bootstrap cluster labels
	bdist<-matrix(0,nrow=p,ncol=k) #holds distance to medoids (recycled)
	props<-matrix(0,nrow=p,ncol=k) #holds the cluster probabilities
	for(i in 1:I){
		samp<-sample(1:n,replace=TRUE)
		for(j in 1:k){
			bdist[,j]<-distancevector(data[,samp],data[medoids[j],samp],d)
		}
       	blabs[,i]<-apply(bdist,1,order)[1,]
	}
	for(i in 1:k)
		props[,i]<-apply(blabs==i,1,mean,na.rm=TRUE)
	props[medoids,]<-diag(k)
	if(hopachlabels)
		labs<-as.character(hopachobj$clust$lab[medoids])
	else
		labs<-paste("Cluster",0:(k-1),sep="")
	dimnames(props)<-list(dimnames(data)[[1]],labs)
	return(props)
}
