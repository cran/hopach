## FIXME:  remove arg I after version 1.2
bootmedoids<-function(data,medoids,d="cosangle",B=1000,I){
	if(inherits(data,"exprSet")) 
		data<-exprs(data)
	data<-as.matrix(data)
	p<-length(data[,1])
	n<-length(data[1,])
	k<-length(medoids)
	if(!missing(I)){
         	if(!missing(B))
          		stop("You can only specify one of I and B")
      	 	else{
        		warning("Argument I is deprecated, please use argument B")
        		B<-I
		}
	}
	blabs<-matrix(0,nrow=p,ncol=B) #holds the bootstrap cluster labels
	bdist<-matrix(0,nrow=p,ncol=k) #holds distance to medoids (recycled)
	props<-matrix(0,nrow=p,ncol=k) #holds the cluster probabilities
	for(b in 1:B){
		samp<-sample(1:n,replace=TRUE)
		for(j in 1:k){
			bdist[,j]<-distancevector(data[,samp],data[medoids[j],samp],d)
		}
       	blabs[,b]<-apply(bdist,1,order)[1,]
	}
	for(i in 1:k)
		props[,i]<-rowMeans(blabs==i,na.rm=TRUE)
	props[medoids,]<-diag(k)
	dimnames(props)<-list(dimnames(data)[[1]],paste("Cluster",0:(k-1),sep=""))
	return(props)
}

boothopach<-function(data,hopachobj,B=1000,I,hopachlabels=FALSE){
	if(inherits(data,"exprSet")) 
		data<-exprs(data)
	data<-as.matrix(data)
	p<-length(data[,1])
	n<-length(data[1,])
	medoids<-hopachobj$clust$medoids
	d<-hopachobj$metric
	k<-length(medoids)
	if(!missing(I)){
         	if(!missing(B))
          		stop("You can only specify one of I and B")
      	 	else{
        		warning("Argument I is deprecated, please use argument B")
        		B<-I
		}
	}
	blabs<-matrix(0,nrow=p,ncol=B) #holds the bootstrap cluster labels
	bdist<-matrix(0,nrow=p,ncol=k) #holds distance to medoids (recycled)
	props<-matrix(0,nrow=p,ncol=k) #holds the cluster probabilities
	for(b in 1:B){
		samp<-sample(1:n,replace=TRUE)
		for(j in 1:k){
			bdist[,j]<-distancevector(data[,samp],data[medoids[j],samp],d)
		}
       	blabs[,b]<-apply(bdist,1,order)[1,]
	}
	for(i in 1:k)
		props[,i]<-rowMeans(blabs==i,na.rm=TRUE)
	props[medoids,]<-diag(k)
	if(hopachlabels)
		labs<-as.character(hopachobj$clust$lab[medoids])
	else
		labs<-paste("Cluster",0:(k-1),sep="")
	dimnames(props)<-list(dimnames(data)[[1]],labs)
	return(props)
}
