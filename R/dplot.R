dplot<-function(dist,hopachobj,ord="final",col=heat.colors(12),main=NULL,xlab=NULL,ylab=NULL,labels=NULL,showclusters=TRUE,...){
	dist<-as.matrix(dist)
	p<-nrow(dist)
# Don't hard code main
#	if(ord!="none")
#		main<-paste(main,"Ordered Distance Matrix",sep="\n")
	if(is.null(ylab))
		ylab=""
	if(is.null(xlab))
		xlab=""
	if(showclusters)
		boundary<-cumsum(hopachobj$clustering$size)+0.5
	ordering<-1:p
	if(ord=="final")
		ordering<-hopachobj$final$ord
	if(ord=="cluster")
			ordering<-hopachobj$clustering$ord
	distplot<-dist[ordering,ordering]
	diag(distplot)<-min(dist[upper.tri(dist)])
	par(mfrow=c(1,1),pty="s")
	image(1:p,1:p,distplot[,p:1],main=main,xlab=xlab,ylab=ylab,axes=FALSE,col=col,...)
	if(!is.null(labels)){
		labels<-as.character(labels)
		axis(2,labels=rev(labels[ordering]),at=1:p,cex.axis=0.75,col.axis=2,las=2,...)
		axis(1,labels=labels[ordering],at=1:p,cex.axis=0.75,col.axis=2,las=2,...)
	}
	box()
	if(showclusters)
		abline(v=boundary,h=p+1-boundary,lty=2)
}
