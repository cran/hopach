bootplot<-function(bootobj,hopachobj,ord="bootp",main=NULL,labels=NULL,showclusters=TRUE,...){
	p<-nrow(bootobj)
	k<-ncol(bootobj)
	shownames<-(p<50)
	ordering<-hopachobj$clustering$ord
	if(ord=="bootp"){
		start<-1
		stop<-hopachobj$clust$sizes[1]
		set<-ordering[start:stop]
		#ordering[start:stop]<-set[rev(order(bootobj[set,1]))]
		ordering[start:stop]<-set[order(bootobj[set,1], decreasing = TRUE)]
		for(i in 2:hopachobj$clust$k){
			start<-stop+1
			stop<-cumsum(hopachobj$clust$sizes)[i]
			set<-ordering[start:stop]
			ordering[start:stop]<-set[order(bootobj[set,i], decreasing = TRUE )]	
		}
	}
	if(ord=="final")
		ordering<-hopachobj$final$ord
	if(ord=="none"){
		ordering<-1:p
		showclusters=FALSE
	}
	bootobj<-bootobj[ordering,]
	colors<-rainbow(k)
	colors<-c(colors[seq(1,k,by=2)],colors[seq(2,k,by=2)])
# don't hard code main
#	main<-paste(main,"Barplot of Bootstrap Reappearance Proportions",sep="\n")
	if(is.null(labels))
		labels<-dimnames(bootobj)[[1]]
	par(oma=c(0,0,0,2))
	barplot(t(bootobj),ylim=c(1,p),border=FALSE,space=0,horiz=TRUE,names.arg=labels[ordering],las=1,main=main,cex.names=0.75,legend.text=FALSE,col=colors,axisnames=shownames,xlab="Proportion",...)
	if(showclusters){
		abline(h=cumsum(hopachobj$clust$sizes))
		mtext(colnames(bootobj),outer=TRUE,side=4,at=cumsum(hopachobj$clust$sizes)/p*0.7+0.16,line=-2,col=colors,las=1,cex=0.6)
	}
}
