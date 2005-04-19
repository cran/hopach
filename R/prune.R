#function to prune a hopach tree to a given level
# can be lower or higher than main clusters
prune<-function(data,hobj,level,dmat=NULL,ord="own"){
	data<-as.matrix(data)
	if(round(level)!=level)
		stop("level must be an integer")
	if(level<1)
		stop("level must be at least 1")
	if(is.null(dmat))
		dmat<-distancematrix(data,d=hobj$metric)
	final<-digits(hobj$final$labels[1])
	main<-digits(hobj$clustering$labels[1])
	if(level>final)
		stop(paste("Tree has",final,"levels, less than",level,"- pick a smaller value for level",sep=" "))
	if(level<main){
		cat("Pruning tree to a level above the main clusters...\n")
		labs<-trunc(hobj$clustering$labels/(10^(main-level)))
		ulabs<-sort(unique(labs))
		meds<-hobj$final$medoids[hobj$final$medoids[,1]%in%ulabs,2]
		hobj$clustering$k<-length(ulabs)
		hobj$clustering$medoids<-meds
		hobj$clustering$sizes<-as.vector(table(labs))
		hobj$clustering$labels<-labs
		cuttree<-list(hobj$clustering$k,meds,hobj$clustering$sizes,labs,0,hobj$final$medoids)
		hobj$clustering$order<-orderelements(cuttree,data,rel=ord,d=hobj$metric,dmat=dmat)[[2]]
	}
	if(level>main){
		cat("Pruning tree to a level below the main clusters...\n")
		labs<-trunc(hobj$final$labels/(10^(final-level)))
		ulabs<-sort(unique(labs))
		meds<-hobj$final$medoids[hobj$final$medoids[,1]%in%ulabs,2]
		hobj$clustering$k<-length(ulabs)
		hobj$clustering$medoids<-meds
		hobj$clustering$sizes<-as.vector(table(labs))
		hobj$clustering$labels<-labs
		cuttree<-list(hobj$clustering$k,meds,hobj$clustering$sizes,labs,0,hobj$final$medoids)
		hobj$clustering$order<-orderelements(cuttree,data,rel=ord,d=hobj$metric,dmat=dmat)[[2]]	
	}
	hobj
}
