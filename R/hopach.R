#Hierarchical Ordered Partitioning and Collapsing Hybrid (HOPACH)#

#1. silhouette related calculations#

#a. silhouettes

#from medoids and distance matrix 
medstosil<-function(medoids,dist){
	if(!is.vector(medoids))
		stop("First arg to medstosil() must be a vector")
	if(!is.matrix(dist))
		stop("Second arg to medstosil() must be a matrix")
	p<-length(dist[1,])
	if(p!=length(dist[,1]))
		stop("Second arg to medstosil() not a square matrix")
	k<-length(medoids)
	if(k<2){
		warning("Less than 2 medoids - can't calculate silhouettes!")
		sil<-rep(NA,p)
		clust<-rep(1,p)
	}
	else{
		#get labels
		clust<-apply(dist[,medoids],1,order)[1,]
		#get clussizes (same order as medoids)
		clussize<-table(clust)
		#get sils
		avgdist<-matrix(0,nrow=p,ncol=k)
		a<-b<-NULL
		for(j in (1:k)){
			for(m in (1:k)){
				subdist<-dist[clust==j,clust==m]
				if(clussize[j]>1) 
					avgdist[clust==j,m]<-rowSums(as.matrix(subdist))
				else
					avgdist[clust==j,m]<-sum(subdist)
				if(m==j){
					if(clussize[m]>1) 
						avgdist[clust==j,m]<-avgdist[clust==j,m]/(clussize[m]-1)
					else 
						avgdist[clust==j,m]<-0
				}
				else 
					avgdist[clust==j,m]<-avgdist[clust==j,m]/clussize[m]
			}
			a[clust==j]<-avgdist[clust==j,j]
			if(clussize[j]>1) 
				b[clust==j]<-apply(as.matrix(avgdist[clust==j,-j]),1,min)
			else 
				b[clust==j]<-min(avgdist[clust==j,-j])
		}
		sil<-(b-a)/pmax(a,b)
	}
	list(clust,sil)
}

#from labels and distance matrix
labelstosil<-function(labels,dist){
	if(!is.vector(labels))
		stop("First arg to labelstosil() must be a vector")
	if(!is.matrix(dist))
		stop("Second arg to labelstosil() must be a matrix")
	p<-length(dist[1,])
	if(p!=length(dist[,1]))
		stop("Second arg to labelstosil() not a square matrix")
	if(length(labels)!=p)
		stop("Distance matrix and labels dimensions do not agree in labelstosil()")
	unlabels<-sort(unique(labels))
	k<-length(unlabels)
	if(k<2){
		warning("Only one label - can't calculate silhouettes!")
		sil<-rep(NA,p)
	}
	else{
		clussize<-table(labels)
		#get sils
		avgdist<-matrix(0,nrow=p,ncol=k)
		a<-b<-NULL
		for(j in (1:k)){
			for(m in (1:k)){
				subdist<-dist[labels==unlabels[j],labels==unlabels[m]]
				if(clussize[j]>1) 
					avgdist[labels==unlabels[j],m]<-rowSums(as.matrix(subdist))
				else 
					avgdist[labels==unlabels[j],m]<-sum(subdist)
				if(m==j){
					if(clussize[m]>1) 
						avgdist[labels==unlabels[j],m]<-avgdist[labels==unlabels[j],m]/(clussize[m]-1)
					else 
						avgdist[labels==unlabels[j],m]<-0
				}
				else 
					avgdist[labels==unlabels[j],m]<-avgdist[labels==unlabels[j],m]/clussize[m]
			}
			a[labels==unlabels[j]]<-avgdist[labels==unlabels[j],j]
			if(clussize[j]>1) 
				b[labels==unlabels[j]]<-apply(as.matrix(avgdist[labels==unlabels[j],-j]),1,min)
			else 
				b[labels==unlabels[j]]<-min(avgdist[labels==unlabels[j],-j])
		}
		sil<-(b-a)/pmax(a,b)
	}
	list(labels,sil)
}


#b. mean/median split silhoutte

labelstomss<-function(labels,dist,khigh=9,within="med",between="med",hierarchical=TRUE){
	if(!is.vector(labels))
		stop("First arg to labelstomss() must be a vector")
	if(!is.matrix(dist))
		stop("Second arg to labelstomss() must be a matrix")
	p<-length(dist[1,])
	if(p!=length(dist[,1]))
		stop("Second arg to labelstomss() not a square matrix")
	if(length(labels)!=p)
		stop("Distance matrix and labels dimensions do not agree in labelstomss()")
	unlabels<-sort(unique(labels))
	k<-length(unlabels)
	ss<-NULL
	for(i in 1:k){
		labs<-(1:p)[labels==unlabels[i]]
		pp<-length(labs)
		if(pp<3)
			ss[i]<-NA
		else{
			dissvec<-dissvector(dist[labs,labs])
			bestk<-silcheck(dissvec,min(khigh,(length(labs)-1),na.rm=TRUE),diss=TRUE)[1]
			if(within=="med") 
				ss[i]<-median(pam(dissvec,bestk,diss=TRUE)$silinfo$widths[,3])
			if(within=="mean") 
				ss[i]<-pam(dissvec,bestk,diss=TRUE)$silinfo$avg.width
		}
	}
	if(sum(is.na(ss))==k)
		out<-NA
	else{
		if(hierarchical==TRUE){
			parentlab<-trunc(labels/10)
			unparentlab<-sort(unique(parentlab))
			pk<-length(unparentlab)
			pss<-NULL
			for(i in 1:pk){
				if(between=="med") 
					pss[i]<-median(ss[trunc(unlabels/10)==unparentlab[i]],na.rm=TRUE)
				if(between=="mean") 
					pss[i]<-mean(ss[trunc(unlabels/10)==unparentlab[i]],na.rm=TRUE)
			}
			if(between=="med") 
				out<-median(pss,na.rm=TRUE)
			if(between=="mean") 
				out<-mean(pss,na.rm=TRUE)
		}
		else{
			if(between=="med") 
				out<-median(ss,na.rm=TRUE)
			if(between=="mean") 
				out<-mean(ss,na.rm=TRUE)
		}
	}
	return(out)
}

#c. optimizing number of clusters with average silhouette or mss

#silcheck (silhouettes)
silcheck<-function(data,kmax=9,diss=FALSE,echo=FALSE,graph=FALSE){
	sil<-NULL
	m<-min(kmax,max((!diss)*(dim(data)[1]-1),(diss)*(0.5*(sqrt(1+8*length(data))-1)),na.rm=TRUE))
	if(m<2)
		out<-c(1,NA)
	else{
		for(i in 1:(m-1))
			sil[i]<-pam(data,k=(i+1),diss=diss)$silinfo$avg.width
		if(echo)
			cat("best k = ",order(sil)[length(sil)]+1,", sil(k) = ",round(max(sil),4),"\n")
		if(graph){
			plot(2:m,sil,type="n",xlab="Number of Clusters",ylab="Average Silhouette")
			text(2:m,sil,2:m)
		}
		out<-c(order(sil)[length(sil)]+1,max(sil))
	}
	return(out)
}

#msscheck (mss)
msscheck<-function(dist,kmax=9,khigh=9,within="med",between="med",force=FALSE,echo=FALSE,graph=FALSE){
	if(!is.matrix(dist))
		stop("First arg to msscheck() must be a matrix")
	p<-length(dist[1,])
	if(p!=length(dist[,1]))
		stop("First arg to msscheck() not a square matrix")
	if(p<3)
		out<-c(1,NA)
	else{
		dvec<-dissvector(dist)
		kmax<-min(kmax,p-1,na.rm=TRUE)
		if(force)
			mss<-0
		else
			mss<-labelstomss(rep(1,p),dist,khigh,within,between)
		for(k in 2:kmax)
			mss[k]<-labelstomss(pam(dvec,k,diss=TRUE)$clust,dist,khigh,within,between)
		shift<-0
		if(force){
			mss<-mss[-1]
			shift<-1
		}
		if(echo)
			cat("best k = ",order(mss)[1]+shift,", mss(k) = ",round(min(mss),4),"\n")
		if(graph){
			kmin<-ifelse(force,2,1)
			plot(kmin:kmax,mss,type="n",xlab="Number of Clusters",ylab="MSS")
			text(kmin:kmax,mss,kmin:kmax)
		}
		out<-c(order(mss)[1]+shift,min(mss))
	}
	return(out)
}

#2. Utility functions

#counts the number of digits in label
digits<-function(label){
	label<-label[1]
	count<-0
	while(label>=1){
		count<-count+1
		label<-label/10
	}
	return(count)
}

#truncates labels to dig digits
cutdigits<-function(labels,dig){
	dl<-NULL
	for(i in 1:length(labels)) 
		dl[i]<-digits(labels[i])
	df<-max(0,dl-dig)
	trunc(labels/(10^df))		
}

#removes trailing zeros from labels
cutzeros<-function(labels){
	for(i in 1:length(labels)){
		while(trunc(labels[i]/10)*10==labels[i]){
			labels[i]<-labels[i]/10
		}
	}
	return(labels)
}

#returns the number of non-zero digits in labels
nonzeros<-function(labels){
	for(i in 1:length(labels)){
		while(trunc(labels[i]/10)*10==labels[i]){
			labels[i]<-labels[i]/10
		}
		labels[i]<-digits(labels[i])
	}
	return(labels)	
}

#3. functions for making the tree#

#a. msssplitcluster: splits a cluster# 
	#clust1 is gene by subjects dataframe for cluster1
	#l1 is an integer label of cluster 1: e.g 11,12,13,22, etc describing its path so far in the tree.
	#id1 is id's of cluster1 indicating row numbers in original data frame subdata
	#kmax is the maximum number of groups
	#khigh  is the maximum number of child groups for each group when computing mss
	#medoid1 is the row number (in full data matrix!) indicating  medoid for cluster 1
	#medtodist is the distance from each gene in clust1 to the neighboring cluster medoid,
	#right is 1 if medoid2 is to the right and is 0 if clust1 is the last cluster so that medoid2
	# is the medoid2 to the left of clust1
	#silh is the silhouette of clust1
	#dist1 is the distance matrix for all genes in clust1 
msssplitcluster<-function(clust1,l1,id1,medoid1,med2dist,right,dist1,kmax=9,khigh=9,within="med",between="med"){
	if(!medoid1) 
		warning("Medoid missing - continue to split cluster")
	else{
		if(sum(medoid1==id1)==0 & medoid1) 
			warning("Medoid not in cluster - continue to split cluster")
	}
	if(is.matrix(clust1)){			
		p1<-length(clust1[,1])
		n<-length(clust1[1,])
	}
	else p1<-1
	if(p1<3) 
		k1<-1				
	else{
		l<-length(clust1[,1])
		dissvec<-dissvector(dist1)
		kmax<-min(p1-1,kmax,na.rm=TRUE)
		khigh<-min(p1-1,khigh,na.rm=TRUE)
		k1<-msscheck(dist1,kmax,khigh,within,between)[1]
		if(k1>1){
			pamobj<-pam(dissvec,k1,diss=TRUE)
			newclussizes<-pamobj$clusinfo[,1]
			newmedoids1<-id1[pamobj$medoids]  
			newlabels1<-pamobj$clustering
			distnewmedoids<-NULL 
			for(j in (1:k1)) 
				distnewmedoids[j]<-mean(med2dist[newlabels1==newlabels1[pamobj$medoids[j]]])  
			if(right==1) 
				ord<-rev(order(distnewmedoids))
			else 
				ord<-order(distnewmedoids)
			newmedoids1<-newmedoids1[ord]
			newclussizes<-newclussizes[ord]
			oldlab<-newlabels1
			for(j in (1:k1)) 
				newlabels1[oldlab==ord[j]]<-j
			newlabels1<-rep(10*l1,l)+newlabels1
		}
	}
	if(k1==1){
		newmedoids1<-medoid1
		newlabels1<-rep(10*l1,p1)
		newclussizes<-p1
	}
	for(a in (1:length(newmedoids1))){
		if(sum(newmedoids1[a]==id1)==0) 
			warning("Problem with new medoids after splitting cluster")
	}
	list(k1,newmedoids1,newlabels1,newclussizes)
}

#b. mssnextlevel: calls mssspltitcluster to produce the next level of the tree#
	#data is the data frame
	#prevlevel is the previous level of the tree
	#dmat is the distance matrix
	#kmax is the maximum number of groups
	#khigh is the maximum number of child groups for each group when computing mss
	#within and between are either "med" for median split silhouette or "mean"
	# for mean split silhouette
mssnextlevel<-function(data,prevlevel,dmat,kmax=9,khigh=9,within="med",between="med"){
	if(!is.matrix(data))
		stop("Frist arg to mssnextlevel() must be a matrix")
	if(!is.matrix(dmat))
		stop("Third arg to mssnextlevel() must be a matrix")
	n<-length(data[1,])
	p<-length(data[,1])
	if(length(dmat[1,])!=p)
		stop("Data and distance matrix dimensions do not agree in mssnextlevel()")
	if(length(dmat[,1])!=p)
		stop("Third arg to mssnextlevel() is not a square matrix")
	id<-1:p
	k<-prevlevel[[1]]
	medoids<-prevlevel[[2]]
	labels<-prevlevel[[4]]
	newk<-0
	newlabels<-newmedoids<-newclussizes<-NULL
	count<-1
	ordlabels<-sort(unique(labels))
	if(length(ordlabels)!=k) 
		warning("Number of unique labels not equal number of clusters in mssnextlevel()")
	if(sum(is.na(medoids))){
		warning("Missing values in medoid vector in nextlevel()")
		medoids[is.na(medoids)]<-FALSE
	}
	if(length(unique(medoids))<k && sum(medoids)) 
		warning("Medoids not unique in mssnextlevel()") 
	checkmeans<-FALSE
	if(length(medoids)==1 && !medoids){
		warning("No medoids provided in mssnxtlevel()")
		usemean<-TRUE
	}
	else{
		if(sum(medoids>1)==k) 
			usemean<-FALSE
		else 
			checkmeans<-TRUE
	}
	for(j in (1:k)){
		clust1<-data[labels==ordlabels[j],]	
		id1<-id[labels==ordlabels[j]]
		if(length(id1)>1) 
			clust1<-as.matrix(clust1)
		l1<-ordlabels[j]
		right<-(j<k)
		medoid1<-ifelse(is.na(medoids[j]),0,medoids[j])
		if (j<k) 
			medoid2<-medoids[j+1]
		else 
			medoid2<-medoids[j-1]
		if(length(id1)>1) 
			med2dist<-rowMeans(as.matrix(dmat[labels==ordlabels[j],labels==labels[medoid2]]))
		else 
			med2dist<-mean(dmat[labels==ordlabels[j],labels==labels[medoid2]])
		splitobj<-msssplitcluster(clust1,l1,id1,medoid1,med2dist,right,dmat[labels==l1,labels==l1],kmax,khigh,within,between)
		newlabels[labels==ordlabels[j]]<-splitobj[[3]]
		k1<-splitobj[[1]]
		start<-count
		end<-count+k1-1
		newmedoids[start:end]<-splitobj[[2]]
		newclussizes[start:end]<-splitobj[[4]]
		count<-count+k1		
	}
	count<-newk<-count-1
	newmedoids<-newmedoids[1:newk]
	newclussizes<-newclussizes[1:newk]
	final<-0
	if(count==k)
		final<-1
	if(max(newclussizes)==3) 
		final<-1
	list(newk,newmedoids,newclussizes,newlabels,final,rbind(prevlevel[[6]],cbind(sort(unique(newlabels)),newmedoids)))
}

#c. orderelements: produces an ordering of elements within a set of clusters 
	#level is a level of the tree
	#data is the data frame
	#rel is an indicator of whether to order elements in each cluster with respect 
	# to their own medoid ("own") or the neighboring medoid to the right ("neighbor")
	# or using improveordering() function ("co"). the default is "own"
	#d is an indicator of which distance function to use
	# choices are: "cosangle" (default),"abscosangle","euclid","abseuclid","cor","abscor".
	#dmat is the distance matrix. if this has already been calculated by the user, it can
	# be passed into the function in order to save calculation time
orderelements<-function(level,data,rel="own",d="cosangle",dmat=NULL){
	if(!is.matrix(data))
		stop("Second arg to orderelements() must be a matrix")
	idn<-1:length(data[,1])
	k<-level[[1]]
	labels<-level[[4]]
	medoids<-level[[2]]
	clussizes<-level[[3]]
	ord<-order(labels)
	idnord<-idn[ord]
	subdataord<-data[ord,]
	if(is.matrix(dmat))
		distord<-dmat[ord,]
	labelsord<-labels[ord]
	count<-1
	for(j in (1:k)){
		start<-count
		end<-count+clussizes[j]-1
		if(clussizes[j]>2){
			tempid<-idnord[start:end]
			if(rel=="co"){
				if(is.matrix(dmat)) 
					distj<-distord[,ord][start:end,start:end]
				else 
					distj<-distancematrix(subdataord[start:end,],d)
				idnord[start:end]<-tempid[improveordering(distj)]
			}
			else{
				if(rel=="neighbor"){
					if(j<k) 
						mednext<-medoids[j+1]
					else 
						mednext<-medoids[j-1]
				}
				else 
					mednext<-medoids[j]
				if(is.matrix(dmat)) 
					dmednext<-distord[start:end,mednext]	
				else 
					dmednext<-distancevector(subdataord[start:end,],as.vector(data[mednext,]),d)
				if(rel=="neighbor"){
					if(j<k) 
						ordtemp<-rev(order(dmednext))
					else 
						ordtemp<-order(dmednext)
				}
				else 
					ordtemp<-order(dmednext)
				idnord[start:end]<-tempid[ordtemp]
			}
		}
		else 
			idnord[start:end]<-idnord[start:end]
		count<-count+clussizes[j]
	}
	list(data[idnord,],idnord)
}

#d. mssinitlevel: creates ordered initial level#
	#data is the data matrix
	#kmax is the maximum number of groups
	#khigh is the maximum number of child groups for each group when computing mss
	#d is an indicator of which distance function to use.
	# choices are: "cosangle" (default),"abscosangle","euclid","abseuclid","cor","abscor"	
	#dmat is the distance matrix. if this has already been calculated by the user, it can
	# be passed into the function in order to save calculation time
	#within and between are either "med" for median split silhouette or "mean"
	# for mean split silhouette
	#ord is an indicator of how to order the clusters. choices are to maximize 
	# correlation ordering ("co") or to build a tree of cluster medoids ("clust")
mssinitlevel<-function(data,kmax=9,khigh=9,d="cosangle",dmat=NULL,within="med",between="med",ord="co"){
	if(!is.matrix(data))
		stop("First arg to mssinitlevel() must be a matrix")
	p<-length(data[,1])
	if(!is.matrix(dmat))
		dmat<-distancematrix(data,d)
	if(length(dmat[1,])!=p)
		stop("Data and distance matrix dimensions do not agree in mssinitlevel()")
	if(length(dmat[,1])!=p)
		stop("Distance matrix must be a square matrix in mssinitlevel()")
	m<-msscheck(dmat,kmax,khigh,within,between)
	if(m[1]==1){
		cat("No strong evidence for clusters in the first level - \n continuing to split root node anyway. \n")
		m<-msscheck(dmat,kmax,khigh,within,between,force=TRUE)
	}
	pamobj<-pam(dissvector(dmat),m[1],diss=TRUE)
	rowmedoids<-pamobj$medoids
	final<-ifelse(max(pamobj$clusinfo[,1])<3,1,0)
	if(m[1]>2){			
		medoidsdata<-as.matrix(data[rowmedoids,])
		l<-length(rowmedoids)
		medoidsdist<-dmat[rowmedoids,rowmedoids]
		if(ord=="co")
			medoidsord<-improveordering(medoidsdist)
		if(ord=="clust"){
			mpamobj<-pam(dissvector(medoidsdist),2,diss=TRUE)
			labelsmed<-mpamobj$clustering
			medmed<-mpamobj$medoids
			clussizes<-mpamobj$clusinfo[,1]
			prevlevel<-mssnextlevel(medoidsdata,list(2,medmed,clussizes,labelsmed,0,cbind(c(1,2),medmed)),dmat=medoidsdist,kmax,khigh,within,between)
			final<-prevlevel[[5]]
			if(final==0){
				depth<-(l-1)
				for(j in (1:depth)){
					if(final==0){
						clustnext<-mssnextlevel(medoidsdata,prevlevel,dmat=medoidsdist,kmax,khigh,within,between)
						final<-clustnext[[5]]
					}
					if(final==1){ 
						j<-depth
						prevlevel<-clustnext
					}
				}
			}					
			medoidsord<-orderelements(prevlevel,medoidsdata,rel="neighbor",d=d,dmat=medoidsdist)[[2]]
		}
		k<-m[1]
		rowmedoids<-rowmedoids[medoidsord]
		labels<-lab2<-pamobj$clustering
		for(j in (1:k)) 
			lab2[labels==medoidsord[j]]<-j
		output<-list(k,rowmedoids,pamobj$clusinfo[,1][medoidsord],lab2,final,cbind(1:k,rowmedoids))
	}
	else
		output<-list(2,pamobj$medoids,pamobj$clusinfo[,1],pamobj$clustering,final,cbind(1:2,pamobj$medoids))
	return(output)
}

#e. collapsing functions#
	#paircoll() collapses a pair of medoids (i,j)
	#collap() calls paircoll() to consider and possibly perform collapsing
	#msscollap() collapses by sequentially calling collap starting with 
	# the closest pair of clusters til there is no more improvement in mss
	#mssmulticollap tries all pairs of clusters and collapses any that improve mss
	###########################################################################################
	#data is the data matrix
	#level is level of the tree
	#d is an indicator of which distance function to use
	# choices are: "cosangle" (default),"abscosangle","euclid","abseuclid","cor","abscor"
	#dmat is the distance matrix. if this has already been calculated by the user, it can
	# be passed into the function in order to save calculation time.
	#newmed is an indicator of which way to find the medoid of the new cluster after collapsing.
	# choices are: "nn" to use the nearest neighbor of the clustersize-weighted 
	# mean of the two medoids as the medoid of a collapsed cluster, "uwnn" to use an unweighted 
	# version of nearest neighbor so that each cluster (rather than each gene) gets equal
	# weight in the mean, "center" to use the cluster center (element with min sum distance 
	# to all others), "medsil" (default) to use the medoid which maximizes the medoid based 
	# silhouette (i.e.: (a-b)/max(a,b), where a=dist(medoid), b=dist(next closest medoid)). 
	#the silhouettes and splits (arguments [[1]] and [[2]] of level) refer to the original
	# splits and loose their meaning if the child cluster(s) are collapsed
	#impr is a margin of improvement required to accept a collapse with msscollap and
	# mssmulticollap. the default is impr=0
paircoll<-function(i,j,data,level,d="cosangle",dmat=NULL,newmed="medsil"){
	if(!is.matrix(data))
		stop("First arg to paircoll() must be a matrix")
	p<-length(data[,1])
	k<-level[[1]]
	labels<-level[[4]]
	medoids<-level[[2]]
	clussizes<-level[[3]]
	N<-length(level[[6]][,1])
	block<-level[[6]][(N-k+1):N,]
	if(N==k)
		prevblock<-NULL
	else
		prevblock<-level[[6]][1:(N-k),]
	if(max(i,j)>k)
		stop("The clusters to collapse do not exist in paircoll()") 
	labeli<-labels[medoids[i]]
	labelj<-labels[medoids[j]]
	oldlabels<-labels
	labels[labels==labelj]<-labeli
	trunclabels<-trunc(oldlabels/10)
	labelparents<-unique(trunclabels)
	parenti<-order(labelparents)[labelparents==trunc(labeli/10)]
	parentj<-order(labelparents)[labelparents==trunc(labelj/10)]
	if(newmed=="nn")
		fakemed<-(data[medoids[i],]*clussizes[i]+data[medoids[j],]*clussizes[j])/(clussizes[i]+clussizes[j])
	if(newmed=="uwnn")
		fakemed<-(data[medoids[i],]+data[medoids[j],])/2
	if(newmed=="nn" || newmed=="uwnn"){ 
		rowsub<-(1:p)[labels==labeli]
		distsfm<-distancevector(data[rowsub,],as.vector(fakemed),d)
		medoids[i]<-rowsub[order(distsfm)[1]]
	}
	else{
		if(is.matrix(dmat)) 
			colldist<-dmat[labels==labeli,labels==labeli]
		else 
			colldist<-distancematrix(data[labels==labeli,],d)
		rowsub<-(1:p)[labels==labeli]
		if(newmed=="center"){
			sumdist<-rowSums(colldist)
			medoids[i]<-rowsub[order(sumdist)==1]
		}
		if(newmed=="medsil"){
			othermed<-medoids[-c(i,j)]
			collp<-length(labels[labels==labeli])
			othern<-length(othermed)
			if(othern==0)
				stop("Not enough medoids to use newmed='medsil' in paircoll()")
			if(is.matrix(dmat)){
			 	if(othern==1) 
					otherdist<-cbind(dmat[labels==labeli,othermed])
				else 
					otherdist<-rbind(dmat[labels==labeli,othermed])
			}
			else{
				if(othern==1)
					otherdist<-distancevector(data[labels==labeli,],data[othermed,],d)
				else{
					othermedmat<-data[othermed,]
					otherdist<-matrix(0,nrow=collp,ncol=othern)
					for(l in 1:othern)
						otherdist[,l]<-distancevector(data[labels==labeli,],othermedmat[l,],d)
				}
			}			
			if(othern==1)
				b<-otherdist
			else
				b<-apply(otherdist,1,min)
			b<-matrix(b,nrow=collp,ncol=collp)
			diag(b)<-0
			b<-abs(b-colldist)/pmax(colldist,b)
			sumdist<-rowSums(b)
			medoids[i]<-rowsub[rev(order(sumdist))==1]
		}
	}	
	k<-k-1	
	clussizes[i]<-clussizes[i]+clussizes[j]
	block[i,2]<-medoids[i]
	if(j<=k){			
		for(l in (j:k)){
			medoids[l]<-medoids[l+1]
 		        clussizes[l]<-clussizes[l+1]
			block[l,]<-block[l+1,]
		}
	}
	medoids<-medoids[1:k]
	clussizes<-clussizes[1:k]
	block<-block[1:k,]
	if(labels[labels==labeli][1]/10==trunclabels[labels==labeli][1]){
		labels[labels==labeli]<-labels[labels==labeli]+1
		labels[labels==labelj]<-labels[labels==labelj]+1
		block[i,1]<-block[i,1]+1
	}
	return(list(k,medoids,clussizes,labels,level[[5]],rbind(prevblock,block)))
}

#note: this version of collap does not have silhbased arg: for use with MSS only (not silhouettes)
collap<-function(data,level,d="cosangle",dmat=NULL,newmed="medsil"){
	if(!is.matrix(data))
		stop("First arg to collap() must be a matrix")
	k<-level[[1]]
	if(k<3){
		warning("Not enough medoids to use newmed='medsil' in collap() - \n using newmed='nn' instead \n") 
		newmed<-"nn"
	}
	medoids<-level[[2]]
	clussizes<-level[[3]]
	if(sum(is.na(clussizes))) 
		warning("NA in clussizes")
	medoidsdata<-data[medoids,]
	if(sum(is.na(medoidsdata))>0) 
		warning("Missing value(s) in medoidsdata in collap()")
	if(is.matrix(dmat)) 
		distmed<-dmat[medoids,medoids]
	else 
		distmed<-distancematrix(medoidsdata,d)
	distv<-dissvector(distmed)
	indexmin<-order(distv)[1]
	best<-vectmatrix(indexmin,k)
	clustfinal<-paircoll(best[1],best[2],data,level,d,dmat,newmed)
	return(clustfinal)
}

msscollap<-function(data,level,khigh,d="cosangle",dmat=NULL,newmed="medsil",within="med",between="med",impr=0){
	if(!is.matrix(data))
		stop("First arg to msscollap() must be a matrix")	
	if(!is.matrix(dmat)) 
		dmat<-distancematrix(data,d)
	newk<-level[[1]]
	mss1<-labelstomss(level[[4]],dmat,khigh,within,between)
	maxncoll<-max(0,newk-2)
        ncoll<-0
	coll<-1
        if(newk<=2) 
		coll<-0
 	while((coll==1) && (ncoll<= maxncoll)){
                levelc<-collap(data,level,d,dmat,newmed)
		mss2<-labelstomss(levelc[[4]],dmat,khigh,within,between)
		if(mss1==0) 
			r<-0
		else 
			r<-(mss1-mss2)/mss1
		if(r<impr) 
			coll<-0 	 
    		else{
			mss1<-mss2
			level<-levelc
                        ncoll<-ncoll+1
                }
	}
	return(level)
} 

mssmulticollap<-function(data,level,khigh,d="cosangle",dmat=NULL,newmed="medsil",within="med",between="med",impr=0){
	if(!is.matrix(data))
		stop("First arg to mssmulticollap() must be a matrix")
	if(!is.matrix(dmat)) 
		dmat<-distancematrix(data,d)
	medoids<-level[[2]]
	medoidsdata<-data[medoids,]
	if(sum(is.na(medoidsdata))>0) 
		warning("Missing value(s) in medoidsdata in mssmulticollap()")
	distmed<-dmat[medoids,medoids]
	k<-level[[1]]
	ord<-order(dissvector(distmed))
	mss1<-labelstomss(level[[4]],dmat,khigh,within,between)
	maxncoll<-max(0,k*(k-1)/2)
	ncoll<-0 
	i<-1
	while(i<=maxncoll){
		clusts<-vectmatrix(ord[i],k)
		if(k<3){
			if(newmed=="medsil") 
				warning("Can't use newmed=medsil with less than 3 clusters. \n Substituting newmed=nn")
			newmed<-"nn"
		}
		levelc<-paircoll(clusts[1],clusts[2],data,level,d,dmat,newmed)
		mss2<-labelstomss(levelc[[4]],dmat,khigh,within,between)
		if(mss1==0) 
			r<-0
		else 
			r<-(mss1-mss2)/mss1
		if(r>=impr){
			mss1<-mss2
			level<-levelc
                        ncoll<-ncoll+1
			k<-level[[1]]
			maxncoll<-max(0,k*(k-1)/2)
			i<-0
			medoids<-level[[2]]
			medoidsdata<-data[medoids,]
			if(sum(is.na(medoidsdata))>0) 
				warning("Missing value(s) in medoidsdata in mssmulticollap()")
			distmed<-dmat[medoids,medoids]
			ord<-order(dissvector(distmed))
	        }
		i<-i+1
	}
	return(level)
}

#f. iterating functions to run down the tree#
	#mssrundown() runs down the tree K levels with a stopping rule to 
	# find the main clusters
	#msscomplete() runs down the tree to the final level from any level
	#############################################################################################
	#data is the data matrix
	#K is the maximum number of levels to compute
	#kmax is the maximum number of groups
	#khigh  is the maximum number of child groups for each group when computing mss
	#d is an indicator of which distance function to use
	# choices are: "cosangle" (default),"abscosangle","euclid","abseuclid","cor","abscor"	
	#dmat is the distance matrix. if this has already been calculated by the user, it can
	# be passed into the function in order to save calculation time
	#coll is an indicator of how to collapse. the choices are to begin with the closest 
	# pair of clusters and collapse til there is no more improvement in mss ("seq") 
	# or to try all pairs of clusters and accept any collapse that improves mss ("all").
	#newmed is an indicator of which way to find the medoid of the new cluster after collapsing.
	# choices are: "nn" to use the nearest neighbor of the clustersize-weighted 
	# mean of the two medoids as the medoid of a collapsed cluster, "uwnn" to use an unweighted 
	# version of nearest neighbor so that each cluster (rather than each gene) gets equal
	# weight in the mean, "center" to use the cluster center (element with min sum distance 
	# to all others), "medsil" (default) to use the medoid which maximizes the medoid based 
	# silhouette (i.e.: (a-b)/max(a,b), where a=dist(medoid), b=dist(next closest medoid)). 
	#stop is an indicator that the tree should stop as soon as there is an increase in 
	# mss moving to the next level
	#finish is an indicator that the tree should compute all K levels and return the
	# one with the minimum mss (when finish==FALSE and stop==FALSE, level K is returned)
	#within and between are either "med" for median split silhouette or "mean"
	# for mean split silhouette.
	#impr is a margin of improvement required to accept a collapse with msscollap and
	# mssmulticollap. the default is impr=0.
mssrundown<-function(data,K=16,kmax=9,khigh=9,d="cosangle",dmat=NULL,initord="co",coll="seq",newmed="medsil",stop=TRUE,finish=FALSE,within="med",between="med",impr=0){
	if(!is.matrix(data))
		stop("First arg to mssrundown() must be a matrix")
	if(!is.matrix(dmat))
                dmat<-distancematrix(data,d)
	bestlevel<-level<-mssinitlevel(data,kmax,khigh,d,dmat,within,between,initord)
	bestmss<-mss<-labelstomss(level[[4]],dmat,khigh,within,between)
	bestl<-l<-1
	ind<-0
	cat("Searching for main clusters... \n")
	if(level[[5]]==1)
		return(level)
	while((l<=K) && (ind==0)){
		cat("Level ",l,"\n")
		if(coll=="seq")	
			levelc<-msscollap(data,level,khigh,d,dmat,newmed,within,between,impr)
		if(coll=="all") 
			levelc<-mssmulticollap(data,level,khigh,d,dmat,newmed,within,between,impr)
		mss<-labelstomss(levelc[[4]],dmat,khigh,within,between)
		if(mss>=bestmss & stop==TRUE)
			ind<-1
		else{
			if(mss<bestmss){
				bestlevel<-levelc
				bestmss<-mss
				bestl<-l
			}
		}
		l<-l+1
		if(l<=K){	
			level<-mssnextlevel(data,levelc,dmat,kmax,khigh,within,between)
			if(finish==TRUE){
				if(sum(trunc(level[[4]]/10)*10==level[[4]])==length(level[[4]]) & l<=K){
					ind<-1
					bestlevel<-levelc
					bestmss<-mss
					bestl<-(l-1)
				}
			}
		}
	}
	cat("Identified",bestlevel[[1]]," main clusters in level",bestl,"with MSS =",bestmss,"\n")
	return(bestlevel)
}

msscomplete<-function(level,data,K=16,khigh=9,d="cosangle",dmat=NULL,within="med",between="med"){
	if(!is.matrix(data))
		stop("First arg to msscomplete() must be a matrix")
	if(!is.matrix(dmat))
                dmat<-distancematrix(data,d)
	count<-digits(level[[4]][1])
	cat("Running down without collapsing from Level",count,"\n")
	while((max(level[[3]])>3) & (count<K)){
		level<-newnextlevel(data,level,dmat,2,khigh)
		count<-count+1
		cat("Level",count,"\n")
	}
	return(level)
}


#newnextlevel and newsplitcluster are needed in msscomplete to rundown 
#completely without collapsing back to the main clusters.
	#data is the data matrix
	#prevlevel is the level from which a next level is produced
	#dmat is the distance matrix, as above
	#klow and khigh are the min and max number of children at each node
	#newnextlevel produces the args to newsplitcluster and calls
	# this function to do the splitting of each node
newnextlevel<-function(data,prevlevel,dmat,klow=2,khigh=6){
	if(!is.matrix(data))
		stop("First arg to newnextlevel() must be a matrix")
	p<-length(data[,1])
	n<-length(data[1,])
	if(length(dmat[1,])!=p)
		stop("Distance and data matrix dimensions do not agree in newnextlevel()")
	if(length(dmat[,1])!=p)
		stop("Third arg to newnextlevel() is not a square matrix")
	id<-1:p
	k<-prevlevel[[1]]
	medoids<-prevlevel[[2]]
	labels<-prevlevel[[4]]
	newk<-0
	newlabels<-newmedoids<-newclussizes<-NULL
	count<-1
	ordlabels<-sort(unique(labels))
	if(length(ordlabels)!=k) 
		warning("Number of unique labels not equal number of clusters in newnextlevel()")
	if(sum(is.na(medoids))){
		warning("Missing value(s) in medoid vector in newnextlevel()")
		medoids[is.na(medoids)]<-FALSE
	}
	if(length(unique(medoids))<k && sum(medoids)) 
		warning("Medoids in newnextlevel() are not unique") 
	checkmeans<-FALSE
	if(length(medoids)==1 && !medoids){
		warning("No medoids provided in newnextlevel()")
		usemean<-TRUE
	}
	else{
		if(sum(medoids>1)==k) 
			usemean<-FALSE
		else 
			checkmeans<-TRUE
	}
	for(j in (1:k)){
		clust1<-data[labels==ordlabels[j],]	
		id1<-id[labels==ordlabels[j]]
		if(length(id1)>1){
			kmax<-min(c(khigh,dim(clust1)[1]-1))
			clust1<-as.matrix(clust1)
		}
		l1<-ordlabels[j]
		right<-(j<k)
		medoid1<-ifelse(is.na(medoids[j]),0,medoids[j])
		if (j<k) 
			medoid2<-medoids[j+1]
		else 
			medoid2<-medoids[j-1]
		if(length(id1)>1) 
			med2dist<-rowMeans(as.matrix(dmat[labels==ordlabels[j],labels==labels[medoid2]]))
		else 
			med2dist<-mean(dmat[labels==ordlabels[j],labels==labels[medoid2]])
		splitobj<-newsplitcluster(clust1,l1,id1,klow,kmax,medoid1,med2dist,right,dmat[labels==l1,labels==l1]) 
		newlabels[labels==ordlabels[j]]<-splitobj[[3]]
		k1<-splitobj[[1]]
		start<-count
		end<-count+k1-1
		newmedoids[start:end]<-splitobj[[2]]
		newclussizes[start:end]<-splitobj[[4]]
		count<-count+k1-1+1		
	}
	count<-count-1
	newk<-count
	newmedoids<-newmedoids[1:newk]
	newclussizes<-newclussizes[1:newk]
	final<-0
	if(count==k)
		final<-1
	if(max(newclussizes)==3) 
		final<-1
	return(list(newk,newmedoids,newclussizes,newlabels,final,rbind(rbind(prevlevel[[6]],cbind(sort(unique(newlabels)),newmedoids)))))
}

newsplitcluster<-function(clust1,l1,id1,klow=2,khigh=2,medoid1,med2dist,right,dist1){
	if(!medoid1) 
		warning("Medoid missing - continue to split cluster")
	else{
		if(sum(medoid1==id1)==0 & medoid1) 
			warning("Medoid not in cluster - continue to split cluster")
	}
	if(is.matrix(clust1)){			
		p1<-length(clust1[,1])
		n<-length(clust1[1,])
	}
	else p1<-1
	if(p1<3){
		k1<-1				
		newmedoids1<-medoid1
		newlabels1<-rep(10*l1,p1)
		newclussizes<-p1
	}
	else{
		l<-length(clust1[,1])
		dissvec<-dissvector(dist1)
		kmax<-min(p1-1,khigh)
		a<-rep(0,(kmax-klow+2))
		best<-2
		for(j in (klow:kmax)){
			a[j]<-pam(dissvec,j,diss=TRUE)$silinfo$avg.width
			if (a[j]>a[best]) best<-j
		}
		k1<-best
		pamobj<-pam(dissvec,k1,diss=TRUE)
		newclussizes<-pamobj$clusinfo[,1]
		newmedoids1<-id1[pamobj$medoids]
		newlabels1<-pamobj$clustering
		distnewmedoids<-NULL
		for(j in (1:k1)) 
			distnewmedoids[j]<-mean(med2dist[newlabels1==newlabels1[pamobj$medoids[j]]])  
		if(right==1) 
			ord<-rev(order(distnewmedoids))
		else 
			ord<-order(distnewmedoids)
		newmedoids1<-newmedoids1[ord]
		newclussizes<-newclussizes[ord]
		oldlab<-newlabels1
		for(j in (1:k1)) 
			newlabels1[oldlab==ord[j]]<-j
		newlabels1<-rep(10*l1,l)+newlabels1
	}
	for(a in (1:length(newmedoids1))){
		if(sum(newmedoids1[a]==id1)==0) 
			warning("Problem with new medoids after splitting cluster")
	}
	return(list(k1,newmedoids1,newlabels1,newclussizes))
}

#g. wrapper function to build the whole tree with clusters#
	#data is the data matrix.
	#clusters (default="best") tells how to identify the main clusters 
	# clusters="greedy" stops at the first level where MSS increases
	# clusters="none" does not identify main clusters
	# clusters="best" identifies the level<=K with best MSS as main clusters
	#K is the maximum number of levels to compute. right now, still have 
	# computational problem with doing more than 16, which is the default.
	#kmax is the maximum number of groups.
	#khigh  is the maximum number of child groups for each group when computing mss.
	#d is an indicator of which distance function to use.
	# choices are: "cosangle" (default),"abscosangle","euclid","abseuclid","cor","abscor".	
	#dmat is the distance matrix. if this has already been calculated by the user, it can
	# be passed into the function in order to save calculation time.
	#coll is an indicator of how to collapse. the choices are to begin with the closest 
	# pair of clusters and collapse til there is no more improvement in mss ("seq") 
	# or to try all pairs of clusters and accept any collapse that improves mss ("all").
	#newmed is an indicator of which way to find the medoid of the new cluster after collapsing.
	# choices are: "nn" to use the nearest neighbor of the clustersize-weighted 
	# mean of the two medoids as the medoid of a collapsed cluster, "uwnn" to use an unweighted 
	# version of nearest neighbor so that each cluster (rather than each gene) gets equal
	# weight in the mean, "center" to use the cluster center (element with min sum distance 
	# to all others), "medsil" (default) to use the medoid which maximizes the medoid based 
	# silhouette (i.e.: (a-b)/max(a,b), where a=dist(medoid), b=dist(next closest medoid)). 
	#mss is either "med" (default) for median split silhouettes or "mean" for mean 
	# split silhouettes.
	#impr is a margin of improvement required to accept a collapse with msscollap and
	# mssmulticollap. the default is impr=0.
	#initord is "co" (default) if improveordering() is used to order the clusters in 
	# the first level or "clust" if clsutering the medoids is used.
	#ord determines how elements are ordered within clusters: "co" is 
	# using improveordering(), "own" is distance to their own medoid, and "nieghbor"
	# is distance to the neighboring medoid (to the right). 
hopach<-function(data,dmat=NULL,d="cosangle",clusters="best",K=15,kmax=9,khigh=9,coll="seq",newmed="medsil",mss="med",impr=0,initord="co",ord="own"){
	if(inherits(data,"exprSet")) 
		data<-exprs(data)
	data<-as.matrix(data)
	if(K>15){
		K<-15
		warning("K set to 15 - can't do more than 15 splits")
	}
	if(K<1){
		K<-1
		warning("K set to 1 - can't do less than 1 level")
	}
	if(clusters!="none"){
		cuttree<-mssrundown(data,K,kmax,khigh,d,dmat,initord,coll,newmed,stop=(clusters=="greedy"),finish=TRUE,within=mss,between=mss,impr)
		if(cuttree[[1]]>1) 
			cutord<-orderelements(cuttree,data,rel=ord,d,dmat)[[2]]
		else 
			cutord<-NULL
		out1<-list(k=cuttree[[1]],medoids=cuttree[[2]],sizes=cuttree[[3]],labels=cuttree[[4]],order=cutord)
		finaltree<-msscomplete(cuttree,data,K,khigh,d,dmat,within=mss,between=mss)
	}
	else{
		out1<-NULL
		finaltree<-msscomplete(mssinitlevel(as.matrix(data),kmax,khigh,d,dmat,within=mss,between=mss,initord),data,K,khigh,d,dmat,within=mss,between=mss)
	}
	dimnames(finaltree[[6]])<-list(NULL,c("label","medoid"))
	out2<-list(labels=finaltree[[4]],order=orderelements(finaltree,data,rel=ord,d,dmat)[[2]],medoids=finaltree[[6]])
	return(list(clustering=out1,final=out2,call=match.call(),metric=d))
}
