#Distance matrix related functions

#a. compute distances between all rows of a matrix#

#cosine-angle
disscosangle<-function(X,na.rm=TRUE){
	if(!is.matrix(X))
		stop("arg to disscosangle() must be a matrix")
	dX<-dim(X)
	p<-dX[1]
	n<-dX[2]
	if(na.rm){
		N<-rowSums(!is.na(X))
		N2<-(!is.na(X))%*%t(!is.na(X))
		X[is.na(X)]<-0
		N<-sqrt(N%*%t(N))/N2
	}
	else
		N<-1
	out<-rowSums(X^2)
	out<-1-N*tcrossprod(X)/sqrt(tcrossprod(out))
	diag(out)<-0
	suppressWarnings(out<-sqrt(out))
	out[out=="NaN"]<-0
	return(out)
}

#absolute cosine-angle
dissabscosangle<-function(X,na.rm=TRUE){
	if(!is.matrix(X))
		stop("arg to dissabscosangle() must be a matrix")
	dX<-dim(X)
	p<-dX[1]
	n<-dX[2]
	if(na.rm){
		N<-rowSums(!is.na(X))
		N2<-(!is.na(X))%*%t(!is.na(X))
		X[is.na(X)]<-0
		N<-sqrt(N%*%t(N))/N2
	}
	else
		N<-1
	out<-rowSums(X^2)
	out<-1-abs(N*tcrossprod(X)/sqrt(tcrossprod(out)))
	diag(out)<-0
	suppressWarnings(out<-sqrt(out))
	out[out=="NaN"]<-0
	return(out)
}

#euclidean
#note: disseuclid(X)=daisy(X)/sqrt(dim(X)[2])
disseuclid<-function(X,na.rm=TRUE){
	if(!is.matrix(X))
		stop("arg to disseuclid() must be a matrix")
	dX<-dim(X)
	p<-dX[1]
	n<-dX[2]
	if(na.rm){
		N<-rowSums(!is.na(X))
		N2<-(!is.na(X))%*%t(!is.na(X))
		X[is.na(X)]<-0
		out<-matrix(rep(rowSums(X^2)/N,p),ncol=p)+t(matrix(rep(rowSums(X^2)/N,p),ncol=p))-2*X%*%t(X)/N2	
	}
        else
		out<-matrix(rep(rowMeans(X^2),p),ncol=p)+t(matrix(rep(rowMeans(X^2),p),ncol=p))-2*X%*%t(X)/n
	diag(out)<-0
	suppressWarnings(out<-sqrt(out))
	out[out=="NaN"]<-0
	return(out)
}

#absolute euclidean
dissabseuclid<-function(X,na.rm=TRUE){    
	if(!is.matrix(X))
		stop("arg to dissabseuclid() must be a matrix")
	dX<-dim(X)
	p<-dX[1]
	n<-dX[2]
	if(na.rm){
		N<-rowSums(!is.na(X))		
		N2<-(!is.na(X))%*%t(!is.na(X))
		X[is.na(X)]<-0
		out1<-matrix(rep(rowSums(X^2)/N,p),ncol=p)+t(matrix(rep(rowSums(X^2)/N,p),ncol=p))-2*X%*%t(X)/N2	
		out2<-matrix(rep(rowSums(X^2)/N,p),ncol=p)+t(matrix(rep(rowSums(X^2)/N,p),ncol=p))+2*X%*%t(X)/N2
	}
	else{
	        out1<-matrix(rep(rowMeans(X^2),p),ncol=p)+t(matrix(rep(rowMeans(X^2),p),ncol=p))-2*X%*%t(X)/n
	        out2<-matrix(rep(rowMeans(X^2),p),ncol=p)+t(matrix(rep(rowMeans(X^2),p),ncol=p))+2*X%*%t(X)/n
	}
        out1<-pmin(out1,out2)
        diag(out1)<-0
        suppressWarnings(out1<-sqrt(out1))
 	out1[out1=="NaN"]<-0
	return(out1)
}

#correlation
disscor<-function(X,na.rm=TRUE){
	if(!is.matrix(X))
		stop("arg to disscor() must be a matrix")
	p<-dim(X)[1]
	if(na.rm)
		na<-"pairwise.complete.obs"
	else
		na<-"all.obs"
	out<-1-cor(t(X),use=na)
	diag(out)<-0
	suppressWarnings(out<-sqrt(out))
	out[out=="NaN"]<-0
	return(out)
}

#absolute correlation
dissabscor<-function(X,na.rm=TRUE){
	if(!is.matrix(X))
		stop("arg to dissabscor() must be a matrix")
	p<-dim(X)[1]
	if(na.rm)
		na<-"pairwise.complete.obs"
	else
		na<-"all.obs"
	out<-1-abs(cor(t(X),use=na))
	diag(out)<-0
	suppressWarnings(out<-sqrt(out))
	out[out=="NaN"]<-0
	return(out)
}

#b. compute distances between rows of a matrix and a vector#

#cosine-angle
vdisscosangle<-function(X,y,na.rm=TRUE){
	if(!is.matrix(X))
		stop("First arg to vdisscosangle() must be a matrix")
	if(!is.vector(y))
		stop("Second arg to vdisscosangle() must be a vector")
	dX<-dim(X)
	p<-dX[1]
	n<-dX[2]
	if(length(y)!=n)
		stop("Matrix and vector dimensions do not agree in vdisscosangle()")
	if(na.rm){
		N<-rowSums(!is.na(X))
		N2<-sum(!is.na(y))
		N3<-(!is.na(X))%*%(!is.na(y))
		X[is.na(X)]<-0
		y[is.na(y)]<-0
		N<-sqrt(N*N2)/N3
	}
	else
		N<-1
	suppressWarnings(out<-sqrt(as.vector(rep(1,p)-N*(X%*%y)/sqrt(rowSums(X^2)*sum(y^2)))))
	out[out=="NaN"]<-0
	return(out)
}

#absolute cosine-angle
vdissabscosangle<-function(X,y,na.rm=TRUE){
	if(!is.matrix(X))
		stop("First arg to vdissabscosangle() must be a matrix")
	if(!is.vector(y))
		stop("Second arg to vdissabscosangle() must be a vector")
	dX<-dim(X)
	p<-dX[1]
	n<-dX[2]
	if(length(y)!=n)
		stop("Matrix and vector dimensions do not agree in vdissabscosangle()")
	if(na.rm){
		N<-rowSums(!is.na(X))
		N2<-sum(!is.na(y))
		N3<-(!is.na(X))%*%(!is.na(y))
		X[is.na(X)]<-0
		y[is.na(y)]<-0
		N<-sqrt(N*N2)/N3
	}
	else
		N<-1
	suppressWarnings(out<-sqrt(as.vector(rep(1,p)-abs(N*(X%*%y)/sqrt(rowSums(X^2)*sum(y^2))))))	
	out[out=="NaN"]<-0
	return(out)
}

#euclidean
vdisseuclid<-function(X,y,na.rm=TRUE){
	if(!is.matrix(X))
		stop("First arg to vdisseuclid() must be a matrix")
	if(!is.vector(y))
		stop("Second arg to vdisseuclid() must be a vector")
	dX<-dim(X)
	p<-dX[1]
	n<-dX[2]
	if(length(y)!=n)
		stop("Matrix and vector dimensions do not agree in vdisseuclid()")
	if(na.rm){
		N1<-rowSums(!is.na(X))
		N2<-sum(!is.na(y))
		N3<-(!is.na(X))%*%(!is.na(y))
		X[is.na(X)]<-0
		y[is.na(y)]<-0
		suppressWarnings(out<-sqrt(as.vector(rowSums(X^2)/N1+sum(y^2)/N2-2*X%*%y/N3)))	
	}
	else
		suppressWarnings(out<-sqrt(as.vector(rowMeans(X^2)+mean(y^2)-2*X%*%y/n)))
	out[out=="NaN"]<-0
	return(out)
}

#absolute euclidean
vdissabseuclid<-function(X,y,na.rm=TRUE){
	if(!is.matrix(X))
		stop("First arg to vdissabseuclid() must be a matrix")
	if(!is.vector(y))
		stop("Second arg to vdissabseuclid() must be a vector")
	dX<-dim(X)
	p<-dX[1]
	n<-dX[2]
	if(length(y)!=n)
		stop("Matrix and vector dimensions do not agree in vdissabseuclid()")
	if(na.rm){
		N1<-rowSums(!is.na(X))
		N2<-sum(!is.na(y))
		N3<-(!is.na(X))%*%(!is.na(y))
		X[is.na(X)]<-0
		y[is.na(y)]<-0
		out1<-rowSums(X^2)/N1+sum(y^2)/N2-2*X%*%y/N3
		out2<-rowSums(X^2)/N1+sum(y^2)/N2+2*X%*%y/N3
	}
	else{
	        out1<-rowMeans(X^2)+mean(y^2)-2*X%*%y/n
       	 	out2<-rowMeans(X^2)+mean(y^2)+2*X%*%y/n
	}
        suppressWarnings(out1<-sqrt(as.vector(pmin(out1,out2))))
	out1[out1=="NaN"]<-0
	return(out1)
}

#correlation
vdisscor<-function(X,y,na.rm=TRUE){
	if(!is.matrix(X))
		stop("First arg to vdisscor() must be a matrix")
	if(!is.vector(y))
		stop("Second arg to vdisscor() must be a vector")
	p<-dim(X)[1]
	if(length(y)!=length(X[1,]))
		stop("Matrix and vector dimensions do not agree in vdisscor()")
	if(na.rm)
		na<-"pairwise.complete.obs"
	else
		na<-"all.obs"
	suppressWarnings(out<-sqrt(as.vector(rep(1,p)-cor(t(X),y,use=na))))
	out[out=="NaN"]<-0
	return(out)
}

#absolute correlation
vdissabscor<-function(X,y,na.rm=TRUE){
	if(!is.matrix(X))
		stop("First arg to vdissabscor() must be a matrix")
	if(!is.vector(y))
		stop("Second arg to vdissabscor() must be a vector")
	p<-dim(X)[1]
	if(length(y)!=length(X[1,]))
		stop("Matrix and vector dimensions do not agree in vdissabscor()")
	if(na.rm)
		na<-"pairwise.complete.obs"
	else
		na<-"all.obs"
	suppressWarnings(out<-sqrt(as.vector(rep(1,p)-abs(cor(t(X),y,use=na)))))
	out[out=="NaN"]<-0
	return(out)
}

#c. wrapper functions#

#makes a distance matrix from X using distance d#
distancematrix<-function(X,d,na.rm=TRUE){
	X<-as.matrix(X)
	if (d=="cosangle") return(disscosangle(X,na.rm))
	if (d=="abscosangle") return(dissabscosangle(X,na.rm))
	if (d=="euclid") return(disseuclid(X,na.rm))
	if (d=="abseuclid") return(dissabseuclid(X,na.rm))
	if (d=="cor") return(disscor(X,na.rm))
	if (d=="abscor") return(dissabscor(X,na.rm))
	#insert your own distance function here
	stop("Distance metric ",d," not available")
}

#makes a distance vector from X and y using distance d#
distancevector<-function(X,y,d,na.rm=TRUE){
	X<-as.matrix(X)
	y<-as.vector(y)
	if (d=="cosangle") return(vdisscosangle(X,y,na.rm))
	if (d=="abscosangle") return(vdissabscosangle(X,y,na.rm))
	if (d=="euclid") return(vdisseuclid(X,y,na.rm))
	if (d=="abseuclid") return(vdissabseuclid(X,y,na.rm))
	if (d=="cor") return(vdisscor(X,y,na.rm))
	if (d=="abscor") return(vdissabscor(X,y,na.rm))
	#insert your own distance function here
	stop("Distance metric ",d," not available")
}

#d. conversions#

#converts distance matrix to a vector#
dissvector<-function(M){
	if(!is.matrix(M))
		stop("arg to dissvector() must be a matrix")
	dM<-dim(M)
	if(dM[1]!=dM[2])
		stop("arg to dissvector() not a sqaure matrix")
	p<-dM[1]
	count<-1
	v<-rep(0,p*(p-1)/2)
	for (i in 1:(p-1)){
		v[count:(count+p-i-1)]<-M[i,(i+1):p]
		count<-count+p-i
	}
	return(v)
}

#converts distance vector to a matrix#
dissmatrix<-function(v){
	if(!is.vector(v))
		stop("arg to dissmatrix() must be a vector")
	p<-(1+sqrt(1+8*length(v)))/2
	M<-matrix(0,nrow=p,ncol=p)
	count<-1
	for (i in 1:(p-1)){
		M[i,(i+1):p]<-v[count:(count+p-i-1)]
		count<-count+p-i
	}
	return(M+t(M))
}

#maps index of distance vector into the [i,j] of corresponding pxp distance matrix#
vectmatrix<-function(index,p){
 	count<-1
 	s<-p-1
 	while((index-s)>0){
		s<-s+(p-1-count)
 		count<-count+1
 	}
 	i<-count
 	j<-p-(s-index)
 	return(c(i,j))
}


#e. correlation ordering#

#computes correlation ordering
correlationordering<-function(dist){
	if(!is.matrix(dist))
		stop("arg to correlationordering() must be a matrix")
	p<-length(dist[1,])
	if(p!=length(dist[,1]))
		stop("arg to correlationordering() not a square matrix")
	a<-dissvector(dist)
	b<-dissvector(abs(matrix(1:p,nrow=p,ncol=p,byrow=TRUE)-matrix(1:p,nrow=p,ncol=p)))
	return(cor(a,b))
}

#optimizes correlation ordering
improveordering<-function(dist,echo=FALSE){
	if(!is.matrix(dist))
		stop("arg to improveordering() must be a matrix")
	p<-length(dist[1,])
	if(p!=length(dist[,1]))
		stop("arg to improveordering() not a square matrix")
	v<-correlationordering(dist)
	if(echo)
		cat("Old order:",v,"\n")
	ord<-neword<-(1:p)
	if(!is.na(v)){
		final<-0
		for(gap in (1:(p-1))){
			while(final==0){
				sum<-0
				for(j in (1:(p-gap))){
					temp<-neword[j]
					neword[j]<-neword[j+gap]
					neword[j+gap]<-temp
					vnew<-correlationordering(dist[neword,neword])
					if(vnew<v) 
						neword<-ord
					else{
				 		sum<-sum+1
						v<-vnew
					}
					ord<-neword
				}
				if(sum==0) 
					final<-1
			}
		}
		if(echo)
			cat("New order:",correlationordering(dist[neword,neword]),"\n")
	}
	return(neword)		
}

