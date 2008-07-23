################################################################################
#                  First Attempt at HDIST object for HOPACH                    #
#                  G.D.Wall 10-3-07                                            #
################################################################################

################################################################################
#                  S4 Class Definition                                         #
################################################################################

setClass( "hdist", representation(Data = "numeric", 
				 Size = "numeric", 
				 Labels = "numeric",
				 Call = "character"),
				 prototype=list(
					Data=numeric(0),
					Size=numeric(0),
					Labels=numeric(0),
					Call=character(0)
	)
)	

"hdist" <- function(Data, Size, Labels = NULL, Call=""){
	new("hdist",Data=Data,Size=Size,Labels=Labels,Call=Call)	
}

################################################################################
#                  Type Checking                                               #
################################################################################

"is.hdist" <- function(x){
	return(class(x) == "hdist")		  
}

################################################################################
#                  HDIST get functions					       #
################################################################################
setMethod("as.vector", signature(x="hdist", mode="missing"), function(x){
	return(x@Data) 
})

setMethod("length", signature(x="hdist"), function(x){
	return(x@Size) 
})

if(!isGeneric("labels")){
	setGeneric("labels", function(object,...) standardGeneric("labels"))
}

setMethod("labels", signature(object="hdist"), function(object){
	return(object@Labels)
})


################################################################################
#                  Class Conversions                                           #
################################################################################

###  Matrix and dist conversions to a HDIST object ###
as.hdist <- function(x){
	if( is.hdist(x) ){
		return(x)
	}else if(is.matrix(x)){
		return(as.hdist(x)) 
	}
	else{
		warning("No coversion function exists for object...")
		return(x)
	}
}

###  Matrix to a HDIST object ###
setMethod("as.hdist", signature(x="matrix"), function(x){
	as(x,"hdist") 
})

setAs("matrix","hdist",function(from)
{
	out<-hdist(Data=dissvector(from), Size=dim(from)[1], Labels=(1:dim(from)[1]), Call="")
	return(out)
})


###  HDIST to a matrix object ###
setMethod("as.matrix", signature(x="hdist"), function(x){
	as(x,"matrix") 
})

setAs("hdist","matrix", function(from)
{
	size <- from@Size
	df <- matrix(0,size,size)
	df[row(df) > col(df)] <- from@Data
	df <- df + t(df)
	labels <- from@Labels
	dimnames(df) <- if(is.null(labels))
			list(1:size,1:size)
		else	
			list(labels,labels)
	df
})

################################################################################
#                  HDIST print method 					       #
################################################################################
setMethod( "show", signature(object="hdist"), function(object){
	.print.hdist(object)
})

.print.hdist <- function(x, diag = NULL, upper = NULL, digits = getOption("digits"),
	justify = "none", right = TRUE, ...)
{
	if(length(x@Data) > 0){
		if(is.null(diag))
			diag <- if(is.null(a <- diag))
				FALSE
			else a
		if(is.null(upper))	
			upper <- if(is.null(a <- upper))
				FALSE
		else a
		m <- as.matrix(x)
		i <- x@Labels
		dimnames(m) <- list(i,i) 
		cf <- format(m, digits = digits, justify = justify)
		if(!upper)
			cf[row(cf)<col(cf)] <- ""
		if(!diag)
			cf[row(cf) == col(cf)] <- ""
		print( if(diag || upper) 
			cf
			else cf[-1,-(dim(x)[1]),drop=FALSE],quote=FALSE,
				right=right,...)
	}
	else{
		cat(data.class(x), "(0)\n", sep="")	
	}	
	invisible(x)
}


################################################################################
#                  HDIST dim method 					       #
################################################################################
setMethod( "dim", "hdist", function(x){
	return(c(x@Size,x@Size))	
})

################################################################################
#                           HDIST Subsetting Method                            # 
################################################################################

setMethod("[", signature="hdist", definition=function(x, i=NA, j=NA, ..., drop=TRUE ){
	ans <- numeric() 

	##### remove all NAs, convert bools to numerics, remove 0 indices #####
	if( is.logical(i) ){
		if( length(i) == x@Size ){
			i <- which(i == TRUE) 
		}else if( length(i) < x@Size ){
			i <- rep(i, length.out=x@Size)		  
			i <- which(i == TRUE) 
		}else{
			stop("Invalid i index")
		}
	}else if( is.numeric(i) ){
		i <- floor(i) 
		i <- i[which(i != 0)]
	}else
		stop("Invalid i index")

	if( is.logical(j) )
		if( length(j) == x@Size ){
			j <- which(j == TRUE) 
		}else if( length(j) < x@Size ){
			j <- rep(j, length.out=x@Size)		  
			j <- which(j == TRUE) 
		}else{
			stop("Invalid j index")
		}
	else if( is.numeric(j) ){
		j <- floor(j) 
		j <- j[which(j != 0)]
	}else
		stop("Invalid j index")
	###### --------------------------------------------------------- ######
		
	##### Expand all empty/default indices                           ######
	if( length(i) == 0 && length(j) == 0 )
		return(as.matrix(x))
	else if( length(i) == 0 )
		i <- 1:x@Size 		  
	else if( length(j) == 0 )
		j <- 1:x@Size		  
	###### --------------------------------------------------------- ######
	
	if( any(i > x@Size) || any(j > x@Size) )
		stop("Invalid index")

	if( length(i) == length(j) ){
		if( all(i == j) ){
			ans <- .Call("R_getSubset",as.numeric(x@Data),as.integer(x@Size),as.integer(i),as.integer(j), as.integer(1))
			ans <- new("hdist",Data=ans, Size=length(i), Labels=i, Call=x@Call)
		}else{
			ans <- .Call("R_getSubset",as.numeric(x@Data),as.integer(x@Size),as.integer(i),as.integer(j), as.integer(0))
			dimnames(ans) <- list(i,j) 
		}
	}
	else{
		ans <- .Call("R_getSubset",as.numeric(x@Data),as.integer(x@Size),as.integer(i),as.integer(j), as.integer(0))
		dimnames(ans) <- list(i,j) 
	}

	return(ans) 
})	 

