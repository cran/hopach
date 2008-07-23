
/***************************************************************************************/
/* C code to expand functionality of HOPACH package. Implementation of an hdist object */
/* subsetting operation. All indices are check in R before passing to these routines.  */
/* Author: Gregory D. Wall                                                             */
/* Updated: 12.10.2007                                                                 */
/***************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <Rinternals.h>
#include <Rdefines.h>

#define max(a,b)(a>b?a:b)

/***************************************************************************************/
/*                                 Internal C Struct                                   */
/***************************************************************************************/

typedef struct{
	int size ;
	int *array ;
} Index ;

/***************************************************************************************/
/*                                    Prototypes                                       */
/***************************************************************************************/

void buildHDIST (const double*, double*, Index *, Index *, const int) ;
void buildMatrix(const double*, double*, Index *, Index *, const int) ;
int  getPos(const int i, const int j, const int dim) ;

/***************************************************************************************/
/*                             C HDIST Subsetting Method                               */
/***************************************************************************************/

SEXP R_getSubset( SEXP dMat, SEXP dim, SEXP subsetRows, SEXP subsetCols, SEXP asHDIST ){
	SEXP ans = R_NilValue ;
	int size=0 ;
	double *data, *ansPtr ;

	Index rows, cols ;

	/* all parameters are checked in R */
	size       = INTEGER_VALUE(dim)          ;
	data       = NUMERIC_POINTER(dMat)       ; 
	rows.array = INTEGER_POINTER(subsetRows) ;
	cols.array = INTEGER_POINTER(subsetCols) ;
	rows.size  = length(subsetRows)          ;
	cols.size  = length(subsetCols)          ;

	if( INTEGER_VALUE(asHDIST) ){
		if( length(subsetRows) == 1 ){
			PROTECT( ans = allocVector(REALSXP,1) );
			ansPtr = REAL( ans ) ;
			ansPtr[0] = 0.0 ;
		}else{
			PROTECT( ans = allocVector( REALSXP,rows.size*(rows.size-1)/2 ) );
			ansPtr = REAL( ans ) ;
			buildHDIST(data, ansPtr, &rows, &cols, size) ;	
		}
		UNPROTECT(1) ;
		return(ans) ;
		
	//}else if( !INTEGER_VALUE(asHDIST) ){
	}else {
		PROTECT( ans = allocVector(REALSXP,rows.size*cols.size) );
		ansPtr = REAL( ans ) ;
		
		PROTECT( dim = allocVector(INTSXP, 2)) ;
		INTEGER(dim)[0] = length(subsetRows)   ;
		INTEGER(dim)[1] = length(subsetCols)   ;
	
		setAttrib( ans, R_DimSymbol, dim ) ;
		buildMatrix(data, ansPtr, &rows, &cols, size) ;	
	
		UNPROTECT(2) ;
		return(ans) ;
	}
	/*}else{
	#	return(ans) ;	
	#}
	*/
}


void buildMatrix(const double *data, double *ansPtr, Index *rows, Index *cols, const int dim){
	int pos=0, A=(dim*(dim-1))/2 ;
	int nrows = rows->size ;
	int ncols = cols->size ;
	
	for(int r=0; r<nrows; r++){
		pos = r ;
		for(int c=0; c<ncols; c++){
			if(rows->array[r] == cols->array[c]) { 
				ansPtr[pos]=0; 
				pos = pos + nrows ;
			}else if( max(rows->array[r],cols->array[c]) == rows->array[r]){
				//ansPtr[pos] = data[getPos(rows->array[r],cols->array[c],dim)] ;
				ansPtr[pos] = data[(A-(((dim-(cols->array[c]))*(dim-(cols->array[c])-1)/2)+(dim-(rows->array[r])+1)))] ;
				pos = pos + nrows;
			}else{
				//ansPtr[pos] = data[getPos(cols->array[c],rows->array[r],dim)] ;
				ansPtr[pos] = data[(A-(((dim-(rows->array[r]))*(dim-(rows->array[r])-1)/2)+(dim-(cols->array[c])+1)))] ;
				pos = pos + nrows; 
			}
		}
	}
}

void buildHDIST(const double *data, double *ansPtr, Index *rows, Index *cols, const int dim){
	int pos=0, A= (dim*(dim-1))/2 ;
	int nrows = rows->size ;
	int ncols = cols->size ;
	for(int r=0; r<nrows; r++){
		for(int c=r+1; c<ncols; c++){
			if( max(rows->array[r],cols->array[c]) == rows->array[r]) 
				//ansPtr[pos] =  data[getPos(rows->array[r],cols->array[c],dim)] ;
				ansPtr[pos] = data[(A-(((dim-(cols->array[c]))*(dim-(cols->array[c])-1)/2)+(dim-(rows->array[r])+1)))] ;
			else if( max(rows->array[r],cols->array[c]) == cols->array[c]) 
				//ansPtr[pos] = data[getPos(cols->array[c],rows->array[r],dim)] ;
				ansPtr[pos] = data[(A-(((dim-(rows->array[r]))*(dim-(rows->array[r])-1)/2)+(dim-(cols->array[c])+1)))] ;
			pos++ ;
		}
	}
}


/* Made this function inline in order to speed up subsetting */
int getPos(const int i, const int j, const int dim){
	return( (dim*(dim-1)/2) - ( ((dim-j)*(dim-j-1)/2) + (dim-i+1) ) ) ;
}

