#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <Rinternals.h>
#include <Rdefines.h>

void get_means(const double *data, double *mu, int row1, int row2, int p, int n) ;

/**************************************************************/
/*                Distance functions Called from R            */
/**************************************************************/
/* NOTES:                                                     */
/* X is a vector of the stacked columns of a data matrix      */
/* nr and nc are the number of rows and columns               */ 
/* all functions return a vector of stacked columns of the    */
/* lower triangle distance matrix just calculated             */
/**************************************************************/


/* ----------------------- correlation ---------------------- */  
SEXP R_disscor( SEXP X, SEXP nr, SEXP nc, SEXP na_rm ){
	SEXP ans ;
	int p, n, index=0, index_i, index_j ;
	double *data, *ansPtr, mu[2], sum1=0, sum2=0, sum3=0, check=0 ;
	p = INTEGER_VALUE(nr) ;
	n = INTEGER_VALUE(nc) ;

	PROTECT( ans = allocVector(REALSXP, p*(p-1)/2) ) ;
	ansPtr = REAL( ans ) ;
	data = NUMERIC_POINTER(X) ; 

	for( int i = 0; i < p-1; i++) {
		for( int j = i+1; j < p; j++ ){
			get_means(data, mu, i, j, p, n) ;

			sum1 = sum2 = sum3 = 0 ;
			index_i = i ; index_j = j ;

			for( int k = 0; k < n; k++) {
				if( !isnan(data[index_i]) && !isnan(data[index_j]) ){
					sum1 = sum1 + (data[index_i]-mu[0])*(data[index_j]-mu[1]) ;
					sum2 = sum2 + pow((data[index_i]-mu[0]),2) ;
					sum3 = sum3 + pow((data[index_j]-mu[1]),2) ;
				}
				index_i = index_i + p ;
				index_j = index_j + p ;
			}
			check = sqrt( sum2 * sum3 ) ;

			/* if check is zero then at least one of the row intensities */
			/* has an error. e.g. containing only one valid intensity    */
			/* -99999 is a temporary indicator of a 'bad' value          */
			if( check ) { ansPtr[index] = 1 - sum1/check ; }
			else { ansPtr[index] = -99999 ; }

			index++ ;
		} 
	}

	UNPROTECT(1) ;
	return(ans) ;
}


/* ------------------- absolute correlation ---------------- */  
SEXP R_dissabscor( SEXP X, SEXP nr, SEXP nc, SEXP na_rm ){
	SEXP ans ;
	int p, n, index=0, index_i, index_j ;
	double *data, *ansPtr, mu[2], sum1=0, sum2=0, sum3=0, check=0 ;
	p = INTEGER_VALUE(nr) ;
	n = INTEGER_VALUE(nc) ;

	PROTECT( ans = allocVector(REALSXP, p*(p-1)/2) ) ;
	ansPtr = REAL( ans ) ;
	data = NUMERIC_POINTER(X) ; 

	for( int i = 0; i < p-1; i++) {
		for( int j = i+1; j < p; j++ ){
			get_means(data, mu, i, j, p, n) ;

			sum1 = sum2 = sum3 = 0 ;
			index_i = i ; index_j = j ;

			for( int k = 0; k < n; k++){
				if( !isnan(data[index_i]) && !isnan(data[index_j]) ){
					sum1 = sum1 + (data[index_i]-mu[0])*(data[index_j]-mu[1]) ;
					sum2 = sum2 + pow((data[index_i]-mu[0]),2) ;
					sum3 = sum3 + pow((data[index_j]-mu[1]),2) ;
				}
				index_i = index_i + p ;
				index_j = index_j + p ;
			}
			check = sqrt( sum2 * sum3 ) ;

			/* if check is zero then at least one of the row intensities */
			/* has an error. e.g. containing only one valid intensity    */
			/* -99999 is a temporary indicator of a 'bad' value          */
			if(check) { ansPtr[index] = 1 - fabs(sum1/check) ; }
			else { ansPtr[index] = -99999 ; }

			index++ ;
		} 
	}

	UNPROTECT(1) ;
	return(ans) ;
}


/* ---------------------- cosine angle  -------------------- */  
SEXP R_disscosangle( SEXP X, SEXP nr, SEXP nc, SEXP na_rm ){
	SEXP ans ;
	int p, n, index=0, index_i, index_j ;
	double *data, *ansPtr, sum1=0, sum2=0, sum3=0, check=0 ;
	p = INTEGER_VALUE(nr) ;
	n = INTEGER_VALUE(nc) ;

	PROTECT( ans = allocVector(REALSXP, p*(p-1)/2) ) ;
	ansPtr = REAL( ans ) ;
	data = NUMERIC_POINTER(X) ; 

	for( int i = 0; i < p-1; i++) {
		for( int j = i+1; j < p; j++ ){
			sum1 = sum2 = sum3 = 0 ;
			index_i = i ; index_j = j ;

			for( int k = 0; k < n; k++){
				if( !isnan(data[index_i]) && !isnan(data[index_j]) ){
					sum1 = sum1 + (data[index_i])*(data[index_j]) ;
					sum2 = sum2 + pow((data[index_i]),2) ;
					sum3 = sum3 + pow((data[index_j]),2) ;
				}
				index_i = index_i + p ;
				index_j = index_j + p ;
			}

			check = sqrt( sum2 * sum3 ) ;

			/* if check is zero then at least one of the row intensities */
			/* has an error. e.g. containing only one valid intensity    */
			/* -99999 is a temporary indicator of a 'bad' value          */

			if( check ) { ansPtr[index] = 1 - sum1/check ; }
			else { ansPtr[index] = -99999 ; }

			index++ ;
		} 
	}

	UNPROTECT(1) ;
	return(ans) ;
}

/* ---------------- absolute cosine angle  ---------------- */  
SEXP R_dissabscosangle( SEXP X, SEXP nr, SEXP nc, SEXP na_rm ){
	SEXP ans ;
	int p, n, index=0, index_i, index_j ;
	double *data, *ansPtr, sum1=0, sum2=0, sum3=0, check=0 ;
	p = INTEGER_VALUE(nr) ;
	n = INTEGER_VALUE(nc) ;

	PROTECT( ans = allocVector(REALSXP, p*(p-1)/2) ) ;
	ansPtr = REAL( ans ) ;
	data = NUMERIC_POINTER(X) ; 

	for( int i = 0; i < p-1; i++) {
		for( int j = i+1; j < p; j++ ){
			sum1 = sum2 = sum3 = 0 ;
			index_i = i ; index_j = j ;

			for( int k = 0; k < n; k++){
				if( !isnan(data[index_i]) && !isnan(data[index_j]) ){
					sum1 = sum1 + (data[index_i])*(data[index_j]) ;
					sum2 = sum2 + pow((data[index_i]),2) ;
					sum3 = sum3 + pow((data[index_j]),2) ;
				}
				index_i = index_i + p ;
				index_j = index_j + p ;
			}
			check = sqrt( sum2 * sum3 ) ;
			
			/* if check is zero then at least one of the row intensities */
			/* has an error. e.g. containing only one valid intensity    */
			/* -99999 is a temporary indicator of a 'bad' value          */
			if( check ) { ansPtr[index] = 1 - fabs(sum1/check) ; }
			else { ansPtr[index] = -99999 ; }

			index++ ;
		} 
	}

	UNPROTECT(1) ;
	return(ans) ;
}

/* ----------------------- Internal Routines -------------------- */
void get_means(const double *data, double *mu, int row1, int row2, int p, int n){
	double sum1 = 0, sum2 = 0, total = 0 ;
	int index_i = row1 ; 
	int index_j = row2 ;

	for( int k = 0; k < n; k++){
		if( !isnan(data[index_i]) && !isnan(data[index_j]) ){
			sum1 = sum1 + data[index_i] ;
			sum2 = sum2 + data[index_j] ;
			total++ ;
		}
		index_i = index_i + p ;
		index_j = index_j + p ;
	}
	mu[0] = sum1/total ;
	mu[1] = sum2/total ;
}




