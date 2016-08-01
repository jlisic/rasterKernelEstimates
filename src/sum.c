
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>


#if defined _OPENMP
    #include <omp.h>
#endif




double sumKernel(
    double * x,    /* naip image */
    double * W,    /* pre computed spatial weights */
    size_t i,      /* current location in columns */
    size_t j,      /* current location in rows */
    size_t dRow,
    size_t dCol,
    size_t nRow,   /* number of Rows */
    size_t nCol    /* number of Columns */
  ) {

  /* adjustment that must be applied for edge effects */
  size_t k, l;

  size_t k_start;
  size_t k_stop;
  size_t l_start;
  size_t l_stop;

  double mu = 0;

  size_t k_local;
  size_t l_local;

  /* the starts */
  if( i < dRow/2 ) {
    k_start = 0; 
  } else {
    k_start = i - dRow/2 ;
  }
  if( j < dCol/2 ) {
    l_start = 0; 
  } else {
    l_start = j - dCol/2 ;
  }
  /* the stops */
  if( i + dRow/2 + 1 > nRow ) {
    k_stop = nRow; 
  } else {
    k_stop = i + dRow/2 + 1;
  }
  if( j + dCol/2 + 1  > nCol ) {
    l_stop = nCol; 
  } else {
    l_stop = j + dCol/2 + 1;
  }
        
  if( x[i*nCol + j] == INFINITY ) return( INFINITY);
  if( x[i*nCol + j] == -INFINITY ) return( -INFINITY);
  if( x[i*nCol + j] == NAN ) return( NAN);

  /* first pass variance */
  for(
      k=k_start, 
      k_local=k_start - i + (dRow/2); 
      k < k_stop; 
      k++, k_local++
      ) {
    for(
        l=l_start, 
        l_local=l_start -j + (dCol/2);
        l < l_stop;
        l++, l_local++
        ) {
  
      if( x[k * nCol + l] == INFINITY ) continue;
      if( x[k * nCol + l] == -INFINITY ) continue;
      if( x[k * nCol + l] == NAN ) continue;
     
        mu += x[k * nCol + l] * W[ k_local*dCol + l_local];

    }
  }
  return( mu ) ;
}


void rSmoothSums( 
    double * x,         /* this is the multi year naip images  */ 
    double * mu,        /* this is the input/returned mu */ 
    double * WMu,      /* weight */
    int * nRowPtr, 
    int * nColPtr,
    int * dRowPtr, 
    int * dColPtr
    ) {
 
  /* move R ints to size_t */

  size_t dRow = *dRowPtr;
  size_t dCol = *dColPtr;
  
  size_t nRow = *nRowPtr;
  size_t nCol = *nColPtr;

  size_t i,j;
  
#pragma omp parallel

#pragma omp parallel for private(j)
  for( i=0; i < nRow; i++) {
    for( j=0; j < nCol; j++) {
      mu[i*nCol + j] = sumKernel( x, WMu, i,j,dRow,dCol,nRow,nCol); 
    }
  }
#pragma omp barrier

  

  return;
}

