#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include "median.h"


#if defined _OPENMP
  #include <omp.h>
  #define CSTACK_DEFNS 7
  #include "Rinterface.h"
#endif





/* integer max */
int intMax ( int x, int y) {
  return( ( x > y) ? x : y ) ; 
}

/* integer min */
int intMin ( int x, int y) {
  return( ( x < y) ? x : y ) ; 
}




/* modal kernel */
double modalKernel(
    int * x,    /* raster image */
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
  
  size_t M = 0;
  size_t m = 0;
  double maxValue = -INFINITY; /* used to determine max weighted value */
  int mu = 0;

  size_t k_start;
  size_t k_stop;
  size_t l_start;
  size_t l_stop;
  
  size_t k_local;
  size_t l_local;


  int *    maxArray      = (int *) calloc( dRow * dCol, sizeof(int) );
  double * maxArrayValue = (double *) calloc( dRow * dCol, sizeof(double) );
  // handle tie breaks 
  double tieBreak;
  double maxTie = runif(0.0,1.0);

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
        
        if( x[k*nCol + l] == NAN ) continue;
      if( x[k*nCol + l] >= 0 ) {  /* only run over non-negative values */
        
        for(m=0; m < M; m++) {
          /* increment found values */
          if( maxArray[m] == x[k*nCol + l]  ) { 
            maxArrayValue[m] += W[ k_local*dCol + l_local];
            break;
          }
        }
        /* if the value is not found add it */
        if( m == M) {
          maxArray[m] = x[k*nCol + l ];
          maxArrayValue[m] = W[ k_local*dCol + l_local];
          M++;
        }
      }
    }
  }
      
  /* handle the all NA case */ 
  if( M == 0 ) {
    free(maxArray);
    free(maxArrayValue);
    return( -1 ) ;
  }
  
  /* determine max value */ 
  for(m=0; m < M ; m++) { 
    if( maxArrayValue[m] > maxValue ) {
      maxValue = maxArrayValue[m];
      mu = maxArray[m];
      // handle ties 
    } else if( maxArrayValue[m] == maxValue ) {
      tieBreak = runif(0.0, 1.0);
      if( tieBreak > maxTie ) {  
        maxValue = maxArrayValue[m];
        mu = maxArray[m];
        maxTie = tieBreak;
      }
    }
  }

  free(maxArray);
  free(maxArrayValue);
  return( mu ) ;
}




/* quantile kernel */
double quantileKernel(
    double * x,    /* naip image */
    double * W,    /* pre computed spatial weights */
    double quantile,      /* quantile */
    size_t i,      /* current location in columns */
    size_t j,      /* current location in rows */
    size_t dRow,
    size_t dCol,
    size_t nRow,   /* number of Rows */
    size_t nCol    /* number of Columns */
  ) {

  /* adjustment that must be applied for edge effects */
  size_t k, l;
  size_t quantile_t; /* size_t quantile */

  size_t k_start;
  size_t k_stop;
  size_t l_start;
  size_t l_stop;


  double * medianArray = (double *) calloc( dRow * dCol, sizeof(double) );
  double ** medianArrayPtr = (double **) calloc( dRow * dCol, sizeof(double *) );
  for(k=0; k <dRow*dCol;k++) medianArrayPtr[k] = &medianArray[k];
  
  double mu;
  int m = 0;
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
  
        if( x[k*nCol + l] == NAN ) continue;
       
        // only consider elements with positive valued weights 
        if( W[ k_local*dCol + l_local] > 0 ) {
          medianArray[m] = x[k*nCol + l];
          m++;
        }
    

    }
  }

  if ( m > 0) {
    quantile_t = (size_t) ( (double) m * quantile ); 
    //Rprintf("m=%d quantile=%f quantile_t=%d\n", (int) m, quantile, (int) quantile_t );
    mu = quantile_quickSelectIndex( medianArrayPtr, quantile_t, m );
  } else {
    mu = NAN;
  }
  

  free(medianArray);
  free(medianArrayPtr);
  return( mu ) ;
}




/* generic kernel */
double meanKernel(
    double * x,    /* naip image */
    double * var,  /*  */
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

  double w = 0;        /* total weight, used to make weight adjustments */
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
        w += W[ k_local*dCol + l_local];

    }
  }

  return( mu/w ) ;
}


/* generic kernel */
double gaussianKernel(
    double * x,    /* naip image */
    double hInv,    /* pre computed spatial weights */
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

  double w = 0;        /* total weight, used to make weight adjustments */
  double w2 = 0;
  double mu = 0;



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
  for( k=k_start; k < k_stop; k++) {
    for( l=l_start; l < l_stop; l++) {
  
      if( x[k * nCol + l] == INFINITY ) continue;
      if( x[k * nCol + l] == -INFINITY ) continue;
      if( x[k * nCol + l] == NAN ) continue;
     
        w = (x[k * nCol + l] - x[i * nCol + j]) *hInv;
        w *= w;
        mu += exp( -0.5 * w ) * 0.3989423 * hInv; 

        w2 += 1.0;
    }
  }

  if( w2 > 0) mu = mu/w2;

  return( mu ) ;
}




/* variance kernel */
double varKernel(
    double * x,    /* naip image */
    double * mu,  /*  */
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

  double w = 0;        /* total weight, used to make weight adjustments */
  double var = 0; /* smoothed x value we are goinng to return */

  double varTmp;

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
       
  /* correctly handle NAN and INF cases */ 
  if( x[i*nCol + j] == INFINITY ) return( INFINITY);
  if( x[i*nCol + j] == -INFINITY ) return( -INFINITY);
  if( x[i*nCol + j] == NAN ) return( NAN);

  /* 
   * k_start creates a link to the original data
   * k_local creates a link to the weights 
   */

  // second pass for variance 
  for(
      k=k_start, 
      k_local=k_start - i + (dRow/2); 
      k < k_stop; 
      k++, k_local++
      ) {
    for(
        l=l_start, 
        l_local=l_start - j + (dCol/2);
        l < l_stop;
        l++, l_local++
        ) {
        
        /* not mathematically correct, but good enough */ 
        if( x[k * nCol + l] == INFINITY ) continue;
        if( x[k * nCol + l] == -INFINITY ) continue;
        if( x[k * nCol + l] == NAN ) continue;

        if( mu[k * nCol + l] == INFINITY ) continue;
        if( mu[k * nCol + l] == -INFINITY ) continue;
        if( mu[k * nCol + l] == NAN ) continue;

        varTmp = x[k * nCol + l] - mu[i* nCol + j];
        var += varTmp * varTmp * W[ k_local*dCol + l_local];
        w += W[ k_local*dCol + l_local] ; 
    }
  }

  return( var/w  ) ;
}


void rSmoothLocalMoments( 
    double * x,         /* this is the multi year naip images  */ 
    double * mu,        /* this is the input/returned mu */ 
    double * var,        /* this is the input/returned Var */ 
    double * WMu,      /* weight */
    double * WVar,      /* weight */
    int * nRowPtr, 
    int * nColPtr,
    int * dRowPtr, 
    int * dColPtr,
    int * momentsPtr
    ) {
 
  /* move R ints to size_t */

  size_t dRow = *dRowPtr;
  size_t dCol = *dColPtr;
  
  size_t nRow = *nRowPtr;
  size_t nCol = *nColPtr;

  size_t i,j;
  

#if defined _OPENMP
  R_CStackLimit=(uintptr_t)-1;
#endif


#pragma omp parallel for private(j)
  for( i=0; i < nRow; i++) {
    for( j=0; j < nCol; j++) {
      mu[i*nCol + j] = meanKernel( x, var,  WMu, i,j,dRow,dCol,nRow,nCol); 
    }
  }
#pragma omp barrier

  if( *momentsPtr > 1) {
#pragma omp parallel for private(j)
    for( i=0; i < nRow; i++) {
      for( j=0; j < nCol; j++) {
        var[i*nCol + j] = varKernel( x, mu, WMu, i,j,dRow,dCol,nRow,nCol); 
      }
    }
#pragma omp barrier
  }
  

  return;
}



void rSmoothCategorical( 
    int * x,         /* this is the multi year naip images  */ 
    int * mu,        /* this is the input/returned mu */ 
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

#if defined _OPENMP
  R_CStackLimit=(uintptr_t)-1;
#endif
  
#pragma omp parallel for private(j)
  for( i=0; i < nRow; i++) {
    for( j=0; j < nCol; j++) {
      if( x[i*nCol+j] >= 0 ) {
        mu[i*nCol + j] = modalKernel( x, WMu, i,j,dRow,dCol,nRow,nCol); 
      } else {
        mu[i*nCol + j] = x[i*nCol + j];
      }
    }
  }
#pragma omp barrier

  return;
}



void rSmoothLocalQuantile( 
    double * x,         /* this is the multi year naip images  */ 
    double * mu,        /* this is the input/returned mu */ 
    double * WMu,      /* weight */
    double * quantile, /* quantile */
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
  
#if defined _OPENMP
  R_CStackLimit=(uintptr_t)-1;
#endif


#pragma omp parallel for private(j)
  for( i=0; i < nRow; i++) {
    for( j=0; j < nCol; j++) {
      mu[i*nCol + j] = quantileKernel( x, WMu, *quantile, i,j,dRow,dCol,nRow,nCol); 
    }
  }
#pragma omp barrier


  return;
}


