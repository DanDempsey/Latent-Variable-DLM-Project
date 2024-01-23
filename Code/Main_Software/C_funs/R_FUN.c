#include "required_libs.h"

//void R_FUN(double *xi, int *ymax, double *R[*ymax+1][*ymax+1]);

void R_FUN(double *xi, int *ymax, double R[][*ymax]) {
  
  // Declarations
  int i, j;
  double term1, term2;
  
  // Construct the matrix of probability distributions
  for( i=2; i<*ymax; i++ ) {
    for( j=1; j<=i; j++ ) {
      //printf( "%d", i );
      term1 = (double) (i - 1) / i;
      term2 = (double) *xi / i;
      R[i][j] = R[i-1][j]*term1 + R[i-1][j-1]*term2;
    }
  }
  
  return;
  
}

