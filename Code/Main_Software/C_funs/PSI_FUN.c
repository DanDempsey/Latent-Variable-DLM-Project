
#include "required_libs.h"

//void PSI_FUN(double *xi, int *y, double R[][*ymax], int *y_seq, int *ymax, int *psi);
int sampleC(int *x, double *p, int len_p);

void PSI_FUN(double *xi, int *y, int *N, int *y_seq, int *ymax, double R[][*ymax], int *psi) {
  
  // Declarations
  int i, j;
  double term1, term2;
  
  // Construct the matrix of probability distributions
  for( i=2; i<*ymax; i++ ) {
    for( j=1; j<=i; j++ ) {
      term1 = (double) (i - 1) / i;
      term2 = (double) *xi / i;
      R[i][j] = R[i-1][j]*term1 + R[i-1][j-1]*term2;
    }
  }
  
  // Sample the latent variable from constructed probability distributions
  for ( i=0; i<*N; i++ ) {
    double *Ri = R[y[i]];
    *psi += sampleC( y_seq, Ri, *ymax );
    //psi[i] = sampleC( y_seq, Ri, *ymax );
  }
  
  return;
  
}

int sampleC(int *x, double *p, int len_p) {
  
  int j;
  int i = 1;
  double u, running_sum = p[0];
  double sum_p = 0.0;
  
  for ( j=0; j<len_p; j++ ) {
    sum_p += p[j];
  }
  
  GetRNGstate();
  
  u = runif( 0, sum_p );
  while ( u > running_sum ) {
    running_sum += p[i];
    i++;
  }
  
  PutRNGstate();  
    
  return x[i-1];
    
}

