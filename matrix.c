#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void print33array(double A[3][3]);
void the_lu_treader(double *A[3][3], double *L[3][3]);

int main()
{
  double A[3][3] = {
    {1, 1, 1},
    {-1, 0, 2},
    {3, 2, -1}
  };

  puts("A =");
  print33array(A);
  
  int i,j,k,f;
  /* The L matrix */
  double L[3][3];
  /*Initialize the diagonal*/
  for(i=0;i<3;++i)
    L[i][i]=1;

  /* The U matrix */
  double U[3][3];

  /*Fill U matrix*/
  for(i=0;i<3;++i)
    for(j=0;j<3;++j)
      {
	U[i][j] = A[i][j];
	for(k=0;k<i-1;++k)
	  {
	    L[i][k] = A[i][k];
	    for(f=0;f<k-1;++f)
	      L[i][k]-=U[f][k]*L[i][f];
	    L[i][k]/=U[k][k];
	    U[i][j]-=U[k][j]*L[i][k];
	  }
      }
  /*Fill L matrix*/
  
  puts("U =");
  print33array(U);
    
  return 0;
}

void print33array(double A[3][3])
{
  int i,j;
  for(i=0;i<3;++i)
    {
      for(j=0;j<3;++j)
	printf("%f ", A[i][j]);
      puts("");
    }
  return;
}

void the_lu_treader(double *A[3][3], double *L[3][3])
{
  *L[3][3] = {};
  int i;
  
  for(i=0;i<3;++i)
    *L[i][i]=1;
  return;
}
