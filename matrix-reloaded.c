#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>

typedef struct
{
  size_t rows;
  size_t columns;
  double **v;
}matrix;

/* Create a matrix dynamically */
matrix create_matrix(size_t rows, size_t columns);

/* Pivot matrix */
matrix pivot(matrix A);

/* Doolittle algorithm. Return U matrix, stores L to pointer */
matrix doolittle_factorize(matrix A, matrix *L);

/* Easily fill a matrix with contents of array */
matrix fill(double *m, size_t rows, size_t columns);

/*Print a matrix*/
void print_matrix(matrix A);

int main()
{
  double m[3][3] = {
    1, 1, 1,
    3, 2, -1
    -1, 0, 2
  };
  
  matrix A =  create_matrix(3,3);
  int i,j;
  A = fill(m,3,3);
  print_matrix(A);
  matrix L = create_matrix(3,3);
  printf("Now I 'm gonna pivot\n");
  /*print_matrix(pivot(A));*/
  A = doolittle_factorize(A, &L);
  puts("U");
  print_matrix(A);
  puts("L");
  print_matrix(L);
  return 0;
}
matrix create_matrix(size_t rows, size_t columns)
{
  matrix m;
  m.v = (double **) malloc(rows*sizeof(double *));

  /*Allocate memory for each row*/
  int i;

  for(i=0;i<rows;++i)
    {
      m.v[i] = (double *) malloc(columns*sizeof(double));
      /* Initialize values to zero:*/
      memset(m.v[i], 0, columns*sizeof(double));
    }
  /* This can be accessed via normal array notation as 
   * v[i][j] <=> *(*(matrix+i)+j)
   */

  m.rows = rows;
  m.columns = columns;
  
  return m;
}

matrix pivot(matrix A)
{
  if(A.rows != A.columns)
    {
      errno = -1;
      return A;
    }

  matrix A_init = A;
  int k,j, i, prow, N = (int) A.rows;
  int pivot[N];
  double c, max;

  for(k=0;k<N;++k)
    {
      /* Find pivot row*/
      pivot[k] = k;
      max=fabs(A.v[k][k]);/*/*/
      for(j=k+1;j<N;++j)
	{
	  /*Compare this to max:*/
	  c=fabs(A.v[j][k]);
	  if(max<c)
	    {
	      max=c;
	      pivot[k]=j; /* Position of the pivot */
	    }
	}

      /*If not current row, interchange current row with pivot*/
      if(pivot[k]!=k)
	for(j=0;j<N;++j)
	  {
	    c=A.v[k][j];
	    A.v[k][j] = A.v[pivot[k]][j];
	    A.v[pivot[k]][j]=c;
	  }      
    }

  return A;
}

matrix doolittle_factorize(matrix A, matrix *L)
{
  int i,j,k;
  if((L->rows != A.rows)||(L->columns != A.columns))
    {
      errno = 5;
      return A;
    }

  /*Convenience*/
  int N = (int) A.rows;

  /* Will need to pivot first if zero in diagonal */
  for(i=0;i<N;++i)
    if(A.v[i][i] == 0)
      {
	A = pivot(A);
	break;
      }
  
  /* The diagonal of L is unit */
  for(i=0;i<N;++i)
    L->v[i][i] =1;
  
  for(i=0;i<N;++i)
    {
      /* Make U elements */
      for(j=i;j<N;++j)
	for(k=0;k<i;++k)
	  A.v[i][j] -= L->v[i][k]*A.v[k][j];
      /* Make L elements */
      for(j=i+1;j<N;++j)
	{
	  L->v[j][i] = A.v[j][i];
	  for(k=0;k<i;++k)      
	    L->v[j][i] -= L->v[j][k]*A.v[k][i];
	  L->v[j][i] /= A.v[i][i];
	}
    }
  /* Empty out remaining values of U from A */
  for(i=0;i<N;++i)
    for(j=0;j<N;++j)
      if(i>j)
	A.v[i][j]=0;
  /* Return U */
  return A;
}

matrix fill(double m[], size_t rows, size_t columns)
{
  matrix A = create_matrix(rows,columns);
  int i,j;
  for(i=0;i<rows;++i)
    for(j=0;j<columns;++j)
      {
	A.v[i][j] = *(m + i * rows + j);
      }
  return A;
}

void print_matrix(matrix A)
{
  int i,j;
  for(i=0;i<A.rows;++i)
    {
      for(j=0;j<A.columns;++j)
	printf("%.2f ", A.v[i][j]);
      puts("");
    }
}
  
