#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* Give a uniform random number in [0,x]*/
#define drand(x) ((random()/((float) RAND_MAX)) *x)

/*The edges of the area of integration*/
#define X_MAX 2
#define Y_MAX 2
#define X_MIN 0
#define Y_MIN 0
int main(int argc, char *argv[])
{
  if((argc<3)||(argc>5))
    {
      printf("Usage: %s iterations seed [points.dat] [output]\n", argv[0]);
      return 1;
    }
  double x,y,k,den,M=0,Mx=0,My=0;
  int i;

  /*For errors:*/
  double dM=0,dMx=0,dMy=0;

  FILE *dat;
  if(argc == 4)
    dat = fopen(argv[3],"w");

  srand(atof(argv[2]));

  double N_ITER = atof(argv[1]);
  
  for(i=0;i<N_ITER;++i)
    {
      x = X_MIN + drand(X_MAX);
      y = Y_MIN + drand(Y_MAX);
      den = x*x + y*y;
      if(den>1 && y<2-x)
	{
	  if(argc == 4)
	    fprintf(dat, "%lf %lf\n",x,y);

	  den += 1;
	  
	  k = x;
	  x = y*den;
	  y = k*den;

	  M += den;
	  Mx += x;
	  My += y;
	  /*errors*/
	  dM += den*den;
	  dMx += x*x;
	  dMy += y*y;
	}
    }

  double V= (Y_MAX-Y_MIN)*(X_MAX-X_MIN);


  M/=N_ITER;  Mx/=N_ITER; My/=N_ITER;
  dM/=N_ITER; dMx/=N_ITER; dMy/=N_ITER;

  /* Errors of the integrals */
  dM = V*sqrt((dM-M*M)/N_ITER);
  dMy = V*sqrt((dMy-My*My)/N_ITER);
  dMx = V*sqrt((dMx-Mx*Mx)/N_ITER);

  /* The values of the integrals */
  M*=V; Mx*=V; My*=V;
  
  /* Store the coordinates here */
  y = My/M;
  x = Mx/M;

  /* Errors of coordinates: square propagation of relative error */
  dM/=M;
  dMy = sqrt(dM*dM+pow(dMy/My,2))*My;
  dMx = sqrt(dM*dM+pow(dMx/Mx,2))*Mx;

  FILE *out;
  if(argc == 5)
    out = fopen(argv[4],"w");
  else
    out = stdout;
  
  fprintf(out,"M = %lf +/- %lf \nx_cm = %lf +/- %lf\ny_cm = %lf +/- %lf\n",M, dM, x, dMx, y, dMy);

  return 0;  
}
