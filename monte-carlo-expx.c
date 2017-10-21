#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define drand(x) ((random()/((float) RAND_MAX)) *x)

int main(int argc, char *argv[])
{
  /* Integrate e^x from bounds given from argv with steps and seed given 
   * from argv 
   */
  if((argc < 5) || (argc > 6) )
    {
      printf("Usage: %s from to steps seed [file.dat]\n", argv[0]);
      return 1;
    }

  FILE *dat;
  if(argc == 6)
    dat = fopen(argv[5], "w");

  printf("Using random seed %f\n", atof(argv[4]));
  srand(atof(argv[4]));

  int N = atoi(argv[3]);
  printf("Will perform %d steps\n", N);

  double l = atof(argv[1]), r = atof(argv[2]);
  printf("Will perform Monte Carlo hit-or-miss integration from %f to %f\n", l, r);

  int i;
  double x,y, integral;

  for(i=0;i<N;++i)
    {
      x = l + drand((r-l));
      y = exp(l) + drand((exp(r) -exp(l)));

      if(y<exp(x))
	{
	  ++integral;
	  if(argc == 6)
	    fprintf(dat, "%f %f\n", x, y);
	}
    }

  double imprecission = (1.0/(sqrt(integral)))*sqrt(1-integral/((float) N));
  integral = (integral/((double) N))*(r - l)*fabs(exp(r)-exp(l));
  if(exp(l)<exp(r))
    integral+=exp(l)*(r-l);
  else
    integral+=exp(r)*(r-l);
  
  printf("Integral of e^x from %f to %f = %.2f +/- %.2f\n", l,r,integral, integral*imprecission);
  if(argc == 6)
    fclose(dat);
  
  return 0;
}
