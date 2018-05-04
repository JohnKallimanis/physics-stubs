#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>

double f(double x);

double Lagrange(double (*f)(),double *p, int order, double x);

int main()
{
  double p[5] = { 1, 1.1, 1.2, 1.3, 1.4 };

  /* Store polynomial factors here:*/
  puts("Will approximate f(x) = e^2x -1 with 4th order lagrange polynomial");
  printf("f(1.25) = %.4f\n", Lagrange(&f, p, 5, 1.25)); 
  return 0;
}

double f(double x)
{
  return exp(2*x) -1;
}

double Lagrange(double (*f)(),double *p, int order, double x)
{
  int j,m;

  double prod, L=0;

  for(j=0;j<order;++j)
    {
      prod=1;
      for(m=0;m<order;++m)
	if(m!=j)
	  prod *= ((x - p[m])/(p[j]-p[m]));
      L += f(p[j])*prod;
    }

  return L;
}
