#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define E 0.00001
double g(double x);

int main()
{
  int i;
  double x = 2;

  for(i=1;fabs(g(x)-x)>E;x=g(x),++i);
  printf("The solution of x^3 -x -1 = 0 in [1,2] is x=%.5f \nNeeded %d iterations to reach this precision.\n",x,i);
  
  return 0;
}

/* The only g(x) = x likely to converge */
double g(double x)
{
  return pow(x+1,1.0/3.0);
}
