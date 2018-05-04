#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

/* Precision asked for */
#define EPSILON 0.0001

/* Solve f(x) = x^2 + 10cosx with fixed-point method */

/* Symmetrical differentiation function */
double deriv(double (*function)(), double x, double h);
/* Fixed point iterator */
double fp_solve(double (*function)(), double x);

/* Most precise possible differentiation h*/
#define H_DIV DBL_MIN 

/*Possible fixed-point functions*/
double g1(double x);
double g2(double x);

int main(int argc, char *argv[])
{
  if(argc != 2)
    {
      printf("Usage: %s start_of_iteration\n", argv[0]);
      return 1;
    }
  double x= atof(argv[1]);

  /* Test each possible function. A third one, g(x)= -10cos(x) will obviously never converge */
  /* This way I avoid analytical work on the functions altogether. */
  if(deriv(&g1,x,H_DIV)<1)
    {
      printf("Using g1(x) = sqrt(10cosx)\n");
      printf("x = %.4f\n",fp_solve(&g1,x));  
    if(deriv(&g2,x,H_DIV)<1)
      {
	printf("Using g2(x) = acos(-x^2/10)\n");
	printf("x = %.4f\n",fp_solve(&g2,x));  
      }
    else
      {
	printf("Will not converge\n");
	return 1;
      }
    }
  /*
   *The function is symmetrical as of y'y. I can find a maximum of two roots with two functions. 
   *The roots are four. Those found, and their opposites
   */
    return 0;
}

double deriv(double (*function)(), double x, double h)
{
  return (function(x+h)-function(x-h))/(2*h);
}

double fp_solve(double (*function)(), double x)
{
  for(;fabs(function(x)-x)>EPSILON;x=function(x));
  return x;
}

double g1(double x)
{
  return sqrt(fabs(10*cos(x)));
}

double g2(double x)
{
  return acos(-pow(x,2)/10.0);
}
