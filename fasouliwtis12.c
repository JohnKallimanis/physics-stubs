#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* Give a uniform random number in [0,x]*/
#define drand(x) ((random()/((float) RAND_MAX)) *x)

double integrable(double r);

struct erval
{
  double value;
  double error;
};

/* Integrator function using Crude Monte Carlo Method. 
 * Make sure to seed random number generator*/
struct erval cr_integrate(double (*function)(),double start, double end, double iterations);

int main(int argc, char *argv[])
{
  if(argc!=4)
    {
      printf("Usage: %s parts iterations seed\n",argv[0]);
      return 1;
    }

  srand(atof(argv[3]));
  int parts = atoi(argv[1]);
  double iterations = atoi(argv[2]);
  if( parts < 1 || iterations < 0)
    return 1;

  printf("Will integrate in %d part(s) doing %d iterations\n", parts, (int) iterations);

  struct erval result;
  result.value=0;
  result.error=0;
  
  int i,j;
  struct erval inter;
  
  double iterations_of_part = iterations/((double) parts);
  for(i=0;i<parts;++i)
    {
      inter=cr_integrate(&integrable,4*i/((double) parts),4*(i+1)/((double) parts),iterations_of_part);
      result.value += inter.value;
      result.error += pow(inter.error,2);
    }


  printf("I = %lf +/- %lf (cgs)\n",result.value, sqrt(result.error/((double) parts)));
  return 0;
}

struct erval cr_integrate(double (*function)(),double start, double end, double iterations)
{
  struct erval result;
  
  result.value=0;
  result.error=0;

  int i;
  double x;
  for(i=0;i<iterations;++i)
    {
      x = start + ((end-start)/fabs(end-start))*drand(fabs(end-start));
      result.value += function(x);
      result.error += pow(function(x),2);
    }
  result.value /= iterations;
  result.error /= iterations;
  result.error = (end-start)*sqrt(result.error-pow(result.value,2));
  result.value *= end-start;

  return result;
}
  
double integrable(double r)
{
  double d = r*r;
  return d*(1+exp(d));
}
