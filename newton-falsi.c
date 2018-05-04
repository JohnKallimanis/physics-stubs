#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

/* h of numerical differentiation*/
#define H_DIV 0.00000000000001

/*Precision*/
#define E 0.00001

#define N_SOLV 3 /*Guess for N-R method*/
double f(double x);
double df(double x);

/* How many iterations. Set by solver functions */
int N_iter;

/*Newton-Raphson Solver*/
double nr_solve(double (*function)(), double (*derivative)(), double x, double precission);

/*Regula Falsi Solver*/
double rf_solve(double (*f)(), double a, double b, double precission);

/*Bisection Solver*/
double bs_solve(double (*f)(), double a, double b, double precission);

/* Symmetrical differentiation function */
double deriv(double (*function)(), double x, double h);

double ndf(double x);

int main()
{
  double x;
  x = nr_solve(&f,&df,N_SOLV,E);
  printf("Newton-Raphson converged after %d iterations to the root x=%.5f\n",N_iter, x);
  /* Try the same with numerical differentiation */
  x = nr_solve(&f,&ndf,N_SOLV,E);
  printf("Newton-Raphson with numerical differntiation converged after %d iterations to the root x=%.5f\n",N_iter, x);
  x = rf_solve(&f,1.5,3,E);
  printf("Regula falsi converged after %d iterations to the root x=%.5f\n",N_iter, x);
  x = bs_solve(&f,1.5,3,E);
  printf("Bisection method converged after %d iterations to the root x=%.5f\n", N_iter,x);
  return 0;
}

double f(double x)
{
  return pow(x,4) + (2*pow(x,3)) - (5*x*x) - (7*x) -5;
}

double df(double x)
{
  return (4*pow(x,3)) + (6*x*x) - (10*x) -7;
}

double nr_solve(double (*function)(), double (*derivative)(), double x, double precission)
{
  for(N_iter=0;fabs(function(x))>precission;x-=(function(x)/derivative(x)), ++N_iter);
  return x;
}

double deriv(double (*function)(), double x, double h)
{
  return (function(x+h)-function(x-h))/(2*h);
}

double ndf(double x)
{
  return deriv(&f, x, H_DIV);
}

double rf_solve(double (*f)(), double a, double b, double precission)
{
  N_iter=0;
  double x=b, x_old;
  do{
    x_old=x;
    x = (a*f(b) -b*f(a))/(f(b)-f(a));
    if(f(x)*f(b)<0)
      a=x;
    else
      b=x;
    
    ++N_iter;
  } while((fabs(x-x_old)>precission*fabs(x)) && fabs(f(x))>precission);
  return x;
}

double bs_solve(double (*f)(), double a, double b, double precission)
{
  N_iter=0;
  if(f(a)*f(b)>0)
    {
      N_iter=-1;
      return 0;
    } 
  if(fabs(f(a))<precission)
    return a;
  if(fabs(f(b))<precission)
    return b;

  double c;
  do{
    c=(a+b)/2.0;
    ++N_iter;
    if(f(a)*f(c)>0)
      a=c;
    else
      b=c;
  }while(fabs(a-b)>precission);

  return c;
}
