#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define H 0.000001

/* Store precision of calculations as side effect of functions */
double PRECISION;

/* Calculate derivative of any order */
double deriv(double (*function)(), double x, double h, unsigned int order);

/* Numerical integrator using trapezoidal method -- !!uses numerical derivation!! -- */
double tr_integrate(double (*f)(), double start, double end, unsigned int steps);
/* This is what you normally use */
double trn_integrate(double (*f)(), double start, double end, unsigned int steps);


/* Numerical integrator using Simpson 1/3 rule */
double sot_integrate(double (*f)(), double start, double end, unsigned int steps);
/* This is what you normally use */
double sotn_integrate(double (*f)(), double start, double end, unsigned int steps);
/* Find the maximum value of a derivative of function (0 order is the function) */
double deriv_maximize(double (*f)(), unsigned int order, double a, double b, double precision);

double g(double x);

int main()
{
  double I;
  int steps = 1;
  do{
    steps++;
    I = trn_integrate(&g, 0, 3, steps);
  }while((PRECISION/I > 0.01) || PRECISION == 0);
  printf("trn_integrate Moment of Inertia from 0 to 3 with %d steps: I = %f +/- %f \%\n", steps, I, 100*PRECISION/I);

  steps=1;
  do{
    steps++;
    I = sotn_integrate(&g, 0, 3, steps);
  }while((PRECISION/I > 0.01) || PRECISION == 0);
  printf("sotn_integrate Moment of Inertia from 0 to 3 with %d steps: I = %f +/- %f \%\n", steps, I, 100*PRECISION/I);  
  
  return 0;
}

double g(double x)
{
  return pow(x,2)*(1+pow(x,3));
}

double deriv(double (*function)(), double x, double h, unsigned int order)
{
  if(order == 0)
    return function(x);
  /* Doubling h for higher makes derivative more precise, as we mitigate subtraction error for values
   * closing together after recursive calls
   */
  return (deriv(function, x+(h*order), h*order, order-1) - deriv(function, x-(h*order), h*order, order-1))/(2.0*h*order);
}

double deriv_maximize(double (*f)(), unsigned int order, double a, double b, double precision)
{
  double c;
  if(b<a)
    {
      c=a;
      a=b;
      b=c;
    }

  int i;
  for(c=deriv(f, a, precision, order),i=1;a<b;++i,a+=i*precision)
    if(f(a)>c)
      c=deriv(f, a, precision, order);

  return c;
}

double tr_integrate(double (*f)(), double start, double end, unsigned int steps)
{
  double h = (end-start)/((double) steps);

  double I=0;

  int i;
  for(i=1;i<steps;++i)
    I += f(start + i*h);

  I*=h;
  I+=0.5*h*(f(start)+f(end));

  /*Estimate precision */
  PRECISION = deriv_maximize(f, 2, start, end, H);
  PRECISION *= pow(h,2)*(start-end)/12.0;
  PRECISION = fabs(PRECISION);
  
  return I;
}

double trn_integrate(double (*f)(), double start, double end, unsigned int steps)
{
  /* Ih has half steps */
  double I=tr_integrate(f, start, end, steps);
  double Ih=tr_integrate(f, start, end, steps/2);

  PRECISION = fabs((Ih-I)/3.0);
  return I;
}

double sot_integrate(double (*f)(), double start, double end, unsigned int steps)
{
  double h = (end-start)/((double) steps);

  double I = 0;
  int i;
  for(i=1;i<steps;++i)
    if(i%2 == 0)
      I+=2*f(start+i*h);
    else
      I+=4*f(start+i*h);

  I += f(start)+f(end);
  I *= h/3.0;

  /* Estimate precision */
  PRECISION = fabs(deriv_maximize(f, 4, start, end, H)*pow(h,4)*(end - start)/2880.0);
  return I;
}

double sotn_integrate(double (*f)(), double start, double end, unsigned int steps)
{
  double I = sot_integrate(f,start,end,steps);
  double Ih = sot_integrate(f,start,end,steps/2);

  /* ... 'cause reasons */
  PRECISION = fabs((Ih-I)/15.0);
  return I;
}
  
  
