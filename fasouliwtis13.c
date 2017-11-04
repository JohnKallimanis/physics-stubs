#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define _USE_MATH_DEFINES
#define g 9.807 /* m/s */
#define dl 0.01
#define dT 0.1

/* Give a uniform random number in [0,x]*/
#define drand(x) ((random()/((float) RAND_MAX)) *x)


/* Return a random number following normal gaussian distribution*/
double rgauss(void);

struct data
{
  double period;
  double length;
};

struct erval
{
  double value;
  double error;
};

int main()
{
  srand(420545);

  /*Calculate the theoritical values we will distort with gaussian error*/
  struct data theoretical[4];
  int i;
  for(i=0;i<4;++i)
    {
      /* For length of 0.2 0.4 ... 1.0...*/
      theoretical[i].length = 0.2*(i+1);
      /* Calculate the period for a given g (see #defines) */
      theoretical[i].period = 2*M_PI*sqrt(theoretical[i].length/g);
    }
  
  struct erval g_exp[4];

  /* For least square fit:*/
  double w,sum_x=0,sum_y=0,sum_2x=0,sum_xy=0;
  
  /*Distort the values with random errors, "experimental values" thus*/
  struct data experimental[4];
  
  for(i=0;i<4;++i)
    {
      /*Error of lenght measurement δl = 0.01m*/
      experimental[i].length = theoretical[i].length + dl*rgauss();
      /*Error of period measurement δt = 0.1s*/
      experimental[i].period = theoretical[i].period + dT*rgauss();

      /* Calculate the g for each experiment */
      g_exp[i].value = 4*M_PI*M_PI*experimental[i].length/pow(experimental[i].period,2);
      /* Error propagation of relative error multiplied by value */
      g_exp[i].error = g_exp[i].value*sqrt(pow(dl/experimental[i].length,2)+4*pow(dT/experimental[i].period,2));

      /* Some calculations for least square fit that is done later x =: T^2 y=: l*/
      w = pow(experimental[i].period,2);
      sum_x += w; sum_2x += w*w;
      sum_y += experimental[i].length;
      sum_xy += experimental[i].length*w;
    }

  puts("Results for the for experiments:");
  puts("l(m) \t\t +/- \t\t  T(s) \t\t +/- \t\t g(m/s^2) \t\t +/-");
  for(i=0;i<4;++i)
    printf("%lf \t %lf \t %lf \t %lf \t %lf \t %lf \n", experimental[i].length, dl, experimental[i].period, dT, g_exp[i].value, g_exp[i].error);

  /* Weighted mean average for g*/
  double sum_wg=0,sum_w=0;
  for(i=0;i<4;++i)
    {
      w=1/pow(g_exp[i].error,2);
      sum_wg +=  w*g_exp[i].value;
      sum_w += w;
    }
  double g_mean = sum_wg/sum_w;
  double g_error_mean = 1/sqrt(sum_w);
  printf("g_mean = %lf +/- %lf m/s\n",g_mean,g_error_mean);

  /* Least square fit y = A +Bx */
  double D = 4*sum_2x-pow(sum_x,2);
  double A = ((sum_2x*sum_y) - (sum_x*sum_xy))/D;
  double B = ((4*sum_xy) - (sum_x*sum_y))/D;
  /* Calculate sigma_y */
  double sigma_y=0;
  for(i=0;i<4;++i)
    sigma_y += pow(experimental[i].length - A - (B*pow(experimental[i].period,2)), 2);
  sigma_y = sqrt(sigma_y/2.0);
  double dA = sigma_y*sqrt(sum_2x/D);
  double dB = sigma_y*sqrt(2.0/D);
  printf("T^2 = %lf + %lfx\n", A, B);
  g_mean = 4*pow(M_PI,2)*B;
  g_error_mean = g_mean*dB/B;
  printf("g_ls = %lf +/- %lf\n",g_mean,g_error_mean);
  
  return 0;
}

double rgauss(void)
{
  static int have_numbers = 0;
  static double ua,ub;
  if(have_numbers)
    {
      have_numbers=0;
      return sqrt(-2*log(ua))*cos(2*M_PI*ub);
    }
  else
    {
      ua = drand(1);
      ub = drand(1);
      have_numbers=1;
      return sqrt(-2*log(ua))*sin(2*M_PI*ub);
    }
}
