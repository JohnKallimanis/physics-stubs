#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define _USE_MATH_DEFINES
/* Numerically Solve throw at an angle*/

const double g = -10.0;

double *make_dynamic_array(size_t slots);
double *init_dynamic_array(double *m, size_t slots);

typedef struct
{
  size_t step;
  size_t size;
  double *t;
  double *y;
  double *v;
} dimension;

typedef struct
{
  dimension x;
  dimension y;
} plane;

dimension make_dimension(size_t size);
dimension dimension_resize(dimension d, size_t new_size);
void init_dimension(size_t size, dimension *d);
void free_dimension(dimension d);

/* Vertical axis */
double ver_dydt(double t, double x, double v);
double ver_dvdt(double t, double x, double v);

/* Horizontal axis */
double hor_dydt(double t, double x, double v);
double hor_dvdt(double t, double x, double v);

/* Solvers need to be fed a properly initialised pointer to initialised plane p */
/* Euler solver */
void euler_solve_free_fall(double angle, double time_step, double in_vel, plane *p);   

/* 4th Order Runge-Kutta solver*/
void rk_solve_free_fall(double angle, double time_step, double in_vel, plane *p);

/* Analytical Solver for free fall*/
void true_solve_free_fall(double angle, double time_step, double in_vel, plane *p);

/* Export a plane to a .dat file for processing */
void export_plane(char *filename, plane p);

/* Export the values of plane p with error values from the difference of each value from corresponding values on plane e*/
void export_compare(char *filename, plane p, plane e);

int main()
{
  plane p,e;
  p.x = make_dimension(100000);
  p.y = make_dimension(100000);
  e.x = make_dimension(100000);
  e.y = make_dimension(100000);
  /* Test multiple time steps with both methods, generate plottable files */
  double time_step;
  char *filename = malloc(20*sizeof(char));
  
  for(time_step=0.01; time_step>0.000001; time_step *= 0.1)
    {
      printf("Generating Euler data files for time step: %f\n", time_step);
      true_solve_free_fall(M_PI/4.0, time_step, 1,  &e);
      euler_solve_free_fall(M_PI/4.0, time_step, 1,  &p);
      sprintf(filename, "Euler_step:%f.dat", time_step);
      export_compare(filename, p, e);
      printf("Generating Runge-Kutta data files for time step: %f\n", time_step);
      rk_solve_free_fall(M_PI/4.0, time_step, 1,  &p);
      sprintf(filename, "RungeKutta_step:%f.dat", time_step);
      export_compare(filename, p, e);
    }
  
  return 0;
}

void export_compare(char *filename, plane p, plane e)
{
  FILE *dat = fopen((const char *) filename, "w");
  fprintf(dat, "t x y xdelta ydelta\n");

  int i;
  for(i=0;i<p.x.step;++i)
    fprintf(dat,"%f %.10f %.10f %.10f %.10f\n", p.x.t[i], p.x.y[i], p.y.y[i], fabs(p.x.y[i] - e.x.y[i]), fabs(p.y.y[i] - e.y.y[i]));
 
  fclose(dat);
  return;
}
