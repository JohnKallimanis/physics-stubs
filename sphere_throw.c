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

int main()
{
  /* Find maximum range with aerodynamic drag */
  double theta;
  double max_range=0;
  double theta_max=0;

  plane p;
  p.x = make_dimension(10000);
  p.y = make_dimension(10000);
  for(theta=0;theta<M_PI/2.0;theta+=M_PI/180.0)
    {
      euler_solve_free_fall(theta, 0.001, 1,  &p);
      if(p.x.y[p.x.step] > max_range)
	{
	  max_range = p.x.y[p.x.step];
	  theta_max = theta;
	}
    }

  printf("Euler Solver: Maximum range at angle %f degrees\n", theta_max*180.0/M_PI);

  for(theta=0;theta<M_PI/2.0;theta+=M_PI/180.0)
    {
      euler_solve_free_fall(theta, 0.001, 1,  &p);
      if(p.x.y[p.x.step] > max_range)
	{
	  max_range = p.x.y[p.x.step];
	  theta_max = theta;
	}
    }
  
  printf("4th Order Runge-Kutta Solver: Maximum range at angle %f degrees\n", theta_max*180.0/M_PI);
  
  return 0;
}
double *init_dynamic_array(double *m, size_t slots)
{
  memset(m,0,slots*sizeof(double));
  return m;
}

double *make_dynamic_array(size_t slots)
{
  double *m = (double *) malloc(slots*sizeof(double));
  m = init_dynamic_array(m, slots);
  
  return m;
}

void init_dimension(size_t size, dimension *d)
{
  d->step = 0;
  d->size = size;
  init_dynamic_array(d->t, size);
  init_dynamic_array(d->y, size);
  init_dynamic_array(d->v, size);
}

dimension make_dimension(size_t size)
{
  dimension d;
  d.step = 0;
  d.size = size;
  d.t = make_dynamic_array(size);
  d.y = make_dynamic_array(size);
  d.v = make_dynamic_array(size);

  return d;
}

dimension dimension_resize(dimension d, size_t new_size)
{
  if(d.size == new_size)
    return d;

  dimension n = make_dimension(new_size);

  int i;

  if(new_size>d.size)
    new_size = d.size;
  
  for(i=0;i<new_size;++i)
    {
      n.t[i] = d.t[i];
      n.y[i] = d.y[i];
      n.v[i] = d.v[i];
    }

  return n;
}

void free_dimension(dimension d)
{
  free(d.t);
  free(d.y);
  free(d.v);

  return;
}

double ver_dydt(double t, double x, double v)
{
  return v;
}
  

double ver_dvdt(double t, double x, double v)
{
  return g - 0.25*fabs(g)*v;
}

double hor_dydt(double t, double x, double v)
{
  return v;
}

double hor_dvdt(double t, double x, double v)
{
  return 0 - 0.25*fabs(g)*v;
}

void euler_solve_free_fall(double angle, double time_step, double in_vel, plane *p)
{
  init_dimension(p->x.size, &(p->x));
  init_dimension(p->y.size, &(p->y));
  
  p->x.v[0] = in_vel*cos(angle);
  p->y.v[0] = in_vel*sin(angle);

  int i = 0;
  do{
    p->y.t[i+1] = p->y.t[i] + time_step;
    p->x.t[i+1] = p->y.t[i] + time_step;
    
    p->y.y[i+1] = p->y.y[i] + time_step*ver_dydt(p->y.t[i], p->y.y[i], p->y.v[i]);
    p->y.v[i+1] = p->y.v[i] + time_step*ver_dvdt(p->y.t[i], p->y.y[i], p->y.v[i]);
    
    p->x.y[i+1] = p->x.y[i] + time_step*hor_dydt(p->x.t[i], p->x.y[i], p->x.v[i]);
    p->x.v[i+1] = p->x.v[i] + time_step*hor_dvdt(p->x.t[i], p->x.y[i], p->x.v[i]);
    
    if(i==p->x.size)
	return;
    
    i++;
    p->y.step = i;
    p->x.step = i;
  }while(p->y.y[i] > 0);
  
  return;
}

void rk_solve_free_fall(double angle, double time_step, double in_vel, plane *p)
{
  init_dimension(p->x.size, &(p->x));
  init_dimension(p->y.size, &(p->y));
  
  p->x.v[0] = in_vel*cos(angle);
  p->y.v[0] = in_vel*sin(angle);

  double k1,k2,k3,k4,l1,l2,l3,l4;
  int i = 0;
  do{
    /* Time */
    p->y.t[i+1] = p->y.t[i] + time_step;
    p->x.t[i+1] = p->y.t[i] + time_step;
    /* Vertical Position and Velocity */
    k1 = time_step*ver_dydt(p->y.t[i], p->y.y[i], p->y.v[i]);
    l1 = time_step*ver_dvdt(p->y.t[i], p->y.y[i], p->y.v[i]);
    
    k2 = time_step*ver_dydt(p->y.t[i] + 0.5*time_step, p->y.y[i] + 0.5*k1, p->y.v[i] + 0.5*l1);
    l2 = time_step*ver_dvdt(p->y.t[i] + 0.5*time_step, p->y.y[i] + 0.5*k1, p->y.v[i] + 0.5*l1);
    
    k3 = time_step*ver_dydt(p->y.t[i] + 0.5*time_step, p->y.y[i] + 0.5*k2, p->y.v[i] + 0.5*l2);
    l3 = time_step*ver_dvdt(p->y.t[i] + 0.5*time_step, p->y.y[i] + 0.5*k2, p->y.v[i] + 0.5*l2);
    
    k4 = time_step*ver_dydt(p->y.t[i] + time_step, p->y.y[i] + k3, p->y.v[i] + l3);
    l4 = time_step*ver_dvdt(p->y.t[i] + time_step, p->y.y[i] + k3, p->y.v[i] + l3);
    
    p->y.y[i+1] = p->y.y[i] + (1.0/6.0)*(k1+2*k2+2*k3+k4);
    p->y.v[i+1] = p->y.v[i] + (1.0/6.0)*(k1+2*k2+2*k3+k4);

    /* Horizontal Position and Velocity */
    k1 = time_step*hor_dydt(p->x.t[i], p->x.y[i], p->x.v[i]);
    l1 = time_step*hor_dvdt(p->x.t[i], p->x.y[i], p->x.v[i]);
    
    k2 = time_step*hor_dydt(p->x.t[i] + 0.5*time_step, p->x.y[i] + 0.5*k1, p->x.v[i] + 0.5*l1);
    l2 = time_step*hor_dvdt(p->x.t[i] + 0.5*time_step, p->x.y[i] + 0.5*k1, p->x.v[i] + 0.5*l1);
    
    k3 = time_step*hor_dydt(p->x.t[i] + 0.5*time_step, p->x.y[i] + 0.5*k2, p->x.v[i] + 0.5*l2);
    l3 = time_step*hor_dvdt(p->x.t[i] + 0.5*time_step, p->x.y[i] + 0.5*k2, p->x.v[i] + 0.5*l2);
    
    k4 = time_step*hor_dydt(p->x.t[i] + time_step, p->x.y[i] + k3, p->x.v[i] + l3);
    l4 = time_step*hor_dvdt(p->x.t[i] + time_step, p->x.y[i] + k3, p->x.v[i] + l3);
    
    p->x.y[i+1] = p->x.y[i] + (1.0/6.0)*(k1+2*k2+2*k3+k4);
    p->x.v[i+1] = p->x.v[i] + (1.0/6.0)*(k1+2*k2+2*k3+k4);

    if(i==p->x.size)
	return;
    
    i++;
    p->y.step = i;
    p->x.step = i;
    
  }while(p->y.y[i] > 0);

  return;
}
