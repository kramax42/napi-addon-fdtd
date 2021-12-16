#ifndef PURE_C_FDTD_2D
#define PURE_C_FDTD_2D

#include <math.h>    // exp , M_PI
#include <stdio.h>   // printf
#include <stdlib.h>  // size_t

typedef struct {
  size_t ticks;

  /* Grid steps. */
  double dx;
  double dt;

  /* const */
  double aa1;   //??
  double Ti;    //??
  double tMax;  //??

  /* Grid size. */
  const size_t Nx;
  const size_t Ny;

  /* Epsilon - the dielectric constant */
  /* Array */
  float *eps;

  /* Magnetic field strength. */
  double *H1;
  double *H2;

  /* Electric field strength */
  double *E1;
  double *E2;

  /* Refractive index */
  const float n1;
  const float lambda;
  const float tau;

} DATA_STRUCT;

void calculate_next_time_layer(DATA_STRUCT *data, double *vector_x,
                               double *vector_y);

void set_params(DATA_STRUCT *data, float lambda, float tau, float n1);

// Getters.
size_t getNx();
float getLambda();
float getTau();
float getN1();

// Setters.
void setLambda(DATA_STRUCT *data, float new_lambda) ;
void setTau(DATA_STRUCT *data, float new_tau) ;
void setN1(DATA_STRUCT *data, float new_n1) ;

void calculations(DATA_STRUCT *data);

// Moor`s 1-st boundary condition
void boundary_conditions_1(DATA_STRUCT *data);

// Moor`s 2-nd boundary condition
void boundary_conditions_2(DATA_STRUCT *data);

#endif