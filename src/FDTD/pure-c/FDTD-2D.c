#include <math.h>    // exp , pi
#include <stdio.h>   // printf
#include <stdlib.h>  // size_t

// gcc FDTD-2D.c -o output -lm

const double PI = M_PI;

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

} DATA_STRUCT;

void set_params(DATA_STRUCT *data) {
  /* Input conditions */
  float lambda = 1;
  float tau = 3;
  float n1 = 1.6;

  /* Initializing const refreactive index. */
  *(float *)&data->n1 = n1;
  *(float *)&data->lambda = lambda;

  data->ticks = 0;

  /* Grid steps. */
  data->dx = 0.05;
  data->dt = 0.025;
  /* Physics params. */
  data->aa1 = lambda * lambda / (0.09 * tau * tau);
  data->tMax = 4 * tau / (lambda / 0.3);

  /* Initializing const grid size. */
  *(size_t *)&data->Nx = 20;  // 2000
  *(size_t *)&data->Ny = 5;   // 500

  data->eps = (float *)malloc(data->Nx * sizeof(float));
  data->H1 = (double *)malloc(data->Nx * sizeof(double));
  data->H2 = (double *)malloc(data->Nx * sizeof(double));
  data->E1 = (double *)malloc(data->Nx * sizeof(double));
  data->E2 = (double *)malloc(data->Nx * sizeof(double));

  for (size_t i = 0; i < data->Nx; ++i) {
    data->E1[i] = 1e-8;
    data->E2[i] = 1e-8;
    data->H1[i] = 1e-8;
    data->H2[i] = 1e-8;
    data->eps[i] = data->n1;
  }
}

void show_data(double *data, size_t size) {
  for (size_t i = 0; i < size; ++i) {
    printf("%0.2f ", data[i]);
  }
}

/*------------------------------------------*/
void calculations(DATA_STRUCT *data) {
  for (int i = 1; i <= data->Ny - 2; i++) {
    data->H2[i] =
        data->H1[i] * data->dt / data->dx - (data->E1[i] - data->E1[i - 1]);

    data->E2[i - 1] = data->E1[i - 1] - (data->H2[i] - data->H2[i - 1]) *
                                            data->dt /
                                            (data->eps[i - 1] * data->dx);
  }

  data->E1[data->Ny - 1] =
      exp(data->aa1 * (data->tMax - data->dt * data->ticks) *
          (data->dt * data->ticks - data->tMax)) *
      sin(2 * PI * data->dt * data->ticks);

  data->H1[data->Ny - 1] = data->eps[data->Ny - 1] * data->E1[data->Ny - 1];

  for (int i = data->Ny; i < data->Nx; i++) {
    data->H2[i] =
        data->H1[i] - (data->E1[i] - data->E1[i - 1]) * data->dt / data->dx;

    data->E2[i - 1] = data->E1[i - 1] - (data->H2[i] - data->H2[i - 1]) *
                                            data->dt /
                                            (data->eps[i - 1] * data->dx);
  }
}

// Moor`s boundary condition.
void boundary_conditions_1(DATA_STRUCT *data) {
  data->H2[0] = data->H1[1] + (data->dt / data->eps[1] - data->dx) /
                                  (data->dt / data->eps[1] + data->dx) *
                                  (data->H2[1] - data->H1[0]);

  data->E2[0] = data->E1[1] + (data->dt / data->eps[1] - data->dx) /
                                  (data->dt / data->eps[1] + data->dx) *
                                  (data->E2[1] - data->E1[0]);
}

// Moor`s boundary condition
void boundary_conditions_2(DATA_STRUCT *data) {
  data->H2[data->Nx - 1] =
      data->H1[data->Nx - 2] +
      (data->dt / data->eps[data->Nx - 2] - data->dx) *
          (data->H2[data->Nx - 2] - data->H1[data->Nx - 1]) /
          (data->dt / data->eps[data->Nx - 2] + data->dx);

  data->E2[data->Nx - 1] =
      data->E1[data->Nx - 2] +
      (data->dt / data->eps[data->Nx - 2] - data->dx) *
          (data->E2[data->Nx - 2] - data->E1[data->Nx - 1]) /
          (data->dt / data->eps[data->Nx - 2] + data->dx);
}

void calculate_next_time_layer(DATA_STRUCT *data, double *vector_x,
                               double *vector_y) {
  calculations(data);

  boundary_conditions_1(data);
  boundary_conditions_2(data);

  for (int i = 0; i < data->Nx; i++) {
    data->H1[i] = data->H2[i];
    data->E1[i] = data->E2[i];
    vector_x[i] = data->dx * data->lambda * (i - 1);
    vector_y[i] = data->E1[i];
  }

  data->ticks++;
}

int main(int argc, char *argv[]) {
  DATA_STRUCT *data_struct = (DATA_STRUCT *)malloc(sizeof(DATA_STRUCT));

  set_params(data_struct);

  // calc time layers
  double *vect_x = (double *)malloc(data_struct->Nx * sizeof(double));
  double *vect_y = (double *)malloc(data_struct->Nx * sizeof(double));

  for (int i = 0; i <= 100; ++i) {
    calculate_next_time_layer(data_struct, vect_x, vect_y);
  }

  show_data(vect_y, data_struct->Nx);
  // unsigned decimal
  printf("%zu", data_struct->ticks);

  /* Free memory */
  free(vect_x);
  free(vect_y);

  free(data_struct->eps);
  free(data_struct->H1);
  free(data_struct->H2);
  free(data_struct->E1);
  free(data_struct->E2);
  free(data_struct);

  return 0;
}