#include "./FDTD-2D.h"


 void set_params(DATA_STRUCT *data, float lambda, float tau, float n1) {
  /* Input conditions */


  /* Initializing const refreactive index. */
  *(float *)&data->n1 = n1;
  *(float *)&data->lambda = lambda;
  *(float *)&data->tau = tau;

  data->ticks = 0;

  /* Grid steps. */
  data->dx = 0.05;
  data->dt = 0.025;
  /* Physics params. */
  data->aa1 = lambda * lambda / (0.09 * tau * tau);
  data->tMax = 4 * tau / (lambda / 0.3);

  /* Initializing const grid size. */
  *(size_t *)&data->Nx = 2000;  // 2000
  *(size_t *)&data->Ny = 500;   // 500

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
      sin(2 * M_PI * data->dt * data->ticks);

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

