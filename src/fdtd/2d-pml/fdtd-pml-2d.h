// 2D FDTD with PML.

// Based on code from book:
// Electromagnetic simulation using the fdtd method with python.
// Chapter 3.

#ifndef FDTD_PML_2D
#define FDTD_PML_2D

#include <math.h>
#include <napi.h>

#include <algorithm>  // std::fill_n
#include <vector>

class FdtdPml2D {
    // Grid sizes.
    static const size_t pml_width = 40;

    static const size_t rows = 220 + pml_width * 2;
    static const size_t cols = 220 + pml_width * 2;

    // Set source position.
    size_t src_row = rows / 2;
    size_t src_col = cols / 2;

    // Light speed.
    const double light_spd = 2.99792458e8;

    const double pi = 3.1415926535897932385;

    // Permeability of free space.
    const double mu0 = 4.0 * pi * 1.0e-7;

    // Permittivity of free space.
    const double epsz = 8.8e-12;

    // Permittivity.
    double epsilon0 = 1;
    double epsilon1 = 1;
    double epsilon2 = 4.2;

    // Conductivity.
    double sigma0 = 5000;
    double sigma1 = 0.00001;
    double sigma2 = 0.001;

    // Field arrays.
    double dz[rows][cols];
    double ez[rows][cols];
    double hx[rows][cols];
    double hy[rows][cols];

    double gaz[rows][cols];

    // Space grid step.
    const double ddx = 1.0e-3;
    const double ddy = ddx;

    // Courant factor.
    const double cfl_factor = 0.99;

    // Time step corressponding Courant factor.
    const double dt = cfl_factor / (light_spd * sqrt(std::pow(1 / ddx, 2) + std::pow(1 / ddy, 2)));

    size_t time_step = 0;

    // Source params.
    // Gaussian beam.
    const size_t t0 = 20;

    // Beam width.
    const size_t tau = 20;

    // PML params.
    double gi2[rows];
    double gi3[rows];
    double fi1[rows];
    double fi2[rows];
    double fi3[rows];

    double gj2[rows];
    double gj3[rows];
    double fj1[rows];
    double fj2[rows];
    double fj3[rows];

   public:
    // Getters.
    static size_t GetRows() {
        return rows - pml_width * 2;
    }

    static size_t GetCols() {
        return cols - pml_width * 2;
    }

    size_t GetStep();
    size_t GetCurrentTimeStep();

    void SetParams(std::vector<std::vector<double>> &eps,
                   std::vector<std::vector<double>> &mu,
                   std::vector<std::vector<double>> &sigma,
                   size_t new_src_position_row,
                   size_t new_src_position_col);

    void Calculation();

    FdtdPml2D();

    FdtdPml2D(std::vector<std::vector<double>> &eps,
              std::vector<std::vector<double>> &mu,
              std::vector<std::vector<double>> &sigma,
              size_t new_src_position_row,
              size_t new_src_position_col);

    void CalcNextLayer(std::vector<double> &vectX,
                       std::vector<double> &vectY,
                       std::vector<double> &vectEz,
                       std::vector<double> &vectHy,
                       std::vector<double> &vectHx,
                       std::vector<double> &vectEnergy,
                       double &max, double &min);
};

#endif