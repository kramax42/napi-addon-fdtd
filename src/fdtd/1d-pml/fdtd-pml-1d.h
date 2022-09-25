#ifndef FDTD_PML_1D
#define FDTD_PML_1D

#include <math.h>
#include <napi.h>

#include <algorithm>  // std::fill
#include <iostream>
#include <vector>

class FdtdPml1D {
    // Speed of light in nm/fs.
    const double light_spd = 299.792458;

    // Courant factor.
    const double cfl_factor = 0.99;

    // Space step in nanometres.
    const double dx = 40.0;

    // Time step in femto seconds.
    const double dt = dx * cfl_factor / light_spd;

    //
    const double t0 = 20.0;

    // Permeability of free space in (V fs^2/e nm).
    const double mu0 = 2.013354451e-4;

    // Permittivity of free space in (e / V nm).
    const double eps0 = 55.26349597e-3;

    // Planck's constant.
    const double h = 0.6582119514;

    // Frequency of light in eV (electron-volts).
    const double omega_ev = 1.5;

    // Frequency of light in (1/fs).
    double omega = omega_ev / h;

    // Width of electric beam.
    double tau = 8.0;

    // Position index of source.
    size_t src_position = 5;

    // Width of PML layer.
    static const size_t pml_width = 50;

    // 1D grid size.
    // static const size_t grid_size = 500 + pml_width * 2;
    static const size_t grid_size = 500;

    // Permitivity array.
    std::array<double, grid_size> eps;

    // Permeability array.
    std::array<double, grid_size> mu;

    // Conductivity array.
    std::array<double, grid_size> sigma;
    std::array<double, grid_size> sigma_star;

    // Optimal polynomial order for grading sigma array (pp 292, Taflove).
    double m = 3;

    // Impedance.
    const double eta = std::sqrt(mu0 / eps0);

    // Required reflection factor.
    const double R = 1e-8;

    // Taflove, page 292, equation 7.61.
    const double sigma_max = -(m + 1) * std::log(R) / (2 * eta * pml_width * dx);

    // PML constants arrays.
    std::array<double, grid_size> A;
    std::array<double, grid_size> B;
    std::array<double, grid_size> C;
    std::array<double, grid_size> D;

    // Electromagnetic field projections in space arrays.
    std::array<double, grid_size> hy;
    std::array<double, grid_size> ex;

    size_t time_step = 0;

   public:
    void SetParams(double new_tau,
                   double new_omega,
                   std::vector<double> &new_eps,
                   std::vector<double> &new_mu,
                   std::vector<double> &new_sigma,
                   size_t new_src_position);

    // Getters.
    size_t GetCurrentTimeStep();
    double GetTau();
    size_t GetSourcePosition();
    static size_t GetGridSize() {
        // return grid_size - pml_width * 2;
        return grid_size;
    }

    void Calculation(std::vector<double> &vect_x,
                     std::vector<double> &vect_ex,
                     std::vector<double> &vect_hy);

    FdtdPml1D();

    FdtdPml1D(double new_tau,
              double new_omega,
              std::vector<double> &new_eps,
              std::vector<double> &new_mu,
              std::vector<double> &new_sigma,
              size_t new_src_position);
};

#endif