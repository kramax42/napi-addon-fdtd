#pragma once

#include <math.h>
#include <napi.h>

#include <vector>
#include <algorithm> // std::fill

class FdtdPml1D
{

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
    size_t src_position = 100;

    // 1D grid size.
    static const size_t grid_size = 500;

    // Width of PML layer.
    const size_t pml_width = 50;

    // Permitivity array.
    // double eps[grid_size];
    std::array<double, grid_size> eps;

    // Permeability array.
    // double mu[grid_size];
    std::array<double, grid_size> mu;

    // Conductivity array.
    // double sigma[grid_size];
    // double sigma_star[grid_size];
    std::array<double, grid_size> sigma;
    std::array<double, grid_size> sigma_star;

    // Optimal polynomial order for grading sigma array (pp 292, Taflove).
    double m = 3;

    // Impedance.
    const double eta = std::sqrt(mu0 / eps0);

    // Required reflection factor.
    const double R = 1e-8;

    // Taflove, pp 292, Eq 7.61.
    const double sigma_max = -(m + 1) * std::log(R) / (2 * eta * pml_width * dx);

    // // PML constants.
    // double A[grid_size];
    // double B[grid_size];
    // double C[grid_size];
    // double D[grid_size];

    // // Electromagnetic field projections in space array.
    // double hy[grid_size];
    // double ex[grid_size];



    // PML constants.
    std::array<double, grid_size> A;
    std::array<double, grid_size> B;
    std::array<double, grid_size> C;
    std::array<double, grid_size> D;

    // Electromagnetic field projections in space array.
    std::array<double, grid_size> hy;
    std::array<double, grid_size> ex;

    

    size_t time_step = 0;

public:
    void SetParams(
        double new_tau,
        double new_omega,
        std::vector<double> &new_eps,
        std::vector<double> &new_mu,
        std::vector<double> &new_sigma,
        size_t new_src_position)
    {
        // tau = new_tau;
        // omega = new_omega;
        src_position = new_src_position;

        // std::fill(&eps, grid_size, eps0);
        // std::fill(&mu, grid_size, mu0);
        // std::fill(&sigma, grid_size, 0.0);

        // // Electromagnetic field projections in space array.
        // std::fill(&hy, grid_size, 0.0);
        // std::fill(&ex, grid_size, 0.0);

        for (int i = 0; i < grid_size; ++i)
        {
            // eps[i] = new_eps[i];
            // mu[i] = new_mu[i];
            // sigma[i] = new_sigma[i];
            eps[i] = eps0;
            mu[i] = mu0;
            sigma[i] = 0.0;

            // Electromagnetic field projections in space array.
            hy[i] = 0;
            ex[i] = 0;
        }

        // Taflove, pp 292, Eq 7.60a.
        for (int i = 0; i < pml_width; ++i)
        {
            double sigma_in_pml = std::pow( (double)i / pml_width, m) * sigma_max;

            // Lossy electric conductivity profile.
            sigma[grid_size - pml_width + i] = sigma_in_pml;
            sigma[pml_width - i - 1] = sigma_in_pml;
        }

        for (int i = 0; i < grid_size; ++i)
        {

            // Eq 7.8 Taflove, pp 275
            // Magnetic conductivity loss.
            sigma_star[i] = sigma[i] * mu0 / eps0;

            // PML constants.
            A[i] = ((mu[i] - 0.5 * dt * sigma_star[i]) / (mu[i] + 0.5 * dt * sigma_star[i]));

            B[i] = (dt / dx) / (mu[i] + 0.5 * dt * sigma_star[i]);

            C[i] = ((eps[i] - 0.5 * dt * sigma[i]) / (eps[i] + 0.5 * dt * sigma[i]));

            D[i] = (dt / dx) / (eps[i] + 0.5 * dt * sigma[i]);
        }

        // printArr(sigma);
        std::cout << src_position;
        // std::cout << (2 * eta * pml_width * dx);

    }

    // Electric field is an Gaussian envelop.

    size_t GetCurrentTimeStep()
    {
        return time_step;
    }

    double GetTau()
    {
        return tau;
    }
    size_t GetSourcePosition()
    {
        return src_position;
    }

    static size_t GetGridSize()
    {
        return grid_size;
    }

    void Calculation(std::vector<double> &vect_x,
                     std::vector<double> &vect_ex,
                     std::vector<double> &vect_hy)
    {
        

        // Insert source in certain space grid.
        double t = (double)time_step * dt;
        double source = std::exp(-(std::pow((t0 - t) / tau, 2))) * std::cos(omega * t);
        ex[src_position] += source;

        for (int i = 0; i < grid_size - 1; ++i)
        {
            hy[i] = A[i] * hy[i] - B[i] * (ex[i + 1] - ex[i]);
        }

        for (int i = 1; i < grid_size - 1; ++i)
        {
            ex[i] = C[i] * ex[i] - D[i] * (hy[i] - hy[i - 1]);
        }
        // ex[grid_size - 1] = ex[grid_size - 2];
        // ex[grid_size-1] = ex[grid_size-2];

        for (int i = 0; i < grid_size; ++i)
        {
            vect_x.push_back(i);
            vect_ex.push_back(ex[i]);
            vect_hy.push_back(hy[i]);
        }

        // ex[grid_size] = ex[grid_size - 1];
        ++time_step;
    }

    FdtdPml1D(double new_tau,
              double new_omega,
              std::vector<double> &new_eps,
              std::vector<double> &new_mu,
              std::vector<double> &new_sigma,
              size_t new_src_position)
    {
        SetParams(new_tau,
                  new_omega,
                  new_eps,
                  new_mu,
                  new_sigma,
                  new_src_position);
    }
};