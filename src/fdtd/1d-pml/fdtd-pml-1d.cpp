#include "fdtd-pml-1d.h"

void FdtdPml1D::SetParams(double new_tau,
                          double new_omega,
                          std::vector<double> &new_eps,
                          std::vector<double> &new_mu,
                          std::vector<double> &new_sigma,
                          size_t new_src_position)
{

    time_step = 0;

    tau = new_tau;
    omega = new_omega;
    src_position = new_src_position;

    for (int i = 0; i < grid_size; ++i)
    {
        eps[i] = new_eps[i] * eps0;
        mu[i] = new_mu[i] * mu0;
        sigma[i] = new_sigma[i];

        // Electromagnetic field projections in space array.
        hy[i] = 0;
        ex[i] = 0;
    }

    // Taflove, page 292, equation 7.60a.
    for (int i = 0; i < pml_width; ++i)
    {

        double sigma_in_pml = std::pow((double)i / pml_width, m) * sigma_max;

        // Lossy electric conductivity profile.
        sigma[grid_size - pml_width + i] = sigma_in_pml;
        sigma[pml_width - i - 1] = sigma_in_pml;
    }

    for (int i = 0; i < grid_size; ++i)
    {

        // Taflove, page 275, equation 7.8.
        // Magnetic conductivity loss.
        sigma_star[i] = sigma[i] * mu0 / eps0;

        // PML constants.
        A[i] = ((mu[i] - 0.5 * dt * sigma_star[i]) / (mu[i] + 0.5 * dt * sigma_star[i]));

        B[i] = (dt / dx) / (mu[i] + 0.5 * dt * sigma_star[i]);

        C[i] = ((eps[i] - 0.5 * dt * sigma[i]) / (eps[i] + 0.5 * dt * sigma[i]));

        D[i] = (dt / dx) / (eps[i] + 0.5 * dt * sigma[i]);
    }
}

size_t FdtdPml1D::GetCurrentTimeStep()
{
    return time_step;
}

double FdtdPml1D::GetTau()
{
    return tau;
}

size_t FdtdPml1D::GetSourcePosition()
{
    return src_position;
}

void FdtdPml1D::Calculation(std::vector<double> &vect_x,
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

    for (int i = 0; i < grid_size; ++i)
    {
        vect_x.push_back(i);
        vect_ex.push_back(ex[i]);
        vect_hy.push_back(hy[i]);
    }

    ++time_step;
}

FdtdPml1D::FdtdPml1D(double new_tau,
                     double new_omega,
                     std::vector<double> &new_eps,
                     std::vector<double> &new_mu,
                     std::vector<double> &new_sigma,
                     size_t new_src_position)
{

    SetParams(new_tau, new_omega, new_eps,
              new_mu, new_sigma, new_src_position);
}