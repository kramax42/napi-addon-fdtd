//http://zfdtd.narod.ru/ -- ������ ������ FDTD

#include "FDTD_2D.h"
#include <iostream>


FDTD_2D::FDTD_2D(double lambda, double tau, double refractive_index)
                : lambda(lambda), tau(tau), refractive_index(refractive_index)
{
    setParams();
}

//Getters
size_t FDTD_2D::GetNx() {
     return Nx;
}

double FDTD_2D::GetTau(){
    return tau;
}

double FDTD_2D::GetLambda()
{
    return lambda;
}

double FDTD_2D::GetRefractiveIndex()
{
    return refractive_index;
}

//double lambda = 1, double tau = 10, double refractive_index = 1
void FDTD_2D::setParams()
{
    ticks = 0;

    // Grid steps.
    // dx = 1 / 20
    dx = 0.05;
    dt = dx / 2;


    // Physics params.
    aa1 = lambda * lambda / (c0 * c0 * tau * tau);
    tMax = 3 * tau / (lambda / c0);

    for (int i = 0; i < Nx; i++)
    {
        Ex_prev[i] = 1e-8;
        Ex[i] = 1e-8;
        Hy_prev[i] = 1e-8;
        Hy[i] = 1e-8;
        eps[i] = refractive_index;
    }
}

// Moor`s boundary condition.
void FDTD_2D::BoundaryConditionsFirst()
{
    Hy[0] = Hy_prev[1] + (dt / eps[1] - dx) / (dt / eps[1] + dx) * (Hy[1] - Hy_prev[0]);

    Ex[0] = Ex_prev[1] + (dt / eps[1] - dx) / (dt / eps[1] + dx) * (Ex[1] - Ex_prev[0]);
}

//Moor`s boundary condition
void FDTD_2D::BoundaryConditionsSecond()
{
    Hy[Nx - 1] =
        Hy_prev[Nx - 2] + (dt / eps[Nx - 2] - dx) * (Hy[Nx - 2] - Hy_prev[Nx - 1]) /
        (dt / eps[Nx - 2] + dx);

    Ex[Nx - 1] =
        Ex_prev[Nx - 2] + (dt / eps[Nx - 2] - dx) * (Ex[Nx - 2] - Ex_prev[Nx - 1]) /
        (dt / eps[Nx - 2] + dx);
}

void FDTD_2D::Calculation()
{
    for (int i = 1; i <= Ny - 2; i++)
    {
        Hy[i] =
            Hy_prev[i] * dt / dx - (Ex_prev[i] - Ex_prev[i - 1]);

        Ex[i - 1] =
            Ex_prev[i - 1] - (Hy[i] - Hy[i - 1]) * dt / (eps[i - 1] * dx);
    }

    Ex_prev[Ny - 1] =
        std::exp(aa1 * (tMax - dt * ticks) * (dt * ticks - tMax)) * std::sin(2 * PI * dt * ticks);

    Hy_prev[Ny - 1] = eps[Ny - 1] * Ex_prev[Ny - 1];

    for (int i = Ny; i < Nx; i++)
    {
        Hy[i] =
            Hy_prev[i] - (Ex_prev[i] - Ex_prev[i - 1]) * dt / dx;

        Ex[i - 1] =
            Ex_prev[i - 1] - (Hy[i] - Hy[i - 1]) * dt / (eps[i - 1] * dx);
    }
}

size_t FDTD_2D::GetCurrentTick()
{
    return ticks;
}

void FDTD_2D::CalcNextLayer( std::vector<double> &vectX,
                        std::vector<double> &vectEx,
                        std::vector<double> &vectHy)
{
    Calculation();
    BoundaryConditionsFirst();
    BoundaryConditionsSecond();

    for (int i = 0; i < Nx; i++)
    {
        Hy_prev[i] = Hy[i];
        Ex_prev[i] = Ex[i];
        vectX.push_back(dx * lambda * (i - 1));
        vectEx.push_back(Ex_prev[i]);
        vectHy.push_back(Hy_prev[i]);
    }

    ticks++;
}

