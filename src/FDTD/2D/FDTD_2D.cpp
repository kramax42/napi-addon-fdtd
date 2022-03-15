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
    dx = 0.05;
    dt = 0.025;


    // Physics params.
    aa1 = lambda * lambda / (0.09 * tau * tau);
    tMax = 4 * tau / (lambda / 0.3);

    for (int i = 0; i < Nx; i++)
    {
        Ex1[i] = 1e-8;
        Ex2[i] = 1e-8;
        Hy1[i] = 1e-8;
        Hy2[i] = 1e-8;
        eps[i] = refractive_index;
    }
}

// Moor`s boundary condition.
void FDTD_2D::BoundaryConditionsFirst()
{
    Hy2[0] = Hy1[1] + (dt / eps[1] - dx) / (dt / eps[1] + dx) * (Hy2[1] - Hy1[0]);

    Ex2[0] = Ex1[1] + (dt / eps[1] - dx) / (dt / eps[1] + dx) * (Ex2[1] - Ex1[0]);
}

//Moor`s boundary condition
void FDTD_2D::BoundaryConditionsSecond()
{
    Hy2[Nx - 1] =
        Hy1[Nx - 2] + (dt / eps[Nx - 2] - dx) * (Hy2[Nx - 2] - Hy1[Nx - 1]) /
        (dt / eps[Nx - 2] + dx);

    Ex2[Nx - 1] =
        Ex1[Nx - 2] + (dt / eps[Nx - 2] - dx) * (Ex2[Nx - 2] - Ex1[Nx - 1]) /
        (dt / eps[Nx - 2] + dx);
}

void FDTD_2D::Calculation()
{
    for (int i = 1; i <= Ny - 2; i++)
    {
        Hy2[i] =
            Hy1[i] * dt / dx - (Ex1[i] - Ex1[i - 1]);

        Ex2[i - 1] =
            Ex1[i - 1] - (Hy2[i] - Hy2[i - 1]) * dt / (eps[i - 1] * dx);
    }

    Ex1[Ny - 1] =
        std::exp(aa1 * (tMax - dt * ticks) * (dt * ticks - tMax)) * std::sin(2 * PI * dt * ticks);

    Hy1[Ny - 1] = eps[Ny - 1] * Ex1[Ny - 1];

    for (int i = Ny; i < Nx; i++)
    {
        Hy2[i] =
            Hy1[i] - (Ex1[i] - Ex1[i - 1]) * dt / dx;

        Ex2[i - 1] =
            Ex1[i - 1] - (Hy2[i] - Hy2[i - 1]) * dt / (eps[i - 1] * dx);
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
        Hy1[i] = Hy2[i];
        Ex1[i] = Ex2[i];
        vectX.push_back(dx * lambda * (i - 1));
        vectEx.push_back(Ex1[i]);
        vectHy.push_back(Hy1[i]);
    }

    ticks++;
}

