
//http://zfdtd.narod.ru/ -- ������ ������ FDTD

#include "FDTD_2D.h"
#include <iostream>


FDTD_2D::FDTD_2D(double lambda, double tau, double n1)
                : lambda(lambda), tau(tau), n1(n1)
{
    setParams();
}

//getters
size_t FDTD_2D::getNx() {
     return Nx;
}

double FDTD_2D::getTau(){
    return tau;
}

double FDTD_2D::getLambda()
{
    return lambda;
}

double FDTD_2D::getN1()
{
    return n1;
}

//double lambda = 1, double tau = 10, double n1 = 1
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
        E1[i] = 1e-8;
        E2[i] = 1e-8;
        H1[i] = 1e-8;
        H2[i] = 1e-8;
        eps[i] = n1;
    }
}

// Moor`s boundary condition.
void FDTD_2D::boundary_conditions_1()
{
    H2[0] = H1[1] + (dt / eps[1] - dx) / (dt / eps[1] + dx) * (H2[1] - H1[0]);

    E2[0] = E1[1] + (dt / eps[1] - dx) / (dt / eps[1] + dx) * (E2[1] - E1[0]);
}

//Moor`s boundary condition
void FDTD_2D::boundary_conditions_2()
{
    H2[Nx - 1] =
        H1[Nx - 2] + (dt / eps[Nx - 2] - dx) * (H2[Nx - 2] - H1[Nx - 1]) /
        (dt / eps[Nx - 2] + dx);

    E2[Nx - 1] =
        E1[Nx - 2] + (dt / eps[Nx - 2] - dx) * (E2[Nx - 2] - E1[Nx - 1]) /
        (dt / eps[Nx - 2] + dx);
}

void FDTD_2D::Calculation()
{
    for (int i = 1; i <= Ny - 2; i++)
    {
        H2[i] =
            H1[i] * dt / dx - (E1[i] - E1[i - 1]);

        E2[i - 1] =
            E1[i - 1] - (H2[i] - H2[i - 1]) * dt / (eps[i - 1] * dx);
    }

    E1[Ny - 1] =
        std::exp(aa1 * (tMax - dt * ticks) * (dt * ticks - tMax)) * std::sin(2 * PI * dt * ticks);

    H1[Ny - 1] = eps[Ny - 1] * E1[Ny - 1];

    for (int i = Ny; i < Nx; i++)
    {

        H2[i] =
            H1[i] - (E1[i] - E1[i - 1]) * dt / dx;

        E2[i - 1] =
            E1[i - 1] - (H2[i] - H2[i - 1]) * dt / (eps[i - 1] * dx);
    }
}

size_t FDTD_2D::getCurrentTick()
{
    return ticks;
}

void FDTD_2D::calcNextLayer( vector<double> &vectX,
                             vector<double> &vectY )
{
    Calculation();
    boundary_conditions_1();
    boundary_conditions_2();

    for (int i = 0; i < Nx; i++)
    {
        H1[i] = H2[i];
        E1[i] = E2[i];
        vectX.push_back(dx * lambda * (i - 1));
        vectY.push_back(E1[i]);
    }

    ticks++;
}

