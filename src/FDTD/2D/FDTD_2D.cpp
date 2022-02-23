
// // http://zfdtd.narod.ru/ --- tutorial FDTD.

// #include "FDTD_2D.h"

// #include <iostream>

// FDTD_2D::FDTD_2D(float lambda, float tau, float refractive_index)
//     : lambda(lambda), tau(tau), refractive_index(refractive_index) {
//   SetParams();
// }

// // Setters.
// void FDTD_2D::SetLambda(float new_lambda) { lambda = new_lambda; }
// void FDTD_2D::SetTau(float new_tau) { tau = new_tau; }
// void FDTD_2D::SetRefractiveIndex(float new_refractive_index) {
//   refractive_index = new_refractive_index;
// }

// // Getters.
// size_t FDTD_2D::GetNx() { return Nx; }

// float FDTD_2D::GetTau() { return tau; }

// float FDTD_2D::GetLambda() { return lambda; }

// float FDTD_2D::GetRefractiveIndex() { return refractive_index; }

// // double lambda = 1,
// // double tau = 10,
// // double refractive_index = 1
// void FDTD_2D::SetParams() {
//   ticks = 0;

//   // Grid steps.
//   dx = 0.05;
//   dt = 0.025;

//   // Physics params.
//   aa1 = lambda * lambda / (0.09 * tau * tau);
//   tMax = 4 * tau / (lambda / 0.3);

//   for (int i = 0; i < Nx; i++) {
//     E1[i] = 1e-8;
//     E2[i] = 1e-8;
//     H1[i] = 1e-8;
//     H2[i] = 1e-8;
//     eps[i] = refractive_index;
//   }
// }

// // Moor`s boundary condition.
// void FDTD_2D::BoundaryConditionsFirst() {
//   H2[0] = H1[1] + (dt / eps[1] - dx) / (dt / eps[1] + dx) * (H2[1] - H1[0]);
//   E2[0] = E1[1] + (dt / eps[1] - dx) / (dt / eps[1] + dx) * (E2[1] - E1[0]);
// }

// // Moor`s boundary condition
// void FDTD_2D::BoundaryConditionsSecond() {
//   H2[Nx - 1] = H1[Nx - 2] + (dt / eps[Nx - 2] - dx) *
//                                 (H2[Nx - 2] - H1[Nx - 1]) /
//                                 (dt / eps[Nx - 2] + dx);

//   E2[Nx - 1] = E1[Nx - 2] + (dt / eps[Nx - 2] - dx) *
//                                 (E2[Nx - 2] - E1[Nx - 1]) /
//                                 (dt / eps[Nx - 2] + dx);
// }

// void FDTD_2D::Calculation() {
//   for (int i = 1; i <= Ny - 2; i++) {
//     H2[i] = H1[i] * dt / dx - (E1[i] - E1[i - 1]);
//     E2[i - 1] = E1[i - 1] - (H2[i] - H2[i - 1]) * dt / (eps[i - 1] * dx);
//   }

//   E1[Ny - 1] = std::exp(aa1 * (tMax - dt * ticks) * (dt * ticks - tMax)) *
//                std::sin(2 * PI * dt * ticks);
//   H1[Ny - 1] = eps[Ny - 1] * E1[Ny - 1];

//   for (int i = Ny; i < Nx; i++) {
//     H2[i] = H1[i] - (E1[i] - E1[i - 1]) * dt / dx;
//     E2[i - 1] = E1[i - 1] - (H2[i] - H2[i - 1]) * dt / (eps[i - 1] * dx);
//   }
// }

// size_t FDTD_2D::GetCurrentTick() { return ticks; }

// void FDTD_2D::CalcNextLayer(std::vector<double> &vectX,
//                             std::vector<double> &vectY) {
//   Calculation();
//   BoundaryConditionsFirst();
//   BoundaryConditionsSecond();

//   for (int i = 0; i < Nx; i++) {
//     H1[i] = H2[i];
//     E1[i] = E2[i];
//     vectX.push_back(dx * lambda * (i - 1));
//     vectY.push_back(E1[i]);
//   }

//   ticks++;
// }



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
        E1[i] = 1e-8;
        E2[i] = 1e-8;
        H1[i] = 1e-8;
        H2[i] = 1e-8;
        eps[i] = refractive_index;
    }
}

// Moor`s boundary condition.
void FDTD_2D::BoundaryConditionsFirst()
{
    H2[0] = H1[1] + (dt / eps[1] - dx) / (dt / eps[1] + dx) * (H2[1] - H1[0]);

    E2[0] = E1[1] + (dt / eps[1] - dx) / (dt / eps[1] + dx) * (E2[1] - E1[0]);
}

//Moor`s boundary condition
void FDTD_2D::BoundaryConditionsSecond()
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
        std::exp(aa1 * (tMax - dt * ticks) * (dt * ticks - tMax)) * std::sin(2 * 3.1416f * dt * ticks);

    H1[Ny - 1] = eps[Ny - 1] * E1[Ny - 1];

    for (int i = Ny; i < Nx; i++)
    {
        H2[i] =
            H1[i] - (E1[i] - E1[i - 1]) * dt / dx;

        E2[i - 1] =
            E1[i - 1] - (H2[i] - H2[i - 1]) * dt / (eps[i - 1] * dx);
    }
}

size_t FDTD_2D::GetCurrentTick()
{
    return ticks;
}

void FDTD_2D::CalcNextLayer( std::vector<double> &vectX,
                        std::vector<double> &vectY)
{
    Calculation();
    BoundaryConditionsFirst();
    BoundaryConditionsSecond();

    for (int i = 0; i < Nx; i++)
    {
        H1[i] = H2[i];
        E1[i] = E2[i];
        vectX.push_back(dx * lambda * (i - 1));
        vectY.push_back(E1[i]);
    }

    ticks++;
}

