//http://zfdtd.narod.ru/ -- ������ ������ FDTD

#include "FDTD_2D_UPDATED.h"
#include <iostream>


FDTD_2D_UPDATED::FDTD_2D_UPDATED(double lambda, double tau, std::vector<double>& Epsilon)
                : lambda(lambda), tau(tau)
{
    setParams(Epsilon);
}

//Getters
size_t FDTD_2D_UPDATED::GetNx() {
     return jmax;
}

double FDTD_2D_UPDATED::GetTau(){
    return tau;
}

double FDTD_2D_UPDATED::GetLambda()
{
    return lambda;
}

// double FDTD_2D_UPDATED::GetRefractiveIndex()
// {
//     return refractive_index;
// }

//double lambda = 1, double tau = 10, double refractive_index = 1
void FDTD_2D_UPDATED::setParams(std::vector<double>& Epsilon)
{
    time_step = 0;

    // Grid steps.
    // dx = 0.05;
    // dt = 0.025;


    // Physics params.
    // aa1 = lambda * lambda / (0.09 * tau * tau);
    // tMax = 4 * tau / (lambda / 0.3);

    for (int i = 0; i < jmax; i++)
    {
        Ex[i] = 0;
        Ex_prev[i] = 0;

        Hy[i] = 0;
        Hy_prev[i] = 0;

        // eps[i] = refractive_index;
        eps[i] = Epsilon[i];
    }

    // for (int i = 200; i < 300; i++){
    //     eps[i] = 4 * eps0;
    // }

    
}

// // Moor`s boundary condition.
// void FDTD_2D_UPDATED::BoundaryConditionsFirst()
// {
//     Hy2[0] = Hy1[1] + (dt / eps[1] - dx) / (dt / eps[1] + dx) * (Hy2[1] - Hy1[0]);

//     Ex2[0] = Ex1[1] + (dt / eps[1] - dx) / (dt / eps[1] + dx) * (Ex2[1] - Ex1[0]);
// }

// //Moor`s boundary condition
// void FDTD_2D_UPDATED::BoundaryConditionsSecond()
// {
//     Hy2[Nx - 1] =
//         Hy1[Nx - 2] + (dt / eps[Nx - 2] - dx) * (Hy2[Nx - 2] - Hy1[Nx - 1]) /
//         (dt / eps[Nx - 2] + dx);

//     Ex2[Nx - 1] =
//         Ex1[Nx - 2] + (dt / eps[Nx - 2] - dx) * (Ex2[Nx - 2] - Ex1[Nx - 1]) /
//         (dt / eps[Nx - 2] + dx);
// }

// void FDTD_2D_UPDATED::Calculation()
// {
//     for (int i = 1; i <= Ny - 2; i++)
//     {
//         Hy2[i] =
//             Hy1[i] * dt / dx - (Ex1[i] - Ex1[i - 1]);

//         Ex2[i - 1] =
//             Ex1[i - 1] - (Hy2[i] - Hy2[i - 1]) * dt / (eps[i - 1] * dx);
//     }

//     Ex1[Ny - 1] =
//         std::exp(aa1 * (tMax - dt * ticks) * (dt * ticks - tMax)) * std::sin(2 * PI * dt * ticks);

//     Hy1[Ny - 1] = eps[Ny - 1] * Ex1[Ny - 1];

//     for (int i = Ny; i < Nx; i++)
//     {
//         Hy2[i] =
//             Hy1[i] - (Ex1[i] - Ex1[i - 1]) * dt / dx;

//         Ex2[i - 1] =
//             Ex1[i - 1] - (Hy2[i] - Hy2[i - 1]) * dt / (eps[i - 1] * dx);
//     }
// }

size_t FDTD_2D_UPDATED::GetCurrentTick()
{
    return time_step;
}

double FDTD_2D_UPDATED::SourceFunction(double time_step) 
{
    double lambda_0 = 500e-9;

    //  Frequency
    double w0 = 2 * PI * c0 / lambda_0;
    double tau = 30;
    double t0 = tau * 6;

    return exp(-pow((time_step - t0),2) / pow(tau,2)) * sin(w0 * time_step * dt);
}

void FDTD_2D_UPDATED::CalcNextLayer( std::vector<double> &vectX,
                        std::vector<double> &vectEx,
                        std::vector<double> &vectHy)
{

    // # Calculate the Ex field.
    // for (int k = 1; k < Nx; ++k) {
    //     Ex[k] = Ex[k] + coeff * (Hy[k - 1] - Hy[k]);
    // }

    //  Electromagnetic "hard" source.
    //  Put a Gaussian pulse in the middle.
    // double tmp  = ((t0 - time_step) / spread);
    // double pulse = exp(-coeff *  tmp * tmp);

    // int source_dist = 0;

    //  Two sources.
    // Ex[source_position - source_dist] = pulse;
    // Ex[source_position + source_dist] = pulse;

    // //  Absorbing Boundary Conditions
    // Ex[0] = boundary_low.front();
    // boundary_low.push_back(Ex[1]);
    // Ex[Nx - 1] = boundary_high.front();
    // boundary_high.push_back(Ex[Nx - 2]);

    // //  Calculate the Hy field.
    // for (int k = 0; k < Nx-1; ++k) {
    //     Hy[k] = Hy[k] + coeff * (Ex[k] - Ex[k + 1]);
    // }
    // Calculation();
    // BoundaryConditionsFirst();
    // BoundaryConditionsSecond();

    //  Update magnetic field boundaries.
    Hy[jmax - 1] = Hy_prev[jmax - 2];

    
    for(int j = 0; j < (jmax - 1); ++j) 
    {
        Hy[j] = Hy_prev[j] + dt / (dx * mu0) * (Ex[j + 1] - Ex[j]);
        Hy_prev[j] = Hy[j];
    }

    //  Magnetic field source.
    Hy[jsource - 1] -= SourceFunction(time_step) / imp0;
    Hy_prev[jsource - 1] = Hy[jsource - 1];

    //  Update electric field boundaries.
    Ex[0] = Ex_prev[1];

    //  Update electric field.
    for(int j = 1; j < jmax; ++j) 
    {
        Ex[j] = Ex_prev[j] + dt / (dx * eps[j]) * (Hy[j] - Hy[j - 1]);
        Ex_prev[j] = Ex[j];
    }

    //  Electric field source.
    Ex[jsource - 1] -= SourceFunction(time_step + 1);
    Ex_prev[jsource] = Ex[jsource];



    for (int i = 0; i < jmax; i++)
    {
      
        vectX.push_back(i);
        vectEx.push_back(Ex[i]);
        vectHy.push_back(Hy[i]);
    }

    time_step++;
}

