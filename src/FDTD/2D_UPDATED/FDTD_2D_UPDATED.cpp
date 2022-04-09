//http://zfdtd.narod.ru/ -- ������ ������ FDTD

#include "FDTD_2D_UPDATED.h"
#include <iostream>


FDTD_2D_UPDATED::FDTD_2D_UPDATED(double lambda, double tau, std::vector<double>& Epsilon, std::vector<double>& Sigma, std::vector<int> src)
                : lambda(lambda), tau(tau)
{
    setParams(Epsilon, Sigma, src);
}

//Getters
size_t FDTD_2D_UPDATED::GetNx() 
{
     return jmax;
}

double FDTD_2D_UPDATED::GetTau()
{
    return tau;
}

double FDTD_2D_UPDATED::GetLambda()
{
    return lambda;
}

int FDTD_2D_UPDATED::GetSourcePosition()
{
    return source_position;
}

//double lambda = 1, double tau = 10, double refractive_index = 1
void FDTD_2D_UPDATED::setParams(std::vector<double>& Epsilon, std::vector<double>& Sigma, std::vector<int> src)
{
    time_step = 0;

    float eaf;

    // for (int i = 0; i < src.size(); i++){
    //     source_position_vector.push_back(src[i]);
    // }
    source_position_vector = src;

    for (int i = 0; i < jmax; i++)
    {
        Ex[i] = 0;
        Hy[i] = 0;
    
        eps[i] = Epsilon[i];
        sigma[i] = Sigma[i];

        // eaf = dt * sigma[i] / (2 * epsz * eps[i]);
        // ca[i] = (1 - eaf) / (1 + eaf);
        // cb[i] = cfl_factor / (eps[i] * (1 + eaf));
        cb[i] = cfl_factor / eps[i];
        
    }

}


size_t FDTD_2D_UPDATED::GetCurrentTick()
{
    return time_step;
}

// double FDTD_2D_UPDATED::SourceFunction(double time_step) 
// {
//     double lambda_0 = 800e-6;

//     //  Frequency
//     double w0 = c0 / lambda_0;
//     double tau = 50;
//     double t0 = tau * 1;

//     return exp(-pow((time_step - t0),2) / pow(tau,2)) * sin(w0 * time_step * dt);
// }

void FDTD_2D_UPDATED::CalcNextLayer( std::vector<double> &vectX,
                        std::vector<double> &vectEx,
                        std::vector<double> &vectHy)
{

    // Calculate the Ex field.
    for (int k = 1; k < jmax; ++k) 
    {
        // Ex[k] = ca[k] * Ex[k] + cb[k] * (Hy[k - 1] - Hy[k]);
        Ex[k] = Ex[k] + cb[k] * (Hy[k - 1] - Hy[k]);
    }

    //  Electromagnetic "hard" source.
    //  Put a Gaussian pulse in the middle.
    double pulse = exp(-cfl_factor * pow((t0 - time_step) / spread, 2));//* sin(freq_in * time_step * dt);

    for (int i = 0; i < source_position_vector.size(); i++){
        Ex[source_position_vector[i]] += pulse;
    }
    


    // Absorbing Boundary Conditions.
    Ex[0] = boundary_low.front();
    boundary_low.erase(boundary_low.begin());
    boundary_low.push_back(Ex[1]);
    Ex[jmax - 1] = boundary_high.front();
    boundary_high.erase(boundary_high.begin());
    boundary_high.push_back(Ex[jmax - 2]);

    //  Calculate the Hy field.
    for (int k = 0; k < jmax-1; ++k) 
    {
        Hy[k] = Hy[k] + cfl_factor * (Ex[k] - Ex[k + 1]);
    }
    

    for (int i = 0; i < jmax; i++) 
    {
        vectX.push_back(i);
        vectEx.push_back(Ex[i]);
        vectHy.push_back(Hy[i]);
    }

    time_step++;
}










































// //http://zfdtd.narod.ru/ -- ������ ������ FDTD

// #include "FDTD_2D_UPDATED.h"
// #include <iostream>


// FDTD_2D_UPDATED::FDTD_2D_UPDATED(double lambda, double tau, std::vector<double>& Epsilon, std::vector<double>& sigma, int source_position)
//                 : lambda(lambda), tau(tau), source_position(source_position)
// {
//     setParams(Epsilon, sigma);
// }

// //Getters
// size_t FDTD_2D_UPDATED::GetNx() {
//      return jmax;
// }

// double FDTD_2D_UPDATED::GetTau(){
//     return tau;
// }

// double FDTD_2D_UPDATED::GetLambda()
// {
//     return lambda;
// }

// int FDTD_2D_UPDATED::GetSourcePosition()
// {
//     return source_position;
// }

// //double lambda = 1, double tau = 10, double refractive_index = 1
// void FDTD_2D_UPDATED::setParams(std::vector<double>& Epsilon, std::vector<double>& sigma)
// {
//     time_step = 0;

//     // Grid steps.
//     // dx = 0.05;
//     // dt = 0.025;


//     // Physics params.
//     // aa1 = lambda * lambda / (0.09 * tau * tau);
//     // tMax = 4 * tau / (lambda / 0.3);

//     float cfl_factor = 0.5;
//     float eaf;

//     for (int i = 0; i < jmax; i++)
//     {
//         Ex[i] = 0;
//         Ex_prev[i] = 0;

//         Hy[i] = 0;
//         Hy_prev[i] = 0;

//         // eps[i] = refractive_index;
//         eps[i] = Epsilon[i];
//         sigma[i] = sigma[i];

//         eaf = dt * sigma[i] / (2 * eps[i]);
//         ca[i] = (1.0 - eaf) / (1.0 + eaf);
//         cb[i] = cfl_factor / (eps[i] * (1 + eaf));
        
//     }

// }


// size_t FDTD_2D_UPDATED::GetCurrentTick()
// {
//     return time_step;
// }

// double FDTD_2D_UPDATED::SourceFunction(double time_step) 
// {
//     double lambda_0 = 500e-9;

//     //  Frequency
//     double w0 = 2 * PI * c0 / lambda_0;
//     double tau = 30;
//     double t0 = tau * 2;

//     return exp(-pow((time_step - t0),2) / pow(tau,2)) * sin(w0 * time_step * dt);
// }

// void FDTD_2D_UPDATED::CalcNextLayer( std::vector<double> &vectX,
//                         std::vector<double> &vectEx,
//                         std::vector<double> &vectHy)
// {

//     // # Calculate the Ex field.
//     // for (int k = 1; k < Nx; ++k) {
//     //     Ex[k] = Ex[k] + coeff * (Hy[k - 1] - Hy[k]);
//     // }

//     //  Electromagnetic "hard" source.
//     //  Put a Gaussian pulse in the middle.
//     // double tmp  = ((t0 - time_step) / spread);
//     // double pulse = exp(-coeff *  tmp * tmp);

//     // int source_dist = 0;

//     //  Two sources.
//     // Ex[source_position - source_dist] = pulse;
//     // Ex[source_position + source_dist] = pulse;

//     // //  Absorbing Boundary Conditions
//     // Ex[0] = boundary_low.front();
//     // boundary_low.push_back(Ex[1]);
//     // Ex[Nx - 1] = boundary_high.front();
//     // boundary_high.push_back(Ex[Nx - 2]);

//     // //  Calculate the Hy field.
//     // for (int k = 0; k < Nx-1; ++k) {
//     //     Hy[k] = Hy[k] + coeff * (Ex[k] - Ex[k + 1]);
//     // }
//     // Calculation();
//     // BoundaryConditionsFirst();
//     // BoundaryConditionsSecond();

//     //  Update magnetic field boundaries.
//     Hy[jmax - 1] = Hy_prev[jmax - 2];

    
//     for(int j = 0; j < (jmax - 1); ++j) 
//     {
//         // Hy[j] = Hy_prev[j] + dt / (dx * mu0) * (Ex[j + 1] - Ex[j]);
//         Hy[j] = ca[j] * Hy_prev[j] + cb[j] * (Ex[j + 1] - Ex[j]);
//         Hy_prev[j] = Hy[j];
//     }

//     //  Magnetic field source.
//     Hy[source_position - 1] -= SourceFunction(time_step) / imp0;
//     Hy_prev[source_position - 1] = Hy[source_position - 1];

//     //  Update electric field boundaries.
//     Ex[0] = Ex_prev[1];

//     //  Update electric field.
//     for(int j = 1; j < jmax; ++j) 
//     {
//         // Ex[j] = Ex_prev[j] + dt / (dx * eps[j]) * (Hy[j] - Hy[j - 1]);
//         Ex[j] = ca[j] * Ex_prev[j] + cb[j] * (Hy[j] - Hy[j - 1]);
        
//         Ex_prev[j] = Ex[j];
//     }

//     //  Electric field source.
//     Ex[source_position - 1] -= SourceFunction(time_step + 1);
//     Ex_prev[source_position] = Ex[source_position];



//     for (int i = 0; i < jmax; i++)
//     {
      
//         vectX.push_back(i);
//         vectEx.push_back(Ex[i]);
//         vectHy.push_back(Hy[i]);
//     }

//     time_step++;
// }

