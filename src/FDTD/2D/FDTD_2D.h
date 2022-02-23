// #pragma once
// #include <vector>
// #include <napi.h>
// #include <math.h>
// #include "../../globals.h"


// class FDTD_2D
// {
//     size_t ticks;
//     const double PI = 3.141592653;  

//     //grid steps
//     double dx;
//     double dt;

//     //const
//     double aa1;  //??
//     double Ti;   //??
//     double tMax; //??

//     // grid size
//     // static const size_t Nx = 2001;
//     // static const size_t Ny = 501;

//     // epsilon - the dielectric constant
//     double eps[Nx];

//     // Magnetic field strength.
//     double H1[Nx];
//     double H2[Nx];

//     // Electric field strength.
//     double E1[Nx];
//     double E2[Nx];

//     // lambda - wave length
//     float lambda;

//     // tau - pulse duration
//     float tau;

//     // refractive index
//     float refractive_index;


//     //Moor`s boundary condition.
//     void BoundaryConditionsFirst();

//     //Moor`s boundary condition.
//     void BoundaryConditionsSecond();

//     // Updating values for new time layer.
//     void Calculation();



// public:
//     FDTD_2D(float lambda, float tau, float refractive_index);

//     //double lambda, double tau, double refractive_index
//     void SetParams();

//     // Getters.
//     size_t GetNx();
//     float GetLambda();
//     float GetTau();
//     float GetRefractiveIndex();

//     // Setters.
//     void SetLambda(float new_lambda);
//     void SetTau(float new_tau);
//     void setRefractiveIndex(float new_refractive_index);

//     size_t GetCurrentTick();


// };



 #include <vector>
#include <napi.h>
#include <math.h>



using namespace std;
class FDTD_2D
{
    size_t ticks;
    const double PI = 3.141592653;  

     // Constants
    double aa1;  //??
    double Ti;   //??
    double tMax; //??

    // Grid size
    static const size_t Nx = 2001;
    static const size_t Ny = 501;

    // Grid steps.
    double dx;
    double dt;

    //epsilon - the dielectric constant
    double eps[Nx];

    //magnetic field strength
    double H1[Nx];
    double H2[Nx];


    //electric field strength
    double E1[Nx];
    double E2[Nx];

    // lambda - wave length
    double lambda;

    // tau - pulse duration
    double tau;

    // refractive index
    double refractive_index;


    // Moor`s boundary condition.
    void BoundaryConditionsFirst();

    // Moor`s boundary condition.
    void BoundaryConditionsSecond();

    // Updating values for new time layer.
    void Calculation();



public:
    FDTD_2D(double lambda, double tau, double refractive_index);

    //double lambda, double tau, double refractive_index
    void setParams();

    // Getters.
    size_t GetNx();
    double GetLambda();
    double GetTau();
    double GetRefractiveIndex();

    // Setters.
    void setLambda(double l) { lambda = l; }
    void setTau(double t) { tau = t; }
    void setRefractiveIndex(double n) { refractive_index = n; }

    size_t GetCurrentTick();

    //start calculation
    void CalcNextLayer( std::vector<double> &vectX,
                        std::vector<double> &vectY);
};

