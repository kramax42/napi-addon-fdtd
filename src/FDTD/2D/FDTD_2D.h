#pragma once
#include <vector>
#include <napi.h>
#include <math.h>
#include "../../globals.h"


class FDTD_2D
{
    size_t ticks;
    const double PI = 3.141592653;

    //grid steps
    double dx;
    double dt;

    //const
    double aa1;  //??
    double Ti;   //??
    double tMax; //??

    // grid size
    // static const size_t Nx = 2001;
    // static const size_t Ny = 501;

    // epsilon - the dielectric constant
    double eps[Nx];

    // Magnetic field strength.
    double H1[Nx];
    double H2[Nx];

    // Electric field strength.
    double E1[Nx];
    double E2[Nx];

    // lambda - wave length
    float lambda;

    // tau - pulse duration
    float tau;

    // refractive index
    float refractive_index;


    //Moor`s boundary condition.
    void BoundaryConditionsFirst();

    //Moor`s boundary condition.
    void BoundaryConditionsSecond();

    // Updating values for new time layer.
    void Calculation();



public:
    FDTD_2D(float lambda, float tau, float refractive_index);

    //double lambda, double tau, double n1
    void SetParams();

    // Getters.
    size_t GetNx();
    float GetLambda();
    float GetTau();
    float GetRefractiveIndex();

    // Setters.
    void SetLambda(float new_lambda);
    void SetTau(float new_tau);
    void SetRefractiveIndex(float new_refractive_index);

    size_t GetCurrentTick();

    //start calculation
    void CalcNextLayer( std::vector<double> &vectX,
                        std::vector<double> &vectY);
};

