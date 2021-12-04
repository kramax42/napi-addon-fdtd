#pragma once
#include <vector>
#include <napi.h>
#include <math.h>



using namespace std;
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

    //grid size
    static const size_t Nx = 2001;
    static const size_t Ny = 501;

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

    // n1 - refractive index
    double n1;


    //Moor`s boundary condition.
    void boundary_conditions_1();

    //Moor`s boundary condition.
    void boundary_conditions_2();

    // Updating values for new time layer.
    void Calculation();



public:
    FDTD_2D(double lambda, double tau, double n1);

    //double lambda, double tau, double n1
    void setParams();

    // Getters.
    size_t getNx();
    double getLambda();
    double getTau();
    double getN1();

    // Setters.
    void setLambda(double l) { lambda = l; }
    void setTau(double t) { tau = t; }
    void setN1(double n) { n1 = n; }

    size_t getCurrentTick();

    //start calculation
    void calcNextLayer( vector<double> &vectX,
                        vector<double> &vectY);
};

