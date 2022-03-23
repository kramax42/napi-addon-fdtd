#include <vector>
#include <napi.h>
#include <math.h>
#include <vector>

using namespace std;
class FDTD_2D_UPDATED
{
    size_t time_step;
    const double PI = 3.141592653;  

     // Constants
    double aa1;  //??
    double Ti;   //??
    double tMax; //??

    // Grid size
    static const size_t jmax = 400;
    // static const size_t Ny = 501;

   

    //epsilon - the dielectric constant
    double eps[jmax];

    //magnetic field strength
    double Hy[jmax];
    double Hy_prev[jmax];

    //electric field strength
    double Ex[jmax];
    double Ex_prev[jmax];

    // lambda - wave length
    double lambda;

    // tau - pulse duration
    double tau;

    //  Constants/
    double eps0 = 4.85418e-12;
    double mu0 = 1.25664e-10;

    // Light speed.
    double c0 = 1 / sqrt(eps0 * mu0);

    double imp0 = sqrt(mu0 / eps0);

    int jsource = 30;
    int nmax = 400;


    // Wave lenght in meters.
    double lambda_min = 350e-9;

    //  Grid steps.
    double dx = lambda_min / 20;
    double dt = dx / c0;



    //  Pulse parameters
    // int source_position = int(Nx / 2);
    // int t0 = 40;
    // int spread = 12;

    //  Coefficient corresponds Courant Condition.
    // double coeff = 0.5;

    // vector <double> boundary_low = {0, 0};
    // vector <double> boundary_high = {0, 0};




    // Moor`s boundary condition.
    // void BoundaryConditionsFirst();

    // Moor`s boundary condition.
    // void BoundaryConditionsSecond();

    // Updating values for new time layer.
    // void Calculation();



public:
    FDTD_2D_UPDATED(double lambda, double tau, std::vector<double>& Epsilon);

    //double lambda, double tau, double refractive_index
    void setParams(std::vector<double>& Epsilon);

    // Getters.
    size_t GetNx();
    double GetLambda();
    double GetTau();
    // double GetRefractiveIndex();

    // Setters.
    double SourceFunction(double time_step);
    void setLambda(double l) { lambda = l; }
    void setTau(double t) { tau = t; }
    // void setRefractiveIndex(double n) { refractive_index = n; }

    size_t GetCurrentTick();

    //start calculation
    void CalcNextLayer( std::vector<double> &vectX,
                        std::vector<double> &vectEx,
                        std::vector<double> &vectHy);
};

