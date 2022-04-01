#include <vector>
#include <napi.h>
#include <math.h>
#include <vector>

using namespace std;
class FDTD_2D_UPDATED
{
    size_t time_step;
    const double PI = 3.141592653;  

 
    // Grid size
    static const size_t jmax = 400;

   
    double ca[jmax];
    double cb[jmax];


    // epsilon - electric permitivitty
    double eps[jmax];

    // omega - conductivity
    double omega[jmax];

    // Magnetic field strength.
    double Hy[jmax];

    // Electric field strength.
    double Ex[jmax];

    // lambda - wave length
    double lambda;

    // tau - pulse duration
    double tau;

    // Light speed.
    double c0 = 3e8;


    int source_position;
    int nmax = 400;

    //  Coefficient corresponds Courant Condition.
    float cfl_factor = 0.5;

   
    vector <double> boundary_low = {0, 0};
    vector <double> boundary_high = {0, 0};




    // Frequency in MHz.
    double freq_in = 700e6;

    double lambda0 = c0 / freq_in;

    double epsz = 8.854e-12;

    //  Grid steps.
    double dx = lambda0 / 25;
    double dt = dx / (2*c0);

    double t0 = 50;
    int spread = 10;



    // Moor`s boundary condition.
    // void BoundaryConditionsFirst();

    // Moor`s boundary condition.
    // void BoundaryConditionsSecond();

    // Updating values for new time layer.
    // void Calculation();

    


public:
    FDTD_2D_UPDATED(double lambda, double tau, std::vector<double>& Epsilon, std::vector<double>& Omega, int source_position);

    void setParams(std::vector<double>& Epsilon, std::vector<double>& Omega);

    // Getters.
    size_t GetNx();
    double GetLambda();
    double GetTau();
    int GetSourcePosition();
    // double GetRefractiveIndex();

    // Setters.
    double SourceFunction(double time_step);
    void setLambda(double l) { lambda = l; }
    void setTau(double t) { tau = t; }
    void setSourcePosition(int new_source_pos) { source_position = new_source_pos; }
    // void setRefractiveIndex(double n) { refractive_index = n; }

    size_t GetCurrentTick();

    //start calculation
    void CalcNextLayer( std::vector<double> &vectX,
                        std::vector<double> &vectEx,
                        std::vector<double> &vectHy);
};

