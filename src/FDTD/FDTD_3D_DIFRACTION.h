#include "FDTD_3D.h"

class FDTD_3D_DIFRACTION : public FDTD_3D
{
    // n2 - refractive index
    double n2;

public:
    FDTD_3D_DIFRACTION(double lambda, double beamsize, double n1, double n2);
};