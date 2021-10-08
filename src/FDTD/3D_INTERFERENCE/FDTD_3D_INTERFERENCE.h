#include "../3D/FDTD_3D.h"

class FDTD_3D_INTERFERENCE : public FDTD_3D
{

    const int P11 = 161;

    const int P12 = 181;

    const int P21 = 241;

    const int P22 = 221;

    const int Nxp = 100;

    virtual void Calculation();

public:
    FDTD_3D_INTERFERENCE(double lambda, double beamsize, double n1);
                    // : FDTD_3D(lambda, beamsize, n1){}

};