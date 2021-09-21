#include "../3D/FDTD_3D.h"

class FDTD_3D_DIFRACTION : public FDTD_3D
{
    // n2 - refractive index
    double n2;

    // Difraction grid sizes.
    const size_t gridWidth = 10;
    const size_t gridGap = 20;
    size_t gridGapCount = static_cast<size_t>(Ny / gridGap);
    const size_t gridBeginX = 50;
    const size_t gridEndX = gridBeginX + gridGap;

public:
    FDTD_3D_DIFRACTION(double lambda, double beamsize, double n1, double n2);

    virtual void setParams();

    void setN2(double n) { n2 = n; }
};