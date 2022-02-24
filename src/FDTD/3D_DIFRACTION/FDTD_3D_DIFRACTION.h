#include "../3D/FDTD_3D.h"

class FDTD_3D_DIFRACTION : public FDTD_3D
{
    // Difraction grid sizes.
    const size_t gridWidth = 10;
    const size_t gridGap = 20;

    // const size_t gridWidth = 5; // temporary value
    // const size_t gridGap = 3;    // temporary value

    size_t gridGapCount = static_cast<size_t>(Ny / gridGap);
    // const size_t gridBeginX = 5;
    
    const size_t gridBeginX = 35;
    const size_t gridEndX = gridBeginX + gridGap;

public:
    FDTD_3D_DIFRACTION(double lambda, double beamsize, std::vector<std::vector<double>> &matrixRefrIndex);
    virtual void setParams(std::vector<std::vector<double>> &matrixRefrIndex);
};