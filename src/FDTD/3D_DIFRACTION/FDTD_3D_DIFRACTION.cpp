#include "FDTD_3D_DIFRACTION.h"

FDTD_3D_DIFRACTION::FDTD_3D_DIFRACTION(double lambda, double beamsize, double n1, double n2)
    : FDTD_3D(lambda, beamsize, n1), n2(n2)
{
    setParams();
}

void FDTD_3D_DIFRACTION::setParams()
{
    ticks = 0;

    dc = 0.3;
    Ti = lambda / dc;
    dx = 1.0 / 20;
    dt = 1.0 / 40;
    yMax = (Ny - 1) / 2;
    gammar = lambda * lambda / (2 * beamsize * beamsize);
    w0 = 1.0;

    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            Ez1[i][j] = 1e-16;
            Ez2[i][j] = 1e-16;
            Hy1[i][j] = 1e-16;
            Hy2[i][j] = 1e-16;
            Hx1[i][j] = 1e-16;
            Hx2[i][j] = 1e-16;
            yy1[i][j] = n1;
        }
    }


    for (size_t i = gridBeginX; i <= gridEndX; i++)
    {
        // Each grid gap.
        for (size_t j = 0; j < gridGapCount; j += 2)
        {
            for (int k = gridGap * j; k < gridGap * (j + 1); k++)
            {
                yy1[i][j] = n2;
            }
        }
    }
}