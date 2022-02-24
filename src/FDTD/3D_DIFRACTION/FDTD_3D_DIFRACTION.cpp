#include "FDTD_3D_DIFRACTION.h"
#include <iostream>

FDTD_3D_DIFRACTION::FDTD_3D_DIFRACTION(double lambda, double beamsize, std::vector<std::vector<double>> &matrixRefrIndex)
    : FDTD_3D(lambda, beamsize, 1.5f)
{
    setParams(matrixRefrIndex);
}

void FDTD_3D_DIFRACTION::setParams(std::vector<std::vector<double>>& matrixRefrIndex)
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
            yy1[i][j] = matrixRefrIndex[i][j];
        }
    }
}