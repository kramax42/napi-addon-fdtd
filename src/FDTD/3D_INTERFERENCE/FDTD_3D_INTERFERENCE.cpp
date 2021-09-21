#include "FDTD_3D_INTERFERENCE.h"

FDTD_3D_INTERFERENCE::FDTD_3D_INTERFERENCE(double lambda, double beamsize, double n1)
                                          : FDTD_3D(lambda, beamsize, n1){}

void FDTD_3D_INTERFERENCE::Calculation()
{

    for (int j = 1; j < Ny - 1; j++)
    {
        int i = 0;
        Ez1[i][j] = sin(2 * PI * dt * ticks) * exp(gammar * dx * dx * (yMax - j) * (j - yMax));
        Hy1[i][j] = Ez1[i][j];
    }

    for (int i = 0; i < Nxp - 1; i++)
    {
        for (int j = 0; j < Ny - 1; j++)
        {
            Hx2[i][j + 1] = Hx1[i][j + 1] - (Ez1[i][j + 1] - Ez1[i][j]) * dt / dx;
            Hy2[i + 1][j] = Hy1[i + 1][j] + (Ez1[i + 1][j] - Ez1[i][j]) * dt / dx;
            Ez2[i][j] = Ez1[i][j] + (Hy2[i + 1][j] - Hy2[i][j] - Hx2[i][j + 1] + Hx2[i][j]) * dt / (yy1[i][j] * dx);
        }
    }


    /// NEED IN REFACTORING.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    int i = Nxp;

    for (int j = 0; j < P11; j++)
    {
        Hx2[i][j + 1] = Hx1[i][j + 1] - (Ez1[i][j + 1] - Ez1[i][j]) * dt / dx;
        Hy2[i + 1][j] = Hy1[i + 1][j] + (Ez1[i + 1][j] - Ez1[i][j]) * dt / dx;
        Ez2[i][j] = Ez1[i][j] + (Hy2[i + 1][j] - Hy2[i][j] - Hx2[i][j + 1] + Hx2[i][j]) * dt / (yy1[i][j] * dx);
    }

    for (int j = P11 + 1; j < P12; j++)
    {
        Hy2[Nxp][j] = Hy1[Nxp - 1][j] + (dt / n1 - dx) * (Hy2[Nxp - 1][j] - Hy1[Nxp][j]) / (dt / n1 + dx);
        Ez2[Nxp][j] = Ez1[Nxp - 1][j] + (dt / n1 - dx) * (Ez2[Nxp - 1][j] - Ez1[Nxp][j]) / (dt / n1 + dx);
        Hx2[Nxp][j] = Hx1[Nxp - 1][j] + (dt / n1 - dx) * (Hx2[Nxp - 1][j] - Hx1[Nxp][j]) / (dt / n1 + dx);
    }

    for (int j = P12 + 1; j < P22; j++)
    {
        Hx2[i][j + 1] = Hx1[i][j + 1] - (Ez1[i][j + 1] - Ez1[i][j]) * dt / dx;
        Hy2[i + 1][j] = Hy1[i + 1][j] + (Ez1[i + 1][j] - Ez1[i][j]) * dt / dx;
        Ez2[i][j] = Ez1[i][j] + (Hy2[i + 1][j] - Hy2[i][j] - Hx2[i][j + 1] + Hx2[i][j]) * dt / (yy1[i][j] * dx);
    }

    for (int j = P22 + 1; j < P21; j++)
    {
        Hy2[Nxp][j] = Hy1[Nxp - 1][j] + (dt / n1 - dx) * (Hy2[Nxp - 1][j] - Hy1[Nxp][j]) / (dt / n1 + dx);
        Ez2[Nxp][j] = Ez1[Nxp - 1][j] + (dt / n1 - dx) * (Ez2[Nxp - 1][j] - Ez1[Nxp][j]) / (dt / n1 + dx);
        Hx2[Nxp][j] = Hx1[Nxp - 1][j] + (dt / n1 - dx) * (Hx2[Nxp - 1][j] - Hx1[Nxp][j]) / (dt / n1 + dx);
    }

    for (int j = P21 + 1; j < Ny; j++)
    {
        Hx2[i][j + 1] = Hx1[i][j + 1] - (Ez1[i][j + 1] - Ez1[i][j]) * dt / dx;
        Hy2[i + 1][j] = Hy1[i + 1][j] + (Ez1[i + 1][j] - Ez1[i][j]) * dt / dx;
        Ez2[i][j] = Ez1[i][j] + (Hy2[i + 1][j] - Hy2[i][j] - Hx2[i][j + 1] + Hx2[i][j]) * dt / (yy1[i][j] * dx);
    }




    for (int i = Nxp; i < Nx - 1; i++)
    {
        for (int j = 0; j < Ny - 1; j++)
        {
            Hx2[i][j + 1] = Hx1[i][j + 1] - (Ez1[i][j + 1] - Ez1[i][j]) * dt / dx;
            Hy2[i + 1][j] = Hy1[i + 1][j] + (Ez1[i + 1][j] - Ez1[i][j]) * dt / dx;
            Ez2[i][j] = Ez1[i][j] + (Hy2[i + 1][j] - Hy2[i][j] - Hx2[i][j + 1] + Hx2[i][j]) * dt / (yy1[i][j] * dx);
        }
    }


    /// NEED IN REFACTORING.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
}

