#include "FDTD_3D.h"

FDTD_3D::FDTD_3D(double lambda, double beamsize, double n1)
			   : lambda(lambda), beamsize(beamsize), n1(n1)
{
	setParams();
}

void FDTD_3D::boundary_conditions_1()
{
	for (int j = 0; j < Ny; j++) {
		Ez2[0][j] = Ez1[1][j] + (dt / n1 - dx) * (Ez2[1][j] - Ez1[0][j]) / (dt / n1 + dx);
		Hy2[0][j] = Hy1[1][j] + (dt / n1 - dx) * (Hy2[1][j] - Hy1[0][j]) / (dt / n1 + dx);
		Hx2[0][j] = Hx1[1][j] + (dt / n1 - dx) * (Hx2[1][j] - Hx1[0][j]) / (dt / n1 + dx);
	}
}


void FDTD_3D::boundary_conditions_3()
{
	for (int i = 0; i < Nx; i++) {
		Ez2[i][0] = Ez1[i][1] + (dx - dt / n1) * (Ez2[i][1] - Ez1[i][0]) / (dt / n1 + dx);
		Hy2[i][0] = Hy1[i][1] + (dx - dt / n1) * (Hy2[i][1] - Hy1[i][0]) / (dt / n1 + dx);
		Hx2[i][0] = Hx1[i][1] + (dx - dt / n1) * (Hx2[i][1] - Hx1[i][0]) / (dt / n1 + dx);
	}
}


void FDTD_3D::boundary_conditions_2()
{
	for (int j = 0; j < Ny; j++) {
		Hy2[Nx - 1][j] = Hy1[Nx - 2][j] + (dt / n1 - dx) * (Hy2[Nx - 2][j] - Hy1[Nx - 1][j]) / (dt / n1 + dx);
		Ez2[Nx - 1][j] = Ez1[Nx - 2][j] + (dt / n1 - dx) * (Ez2[Nx - 2][j] - Ez1[Nx - 1][j]) / (dt / n1 + dx);
		Hx2[Nx - 1][j] = Hx1[Nx - 2][j] + (dt / n1 - dx) * (Hx2[Nx - 2][j] - Hx1[Nx - 1][j]) / (dt / n1 + dx);
	}
}

void FDTD_3D::boundary_conditions_4()
{
	for (int i = 0; i < Nx; i++) {
		Hy2[i][Ny - 1] = Hy1[i][Ny - 2] + (dt / n1 - dx) * (Hy2[i][Ny - 2] - Hy1[i][Ny - 1]) / (dt / n1 + dx);
		Ez2[i][Ny - 1] = Ez1[i][Ny - 2] + (dt / n1 - dx) * (Ez2[i][Ny - 2] - Ez1[i][Ny - 1]) / (dt / n1 + dx);
		Hx2[i][Ny - 1] = Hx1[i][Ny - 2] + (dt / n1 - dx) * (Hx2[i][Ny - 2] - Hx1[i][Ny - 1]) / (dt / n1 + dx);
	}
}

void FDTD_3D::Calculation()
{

	for (int j = 1; j < Ny - 1; j++) {
		int i = 0;
		Ez1[i][j] = sin(2 * PI * dt * ticks) * std::exp(gammar * dx * dx * (yMax - j) * (j - yMax));
		Hy1[i][j] = Ez1[i][j];
	}

	for (int i = 0; i < Nx - 1; i++) {
		for (int j = 0; j < Ny - 1; j++) {
			Hx2[i][j + 1] = Hx1[i][j + 1] - (Ez1[i][j + 1] - Ez1[i][j]) * dt / dx;
			Hy2[i + 1][j] = Hy1[i + 1][j] + (Ez1[i + 1][j] - Ez1[i][j]) * dt / dx;
			Ez2[i][j] = Ez1[i][j] + (Hy2[i + 1][j] - Hy2[i][j] - Hx2[i][j + 1] + Hx2[i][j]) * dt / (yy1[i][j] * dx);
		}
	}
}

size_t FDTD_3D::getCurrentTick()
{
	return ticks;
}

void FDTD_3D::setParams()
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
}
//---------------------------------------------------------------------------

void FDTD_3D::calcNextLayer(vector<double> &vectX,
							vector<double> &vectY,
							vector<double> &vectEz,
							vector<double> &vectHy,
							vector<double> &vectHx,
							vector<double> &vectEnergy)
{
	//setParams();
	//std::fstream out("out.txt");
	//out << "Ez" << "\t" << "\t" << "Hy" << "\t" << "\t" << "Hx" << "\t" << "Energy: " << "\n";

		Calculation();
		boundary_conditions_1();
		boundary_conditions_3();
		boundary_conditions_2();
		boundary_conditions_4();

		for (int i = 0; i < Nx; i++) {
			for (int j = 0; j < Ny; j++) {
				Hx1[i][j] = Hx2[i][j];
				Hy1[i][j] = Hy2[i][j];
				Ez1[i][j] = Ez2[i][j];
			}
		}

		size_t step = getStep();

		for (int xx = 1; xx < Nx - 1; xx += step)
		{
			// vectX.push_back(xx);
			// vectY.push_back(xx);
			for (int yy = 1; yy < Ny - 1; yy += step)
			{
				// Wem
				double energy = yy1[xx][yy] * yy1[xx][yy] * Ez1[xx][yy] * Ez1[xx][yy] +
					Hy1[xx][yy] * Hy1[xx][yy] + Hx1[xx][yy] * Hx1[xx][yy];

				vectX.push_back(xx);
				vectY.push_back(yy);
				vectEz.push_back(Ez1[xx][yy]);
				vectHy.push_back(Hy1[xx][yy]);
				vectHx.push_back(Hx1[xx][yy]);
				vectEnergy.push_back(energy);
				}
			}

			ticks++;
}


