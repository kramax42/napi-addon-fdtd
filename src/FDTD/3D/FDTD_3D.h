#pragma once
#include <vector>
#include <math.h>
#include "../../globals.h"

using namespace std;

class FDTD_3D {

protected:
	size_t ticks;
	const double PI = 3.141592653;

	// static const size_t Nx = 400;
	// static const size_t Ny = 400;

	double yy1[Nx][Ny];

	//Ez
	double Ez1[Nx][Ny], Ez2[Nx][Ny];
	// Hy
	double Hy1[Nx][Ny], Hy2[Nx][Ny];
	//Hx
	double Hx1[Nx][Ny], Hx2[Nx][Ny];

	// lambda - wave length, mkm
	double lambda;

	// n1 - refractive index
	double n1;

	// radius of beam  mkm
	double beamsize;

	double Ti;
	double dc;
	double gammar;
	double w0;
	double yMax;
	double dt;
	double dx;

	void boundary_conditions_1();
	void boundary_conditions_2();
	void boundary_conditions_3();
	void boundary_conditions_4();

	virtual void Calculation();

public:
	FDTD_3D(double lambda, double beamsize, double n1);


	virtual void setParams();

	// Getters.
	static size_t getNx() { return Nx; }
	static size_t getNy() { return Ny; }
	size_t getStep() { return 4; }

	double getLambda() { return lambda; }
	double getBeamsize() { return beamsize; }
	double getN1() { return n1; }

	// Setters.
	void setLambda(double l) { lambda = l; }
	void setBeamsize(double b) { beamsize = b; }
	void setN1(double n) { n1 = n; }

	size_t getCurrentTick();


	//start calculation
	void calcNextLayer(vector<double> &vectX,
					   vector<double> &vectY,
					   vector<double> &vectEz,
					   vector<double> &vectHy,
					   vector<double> &vectHx,
					   vector<double> &vectEnergy);
};

