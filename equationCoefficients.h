#pragma once
#include "potentialSignal.h"
#include "mesh.h"
#include "thermodynamics.h"
#include "IonSystem.h"
#include <string.h>

class equationCoefficients
{
public:
	equationCoefficients(const mesh& fMesh, const potentialSignal& fSignal, const thermoDynamics& fThermo, const IonSystem& membraneIons, const IonSystem& solutionIons);
	
	Eigen::VectorXd m1, m2, m3, m3_, m4, m5, m6, m7, m8, m9, m10, m11; // equation coefficients
	double m12m, m12s; // permeability coefficient for membrane and solution phases

private:
	void calculateCoeff();
	void printParams(std::string paramName);

	const double dt;
	const double F_RT;
	const double re;
	const double dy;
	const double A;
	const Eigen::VectorXd& y;
	const double _EpsilonrEpsilon0s;
	const double _EpsilonrEpsilon0m;
};