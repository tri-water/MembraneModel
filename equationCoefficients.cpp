#include "equationCoefficients.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include "ErrorTypes.h"


equationCoefficients::equationCoefficients(const mesh& fMesh, const potentialSignal& fSignal, const thermoDynamics& fThermo, const IonSystem& membraneIons, const IonSystem& solutionIons) :
	m1(fMesh.NumberOfNodes), m2(fMesh.NumberOfNodes), m3(fMesh.NumberOfNodes), m3_(fMesh.NumberOfNodes), m4(fMesh.NumberOfNodes), 
	m5(fMesh.NumberOfNodes), m6(fMesh.NumberOfNodes), m7(fMesh.NumberOfNodes), m8(fMesh.NumberOfNodes), m9(fMesh.NumberOfNodes), 
	m10(fMesh.NumberOfNodes), m11(fMesh.NumberOfNodes), dt(fSignal.dt), F_RT(fThermo.F_R_T), re(fMesh.re), dy(fMesh.dy), y(fMesh.y), A(fMesh.A),
	_EpsilonrEpsilon0m(membraneIons.ReciprocalEpsilon_rEpsilon_0), _EpsilonrEpsilon0s(solutionIons.ReciprocalEpsilon_rEpsilon_0)
{
	calculateCoeff();
}

/*
void equationCoefficients::calculateCoeff()
{
	double dydt = dy*dt;
	double rere = re*re;
	double dydy = dy*dy;
	double reredy = re*re*dy;
	double reredydy = re*re*dy*dy;

	double dydt_rere = dydt / rere;
	double dt_rere = dt / rere;

	double int4dydy = 4 * dydy;
	double int4F_RT = 4 * F_RT;
	double int8F_RT = 8 * F_RT;

#pragma omp parallel for
	for (int i = 0; i < static_cast<int>(m1.rows()); ++i) {
		m1[i] = 4 * pow(1 - y[i], 3)*dydt_rere + 4 * pow(1 - y[i], 4)*dt_rere;
		m2[i] = 4 * pow(1 - y[i], 4)*dt_rere - 4 * pow(1 - y[i], 3)*dydt_rere;
		m3[i] = -int4dydy;
		m3_[i] = -8 * pow(1 - y[i], 4)*dt_rere;
		m4[i] = F_RT*pow(1 - y[i], 4)*dt_rere;
		m5[i] = int4F_RT*pow(1 - y[i], 3)*dydt_rere + int4F_RT*pow(1 - y[i], 4)*dt_rere;
		m6[i] = int4F_RT*pow(1 - y[i], 4)*dt_rere - 4*F_RT*pow(1 - y[i], 3)*dydt_rere;
		m7[i] = -int8F_RT*pow(1 - y[i], 4)*dt_rere;
		m8[i] = int4dydy;
		m9[i] = pow(1 - y[i], 3) / (reredy) + pow(1 - y[i], 4) / (reredydy);
		m10[i] = pow(1 - y[i], 4) / (reredydy) - pow(1 - y[i], 3) / (reredy);
		m11[i] = -2 * pow(1 - y[i], 4) / (reredydy);
	}

	m12m = _EpsilonrEpsilon0m*96485.33289*1000; // faraday constant
	m12s = _EpsilonrEpsilon0s*96485.33289*1000;
}
*/

void equationCoefficients::calculateCoeff()
{
	double dydt = dy*dt;
	double dydy = dy*dy;
	double AA = A*A;
	double re_A = re - A;
	double dt_AA = dt / AA;
	double dydt_A = dydt / A;
	double AAdydy = AA*dydy;
	double Ady = A*dy;

	double int4F_RT = 4 * F_RT;
	double int8F_RT = 8 * F_RT;

	double int1_y3 = 0;
	double int1_y4 = 0;
	double rePre_ATy = 0;

#pragma omp parallel for private(int1_y3, int1_y4, rePre_ATy)
	for (int i = 0; i < static_cast<int>(m1.rows()); ++i) {
		int1_y3 = pow(1 - y[i], 3); // (1 - y)^3
		int1_y4 = pow(1 - y[i], 4); // (1 - y)^4
		rePre_ATy = re + re_A*y[i]; // re + re_A*y

		m1[i] = 4 * int1_y3 / rePre_ATy * dydt_A + 4 * int1_y4 * dt_AA;
		m2[i] = 4 * int1_y4 * dt_AA - 4 * int1_y3 / rePre_ATy * dydt_A;
		m3[i] = -4 * dydy;
		m3_[i] = -8 * int1_y4 * dt_AA;
		m4[i] = F_RT*int1_y4 * dt_AA;
		m5[i] = int4F_RT*int1_y3 / rePre_ATy*dydt_A + int4F_RT*int1_y4 * dt_AA;
		m6[i] = int4F_RT*int1_y4 * dt_AA - int4F_RT*int1_y3 / rePre_ATy * dydt_A;
		m7[i] = -int8F_RT*int1_y4 * dt_AA;
		m8[i] = 4 * dydy;
		m9[i] = int1_y3 / rePre_ATy / Ady + int1_y4 / AAdydy;
		m10[i] = int1_y4 / AAdydy - int1_y3 / rePre_ATy / Ady;
		m11[i] = -2 * int1_y4 / AAdydy;
	}

	m12m = _EpsilonrEpsilon0m*96485.33289 * 1000; // faraday constant
	m12s = _EpsilonrEpsilon0s*96485.33289 * 1000;
}

void equationCoefficients::printParams(std::string paramName) {
	auto data = m1.head(10);
	if (paramName == "m1") {
		data = m1.head(10);
	}
	else if (paramName == "m2") {
		data = m2.head(10);
	}
	else if (paramName == "m3") {
		data = m3.head(10);
	}
	else if (paramName == "m4") {
		data = m4.head(10);
	}
	else if (paramName == "m5") {
		data = m5.head(10);
	}
	else if (paramName == "m6") {
		data = m6.head(10);
	}
	else if (paramName == "m7") {
		data = m7.head(10);
	}
	else if (paramName == "m8") {
		data = m8.head(10);
	}
	else if (paramName == "m9") {
		data = m9.head(10);
	}
	else if (paramName == "m10") {
		data = m10.head(10);
	}
	else if (paramName == "m11") {
		data = m11.head(10);
	}
	else {
		std::cout << "Please enter a right parameter name";
		throw(ErrorParamName());
	}

	std::ofstream printfile;
	printfile.open(paramName + ".txt");
	printfile << data;
	printfile.close();
}