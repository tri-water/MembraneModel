#include "thermodynamics.h"
#include <iostream>

thermoDynamics::thermoDynamics(double fE_formal, int fn, double fk0, double falfa):
	E_formal(fE_formal), R(8.3144621), T(293), F(9.64853399 * 10000), n(fn), k0(fk0), alfa(falfa), 
	F_R_T(9.64853399 * 10000 / 8.3144621 / 293), nF_R_T(fn*9.64853399 * 10000 / 8.3144621 / 293),
	minusAlfaNF_R_T(-falfa*fn*9.64853399 * 10000 / 8.3144621 / 293), OneMinusAlfaNF_R_T((1 - falfa)*fn*9.64853399 * 10000 / 8.3144621 / 293)
{
	ReciprocalThermConst = 1 / nF_R_T;
}

double thermoDynamics::ratio_ox2red(double E) const {
	double k = exp(nF_R_T*(E - E_formal));
	// return the Cox/Cred ratio
	return k;
}

electrodeReaction::electrodeReaction(double fEpsilon_d, double fEpsilon_oc, double fEpsilon_ic, double fmu_i, thermoDynamics& fInnerThermo, mesh& fMesh) :
	Epsilon_d(fEpsilon_d), Epsilon_oc(fEpsilon_oc), Epsilon_ic(fEpsilon_ic), mu_i(fmu_i), mu(fMesh.mu), InnerThermo(fInnerThermo),
	DrivingPotentialCoeff(fEpsilon_d/fEpsilon_oc*((fMesh.mu -fmu_i)/(fMesh.re+ fMesh.mu)/(fMesh.re+fmu_i) + fEpsilon_d/fEpsilon_ic*fMesh.mu/(fMesh.re+fmu_i)/fMesh.re)*(fMesh.re+ fMesh.mu)*(fMesh.re+ fMesh.mu)),
	nF_R_T(fInnerThermo.nF_R_T), minusAlfaNF_R_T(fInnerThermo.minusAlfaNF_R_T), OneMinusAlfaNF_R_T(fInnerThermo.OneMinusAlfaNF_R_T), E_formal(fInnerThermo.E_formal),
	F_RT(fInnerThermo.F_R_T)
{
	std::cout << "DrivingPotentialCoeff: " << DrivingPotentialCoeff << std::endl << "mu: " << mu << std::endl;

}
/*
interfaceReaction::interfaceReaction(double fkf, double fkb) :
	kf(fkf*1000), kb(fkb*1000)
{

}
*/

interfaceReaction::interfaceReaction(thermoDynamics& fInnerThermo) :
	InnerThermo(fInnerThermo), nF_R_T(fInnerThermo.nF_R_T), minusAlfaNF_R_T(fInnerThermo.minusAlfaNF_R_T), 
	OneMinusAlfaNF_R_T(fInnerThermo.OneMinusAlfaNF_R_T), F_RT(fInnerThermo.F_R_T), dE_formal(fInnerThermo.E_formal)
{

}

nonChargedInterfaceReaction::nonChargedInterfaceReaction(double fkf, double fkb) :
	innerkf(fkf), innerkb(fkb)
{

}