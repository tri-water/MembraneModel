#pragma once
#include <cmath>
#include "mesh.h"

class thermoDynamics
{
	friend class electrodeReaction;
	friend class interfaceReaction;
public:
	thermoDynamics(double fE_formal, int fn, double fk0, double falfa);
	~thermoDynamics() {};

	const double F; // Faraday constant, C mol^-1
	const int n; // the number of electron transfered
	const double E_formal; // standard potential, V
	const double F_R_T; // F/(RT)
	const double nF_R_T; // nF/(RT)
	const double minusAlfaNF_R_T; // -alfa*nF/(RT)
	const double OneMinusAlfaNF_R_T; // (1-alfa)*nF/(RT)
	const double k0; // standard rate constant, convert the unit of k0 from m s^-1 to dm^3 m^-2 s^-1

	double ratio_ox2red(double E) const; 	// calculate the Red concentration under thermodynamic equilibarium
											// E: the potential subjected to the reaction
											// C_total: the total concentration of the redox pair
											// C_s/C_m equivalent

private:
	const double R; // gas constant, J K^-1 mol^-1
	const double T; // temperature, K
	double ReciprocalThermConst;
	const double alfa; // reduction charge transfer coefficient
};

class electrodeReaction
{
public:
	electrodeReaction(double fEpsilon_d, double fEpsilon_oc, double fEpsilon_ic, double fmu_i, thermoDynamics& fInnerThermo, mesh& fMesh);
	double DrivingPotential(double appliedPotential, double PotentialOHP) const { return appliedPotential- PotentialOHP - E_formal; }
	double potentialOHP(double appliedPotential, double potentialGrad) const { return appliedPotential + DrivingPotentialCoeff*potentialGrad; }
	double kf(double drivingPotnetial) const { return InnerThermo.k0*exp(InnerThermo.minusAlfaNF_R_T*drivingPotnetial); }; // forward reaction rate
	double kb(double drivingPotential) const { return InnerThermo.k0*exp(InnerThermo.OneMinusAlfaNF_R_T*drivingPotential); }; // backward reaction rate
	const double nF_R_T; // nF/(RT)
	const double F_RT;
	const double minusAlfaNF_R_T; // -alfa*nF/(RT)
	const double OneMinusAlfaNF_R_T; // (1-alfa)*nF/(RT)
	const double E_formal; // standard potential, V
	const double DrivingPotentialCoeff; // Electrode potential minus Phi_OHP

private:
	const double Epsilon_d; // the effective dielectric constants of the diffuse double-layer
	const double Epsilon_oc; // the effective dielectric constants of the outer part of electric double-layer
	const double Epsilon_ic; // the effective dielectric onstants of the inner part of electric double layer
	const double mu_i; // the thickness of the inner part of the electric double-layer
	const double mu; // the thickness of the electric double-layer
	

	const thermoDynamics& InnerThermo;
};

/*
class interfaceReaction
{
public:
	interfaceReaction(double fkf, double fkb);
	const double kf; //convert the unit of kf from m s^-1 to dm^3 m^-2 s^-1
	const double kb; // convert the unit of kb from m s^-1 to dm^3 m^-2 s^-1
};
*/

class interfaceReaction
{
public:
	interfaceReaction(thermoDynamics& fInnerThermo);
	interfaceReaction() : F_RT(0), nF_R_T(0), dE_formal(0), minusAlfaNF_R_T(0), OneMinusAlfaNF_R_T(0), InnerThermo(thermoDynamics(0, 0, 0, 0)) {};

	double drivingPotential(double dE) const { return dE - dE_formal; }
	const double F_RT;
	const double nF_R_T; // nF/(RT)
	const double dE_formal; // standard inner electrical potential
	const double minusAlfaNF_R_T; // -alfa*nF/(RT)
	const double OneMinusAlfaNF_R_T; // (1-alfa)*nF/(RT)
	virtual double kf(double drivingPotnetial) const { return InnerThermo.k0*exp(InnerThermo.OneMinusAlfaNF_R_T*drivingPotnetial); }; // forward reaction rate
	virtual double kb(double drivingPotential) const { return InnerThermo.k0*exp(InnerThermo.minusAlfaNF_R_T*drivingPotential); }; // backward reaction rate

private:
	thermoDynamics& InnerThermo;
};

class nonChargedInterfaceReaction : public interfaceReaction
{
public:
	nonChargedInterfaceReaction(double fkf, double fkb);

	virtual double kf(double drivingPotnetial) const { return innerkf; }; // forward reaction rate
	virtual double kb(double drivingPotential) const { return innerkb; }; // backward reaction rate
private:
	const double innerkf;
	const double innerkb;
};