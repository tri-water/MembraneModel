#pragma once
#include <Eigen/Sparse>
#include "mesh.h"
#include <map>
#include "EnumStruct.h"

typedef Eigen::SparseVector<double> SpVectorXd;
typedef Eigen::SparseMatrix<double> SpMatrixXd;
using std::map;

class Ion
{
	friend class solver;
public:
	Ion(double fD, int fZ, double fCinitial, int startIndex, int endIndex);
	const double D; // the diffusion coefficient in membrane phase
	const int Z;
	const double Cinitial;
	const double ZxD;

	void printDense(std::string fileName) const;
	void printLastPeriodDense(std::string fileName) const;
	double DxCi(int i) {return D*((*this)[i]);}
	double ZxDxCi(int i) {return static_cast<double>(Z)*D*((*this)[i]);}
	double& operator[](int i) {return DensityN[i - startIndex];} // get density using mesh index

private:
	const int startIndex;

	Eigen::VectorXd DensityLP; // concentration value at the end of last positive pulse
	Eigen::VectorXd DensityLN; // concentration value at the end of last negative pulse
	Eigen::VectorXd DensityCP; // concentration value at the end of current positive pulse
	Eigen::VectorXd DensityCN; // concentration value at the end of current negative pulse
	Eigen::VectorXd DensityN; // for ions species, it stores the concentration values; for potential, it stores the potential values
};

class ImmobileCharge
{
public:
	ImmobileCharge(int fz, double fC);
	const double zxC;
private:
	const int z;
	const double C;
};

class IonSystem
{
public:
	IonSystem(double fEpsilon_r, map<Species, Ion*>& fIons, Ion& fPotential, ImmobileCharge& fImmobileCharge = ImmobileCharge(0, 0.0));
	Ion& Potential;
	const double ReciprocalEpsilon_rEpsilon_0;
	const double epsilon_rEpsilon_0;
	const double CxZImmobileCharge;
	int IonNum;

	Ion& operator[](Species i) const { return *Ions[i]; }
private:
	const double Epsilon_r;
	const double Epsilon_0; // 8.854187817620 × 10−12 F⋅m−1 -> 8.854187817620 × 10−14 F⋅cm−1
	map<Species, Ion*>& Ions;
};