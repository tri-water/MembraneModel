#include "IonSystem.h"
#include <fstream>
#include <cmath>
#include <iomanip>

Ion::Ion(double fD, int fZ, double fCinitial, int fstartIndex, int fendIndex) :
	D(fD), Z(fZ), Cinitial(fCinitial), DensityN(fendIndex - fstartIndex + 1), 
	DensityLP(fendIndex - fstartIndex + 1), DensityLN(fendIndex - fstartIndex + 1),
	DensityCP(fendIndex - fstartIndex + 1), DensityCN(fendIndex - fstartIndex + 1),
	startIndex(fstartIndex), ZxD(static_cast<double>(fZ)*fD)
{
	for (int i = 0; i < static_cast<int>(DensityN.rows()); ++i) {
		DensityN[i] = Cinitial;
		DensityLP[i] = Cinitial;
		DensityLN[i] = Cinitial;
		DensityCP[i] = Cinitial;
		DensityCN[i] = Cinitial;
	}
}

void Ion::printDense(std::string fileName) const
{
	std::ofstream fout;
	fout.open(fileName);

	if (fout.is_open()) fout << std::setprecision(20) << DensityN;

	fout.close();
}

void Ion::printLastPeriodDense(std::string fileName) const
{
	std::ofstream fout;
	fout.open("Positive " + fileName);

	if (fout.is_open()) fout << std::setprecision(20) << DensityLP;

	fout.close();

	fout.open("Negative " + fileName);

	if (fout.is_open()) fout << std::setprecision(20) << DensityLN;

	fout.close();
}

ImmobileCharge::ImmobileCharge(int fz, double fC) :
	z(fz), C(fC), zxC(fz*fC)
{

}

IonSystem::IonSystem(double fEpsilon_r, map<Species, Ion*>& fIons, Ion& fPotential, ImmobileCharge& fImmobileCharge) :
	Epsilon_r(fEpsilon_r), Epsilon_0(8.854187817620e-12), Ions(fIons), Potential(fPotential), ReciprocalEpsilon_rEpsilon_0(1 / 8.854187817620e-12 / fEpsilon_r),
	epsilon_rEpsilon_0(8.854187817620e-12*fEpsilon_r), CxZImmobileCharge(fImmobileCharge.zxC), IonNum(static_cast<int>(fIons.size()))
{

}