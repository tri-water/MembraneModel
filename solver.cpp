#include "solver.h"
#include <iostream>
#include <fstream>
#include "ErrorTypes.h"

solver::solver(IonSystem& fmIons, IonSystem& fsIons, mesh& fMesh, equationCoefficients& fECoeff, electrodeReaction& feReaction,
	interfaceReaction& fproReaction, interfaceReaction& freaReaction, interfaceReaction& fcatReaction, potentialSignal& finput) :
	mIons(fmIons), sIons(fsIons), Mesh(fMesh), ECoeff(fECoeff), NumOfVars(fMesh.IndexOfBoundary * 4 + (fMesh.NumberOfNodes - fMesh.IndexOfBoundary) * 5),
	sBottomIndex(fMesh.IndexOfBoundary), mUpperIndex(fMesh.IndexOfBoundary - 1), sUpperIndex(fMesh.NumberOfNodes - 1), index1d(fMesh),
	eReaction(feReaction), proReaction(fproReaction), reaReaction(freaReaction), catReaction(fcatReaction), input(finput), sfEfield0(0), sfEfield1(0)
{
	J = SpMatrixXd(NumOfVars, NumOfVars);
	dX = Eigen::VectorXd(NumOfVars);
	X = Eigen::VectorXd(NumOfVars);
	F = Eigen::VectorXd(NumOfVars);

	JList.reserve(NumOfVars * 6);
	JIndex = 0;
}

void solver::solve()
{
	/*
	Newto-Raphson method.
	1. F(Xn)
	2. F'(Xn)dX = -F(Xn)
	3. Xn+1 = Xn + dX
	*/
	//Eigen::SparseLU<SpMatrixXd> dXSolver;
	size_t itercounter = 0;
	
	do {
		calculateF();
		J.makeCompressed();
		dXSolver.analyzePattern(J);
		dXSolver.factorize(J);
		
		if (dXSolver.info() != Eigen::Success) {
			//decomposition failed
			std::cout << "Decomposition of J Failed" << std::endl;
			std::ofstream fout;
			fout.open("errorX.txt");
			fout << X;
			fout.close();
			fout.open("errorF.txt");
			fout << F;
			fout.close();
			std::cout << dXSolver.info() << std::endl;
			throw(DecompositionFailed());
		}
		dX = dXSolver.solve(-F);
		X += dX;
		if (dXSolver.info() != Eigen::Success) {
			// Solving failed
			std::cout << "Solving of derivative matrix Failed" << std::endl;
			throw(SolvingFailed());
		}
		updateJ();
	} while ((dX.array()/X.array()).abs().maxCoeff() > 1e-2 && itercounter++ < 50);
	saveDensity();
	
	sfEfield0 = sfEfield1;
	int ei = index1d(Species::mPotential, 0);
	sfEfield1 = (X[ei + 1] - X[ei]) / (Mesh.R(1) - Mesh.R(0));
}

void solver::initialise()
{
	initialiseX();
	initialiseJ(&solver::LockedPushBack);
	JList.shrink_to_fit();
	int ei = index1d(Species::mPotential, 0);
	sfEfield1 = (X[ei + 1] - X[ei]) / (Mesh.R(1) - Mesh.R(0));
}

void solver::initialiseX()
{
	int upperIndexboundary = sUpperIndex + 1;

	auto& mReactant = mIons[Species::mReactant];
	auto& mProduct = mIons[Species::mProduct];
	auto& mCation = mIons[Species::mCation];
	auto& mPotential = mIons[Species::mPotential];

	auto& sReactant = sIons[Species::sReactant];
	auto& sProduct = sIons[Species::sProduct];
	auto& sCation = sIons[Species::sCation];
	auto& sAnion = sIons[Species::sAnion];
	auto& sPotential = sIons[Species::sPotential];

	for (int i = 0; i < sBottomIndex; ++i) {
		X(index1d(Species::mReactant, i)) = mReactant[i];
		X(index1d(Species::mProduct, i)) = mProduct[i];
		X(index1d(Species::mCation, i)) = mCation[i];
		X(index1d(Species::mPotential, i)) = mPotential[i];
	}

	for (int i = sBottomIndex; i < sUpperIndex + 1; ++i) {
		X(index1d(Species::sReactant, i)) = sReactant[i];
		X(index1d(Species::sProduct, i)) = sProduct[i];
		X(index1d(Species::sCation, i)) = sCation[i];
		X(index1d(Species::sAnion, i)) = sAnion[i];
		X(index1d(Species::sPotential, i)) = sPotential[i];
	}
}

void solver::initialiseJ(void (solver::*Assign)(Tt))
{
	int mUpperIndex = Mesh.IndexOfBoundary - 1;
	int sBottomIndex = Mesh.IndexOfBoundary;

	omp_init_lock(&writeLock);
	// membrane bulk
#pragma omp parallel for
	for (int i = 1; i < mUpperIndex; ++i) {
		bulkMTDerivative(i, mIons, Species::mCation, Species::mPotential, Assign);
		bulkMTDerivative(i, mIons, Species::mProduct, Species::mPotential, Assign);
		bulkMTDerivative(i, mIons, Species::mReactant, Species::mPotential, Assign);
		mBulkPsDerivative(i, mIons, Assign);
	}
	// membrane bottom
	mBottomDerivative(Assign);
	// membrane top
	mTopDerivative(Assign);
	// solution bulk
#pragma omp parallel for
	for (int i = sBottomIndex + 1; i < sUpperIndex; ++i) {
		bulkMTDerivative(i, sIons, Species::sAnion, Species::sPotential, Assign);
		bulkMTDerivative(i, sIons, Species::sCation, Species::sPotential, Assign);
		bulkMTDerivative(i, sIons, Species::sProduct, Species::sPotential, Assign);
		bulkMTDerivative(i, sIons, Species::sReactant, Species::sPotential, Assign);
		sBulkPsDerivative(i, sIons, Assign);
	}
	// solution bottom
	sBottomDerivative(Assign);
	// solution top
	sTopDerivative(Assign);

	omp_destroy_lock(&writeLock);

	// initialise J
	J.setFromTriplets(JList.begin(), JList.end());
}

void solver::updateJ()
{
	JIndex = 0;
	initialiseJ(&solver::LockedIndexAssign);
}

void solver::calculateF()
{
	// membrane bulk
#pragma omp parallel for
	for (int i = 1; i < mUpperIndex; ++i) {
		mBulkF(i);
	}
	// membrane bottom
	mBottomF();
	// membrane top
	mTopF();

	// solution bulk
#pragma omp parallel for
	for (int i = sBottomIndex + 1; i < sUpperIndex; ++i) {
		sBulkF(i);
	}
	// solution bottom
	sBottomF();
	// soluiton top
	sTopF();
}

void solver::bulkMTDerivative(int i, IonSystem& ions, Species k, Species potential, void (solver::*Assign)(Tt))
{
	int index = index1d(k, i);
	auto& ion = ions[k];
	auto& eField = ions.Potential;

	double Cp1 = ECoeff.m1[i] * ion.D + ECoeff.m4[i] * ion.ZxD*(X(index1d(potential, i + 1)) - X(index1d(potential, i - 1)));
	double Cm1 = ECoeff.m2[i] * ion.D - ECoeff.m4[i] * ion.ZxD*(X(index1d(potential, i + 1)) - X(index1d(potential, i - 1)));
	double C = (ECoeff.m3[i] + ECoeff.m3_[i]*ion.D) + ECoeff.m5[i]*ion.ZxD*X(index1d(potential, i + 1)) + ECoeff.m6[i]*ion.ZxD*X(index1d(potential, i - 1)) + ECoeff.m7[i]*ion.ZxD*X(index1d(potential, i));
	double Pp1 = ECoeff.m4[i] * ion.ZxD*(X(index1d(k, i + 1))- X(index1d(k, i - 1))) + X(index1d(k, i)) * ion.ZxD*ECoeff.m5[i];
	double Pm1 = -ECoeff.m4[i] * ion.ZxD*(X(index1d(k, i + 1)) - X(index1d(k, i - 1))) + X(index1d(k, i)) * ion.ZxD*ECoeff.m6[i];
	double P = X(index1d(k, i)) * ECoeff.m7[i] * ion.ZxD;

	(this->*Assign)(Tt(index, index1d(k, i + 1), Cp1));
	(this->*Assign)(Tt(index, index1d(k, i - 1), Cm1));
	(this->*Assign)(Tt(index, index1d(k, i), C));
	(this->*Assign)(Tt(index, index1d(potential, i + 1), Pp1));
	(this->*Assign)(Tt(index, index1d(potential, i - 1), Pm1));
	(this->*Assign)(Tt(index, index1d(potential, i), P));
}

void solver::mBulkPsDerivative(int i, IonSystem& ions, void(solver::*Assign)(Tt))
{
	int index = index1d(Species::mPotential, i);

	
	double Pp1 = ECoeff.m9[i];
	double Pm1 = ECoeff.m10[i];
	double P = ECoeff.m11[i];
	double Cp = ECoeff.m12m*ions[Species::mProduct].Z;
	double Cc = ECoeff.m12m*ions[Species::mCation].Z;

	(this->*Assign)(Tt(index, index + 1, Pp1));
	(this->*Assign)(Tt(index, index - 1, Pm1));
	(this->*Assign)(Tt(index, index, P));
	(this->*Assign)(Tt(index, index1d(Species::mProduct, i), Cp));
	(this->*Assign)(Tt(index, index1d(Species::mCation, i), Cc));
}

void solver::mBottomDerivative(void(solver::*Assign)(Tt))
{
	int ei = index1d(Species::mPotential, 0);
	double y0 = Mesh.y(0);
	double adjCoeff = -pow(1 - y0, 2) / Mesh.A / Mesh.dy;

	auto J_ = [&](Species species) {
		Ion& ion = mIons[species];
		int ci = index1d(species, 0);

		double X1 = adjCoeff*ion.D;
		double X0 = adjCoeff*(-ion.D +
			ion.ZxD*eReaction.F_RT*(X(ei + 1) - X(ei)));
		double P1 = adjCoeff*ion.ZxD*X(ci)*eReaction.F_RT;
		double P0 = -P1;

		if (species == Species::mReactant) {
			double dR = Mesh.A / pow(1 - Mesh.y(0), 2)*Mesh.dy;
			double drivingPotential = eReaction.DrivingPotential(input.AppliedPotential, X(ei));

			double kf = eReaction.kf(drivingPotential);
			double kb = eReaction.kb(drivingPotential);

			double potentialTerm = -kf*eReaction.minusAlfaNF_R_T*X(index1d(Species::mProduct, 0)) + kb*eReaction.OneMinusAlfaNF_R_T*X(index1d(Species::mReactant, 0));

			P0 -= potentialTerm;
			X0 += kb;
			double Xpro0 = -kf;
			(this->*Assign)(Tt(ci, index1d(Species::mProduct, 0), Xpro0));
		}
		else if (species == Species::mProduct) {
			double dR = Mesh.A / pow(1 - Mesh.y(0), 2)*Mesh.dy;
			double drivingPotential = eReaction.DrivingPotential(input.AppliedPotential, X(ei));

			double kf = eReaction.kf(drivingPotential);
			double kb = eReaction.kb(drivingPotential);

			double potentialTerm = -kf*eReaction.minusAlfaNF_R_T*X(index1d(Species::mProduct, 0)) + kb*eReaction.OneMinusAlfaNF_R_T*X(index1d(Species::mReactant, 0));

			P0 += potentialTerm;
			double Xrea0 = -kb;
			X0 += kf;
			(this->*Assign)(Tt(ci, index1d(Species::mReactant, 0), Xrea0));
		}

		(this->*Assign)(Tt(ci, ci + 1, X1));
		(this->*Assign)(Tt(ci, ci, X0));
		(this->*Assign)(Tt(ci, ei + 1, P1));
		(this->*Assign)(Tt(ci, ei, P0));
	};

	// Cation
	J_(Species::mCation);
	// Reactant
	J_(Species::mReactant);
	// Product
	J_(Species::mProduct);

	
	// potential derivative
	double P1 = -eReaction.DrivingPotentialCoeff*adjCoeff;
	double P0 = eReaction.DrivingPotentialCoeff*adjCoeff - 1;

	(this->*Assign)(Tt(ei, ei + 1, P1));
	(this->*Assign)(Tt(ei, ei, P0));
}

void solver::mTopDerivative(void(solver::*Assign)(Tt))
{
	double adjCoeff = -pow(1 - Mesh.y(mUpperIndex), 2) / Mesh.A / Mesh.dy;
	int ei = index1d(Species::mPotential, mUpperIndex);
	double dR = Mesh.A/pow(1 - Mesh.y(mUpperIndex), 2)*Mesh.dy;

	auto F_ = [&](Ion& ion, interfaceReaction& transReaction, int ci) { //This lambda is for non charged species, particularly reactant hear

		double Ci = adjCoeff*(ion.D + ion.ZxD*eReaction.F_RT*(X(ei) - X(ei - 1))) - transReaction.kf(0);
		double Cim1 = -adjCoeff*ion.D;
		double Cip1 = transReaction.kb(0);
		double Pi = adjCoeff*ion.ZxD*X(ci)*eReaction.F_RT;
		double Pim1 = -Pi;

		(this->*Assign)(Tt(ci, ci, Ci));
		(this->*Assign)(Tt(ci, ci - 1, Cim1));
		(this->*Assign)(Tt(ci, ci + 1, Cip1));
		(this->*Assign)(Tt(ci, ei, Pi));
		(this->*Assign)(Tt(ci, ei - 1, Pim1));
	};

	auto F_charged = [&](Ion& ion, interfaceReaction& transReaction, int ci) {// This lambda is for charged species, particularly product and cation
		double drivingPotential = transReaction.drivingPotential(X(ei) - X(ei + 1));
		double kf = transReaction.kf(drivingPotential);
		double kb = transReaction.kb(drivingPotential);

		double Ci = adjCoeff*(ion.D + ion.ZxD*transReaction.F_RT*(X(ei) - X(ei - 1))) - kf;
		double Cim1 = -adjCoeff*ion.D;
		double Cip1 = kb;
		double Pi = adjCoeff*(ion.ZxD*X(ci)*transReaction.F_RT) - (X(ci)*transReaction.OneMinusAlfaNF_R_T*kf - X(ci + 1)*transReaction.minusAlfaNF_R_T*kb);
		double Pim1 = adjCoeff*(-ion.ZxD*X(ci)*transReaction.F_RT);
		double Pip1 = X(ci)*transReaction.OneMinusAlfaNF_R_T*kf - X(ci + 1)*transReaction.minusAlfaNF_R_T*kb;

		(this->*Assign)(Tt(ci, ci, Ci));
		(this->*Assign)(Tt(ci, ci-1, Cim1));
		(this->*Assign)(Tt(ci, ci+1, Cip1));
		(this->*Assign)(Tt(ci, ei, Pi));
		(this->*Assign)(Tt(ci, ei-1, Pim1));
		(this->*Assign)(Tt(ci, ei+1, Pip1));
	};

	// Reactant
	int ci = index1d(Species::mReactant, mUpperIndex);
	F_(mIons[Species::mReactant], reaReaction, ci);
	// Product
	ci = index1d(Species::mProduct, mUpperIndex);
	F_charged(mIons[Species::mProduct], proReaction, ci);
	// Cation
	ci = index1d(Species::mCation, mUpperIndex);
	F_charged(mIons[Species::mCation], catReaction, ci);

	// Potential
	/*
	(this->*Assign)(Tt(ei, ei, -1));
	(this->*Assign)(Tt(ei, ei + 1, 1));
	*/

	double Cr = mIons[Species::mReactant].Z;
	double Cc = mIons[Species::mCation].Z;
	double Cp = mIons[Species::mProduct].Z;

	(this->*Assign)(Tt(ei, index1d(Species::mReactant, mUpperIndex), Cr));
	(this->*Assign)(Tt(ei, index1d(Species::mCation, mUpperIndex), Cc));
	(this->*Assign)(Tt(ei, index1d(Species::mProduct, mUpperIndex), Cp));

}

void solver::sBulkPsDerivative(int i, IonSystem& ions, void(solver::*Assign)(Tt))
{
	int index = index1d(Species::sPotential, i);
	
	double Pp1 = ECoeff.m9[i];
	double Pm1 = ECoeff.m10[i];
	double P = ECoeff.m11[i];
	double Cp = ECoeff.m12s*ions[Species::sProduct].Z;
	double Cc = ECoeff.m12s*ions[Species::sCation].Z;
	double Ca = ECoeff.m12s*ions[Species::sAnion].Z;

	(this->*Assign)(Tt(index, index + 1, Pp1));
	(this->*Assign)(Tt(index, index - 1, Pm1));
	(this->*Assign)(Tt(index, index, P));
	(this->*Assign)(Tt(index, index1d(Species::sProduct, i), Cp));
	(this->*Assign)(Tt(index, index1d(Species::sCation, i), Cc));
	(this->*Assign)(Tt(index, index1d(Species::sAnion, i), Ca));
}

void solver::mBulkF(int i)
{
	int ei = index1d(Species::mPotential, i); // current index of the potential

	// Mass Transfer equation
	auto ionF = [&](Ion& ion, int ci) {
		F(ci) = ion.D*(ECoeff.m1[i] * X(ci + 1) + ECoeff.m2[i] * X(ci - 1)) + (ECoeff.m3[i] + ECoeff.m3_[i] * ion.D)*X(ci)
			+ ion.ZxD*(ECoeff.m4[i] * (X(ci + 1) - X(ci - 1))*(X(ei + 1) - X(ei - 1)) + X(ci)*(ECoeff.m5[i] * X(ei + 1) + ECoeff.m6[i] * X(ei - 1) + ECoeff.m7[i] * X(ei)))
			+ ECoeff.m8[i] * ion[i];

		//if (abs(F(ci)) < 1e-15) F(ci) = 0;
	};

	ionF(mIons[Species::mCation], index1d(Species::mCation, i)); 
	ionF(mIons[Species::mProduct], index1d(Species::mProduct, i));
	ionF(mIons[Species::mReactant], index1d(Species::mReactant, i));
	
	// Poisson equation
	F(ei) = ECoeff.m9[i] * X(ei + 1) + ECoeff.m10[i] * X(ei - 1) + ECoeff.m11[i] * X(ei)
			+ ECoeff.m12m*(mIons[Species::mCation].Z*X(index1d(Species::mCation, i)) 
			+ mIons[Species::mProduct].Z*X(index1d(Species::mProduct, i)) 
			+ mIons[Species::mReactant].Z*X(index1d(Species::mReactant, i)) 
			+ mIons.CxZImmobileCharge);
}

void solver::mBottomF()
{
	// potential equation
	int ei = index1d(Species::mPotential, 0);
	double adjCoeff = -pow(1 - Mesh.y(0), 2) / Mesh.A / Mesh.dy;
	double potentialGrad = -adjCoeff*(X(ei + 1) - X(ei));

	F(ei) = eReaction.potentialOHP(input.AppliedPotential, potentialGrad) - X(ei);
	
	// Nernst-Planck expression
	auto J = [&](Ion& ion, int ci) {
		F(ci) = adjCoeff*(ion.D*(X(ci + 1) - X(ci)) + 
			ion.ZxD*X(ci) * eReaction.F_RT*(X(ei + 1) - X(ei)));
	};

	double drivingPotential = eReaction.DrivingPotential(input.AppliedPotential, X(ei));
	double dR = Mesh.A / pow(1 - Mesh.y(0), 2)*Mesh.dy;
	double kf = eReaction.kf(drivingPotential);
	double kb = eReaction.kb(drivingPotential);

	double reactionRate = (kf*X(index1d(Species::mProduct, 0)) - kb*X(index1d(Species::mReactant, 0)));

	int ci = index1d(Species::mCation, 0);
	J(mIons[Species::mCation], ci);

	ci = index1d(Species::mReactant, 0);
	J(mIons[Species::mReactant], ci);
	F[ci] -= reactionRate;

	ci = index1d(Species::mProduct, 0);
	J(mIons[Species::mProduct], ci);
	F[ci] += reactionRate;
}

void solver::mTopF()
{
	// potential equation
	int ei = index1d(Species::mPotential, mUpperIndex);
	F(ei) = X(ei + 1) - X(ei);
	F(ei) = mIons[Species::mCation].Z*X(index1d(Species::mCation, mUpperIndex))
		+ mIons[Species::mProduct].Z*X(index1d(Species::mProduct, mUpperIndex))
		+ mIons[Species::mReactant].Z*X(index1d(Species::mReactant, mUpperIndex))
		+ mIons.CxZImmobileCharge;
	
	double adjCoeff = -pow(1 - Mesh.y(mUpperIndex), 2) / Mesh.A / Mesh.dy;

	auto J = [&](Ion& ion, int ci) {
		F(ci) = adjCoeff*(ion.D*(X(ci) - X(ci - 1)) +
			ion.ZxD*X(ci)*eReaction.F_RT*(X(ei) - X(ei - 1)));
	};

	double dR = Mesh.A / pow(1 - Mesh.y(mUpperIndex), 2)*Mesh.dy;
	double dE = X(ei) - X(ei + 1);
	// product interphase reaction
	double drivingPotential = proReaction.drivingPotential(dE);
	double kf = proReaction.kf(drivingPotential);
	double kb = proReaction.kb(drivingPotential);
	int ci = index1d(Species::mProduct, mUpperIndex);
	double reactionRate = (kf*X(ci) - kb*X(ci + 1));
	J(mIons[Species::mProduct], ci);
	F(ci) -= reactionRate;
	// reactant interphase reaction
	kf = reaReaction.kf(0);
	kb = reaReaction.kb(0);
	ci = index1d(Species::mReactant, mUpperIndex);
	reactionRate = (kf*X(ci) - kb*X(ci + 1));
	J(mIons[Species::mReactant], ci);
	F(ci) -= reactionRate;
	// cation interphase reaction
	drivingPotential = catReaction.drivingPotential(dE);
	kf = catReaction.kf(drivingPotential);
	kb = catReaction.kb(drivingPotential);
	ci = index1d(Species::mCation, mUpperIndex);
	reactionRate = (kf*X(ci) - kb*X(ci + 1));
	J(mIons[Species::mCation], ci);
	F(ci) -= reactionRate;
}

void solver::sBulkF(int i)
{
	int ei = index1d(Species::sPotential, i); // current index of the potential

	// Mass Transfer equation
	auto ionF = [&](Ion& ion, int ci) {
		F(ci) = ion.D*(ECoeff.m1[i] * X(ci + 1) + ECoeff.m2[i] * X(ci - 1)) + (ECoeff.m3[i] + ECoeff.m3_[i] * ion.D)*X(ci)
			+ ion.ZxD*(ECoeff.m4[i] * (X(ci + 1) - X(ci - 1))*(X(ei + 1) - X(ei - 1)) + X(ci)*(ECoeff.m5[i] * X(ei + 1) + ECoeff.m6[i] * X(ei - 1) + ECoeff.m7[i] * X(ei)))
			+ ECoeff.m8[i] * ion[i];

	};

	ionF(sIons[Species::sAnion], index1d(Species::sAnion, i));
	ionF(sIons[Species::sCation], index1d(Species::sCation, i));
	ionF(sIons[Species::sProduct], index1d(Species::sProduct, i));
	ionF(sIons[Species::sReactant], index1d(Species::sReactant, i));

	// Poisson 
	F(ei) = ECoeff.m9[i] * X(ei + 1) + ECoeff.m10[i] * X(ei - 1) + ECoeff.m11[i] * X(ei)
			+ ECoeff.m12s*(sIons[Species::sCation].Z*X(index1d(Species::sCation, i)) 
			+ sIons[Species::sProduct].Z*X(index1d(Species::sProduct, i)) 
			+ sIons[Species::sReactant].Z*X(index1d(Species::sReactant, i)) 
			+ sIons[Species::sAnion].Z*X(index1d(Species::sAnion, i)));
}

void solver::sBottomF()
{
	int ei = index1d(Species::sPotential, sBottomIndex);
	
	F(ei) = sIons[Species::sCation].Z*X(index1d(Species::sCation, sBottomIndex)) 
			+ sIons[Species::sProduct].Z*X(index1d(Species::sProduct, sBottomIndex)) 
			+ sIons[Species::sReactant].Z*X(index1d(Species::sReactant, sBottomIndex)) 
			+ sIons[Species::sAnion].Z*X(index1d(Species::sAnion, sBottomIndex));

	double adjCoeff = -pow(1 - Mesh.y(sBottomIndex), 2) / Mesh.A / Mesh.dy;

	auto J = [&](Ion& ion, int ci) {
		F(ci) = adjCoeff*(ion.D*(X(ci + 1) - X(ci)) +
			ion.ZxD*X(ci)*eReaction.F_RT*(X(ei + 1) - X(ei)));
	};

	double dR = Mesh.A / pow(1 - Mesh.y(sBottomIndex), 2)*Mesh.dy;
	double dE = X(ei - 1) - X(ei);
	// product interphase reaction
	double drivingPotential = proReaction.drivingPotential(dE);
	double kf = proReaction.kf(drivingPotential);
	double kb = proReaction.kb(drivingPotential);
	int ci = index1d(Species::sProduct, sBottomIndex);
	double reactionRate = (kf*X(ci - 1) - kb*X(ci));
	J(sIons[Species::sProduct], ci);
	F(ci) -= reactionRate;
	// reactant interphase reaction
	kf = reaReaction.kf(0);
	kb = reaReaction.kb(0);
	ci = index1d(Species::sReactant, sBottomIndex);
	reactionRate = (kf*X(ci - 1) - kb*X(ci));
	J(sIons[Species::sReactant], ci);
	F(ci) -= reactionRate;
	//Cation interphase reaction
	drivingPotential = catReaction.drivingPotential(dE);
	kf = catReaction.kf(drivingPotential);
	kb = catReaction.kb(drivingPotential);
	ci = index1d(Species::sCation, sBottomIndex);
	reactionRate = (kf*X(ci - 1) - kb*X(ci));
	J(sIons[Species::sCation], ci);
	F(ci) -= reactionRate;
	//Anion
	ci = index1d(Species::sAnion, sBottomIndex);
	J(sIons[Species::sAnion], ci);
}

void solver::sBottomDerivative(void(solver::*Assign)(Tt))
{
	int ei = index1d(Species::sPotential, sBottomIndex);
	double adjCoeff = -pow(1 - Mesh.y(sBottomIndex), 2) / Mesh.A / Mesh.dy;
	double dR = Mesh.A / pow(1 - Mesh.y(sBottomIndex), 2)*Mesh.dy;

	auto F_ = [&](Ion& ion, interfaceReaction& transReaction, int ci) {// this lambda is for noncharged species, particularly reactant
		double Cip1 = adjCoeff*ion.D;
		double Ci = adjCoeff*(-ion.D + ion.ZxD*eReaction.F_RT*(X(ei + 1) - X(ei))) + transReaction.kb(0);
		double Cim1 = -transReaction.kf(0);
		double Pip1 = adjCoeff*ion.ZxD*X(ci)*eReaction.F_RT;
		double Pi = -Pip1;

		(this->*Assign)(Tt(ci, ci + 1, Cip1));
		(this->*Assign)(Tt(ci, ci, Ci));
		(this->*Assign)(Tt(ci, ci - 1, Cim1));
		(this->*Assign)(Tt(ci, ei + 1, Pip1));
		(this->*Assign)(Tt(ci, ei, Pi));
	};

	auto F_charged = [&](Ion& ion, interfaceReaction& transReaction, int ci) {// this lambda is for charged species, particularly product and cation
		double drivingPotential = transReaction.drivingPotential(X(ei - 1) - X(ei));
		double kf = transReaction.kf(drivingPotential);
		double kb = transReaction.kb(drivingPotential);

		double Ci = adjCoeff*(-ion.D + ion.ZxD*transReaction.F_RT*(X(ei + 1) - X(ei))) + kb;
		double Cim1 = -kf;
		double Cip1 = adjCoeff*ion.D;
		double Pi = adjCoeff*(-ion.ZxD*X(ci)*transReaction.F_RT) + (X(ci - 1)*transReaction.OneMinusAlfaNF_R_T*kf - X(ci)*transReaction.minusAlfaNF_R_T*kb);
		double Pip1 = adjCoeff*ion.ZxD*X(ci)*transReaction.F_RT;
		double Pim1 = -X(ci - 1)*transReaction.OneMinusAlfaNF_R_T*kf - X(ci)*transReaction.minusAlfaNF_R_T*kb;

		(this->*Assign)(Tt(ci, ci + 1, Cip1));
		(this->*Assign)(Tt(ci, ci, Ci));
		(this->*Assign)(Tt(ci, ci - 1, Cim1));
		(this->*Assign)(Tt(ci, ei + 1, Pip1));
		(this->*Assign)(Tt(ci, ei, Pi));
		(this->*Assign)(Tt(ci, ei -1, Pim1));
	};

	// Reactant interphase reaction
	int ci = index1d(Species::sReactant, sBottomIndex);
	F_(sIons[Species::sReactant], reaReaction, ci);
	// Product interphase reaction
	ci = index1d(Species::sProduct, sBottomIndex);
	F_charged(sIons[Species::sProduct], proReaction, ci);
	// Cation interphase reaction
	ci = index1d(Species::sCation, sBottomIndex);
	F_charged(sIons[Species::sCation], catReaction, ci);

	// Anion
	auto ion = sIons[Species::sAnion];
	ci = index1d(Species::sAnion, sBottomIndex);

	(this->*Assign)(Tt(ci, ci + 1, adjCoeff*ion.D));
	(this->*Assign)(Tt(ci, ci, adjCoeff*(-ion.D + ion.ZxD*eReaction.F_RT*(X(ei + 1) - X(ei)))));
	(this->*Assign)(Tt(ci, ei + 1, adjCoeff*ion.ZxD*X(ci)*eReaction.F_RT));
	(this->*Assign)(Tt(ci, ei, -adjCoeff*ion.ZxD*X(ci)*eReaction.F_RT));

	
	// Potential
	//(this->*Assign)(Tt(ei, ei - 1, -1));
	//(this->*Assign)(Tt(ei, ei, 1));
	
	double Cr = sIons[Species::sReactant].Z;
	double Cc = sIons[Species::sCation].Z;
	double Cp = sIons[Species::sProduct].Z;
	double Ca = sIons[Species::sAnion].Z;

	(this->*Assign)(Tt(ei, index1d(Species::sReactant, sBottomIndex), Cr));
	(this->*Assign)(Tt(ei, index1d(Species::sCation, sBottomIndex), Cc));
	(this->*Assign)(Tt(ei, index1d(Species::sProduct, sBottomIndex), Cp));
	(this->*Assign)(Tt(ei, index1d(Species::sAnion, sBottomIndex), Ca));
	
}

void solver::sTopF()
{
	int ei = index1d(Species::sPotential, sUpperIndex);
	F(ei) = X(ei) - sIons[Species::sPotential].Cinitial;

	int ci = index1d(Species::sAnion, sUpperIndex);
	F(ci) = X(ci) - sIons[Species::sAnion].Cinitial;

	ci = index1d(Species::sCation, sUpperIndex);
	F(ci) = X(ci) - sIons[Species::sCation].Cinitial;

	ci = index1d(Species::sProduct, sUpperIndex);
	F(ci) = X(ci) - sIons[Species::sProduct].Cinitial;

	ci = index1d(Species::sReactant, sUpperIndex);
	F(ci) = X(ci) - sIons[Species::sReactant].Cinitial;
}

void solver::sTopDerivative(void(solver::*Assign)(Tt))
{
	int ei = index1d(Species::sPotential, sUpperIndex);
	(this->*Assign)(Tt(ei, ei, 1));

	int ci = index1d(Species::sAnion, sUpperIndex);
	(this->*Assign)(Tt(ci, ci, 1));

	ci = index1d(Species::sCation, sUpperIndex);
	(this->*Assign)(Tt(ci, ci, 1));

	ci = index1d(Species::sProduct, sUpperIndex);
	(this->*Assign)(Tt(ci, ci, 1));

	ci = index1d(Species::sReactant, sUpperIndex);
	(this->*Assign)(Tt(ci, ci, 1));
}

void solver::saveDensity() const
{
	auto SD = [&](Species species, Ion& ion, const int m, const int n) {
#pragma omp parallel for 
		for (int i = m; i < n; ++i) {
			ion[i] = X(index1d(species, i));
		}
	};

	for (int i = 0; i < 4; ++i) {
		SD(Species(i), mIons[Species(i)], 0, mUpperIndex + 1);
	}
	for (int i = 4; i < 9; ++i) {
		SD(Species(i), sIons[Species(i)], sBottomIndex, sUpperIndex + 1);
	}

	if (input.SignalState() == SignalStateSwitch::endOfPositive) {
		for (int i = 0; i < 4; ++i) {
			mIons[Species(i)].DensityLP = mIons[Species(i)].DensityCP;
			mIons[Species(i)].DensityCP = mIons[Species(i)].DensityN;
		}
		for (int i = 4; i < 9; ++i) {
			sIons[Species(i)].DensityLP = sIons[Species(i)].DensityCP;
			sIons[Species(i)].DensityCP = sIons[Species(i)].DensityN;
		}
	}
	else if (input.SignalState() == SignalStateSwitch::endOfNegative) {
		for (int i = 0; i < 4; ++i) {
			mIons[Species(i)].DensityLN = mIons[Species(i)].DensityCN;
			mIons[Species(i)].DensityCN = mIons[Species(i)].DensityN;
		}
		for (int i = 4; i < 9; ++i) {
			sIons[Species(i)].DensityLN = sIons[Species(i)].DensityCN;
			sIons[Species(i)].DensityCN = sIons[Species(i)].DensityN;
		}
	}
}

double solver::faradaicCurrent() const
{
	double adjCoeff = -(1 - Mesh.y(0))*(1 - Mesh.y(0)) / Mesh.A / Mesh.dy;
	const auto& ion = mIons[Species::mProduct];
	int ei = index1d(Species::mPotential, 0);
	int ci = index1d(Species::mProduct, 0);

	return adjCoeff*(ion.D*(X(ci + 1) - X(ci)) + ion.ZxD*X(ci)*eReaction.F_RT*(X(ei + 1) - X(ei)))*ion.Z* 96485.3329*2*3.14*Mesh.re*Mesh.re;
}

double solver::totalCurrent() const
{
	double I = -mIons.epsilon_rEpsilon_0*(sfEfield1 - sfEfield0) / input.dt * 2 * 3.14*Mesh.re*Mesh.re + faradaicCurrent();

	return I;

}



inline void solver::LockedIndexAssign(Tt triplet)
{
	omp_set_lock(&writeLock);
	int index = JIndex++;
	omp_unset_lock(&writeLock);
	JList[index] = triplet;
}

inline void solver::LockedPushBack(Tt triplet)
{
	omp_set_lock(&writeLock);
	JList.push_back(triplet);
	omp_unset_lock(&writeLock);
}