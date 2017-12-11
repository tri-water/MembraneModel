#include "mesh.h"
#include <iostream>
#include "IonSystem.h"
#include "thermodynamics.h"
#include <map>
#include "potentialSignal.h"
#include "solver.h"
#include "equationCoefficients.h"
#include "utility.h"
#include <string>
#include "ErrorTypes.h"
#include <sstream>

using std::cin; using std::cout;
typedef std::pair<Species, Ion*> ionPair;

int main(int argc, char* argv[])
{
	/*
	unit of length: m
	unit of diffusion coefficient: m^2 s^-1
	unit of concentration: mol L^-3
	unit of standard rate const: m s^-1
	*/
	double solutionC(0); // solution phase concentration
	if (argc == 2) {
		std::stringstream ss(argv[1]);
		try {
			ss >> solutionC;
			cout << "Solution electrolyte concentration: " << solutionC << std::endl;
		}
		catch (const std::exception& e) {
			std::cerr << e.what() << '\n' << "Check command line.";
			return 1;
		}
	}
	else {
		std::cerr << "Number of Commandline arguments not right.\n";
		return 1;
	}

	mesh Mesh(2e-7, 1e-8, 5e-5, 100e-2, 5e-3, 100000);
	int meshSize = static_cast<int>(Mesh.R.size());
	int boundaryIndex = Mesh.IndexOfBoundary;
	Mesh.printMesh("mesh ");

	thermoDynamics ElectroThermo(0, 1, 0.01, 0.5);
	electrodeReaction EletroR(12.0, 10.0, 5.0, 3e-9, ElectroThermo, Mesh);
	thermoDynamics ProThermo(-0.1, 1, 0.01, 0.5);
	interfaceReaction ProductIR(ProThermo);
	thermoDynamics CatThermo(-0.1, 1, 0.01, 0.5);
	interfaceReaction CationIR(CatThermo);
	nonChargedInterfaceReaction ReactantIR(0.0001, 0.001);

	map<Species, Ion*> mIons, sIons;
	Ion mReactant(1e-13, 0, 1, 0, boundaryIndex-1); mIons.insert(ionPair(Species::mReactant, &mReactant));
	Ion mProduct(1e-12, 1, 1*ElectroThermo.ratio_ox2red(-0.6), 0, boundaryIndex - 1); mIons.insert(ionPair(Species::mProduct, &mProduct));
	Ion mCation(1e-12, 1, 1 - 1*ElectroThermo.ratio_ox2red(-0.6), 0, boundaryIndex - 1); mIons.insert(ionPair(Species::mCation, &mCation));
	Ion mPotential(0, 0, 0, 0, boundaryIndex - 1); mIons.insert(ionPair(Species::mPotential, &mPotential));
	ImmobileCharge immobileCharge(-1, 1);

	Ion sReactant(1e-9, 0, 0, boundaryIndex, meshSize - 1); sIons.insert(ionPair(Species::sReactant, &sReactant));
	Ion sProduct(1e-9, 1, 0, boundaryIndex, meshSize - 1); sIons.insert(ionPair(Species::sProduct, &sProduct));
	Ion sCation(1e-9, 1, solutionC, boundaryIndex, meshSize - 1); sIons.insert(ionPair(Species::sCation, &sCation));
	Ion sAnion(1e-9, -1, solutionC, boundaryIndex, meshSize - 1); sIons.insert(ionPair(Species::sAnion, &sAnion));
	Ion sPotential(0, 0, 0, boundaryIndex, meshSize - 1); sIons.insert(ionPair(Species::sPotential, &sPotential));

	IonSystem mIonSys(12.0, mIons, mPotential, immobileCharge);
	IonSystem sIonSys(80.0, sIons, sPotential);

	// the time step of each input potential signal must be identical
	SquareWave Esignal(-0.6, 0.3, 0.002, 0.0005, 25, 0.04); // to initialise swSolver and record faradaic current
	SquareWave Esignal_total(-0.6, 0.3, 0.002, 0.0005, 25, 0.04); // to record total current
	LinearPotential LinearEsignal(-0.6, 0, 0.1, 0.0005);
	ConstantPotential statSignal(-0.6, 0.0005, 1);

	equationCoefficients eqCoeff(Mesh, statSignal, ElectroThermo, mIonSys, sIonSys);
	solver SWSolver(mIonSys, sIonSys, Mesh, eqCoeff, EletroR, ProductIR, ReactantIR, CationIR, Esignal);
	solver LinearSolver(mIonSys, sIonSys, Mesh, eqCoeff, EletroR, ProductIR, ReactantIR, CationIR, LinearEsignal);
	solver StatSolver(mIonSys, sIonSys, Mesh, eqCoeff, EletroR, ProductIR, ReactantIR, CationIR, statSignal);

	oneDIndex index1d(Mesh);
	twoDIndex index2d(Mesh, index1d);
	
	LinearSolver.initialise();
	StatSolver.initialise();
	
	for (int i = 0; i < statSignal.GetPeriodNumber(); ++i) {
		statSignal.CalculateAppliedPotential(i);
		try {
			StatSolver.solve();
		}
		catch (DecompositionFailed&) {
			std::cout << "DecompositionFailed at stat\n";
			break;
		}
		
		statSignal.RecordCurrent(StatSolver.totalCurrent());
	}
	
	SWSolver.initialise();
	for (int i = 0; i < Esignal.GetPeriodNumber(); ++i) {
		Esignal.CalculateAppliedPotential(i);
		Esignal_total.CalculateAppliedPotential(i);
		try {
			SWSolver.solve();
		}
		catch (DecompositionFailed&) {
			std::cout << "Decomposition Failed at SWV\n";
			break;
		}
		
		Esignal.RecordCurrent(SWSolver.faradaicCurrent());
		Esignal_total.RecordCurrent(SWSolver.totalCurrent());

		if (Esignal.IsPeak()) {
			std::cout << "peak found\n";
			std::vector<std::string> IonNames =
			{ "mReactant", "mProduct", "mCation", "mPotential", "sReactant", "sProduct", "sCation", "sPotential", "sAnion" };
			
			for (int i = 0; i < 4; ++i) {
				mIonSys[Species(i)].printLastPeriodDense("Peak Current Density - " + IonNames[i] + ".txt");
			}
			for (int i = 4; i < 9; ++i) {
				sIonSys[Species(i)].printLastPeriodDense("Peak Current Density - " + IonNames[i] + ".txt");
			}
		}

		if (Esignal_total.IsPeak()) {
			std::cout << "total peak found\n";
			std::vector<std::string> IonNames =
			{ "mReactant", "mProduct", "mCation", "mPotential", "sReactant", "sProduct", "sCation", "sPotential", "sAnion" };

			for (int i = 0; i < 4; ++i) {
				mIonSys[Species(i)].printLastPeriodDense("Total Peak Current Density - " + IonNames[i] + ".txt");
			}
			for (int i = 4; i < 9; ++i) {
				sIonSys[Species(i)].printLastPeriodDense("Total Peak Current Density - " + IonNames[i] + ".txt");
			}
		}
	}

	Esignal.ExportCurrent("SW Faradaic Current.txt");
	Esignal_total.ExportCurrent("SW Total Current.txt");
	statSignal.ExportCurrent("Constant Potential Total Current.txt");
	utility::exportConstantPotentialParameters("setup parameters.txt", statSignal);
	utility::exportSquareWaveParameters("setup parameters.txt", Esignal);
	utility::exportMeshParameters("setup parameters.txt", Mesh);
	utility::exportIonSystemParameters("setup parameters.txt", mIonSys, sIonSys);
	utility::exportReactionParameters("setup parameters.txt", ElectroThermo, ProductIR, CationIR, ReactantIR);
	

	return 0;
}