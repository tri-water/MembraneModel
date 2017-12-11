#include "utility.h"
#include <fstream>
#include "EnumStruct.h"
#include <vector>
#include <map>

namespace utility
{
	void exportIonSystemParameters(std::string fileName, const IonSystem& mIonSys, const IonSystem& sIonSys)
	{
		std::ofstream file;
		if (!std::ifstream(fileName)) {
			file.open(fileName);
		}
		else {
			file.open(fileName, std::ios_base::app);
		}

		file << "membrane ion system parameters\n";

		std::vector<std::string> allIonNames = { "mReactant", "mProduct", "mCation", "mPotential", "sReactant", 
			"sProduct", "sCation", "sPotential", "sAnion" };
		
		std::vector<Species> mIonEnum = { Species::mReactant, Species::mProduct, Species::mCation };

		for (auto i : mIonEnum) {
			auto& ion = mIonSys[i];
			file << allIonNames[static_cast<int>(i)] << ":\n";
			file << "diffusion coefficient: " << ion.D << " m^2 s^-1\n"
				<< "number of charges: " << ion.Z << "\n"
				<< "inital concentration: " << ion.Cinitial << " mol L^-1\n\n";
		}

		std::vector<Species> sIonEnum = { Species::sReactant, Species::sProduct, Species::sCation, Species::sAnion };

		for (auto i : sIonEnum) {
			auto& ion = sIonSys[i];
			file << allIonNames[static_cast<int>(i)] << ":\n";
			file << "diffusion coefficient: " << ion.D << " m^2 s^-1\n"
				<< "number of charges: " << ion.Z << "\n"
				<< "inital concentration: " << ion.Cinitial << " mol L^-1\n\n";
		}

		file << "initial potential of the membrane phase: " << mIonSys[Species::mPotential].Cinitial << " /V\n"
			<< "initial potential of the solution phase: " << sIonSys[Species::sPotential].Cinitial << " /V\n\n";

		file << "immobile charge (C*Z) of the membrane phase: " << mIonSys.CxZImmobileCharge << " mol L^-1\n\n";

		file << "relative permeability of the membrane phase: " << mIonSys.epsilon_rEpsilon_0 / 8.854187817620e-12 << "\n"
			<< "relative permeability of the solution phase : " << sIonSys.epsilon_rEpsilon_0 / 8.854187817620e-12 << "\n\n";

		file << "-------------------------------------------------------\n";
		file.close();
	}

	void exportReactionParameters(std::string fileName, const thermoDynamics& ElectroRThermo, const interfaceReaction& ProductIR,
		const interfaceReaction& CationIR, const interfaceReaction& ReactantIR)
	{
		std::ofstream file;
		if (!std::ifstream(fileName)) {
			file.open(fileName);
		}
		else {
			file.open(fileName, std::ios_base::app);
		}

		file << "reaction parameters\n\n";

		file << "electrode reaction parameters:\n";
		file << "the number of elecrons transferred: " << ElectroRThermo.n << "\n"
			<< "formal electrode potential: " << ElectroRThermo.E_formal << " /V\n"
			<< "reaction rate constant: " << ElectroRThermo.k0 << " /m s^-1\n"
			<< "alpha: " << -ElectroRThermo.minusAlfaNF_R_T / ElectroRThermo.nF_R_T << "\n\n";

		file << "product inter phase reaction parameters:\n";
		file << "reaction constant: " << ProductIR.kf(0) << " /m s^-1\n"
			<< "formal electrode potential: " << ProductIR.dE_formal << " /E\n\n";

		file << "cation inter phase reaction parameters:\n";
		file << "reaction constant: " << CationIR.kf(0) << " /m s^-1\n"
			<< "frmal electrode potential: " << CationIR.dE_formal << " /V\n\n";

		file << "reactant inter phase reaction parameters:\n";
		file << "forward reaction constant: " << ReactantIR.kf(0) << " /m s^-1\n"
			<< "reverse reaction constant: " << ReactantIR.kb(0) << " /m s^-1\n\n";

		file << "-------------------------------------------------------\n";
		file.close();
	}

	void exportMeshParameters(std::string fileName, const mesh& Mesh)
	{
		std::ofstream file;
		if (!std::ifstream(fileName)) {
			file.open(fileName);
		}
		else {
			file.open(fileName, std::ios_base::app);
		}

		file << "mesh parameters\n\n";
		file << "number of nodes: " << Mesh.NumberOfNodes << " \n"
			<< "the thickness of Helmholtz layer: " << Mesh.mu << " /m \n"
			<< "the farther end of the simulated area: " << Mesh.end << " /m \n"
			<< "the thickness of the membrane: " << Mesh.mThickness << " /m \n"
			<< "the radius of the electrode: " << Mesh.re << " /m \n";

		file << "-------------------------------------------------------\n";
		file.close();
	}

	void exportSquareWaveParameters(std::string fileName, const SquareWave& SWSignal)
	{
		std::ofstream file;
		if (!std::ifstream(fileName)) {
			file.open(fileName);
		}
		else {
			file.open(fileName, std::ios_base::app);
		}

		auto parameters = SWSignal.exportParams();

		file << "square wave parameters\n";
		file << "initial potential: " << parameters[0] << " /V\n"
			<< "end potential: " << parameters[1] << " /V\n"
			<< "cycle potential step: " << parameters[2] << " /V\n"
			<< "time node: " << parameters[3] << " /s\n"
			<< "frequency: " << static_cast<int>(parameters[4]) << " Hz\n"
			<< "amplitude: " << parameters[5] << " /V\n";

		file << "-------------------------------------------------------\n";
		file.close();
	}

	void exportLinearScanParameters(std::string fileName, const LinearPotential& LinearSignal)
	{
		std::ofstream file;
		if (!std::ifstream(fileName)) {
			file.open(fileName);
		}
		else {
			file.open(fileName, std::ios_base::app);
		}

		auto parameters = LinearSignal.exportParams();

		file << "linear scan\n";
		file << "initial potential: " << parameters[0] << " /V\n"
			<< "end potential: " << parameters[1] << " /V\n"
			<< "scan rate: " << parameters[2] << " / Vs^-1\n"
			<< "time node: " << parameters[3] << " /s\n";

		file << "-------------------------------------------------------\n";
		file.close();
	}

	void exportConstantPotentialParameters(std::string fileName, const ConstantPotential& ConstPotential)
	{
		std::ofstream file;
		if (!std::ifstream(fileName)) {
			file.open(fileName);
		}
		else {
			file.open(fileName, std::ios_base::app);
		}

		auto parameters = ConstPotential.exportParams();

		file << "constant potential\n";
		file << "initial potential: " << parameters[0] << " /V\n"
			<< "time node: " << parameters[1] << " /s\n"
			<< "Durance: " << parameters[2] << " /s\n";

		file << "-------------------------------------------------------\n";
		file.close();
	}
}