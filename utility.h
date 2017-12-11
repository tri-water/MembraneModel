#pragma once
#include <string>
#include "IonSystem.h"
#include "mesh.h"
#include "thermodynamics.h"
#include "potentialSignal.h"
namespace utility
{
	void exportIonSystemParameters(std::string fileName, const IonSystem& mIonSys, const IonSystem& sIonSys);
	void exportReactionParameters(std::string fileName, const thermoDynamics& ElectroRThermo, const interfaceReaction& ProductIR,
		const interfaceReaction& CationIR, const interfaceReaction& ReactantIR);
	void exportMeshParameters(std::string fileName, const mesh& Mesh);
	void exportSquareWaveParameters(std::string fileName, const SquareWave& SWSignal);
	void exportLinearScanParameters(std::string fileName, const LinearPotential& LinearSignal);
	void exportConstantPotentialParameters(std::string fileName, const ConstantPotential& ConstPotential);
}