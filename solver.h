#pragma once
#include "IonSystem.h"
#include <Eigen/Sparse>
#include "equationCoefficients.h"
#include <vector>
#include <omp.h>
#include "solverIndex.h"
#include <fstream>

typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SpMatrixXd;
typedef Eigen::SparseVector<double> SpVectorXd;
typedef SpMatrixXd::InnerIterator InnerIterator;
typedef Eigen::Triplet<double> Tt;

class solver
{
public:
	solver(IonSystem& fmIons, IonSystem& fsIons, mesh& fMesh, equationCoefficients& fECoeff, electrodeReaction& feReaction,
		interfaceReaction& fproReaction, interfaceReaction& freaReaction, interfaceReaction& fcatReaction, potentialSignal& finput);
	void initialise();
	void solve();
	double faradaicCurrent() const;
	double solver::totalCurrent() const;

private:
	SpMatrixXd J; // JdX = F for Fc, each element inside is initialised as 0, Jocobian Matrix
	Eigen::VectorXd dX; // the change of x in an iteration
	Eigen::VectorXd X; // the calculated results;
	Eigen::VectorXd F; // values of Newton-Raphson functions 
	Eigen::SparseLU<SpMatrixXd> dXSolver;

	IonSystem& mIons;
	IonSystem& sIons;
	mesh& Mesh;
	const equationCoefficients& ECoeff;
	const oneDIndex index1d;
	electrodeReaction& eReaction;
	interfaceReaction& proReaction;
	interfaceReaction& reaReaction;
	interfaceReaction& catReaction;
	potentialSignal& input;
	double sfEfield0; // electrode surface electrical field of the last time
	double sfEfield1; // electrode surface electrical field of the current time
	
	int NumOfVars; // the number of variables for the system
	std::vector<Tt> JList;
	omp_lock_t writeLock;
	int JIndex;
	int sBottomIndex;
	int mUpperIndex;
	int sUpperIndex;

	void initialiseJ(void (solver::*Assign)(Tt));
	void updateJ();
	void LockedPushBack(Tt triplet);
	void LockedIndexAssign(Tt triplet);
	void bulkMTDerivative(int i, IonSystem& ions, Species k, Species potential, void (solver::*Assign)(Tt)); // mass transfer derivatives
	void mBulkPsDerivative(int i, IonSystem& ions, void(solver::*Assign)(Tt)); // Poisson equation derivatives
	void mBottomDerivative(void(solver::*Assign)(Tt)); // bottom derivative
	void mTopDerivative(void(solver::*Assign)(Tt));
	void sBulkPsDerivative(int i, IonSystem& ions, void(solver::*Assign)(Tt));
	void sBottomDerivative(void(solver::*Assign)(Tt));
	void sTopDerivative(void(solver::*Assign)(Tt));


	void initialiseX();

	void calculateF();
	void mBulkF(int i);
	void mBottomF();
	void mTopF();
	void sBulkF(int i);
	void sBottomF();
	void sTopF();

	void saveDensity() const;

	template <typename T>
	void exportEigen(T eigen, std::string title);
};

template <typename T>
void solver::exportEigen(T eigen, std::string title)
{
	std::ofstream fout(title + ".txt");

	if (fout.is_open()) fout << eigen;
	else std::cout << "Failed to exportEigen";
}