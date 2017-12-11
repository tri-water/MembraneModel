#pragma once
#pragma once
#include <string>
#include <Eigen/Dense>

class mesh
{
public:
	mesh(double fA, double fmu, double mThickness, double fend, double fre, int fn);
	int NumberOfNodes; // number of nodes
	const double A; // parameter controlling grid density
	Eigen::VectorXd R; // the real distance from original point
	Eigen::VectorXd y; // the transformed coordinates
	const double mu; // the thickness of Helmhotz layer
	const double end; // the farther end of the simulated area
	const double mThickness; // the thickness of the membrane phase
	const double dy;
	const double Ady; // A*dy
	const double re; // the radius of the electrode

	int IndexOfBoundary; // return the index of a specific distance
	void printMesh(std::string ex_file_name);

private:
	int calculateIndexOfBoundary();
};