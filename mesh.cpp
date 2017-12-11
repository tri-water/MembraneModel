#include "mesh.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>

using std::cin; using std::cout;

mesh::mesh(double fA, double fmu, double fmThickness, double fend, double fre, int fn) : 
	A(fA), mu(fmu),  end(fend), mThickness(fmThickness), R(10), y(10), dy(1.0/(fn + 1)),
	Ady(fA/ (fn + 1)), re(fre)
{
	double ymu = 1.0 - A / (mu + A);
	double yend = 1.0 - A / (end + A);

	NumberOfNodes = ceil((yend - ymu) / dy);
	R.resize(NumberOfNodes);
	y.resize(NumberOfNodes);

	for (int i = 0; i < NumberOfNodes; ++i) {
		y[i] = i*dy + ymu;
		R[i] = A / (1.0 - y[i]) - A + re;
	}

	IndexOfBoundary = calculateIndexOfBoundary();

	// Print the minimum mesh size
	cout << "Minimum Mesh Size: " << R[1] - R[0] << "\n"
		<< "Total Number of Nodes: " << NumberOfNodes << "\n";
	//cout << "Continue Yes(y)/No(n):";
	//char control;
	//cin >> control;
	//if (control == 'n' || control == 'N') exit(0);
}

int mesh::calculateIndexOfBoundary()
{
	int index = 0;
	double boundaryR = mThickness + re + mu;
	for (int i = 0; i < NumberOfNodes; ++i) {
		if (R[i] >= boundaryR) {
			index = i;
			break;
		}
	}
	return index;
}

void mesh::printMesh(std::string ex_file_name)
{
	int boundaryIndex = IndexOfBoundary;

	std::ofstream myfile1(ex_file_name + "membrane mesh-R.txt");
	std::ofstream myfile2(ex_file_name + "membrane mesh-y.txt");
	std::ofstream myfile3(ex_file_name + "solution mesh-R.txt");
	std::ofstream myfile4(ex_file_name + "solution mesh-y.txt");

	if (myfile1.is_open() && myfile2.is_open() && myfile3.is_open() && myfile4.is_open()) {
		myfile1 << R.topRows(boundaryIndex);
		myfile3 << R.bottomRows(NumberOfNodes - boundaryIndex);
		myfile2 << y.topRows(boundaryIndex);
		myfile4 << y.bottomRows(NumberOfNodes - boundaryIndex);

		myfile1.close();
		myfile2.close();
		myfile3.close();
		myfile4.close();
	}
	else {
		std::cerr << " Unable to open mesh data file" << std::endl;
		exit(EXIT_FAILURE);
	}
}