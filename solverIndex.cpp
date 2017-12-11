#include "solverIndex.h"


oneDIndex::oneDIndex(const mesh& Mesh):
	maxIndex(0), boundaryIndex(Mesh.IndexOfBoundary)
{
	InnerIndex.resize(static_cast<int> (Species::Count));

	for (int k = static_cast<int>(Species::mReactant); k < static_cast<int> (Species::sReactant); ++k) {
		InnerIndex[k].resize(Mesh.IndexOfBoundary);
		for (int i = 0; i < Mesh.IndexOfBoundary; ++i) {
			InnerIndex[k][i] = i + k*Mesh.NumberOfNodes;
			++maxIndex;
		}
	}

	for (int k = static_cast<int>(Species::sReactant); k < static_cast<int> (Species::sAnion); ++k) {
		InnerIndex[k].resize(Mesh.NumberOfNodes - Mesh.IndexOfBoundary);
		for (int i = 0; i < Mesh.NumberOfNodes - Mesh.IndexOfBoundary; ++i) {
			InnerIndex[k][i] = i + (k - 4)*Mesh.NumberOfNodes + Mesh.IndexOfBoundary;
			++maxIndex;
		}
	}

	int k = static_cast<int>(Species::sAnion);
	InnerIndex[k].resize(Mesh.NumberOfNodes - Mesh.IndexOfBoundary);
	for (int i = 0; i < Mesh.NumberOfNodes - Mesh.IndexOfBoundary; ++i) {
		InnerIndex[k][i] = i + (k - 4)*Mesh.NumberOfNodes;
		++maxIndex;
	}

}

twoDIndex::twoDIndex(const mesh& Mesh, const oneDIndex index1d)
{
	int lenth = index1d.maxIndex + 1;
	InnerIndex.resize(lenth);
	auto& indexContainer = index1d.InnerIndex;
	
	int k = 0;
	for (auto kt = indexContainer.begin(); kt != indexContainer.end(); ++kt, ++k) {
		int i = 0;
		for (auto it = (*kt).begin(); it != (*kt).end(); ++it, ++i) {
			if (k < static_cast<int>(Species::sReactant)) {
				InnerIndex[*it] = std::make_pair(static_cast<Species>(k), i);
			}
			else {
				InnerIndex[*it] = std::make_pair(static_cast<Species>(k), i + Mesh.IndexOfBoundary);
			}
		}
	}

}