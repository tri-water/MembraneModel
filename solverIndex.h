#pragma once
#include <vector>
#include "EnumStruct.h"
#include "mesh.h"
#include <utility>
#include <iostream>

class oneDIndex
{
	friend class twoDIndex;
public:
	oneDIndex(const mesh& Mesh);
	int operator () (Species species, int i) const;
	int maxIndex;
private:
	std::vector<std::vector<int>> InnerIndex;
	int boundaryIndex;
};

class twoDIndex
{
	typedef std::pair<Species, int> indexPair;
public:
	twoDIndex(const mesh& Mesh, const oneDIndex index1d);
	indexPair operator () (int oneDIndex) const { return InnerIndex[oneDIndex]; }
private:
	std::vector<indexPair> InnerIndex;
};

inline int oneDIndex::operator()(Species species, int i) const
{
	if (static_cast<int>(species) < static_cast<int>(Species::sReactant)) {
		return InnerIndex[static_cast<int>(species)][i];
	}
	else {
		return InnerIndex[static_cast<int>(species)][i - boundaryIndex];
	}
	
}