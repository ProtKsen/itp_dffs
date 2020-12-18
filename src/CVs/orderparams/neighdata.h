#ifndef NEIGHDATA_H
#define NEIGHDATA_H

#include <vector>
#include "particlesystem.h"
#include "typedefs.h"
#include "constants.h"
#include "Snapshot.h"
#include "lattice.h"

// NeighData is a struct to store data of neighbours

struct NeighData
{
public:
   NeighData(const ParticleSystem& psystem, const Lattice& lattice, const SSAGES::Snapshot& snapshot);

   // we store the number of crystal neighbours for each particle 
   std::vector<int> numneigh;
   std::vector<int> numneigh_local;

};

std::vector<VCCLASS> classifynodes(const Lattice&, const NeighData&);
std::vector<int> largestnodescluster(const Lattice&, const vector<VCCLASS>&);
graph getnodegraph(const Lattice&, const vector<int>&);
bool isnodesneigh(const Node&, const Node&, const Lattice&, double&);
inline void sepnodes(const Node&, const Node&, const int, 
                  const int, const int, const bool,  double*);

#endif
