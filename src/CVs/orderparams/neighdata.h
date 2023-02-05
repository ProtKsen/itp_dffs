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
   NeighData(const ParticleSystem& psystem, const Lattice& lattice, const SSAGES::Snapshot& snapshot, const vector<CVPCLASS>& cvpclass);

   // we store the number of crystal neighbours for each particle 
   std::vector<int> numneigh;
   std::vector<int> numneigh_local;

};

std::vector<CVNCLASS> classifynodes(const Lattice&, const NeighData&);
std::vector<CVPCLASS> classifypars(const ParticleSystem&);
std::vector<int> largestnodescluster(const Lattice&, const vector<CVNCLASS>&);
graph getnodegraph(const Lattice&, const vector<int>&);
bool isnodesneigh(const Node&, const Node&, const Lattice&, double&);
double sepnodes(const Node&, const Node&, const int, 
                  const int, const int, const bool);

#endif
