#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include "boost/multi_array.hpp"
#include "qdata.h"
#include "neighdata.h"
#include "box.h"
#include "particle.h"
#include "qlmfunctions.h"
#include "constants.h"
#include "conncomponents.h"
#include "utility.h"
#include "typedefs.h"
#include "mpi.h"
#include "Snapshot.h"
#include "particlesystem.h"
#include "lattice.h"
#include "getneigh.h"
#include "typedefs.h"
#include <stdio.h> 
#include <time.h> 

using std::vector;
using std::complex;
using std::norm;

// Constructor for NieghData object

NeighData::NeighData(const ParticleSystem& psystem, const Lattice& lattice, const SSAGES::Snapshot& snapshot)
{
   // store number of neighbours and neighbour list
   int Ntot = lattice.allnodes.size();

   numneigh_local.resize(Ntot, 0); // num neighbours for each node in proc
   numneigh.resize(Ntot, 0); // num neighbours for each node
   //numneigh = getneigh_fast(lattice, psystem.allpars, lattice.nsep);

   numneigh_local = getneigh_fast(lattice, psystem.pars, lattice.nsep);
   MPI_Barrier(snapshot.GetCommunicator()); 
   
	int temp_numneigh_local[Ntot];
	int temp_numneigh[Ntot];
	for (int i=0; i < Ntot; ++i) 
   {
	   temp_numneigh_local[i] = numneigh_local[i]; 
	}
   MPI_Allreduce(&temp_numneigh_local, &temp_numneigh, Ntot, MPI_INT, MPI_SUM, snapshot.GetCommunicator());
   for (int i=0; i < Ntot; ++i) 
   {
	   numneigh[i] = temp_numneigh[i];
	}

}


vector<VCCLASS> classifynodes(const Lattice& lattice, const NeighData& neighdata)
{
   int n = lattice.allnodes.size();
   vector<VCCLASS> nodeclass(n, VAP);

   for (int i = 0; i < n; ++i)
   {
      if (neighdata.numneigh[i] >= lattice.ncrneigh)
      {
         nodeclass[i] = CRYST;
      }
   }
   return nodeclass;
}

inline void sepnodes(const Node& p1, const Node& p2, const int lboxx, 
                  const int lboxy, const int lboxz, const bool periodicz,  double* s)
{
   double sepx,sepy,sepz;
   sepx = p1.pos[0] - p2.pos[0];
   sepy = p1.pos[1] - p2.pos[1];
   sepz = p1.pos[2] - p2.pos[2];
     
   if (sepx > 0.5 * lboxx) {
      sepx = sepx - lboxx;
   }
   else if (sepx < -0.5 * lboxx) {
      sepx = sepx + lboxx;
   }
   if (sepy > 0.5 * lboxy) {
      sepy = sepy - lboxy;
   }
   else if (sepy < -0.5 * lboxy) {
      sepy = sepy + lboxy;
   }
   if (periodicz) {               
      if (sepz > 0.5 * lboxz) {
         sepz = sepz - lboxz;
      }
      else if (sepz < -0.5 * lboxz) {
         sepz = sepz + lboxz;
      }
   }

   s[0] = sepx;
   s[1] = sepy;
   s[2] = sepz;
   return;
}

bool isnodesneigh(const Node& n1, const Node& n2, const Lattice& lattice, double& rsq)
{
   double s[3];

   // compute separation between nodes (store in s)
   sepnodes(n1, n2, lattice.lboxx, lattice.lboxy, lattice.lboxz, lattice.zperiodic, s);

   if ((std::fabs(s[0]) < lattice.lcellx * 2) && (std::fabs(s[1]) < lattice.lcelly * 2)
       && (std::fabs(s[2]) < lattice.lcellz * 2)) {
      rsq = s[0] * s[0] + s[1] * s[1] + s[2] * s[2];
      double dist = lattice.lcellx * lattice.lcellx + 
                  lattice.lcelly * lattice.lcelly + lattice.lcellz * lattice.lcellz;
      if (rsq <= dist ) {
         return true;
      }
   }
   return false;
}

graph getnodegraph(const Lattice& lattice,
                const vector<int>& xpars)
{
   graph G;
   vector<int>::size_type nxtal = xpars.size();
   vector<int>::size_type i,j;
   double sep;
   
   for (i = 0; i != nxtal; ++i) {
      for (j = i + 1; j != nxtal; ++j) {
         if (isnodesneigh(lattice.allnodes[xpars[i]], lattice.allnodes[xpars[j]], lattice, sep)) {
            add_edge(i, j, G);
         }
      }
   }
   return G;
}


vector<int> largestnodescluster(const Lattice& lattice, const vector<VCCLASS>& vcclass)
{
   // get vector with indices that are all vapor nodes

   vector<int> xps;
   for (int i = 0; i < lattice.allnodes.size(); ++i) {
      if (vcclass[i] == VAP) {
         xps.push_back(i);
      }
   }
   
   // graph of xtal particles, with each particle a vertex and each
   // link an edge
   graph xgraph = getnodegraph(lattice, xps);

   // largest cluster is the largest connected component of graph
   vector<int> cnums = largestcomponent(xgraph);

   // now largest component returns indexes into array xps, we need
   // to reindex so that it contains indices into psystem.allpars
   // (see utility.cpp)
   reindex(cnums, xps);
   return cnums;    
}
