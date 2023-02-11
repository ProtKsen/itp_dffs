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

NeighData::NeighData(const ParticleSystem& psystem, const Lattice& lattice, 
                     const SSAGES::Snapshot& snapshot, const vector<CVPCLASS>& cvpclass)
{
   // store number of liquid neighbours and neighbour list
   int Ntot = lattice.allnodes.size();

   numneigh_local.resize(Ntot, 0); // num liquid neighbours for each node in proc
   numneigh.resize(Ntot, 0); // num liquid neighbours for each node

   numneigh_local = getneigh_fast(lattice, psystem.allpars, lattice.r_sep, cvpclass);
   
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


vector<CVNCLASS> classifynodes(const Lattice& lattice, const NeighData& neighdata)
{
   int n = lattice.allnodes.size();
   vector<CVNCLASS> nodeclass(n, VAP);

   for (int i = 0; i < n; ++i)
   {
      if (neighdata.numneigh[i] >= 1)
      {
        nodeclass[i] = FLUID;
      }
   }
   return nodeclass;
}

vector<CVPCLASS> classifypars(const ParticleSystem& psystem)
{
   int n = psystem.allpars.size();
   vector<CVPCLASS> parclass(n, VAPOUR);
   
   vector<double> neighs(n, 0);
   for (int i = 0; i < n; ++i)
   {    
        for (int j = 0; j < n; ++j)
        {
            if (i != j)
            {
               double dist;
               dist = psystem.simbox.sepsq(psystem.allpars[i], psystem.allpars[j]);
               if (dist <= psystem.r_near_neigh)
               {
                  neighs[i]++;
               }
               if (neighs[i] >= psystem.num_neighbours_liquid)
               {
                  parclass[i] = FLU;
                  break;
               }
            }
        }
   }
   return parclass;
}

double sepnodes(const Node& p1, const Node& p2, const int lboxx, 
                  const int lboxy, const int lboxz, const bool periodicz)
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

   double dist =  sepx * sepx + sepy * sepy + sepz * sepz;
   return dist;
}

bool isnodesneigh(const Node& n1, const Node& n2, const Lattice& lattice, double& rsq)
{
   double s[3];

   // compute separation between nodes (store in s)
   double dist = sepnodes(n1, n2, lattice.lboxx, lattice.lboxy, lattice.lboxz, lattice.zperiodic);

   if (dist <= lattice.lcellx * lattice.lcellx + lattice.lcelly * lattice.lcelly + lattice.lcellz * lattice.lcellz)
        return true;
   else
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


vector<int> largestnodescluster(const Lattice& lattice, const vector<CVNCLASS>& cvnclass)
{
   // get vector with indices that are all vapor nodes

   vector<int> xps;
   for (int i = 0; i < lattice.allnodes.size(); ++i) {
      if (cvnclass[i] == VAP) {
         xps.push_back(i);
      }
   }

   // graph of xtal particles, with each particle a vertex and each
   // link an edge
   graph xgraph = getnodegraph(lattice, xps);

   // largest cluster is the largest connected component of graph
   vector<int> cnums = largestcomponent(xgraph);

   //std::string FileTest="Test.txt";
   //std::ofstream fout_test(FileTest, std::ios_base::out | std::ios_base::app);
   //for (int i = 0; i < cnums.size(); ++i) {
   // fout_test << cnums[i] << " ";
   //}
   //fout_test << std::endl;
   //fout_test.close();

   // now largest component returns indexes into array xps, we need
   // to reindex so that it contains indices into psystem.allpars
   // (see utility.cpp)
   reindex(cnums, xps);
   return cnums;    
}
