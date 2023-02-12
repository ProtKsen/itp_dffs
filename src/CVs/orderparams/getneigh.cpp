#include "boost/multi_array.hpp"
#include <complex>
#include <iostream>
#include <fstream>
#include <math.h>
#include "constants.h"
#include "particle.h"
#include "node.h"
#include "box.h"
#include "opfunctions.h"
#include "getneigh.h"
#include <string>

using std::complex;
using std::vector;

float sepnp(const Node& p1, const Particle& p2, const int lboxx, 
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

   double dist = sepx * sepx + sepy * sepy + sepz * sepz;
   return dist;
}

vector<int> getneigh(const Lattice& lattice, const vector<Particle>& allpars, const int nsep)
{
   int n = lattice.nodes.size(); // the number of nodes in proc
   int Ntot = allpars.size(); // number of all particles
   std::vector<int> numneigh;
   numneigh.resize(n, 0);

   double s[3];
   // compute separation between particles (store in s)
   for (int i = 0; i < n; ++i)
   {
      numneigh[i] = 0;
      for (int j = 0; j < Ntot; ++j)
      {
        double dist = sepnp(lattice.nodes[i], allpars[j], lattice.lboxx, 
               lattice.lboxy, lattice.lboxz, lattice.zperiodic);
         if (dist <= nsep * nsep)
         {
            ++numneigh[i];
         }
      }
   }
   return numneigh;
}

int jlk_to_number2(int j, int l, int k, int num_cells)
{
   if (j < 0) {j = num_cells + j;}
   if (l < 0) {l = num_cells + l;}
   if (k < 0) {k = num_cells + k;}

   if (j > num_cells) {j = j - num_cells;}
   if (l > num_cells) {l = l - num_cells;}
   if (k > num_cells) {k = k - num_cells;}

   return round(j%num_cells)*num_cells*num_cells+round(l%num_cells)*num_cells+round(k%num_cells);
}

vector<int> getneigh_fast(const Lattice& lattice, const vector<Particle>& pars, const int nsep, const vector<CVPCLASS> cvpclass)
{
   int n = lattice.allnodes.size(); // the number of nodes 
   int np = pars.size(); // number of all particles
   std::vector<int> numneigh;
   numneigh.resize(n, 0);

   for (int i = 0; i < np; ++i)
   {
      if (cvpclass[i] == FLU)
      {
            // find closest node
            int closest_node_j = round(pars[i].pos[0] / lattice.lcellx);
            int closest_node_l = round(pars[i].pos[1] / lattice.lcelly);
            int closest_node_k = round(pars[i].pos[2] / lattice.lcellz);
            int delta = round(nsep / lattice.lcellx) * 4;
            
            for (int j = closest_node_j - delta ; j <= closest_node_j + delta; ++j)
            {
               for (int l = closest_node_l - delta ; l <= closest_node_l + delta; ++l)
               {
                  for (int k = closest_node_k - delta ; k <= closest_node_k + delta; ++k)
                  {
                     int node = jlk_to_number2(j, l, k, lattice.num_cells_1d);       
                     double dist = sepnp(lattice.allnodes[node], pars[i], lattice.lboxx, 
                     lattice.lboxy, lattice.lboxz, lattice.zperiodic);
                     if (dist <= nsep * nsep)
                     {
                       ++numneigh[node];
                     }
                  }
               }
            }
      }
      
   }
   return numneigh;
}