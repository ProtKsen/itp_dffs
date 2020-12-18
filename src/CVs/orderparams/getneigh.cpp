#include "boost/multi_array.hpp"
#include <complex>
#include <iostream>
#include <math.h>
#include "constants.h"
#include "particle.h"
#include "node.h"
#include "box.h"
#include "opfunctions.h"
#include "getneigh.h"

using std::complex;
using std::vector;

inline void sepnp(const Node& p1, const Particle& p2, const int lboxx, 
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

vector<int> getneigh(const Lattice& lattice, const vector<Particle>& allpars, const int nsep)
{
   int n = lattice.nodes.size(); // the number of nodes in proc
   int Ntot = allpars.size(); // number of all particles
   std::vector<int> numneigh;
   numneigh.resize(n);

   double s[3];
   // compute separation between particles (store in s)
   for (int i = 0; i < n; ++i)
   {
      for (int j = 0; j < Ntot; ++j)
      {
         sepnp(lattice.nodes[i], allpars[j], lattice.lboxx, 
               lattice.lboxy, lattice.lboxz, lattice.zperiodic, s);
         if ((s[0] * s[0] + s[1] * s[1] + s[2] * s[2]) <= nsep * nsep)
         {
            ++numneigh[i];
         }
      }
   }
   return numneigh;
}