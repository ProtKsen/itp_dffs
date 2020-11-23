#include <iostream>
#include <vector>
#include "particlesystem.h"
#include "particle.h"
#include "box.h"
#include "conncomponents.h"
#include "utility.h"
#include "diagonalize.h"
#include "gtensor.h"

using std::vector;

// Constructor for gyration tensor (see gtensor.h).

GTensor::GTensor(const ParticleSystem& psystem, const vector<int>& cnums)
{
   // The gyration tensor itself
   gtensor.resize(boost::extents[3][3]);
   gtensor = getgytensor(psystem, cnums);

   // Diagonalise complete tensor
   double fullg[] = {gtensor[0][0], gtensor[0][1], gtensor[0][2],
                     gtensor[1][0], gtensor[1][1], gtensor[1][2],
                     gtensor[2][0], gtensor[2][1], gtensor[2][2]};
   double fullres[9];
   
   // After this call, fulleig stores the 3 eigenvalues
   diagonalize(fullg, 3, fullres, fulleig);

   // Diagonalise 2*2 tensor (x-y)
   double g[] = {gtensor[0][0], gtensor[0][1],
                 gtensor[1][0], gtensor[1][1]};
   double res[4];
   
   // After this call, topeig stores the 2 eigenvalues of top part
   diagonalize(g, 2, res, topeig);
}    

// Take positions of particles in largest cluster.  Return particles
// in the largest cluster, but without periodic BCS.  The trick here
// is to replicate the system in x and y directions, then to find the
// largest cluster in this large system (note the system is assumed
// not to be periodic in z).

vector<Particle> posnoperiodic(const vector<Particle>& cpars, const Box& simbox)
{
   const vector<Particle>::size_type ncl = cpars.size();
   const vector<Particle> repcpars = replicate(cpars, simbox);

   // hack to ignore periodic bcs
   Box bigbox = simbox;
   bigbox.setdims(simbox.lboxx * 20, simbox.lboxy * 20, simbox.lboxz);

   graph xgraph = getxgraph(repcpars, range(0, repcpars.size()), bigbox);
   vector<int> cluspars = largestcomponent(xgraph);

   // create vector of cluster positions
   vector<Particle> ret;
   ret.resize(ncl);
   
   for (vector<Particle>::size_type i = 0; i != ncl; ++i) {
      ret[i] = repcpars[cluspars[i]];
   }

   return ret;
}

// Replicate particles in x and y directions.

vector<Particle> replicate(const vector<Particle>& pars, const Box& simbox)
{
   const vector<Particle>::size_type npar = pars.size();
   vector<Particle> newpars;
   newpars.resize(9 * npar);

   for (vector<Particle>::size_type i = 0; i != npar; ++i) {
      // set the position in the replica as the original position 
      for (int j = 0; j != 9; ++j) {
         newpars[i + j*npar] = pars[i];
      }
      // now add the required offset
      // top
      newpars[i + npar].pos[1] += simbox.lboxy;
      // top right
      newpars[i + 2 * npar].pos[0] += simbox.lboxx;
      newpars[i + 2 * npar].pos[1] += simbox.lboxy;          
      // right
      newpars[i + 3 * npar].pos[0] += simbox.lboxx;
      // bottom right
      newpars[i + 4 * npar].pos[0] += simbox.lboxx;
      newpars[i + 4 * npar].pos[1] -= simbox.lboxy;                   
      // bottom
      newpars[i + 5 * npar].pos[1] -= simbox.lboxy;                   
      // bottom left
      newpars[i + 6 * npar].pos[0] -= simbox.lboxx;
      newpars[i + 6 * npar].pos[1] -= simbox.lboxy;                   
      // left
      newpars[i + 7 * npar].pos[0] -= simbox.lboxx;
      // top left
      newpars[i + 8 * npar].pos[0] -= simbox.lboxx;
      newpars[i + 8 * npar].pos[1] += simbox.lboxy;
   }

   return newpars;
}

// Center of mass of particles.

vector<double> cofmass(const vector<Particle>& particles)
{
   vector<double> cm(3, 0.0);
     
   vector<Particle>::size_type i;
   vector<Particle>::size_type npar = particles.size();
   for (i = 0; i != npar; ++i) {
      cm[0] += particles[i].pos[0];
      cm[1] += particles[i].pos[1];
      cm[2] += particles[i].pos[2];
   }
   cm[0] /= npar;
   cm[1] /= npar;
   cm[2] /= npar;
     
   return cm;
}

// Return radius of gyration tensor (not diagonalised).  The particles
// should have periodic bounary conditions removed.

tensor gytensor(const vector<Particle>& particles)
{
   vector<double> cmass = cofmass(particles);
   tensor gyt(boost::extents[3][3]);
   std::fill(gyt.origin(), gyt.origin() + gyt.size(), 0.0);   

   double rcm[3];
   vector<Particle>::size_type i;
   vector<Particle>::size_type npar = particles.size();
   
   for (i = 0; i != npar; ++i) {
      // distance of particle from center of mass
      rcm[0] = particles[i].pos[0] - cmass[0];
      rcm[1] = particles[i].pos[1] - cmass[1];
      rcm[2] = particles[i].pos[2] - cmass[2];
      // contribution to gyration tensor (top half)
      for (int j = 0; j != 3; ++j) {
         for (int k = j; k != 3; ++k) {
            gyt[j][k] += rcm[j] * rcm[k];
         }
      }
   }
   // compute remaining entries (since symmetric) and normalise
   gyt[1][0] = gyt[0][1];
   gyt[2][0] = gyt[0][2];
   gyt[2][1] = gyt[1][2];
   for (int j = 0; j != 3; ++j) {
      for (int k = 0; k != 3; ++k) {
         gyt[j][k] /= npar;
      }
   }
   return gyt;
}

// Return radius of gyration tensor.  Note that periodic boundary
// conditions need to be removed before computing this. Here this is
// achieved via the function replicate (defined above). note cnums
// gives indices into psystem.allpars of particles in largest cluster.

tensor getgytensor(const ParticleSystem& psystem, const vector<int>& cnums)
{
   // build up vector of particles which are largest cluster only
   vector<Particle>::size_type ncl = cnums.size();
   vector<Particle> clusterpars(ncl);
   
   for (vector<Particle>::size_type i = 0; i != ncl; ++i) {
      clusterpars[i].pos[0] = psystem.allpars[cnums[i]].pos[0];
      clusterpars[i].pos[1] = psystem.allpars[cnums[i]].pos[1];
      clusterpars[i].pos[2] = psystem.allpars[cnums[i]].pos[2];
      // setting the symbol is not really necessary, do it anyway
      clusterpars[i].symbol = 'S';
   }

   // take away periodic bcs
   vector<Particle> cparsnop = posnoperiodic(clusterpars, psystem.simbox);

   return gytensor(cparsnop);
}
