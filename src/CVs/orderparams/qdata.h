#ifndef QDATA_H
#define QDATA_H

#include <vector>
#include "particlesystem.h"
#include "typedefs.h"
#include "constants.h"

// QData is a struct to store data used for the Steindhardt bond order
// parameters. The struct exists since we don't want to recompute
// spherical harmonics many times.

struct QData
{
public:
   QData(const ParticleSystem& psystem, int lval);

   // store the l value, usually either 4 or 6
   int lval;
     
   // we store the number of neighbours for each particle as well as
   // a vector of neighbour particles for each particle in this
   // structure. This is a slightly annoying place to store this
   // information.  It is here since the neighbours are computed when
   // calculating the qlm matrix (see functions qlms).
   vector<int> numneigh;
   vector<vector<int> > lneigh;

   // the complete qlm matrix
   array2d qlm;
     
   // ql, \bar{ql}, wl, \bar{wl} for each particle i
   std::vector<double> ql;
   std::vector<double> qlbar;
   std::vector<double> wl;
   std::vector<double> wlbar;

   // the number of crystalline 'links' that each particle has
   // warning: this only really makes sense for the case l = 6
   vector<int> numlinks;
};

std::vector<TFCLASS> classifyparticlestf(const ParticleSystem&, const QData&);
std::vector<LDCLASS> classifyparticlesld(const ParticleSystem&, const QData&, const QData&);
std::vector<int> largestclusterld(const ParticleSystem&, const std::vector<LDCLASS>&);
std::vector<int> largestclustertf(const ParticleSystem&, const std::vector<TFCLASS>&);

#endif
