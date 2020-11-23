#include <iostream>
#include <vector>
#include <complex>
#include "boost/multi_array.hpp"
#include "qdata.h"
#include "box.h"
#include "particle.h"
#include "qlmfunctions.h"
#include "constants.h"
#include "conncomponents.h"
#include "utility.h"
#include "typedefs.h"

using std::vector;
using std::complex;
using std::norm;

// QData stores qlm in its various forms for each particle in the
// system. Usually l=6 or 4 but the code will accept any value.
//
// The first thing to do is to compute the 2l + 1 dimensional vector
// qlm(i) (m = -l,..,0,..,l) for each particle i in the system.

// This is defined as:
// qlm(i) = 1/N_b(i) \sum_{j=1}^{N_b(i)} Ylm(r_ij)
// Here N_b(i) is the number of neighbours belonging to particle i.
// Two particles are considered neighbours if they are within some cut
// off distance r_cut (This needs to be specified).  The sum is over
// all neighbouring particles.  The functions Ylm are the spherical
// harmonics, and r_ij is a vector from particle i to particle j.
//
// We can use the information qlm(i) in a large number of ways:
//
// i)   To compute the Ql value (e.g. Q6) of a group of particles.
//      (usually all the particles in the system)
//      This is usually called a 'global order parameter'
//      This is defined as:
//      Ql = sqrt( 4pi/(2l + 1) * sum_{m=-l}^{l} |<qlm(i)>|^2 )
//      Where <..> denotes an average over all particles, and |..|
//      denotes absolute value.
// 
// ii)  To compute the ql value of a single particle.
//      This is the same as the formula for Ql, except the < > is
//      removed:
//      ql(i) = sqrt( 4pi/(2l + 1) * sum_{m=-l}^{l} |qlm(i)|^2 )
//      See Lechner and Dellago JCP 129 (2008), equation (3).
//
// iii) To compute the Wl value (e.g. W6) of a group of particles
//      See Steinhart, Nelson, Ronchetti PRB 28 784 (1983)
//      equation (1.5) for the formula.
//
// iv)  To compute the wl value of a single particle
//      See Lechner and Dellago JCP 129 (2008), equation (4).
//      The analogy between i) and ii) is the same as between
//      iii) and iv)
//
// v)   We can normalise the qlm(i).
//      I call the normalised versions \tilde{qlm(i)} or qlmtilde/qlmt
//      in the code below.
//      Then, we can take the dot product (qlmtilde(i) dot qlmtilde(j))
//      between two neighbouring particles i and j
//      Sometimes this dot product is referred to as S_ij in the
//      literature.
//      if S_ij > 0.65 , particles i and j are said to form a 'link'
//      The threshold value 0.65 is arbitrary, others use values
//      between 0.5 and 0.8, the threshold can be specified by the
//      parameter 'linval' below.
//      If particle i has at least 6 links (others use 5-8), it is said
//      to be in a crystalline environment. Again the threshold can
//      be specified by using the parameter 'nlin' in the code.
//      So using this criterion, we can go through every particle in
//      the system, and say whether or not it is in a crystalline
//      (xtal) environment.
//      The size of the largest cluster for which all particles are
//      in a crystalline environment is the familiar Frenkel/ ten Wolde
//      order parameter (which I call N_cl in my papers).

// Constructor for QData object.

QData::QData(const ParticleSystem& psystem, const int _lval) : lval(_lval)
{
   // store number of neighbours and neighbour list
   vector<Particle>::size_type npar = psystem.allpars.size();
   numneigh.resize(npar, 0); // num neighbours for each particle
   lneigh.resize(npar); // neighbour particle nums for each particle

   // matrix of qlm values
   qlm.resize(boost::extents[npar][2 * lval + 1]);
   qlm = qlms(psystem.allpars, psystem.simbox, numneigh, lneigh, lval);

   // Lechner dellago eq 6
   array2d qlmb = qlmbars(qlm, lneigh, lval);

   // get qls and wls
   ql = qls(qlm);
   wl = wls(qlm);
   
   // lechner dellago eq 5
   qlbar = qls(qlmb);
   wlbar = wls(qlmb);

   // compute number of crystalline 'links'
   // first get normalised vectors qlm (-l <= m <= l) for computing
   // dot product Sij
   array2d qlmt = qlmtildes(qlm, numneigh, lval);

   // do dot products Sij to get number of links
   numlinks = getnlinks(qlmt, numneigh, lneigh, psystem.nsurf,
                        psystem.nlinks, psystem.linval, lval);
}

// Classify particles as either Liquid-like or crystalline according
// to the Ten-Wolde Frenkel (TF) method.

vector<TFCLASS> classifyparticlestf(const ParticleSystem& psystem, const QData& q6data)
{
   int npar = q6data.ql.size();
   vector<TFCLASS> parclass(npar, LIQ);

   // from nlinks, work out which particles are xtal
   vector<int> xps = xtalpars(q6data.numlinks, psystem.nlinks);
     
   for (vector<LDCLASS>::size_type i = 0; i != psystem.nsurf; ++i) {
      parclass[i] = SURF;
   }

   for (vector<int>::size_type i = 0; i != xps.size(); ++i) {
      parclass[xps[i]] = XTAL;
   }

   return parclass;
}

// Classify particles as FCC, HCP, BCC, LIQUID, ICOSAHEDRAL or SURFACE
// according to the Lechner Dellago (LD) method.

vector<LDCLASS> classifyparticlesld(const ParticleSystem& psystem, const QData& q4data,
                                    const QData& q6data)
{
   unsigned int npar = q6data.ql.size();
   vector<LDCLASS> parclass(npar);

   for (unsigned int i = 0; i != npar; ++i) {
      if (i < psystem.nsurf) {
         parclass[i] = SURFACE;
      }
      else {
         if (q6data.qlbar[i] < 0.32) {
            parclass[i] = LIQUID;
         }
         else { // particle is solid
            if (abs(q6data.wlbar[i]) > 0.05) {
               parclass[i] = ICOS;
            }
            else if (q6data.wlbar[i] > 0.0) {
               parclass[i] = BCC;
            }
            else { // either HCP or FCC
               if (q4data.wlbar[i] > 0.0) {
                  parclass[i] = HCP;
               }
               else {
                  parclass[i] = FCC;
               }
            }
         }
      }
   }

   return parclass;
}

// Largest cluster using LD classifications.

vector<int> largestclusterld(const ParticleSystem& psystem, const vector<LDCLASS>& ldclass)
{
   // get vector with indices that are all crystal particles
   vector<int> xps;
   for (vector<LDCLASS>::size_type i = 0; i != ldclass.size(); ++i) {
      if ((ldclass[i] == FCC) or (ldclass[i] == HCP) or
          (ldclass[i] == BCC) or (ldclass[i] == ICOS)) {
         xps.push_back(i);
      }
   }

   // graph of xtal particles, with each particle a vertex and each
   // link an edge
   graph xgraph = getxgraph(psystem.allpars, xps, psystem.simbox);

   // largest cluster is the largest connected component of graph
   vector<int> cnums = largestcomponent(xgraph);

   // now largest component returns indexes into array xps, we need
   // to reindex so that it contains indices into psystem.allpars
   // (see utility.cpp)
   reindex(cnums, xps);
   return cnums;
}

// Largest cluster using TF classifications.

vector<int> largestclustertf(const ParticleSystem& psystem, const vector<TFCLASS>& tfclass)
{
   // get vector with indices that are all crystal particles
   vector<int> xps;
   for (vector<LDCLASS>::size_type i = 0; i != tfclass.size(); ++i) {
      if (tfclass[i] == XTAL) {
         xps.push_back(i);
      }
   }

   // graph of xtal particles, with each particle a vertex and each
   // link an edge
   graph xgraph = getxgraph(psystem.allpars, xps, psystem.simbox);

   // largest cluster is the largest connected component of graph
   vector<int> cnums = largestcomponent(xgraph);

   // now largest component returns indexes into array xps, we need
   // to reindex so that it contains indices into psystem.allpars
   // (see utility.cpp)
   reindex(cnums, xps);
   return cnums;    
}
