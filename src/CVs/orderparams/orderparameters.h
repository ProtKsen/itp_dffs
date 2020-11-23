#ifndef ORDERPARAMETERS_H
#define ORDERPARAMETERS_H

#include <iostream>
#include "constants.h"
#include "qdata.h"
#include "gtensor.h"

int csizeld(const std::vector<int>&);
int csizetf(const std::vector<int>&);
double qavgroup(const QData&, const std::vector<int>&);
double eiglarge(const GTensor&);
double eigmid(const GTensor&);
double eigsmall(const GTensor&);
double rogsquared(const GTensor&);
double element33(const GTensor&);
double eiglargetop(const GTensor&);
double eigsmalltop(const GTensor&);
int numconnections(const QData&, const std::vector<int>&);

// template for finding fraction of a particular type of particle in
// a list of particles. Used for n_fcc etc.

template <class etype>
double parfrac(const std::vector<etype>& pclass, const std::vector<int>& cnums, const etype plabel)
{
   int num = 0;
   for (std::vector<int>::size_type i = 0; i != cnums.size(); ++i) {
      if (pclass[cnums[i]] == plabel) {
         ++num;
      }
   }

   return static_cast<double>(num) / cnums.size();
}

// template for finding list of particles of a particular type
// which have at least one neighbour in the list.  Used for N_s, which
// is the number of liquid particles with at least 1 neighbour in
// the cluster.  (It could conversely be used for finding indices of
// particles in the cluster that have at least one liquid neighbour.

template <class etype>
std::vector<int> nparatleastone(const std::vector<etype>& pclass, const std::vector<int>& cnums,
                                const etype plabel, const std::vector<std::vector<int> >& lneigh)
{
   std::vector<int> indexes;
   for (typename std::vector<etype>::size_type i = 0; i != pclass.size(); ++i) {
      if (pclass[i] == plabel) { // e.g. we are a liquid particle
         // go through all neighbours
         for (vector<int>::size_type j = 0; j != lneigh[i].size(); ++j) {
            // is the neighbour in the list (usually cluster)?
            if (find(cnums.begin(), cnums.end(), lneigh[i][j]) != cnums.end()) {
               indexes.push_back(i);
               break;
            }
         }
      }
   }

   return indexes;
}

#endif
