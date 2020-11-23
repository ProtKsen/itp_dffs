#include <vector>
#include "constants.h"
#include "qlmfunctions.h"
#include "qdata.h"
#include "gtensor.h"
#include "orderparameters.h"

using std::vector;

// Size of cluster according to Lechner Dellago (LD) method.

int csizeld(const vector<int>& ldcnums)
{
   return ldcnums.size();
}

// Size of cluster according to Ten-Wolde Frenkel (TF) method.

int csizetf(const vector<int>& tfcnums)
{
   return tfcnums.size();
}

// Average Q value of a group of particles, e.g. those in cluster.

double qavgroup(const QData& qdata, const vector<int>& pnums)
{
   return Qpars(qdata.qlm, pnums, qdata.lval);
}

// Largest eigenvalue of complete gyration tensor.

double eiglarge(const GTensor& gt)
{
   return gt.fulleig[2];
}

// Middle eigenvalue of complete gyration tensor.

double eigmid(const GTensor& gt)
{
   return gt.fulleig[1];
}

// Smallest eigenvalue of complete gyration tensor.

double eigsmall(const GTensor& gt)
{
   return gt.fulleig[0];
}

// squared radius of gyration.  note this is the square of R_g as
// defined by Jungblutt & Dellago.

double rogsquared(const GTensor& gt)
{
   return gt.fulleig[0] + gt.fulleig[1] + gt.fulleig[2];
}

// (3,3) element of non-diagonalise gyration tensor.

double element33(const GTensor& gt)
{
   // indices of course start at zero!
   return gt.gtensor[2][2];
}

// largest eigenvalue of 'top-diagonised' gyration tensor.

double eiglargetop(const GTensor& gt)
{
   return gt.topeig[1];
}

// smallest eigenvalue of 'top-diagonalised' gyration tensor.

double eigsmalltop(const GTensor& gt)
{
   return gt.topeig[0];
}

// total number of crystalline 'connections' for a group of particles
// Note crystalline 'connections' means the same as crystalline
// 'links'.  Note that this only makes sense for q6 (don't use for
// q4).

int numconnections(const QData& q6data, const vector<int>& cnums)
{
   int num = 0;
   for (vector<int>::size_type i = 0; i != cnums.size(); ++i) {
      num += q6data.numlinks[cnums[i]];
   }

   return num;
}
