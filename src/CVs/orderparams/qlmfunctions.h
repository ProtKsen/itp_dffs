#ifndef QLMFUNCTIONS_H
#define QLMFUNCTIONS_H

#include <vector>
#include <complex>
#include "particle.h"
#include "box.h"
#include "typedefs.h"

std::vector<int> xtalpars(const std::vector<int>&, const int);
std::vector<int> getnlinks(const array2d&, const array2d&, const std::vector<int>&,
                           const std::vector<std::vector<int> >&,
                           const int, const int, const double,
                           const int, const int);

array2d qlmtildes(const array2d&, const std::vector<int>&, const int);
array2d qlmbars(const array2d&, const std::vector<std::vector<int> >&, const int, const int);
array2d qlms(const std::vector<Particle>&, const std::vector<Particle>&, const Box&, std::vector<int>&,
             std::vector<std::vector<int> >&, const int, const int);
double Qpars(const array2d&, const std::vector<int>&, const int);

std::vector<double> qls(const array2d&);
std::vector<double> wls(const array2d&);

#endif
