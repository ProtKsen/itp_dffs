#ifndef GETNEIGH_H
#define GETNEIGH_H

#include <vector>
#include <complex>
#include "particle.h"
#include "lattice.h"
#include "node.h"
#include "box.h"
#include "typedefs.h"

//std::vector<int> xtalpars(const std::vector<int>&, const int);
std::vector<int> getneigh(const Lattice&, const std::vector<Particle>&, const int);
std::vector<int> getneigh_fast(const Lattice&, const std::vector<Particle>&, const int, const vector<CVPCLASS> cvpclass);

#endif
