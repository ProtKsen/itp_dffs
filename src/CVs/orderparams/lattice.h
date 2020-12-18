#ifndef LATTICE_H
#define LATTICE_H

#include<vector>
#include<string>
#include "box.h"
#include "node.h"
#include "Snapshot.h"

using std::vector;
using std::string;

// lattice object

struct Lattice
{
public:
   // constructor: this will set the correct values for all of the
   // variables defined below.
   Lattice(string pfile, const SSAGES::Snapshot& snapshot);

   // particle positions
   vector<Node> nodes;
   vector<Node> allnodes;
   // simulation box
   Box simbox;
   // num crystalline neighbours for particle to be in crystalline environment
   unsigned int ncrneigh;
   // the number of nodes in one dimension
   unsigned int numcells_1d;
   // neighbour separation, if r < nsep particles are neighbours
   double nsep;
   double lboxx;
   double lboxy;
   double lboxz;
   double lcellx;
   double lcelly;
   double lcellz; 
   bool zperiodic;
};

#endif
