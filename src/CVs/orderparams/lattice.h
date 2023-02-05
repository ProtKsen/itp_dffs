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
   
   // the number of nodes in one dimension
   unsigned int num_cells_1d;
   
   double r_sep;
   double lboxx;
   double lboxy;
   double lboxz;
   double lcellx;
   double lcelly;
   double lcellz; 
   bool zperiodic;
};

#endif
