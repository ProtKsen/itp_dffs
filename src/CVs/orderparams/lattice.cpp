#include <map>
#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include "lattice.h"
#include "readwrite.h"
#include "box.h"
#include "compile.h"
#include "Snapshot.h"
#include "mpi.h"
#include "math.h"

using std::map;
using std::string;
using std::cout;
using std::endl;

// Constructor for Lattice object.
int jlk_to_number(int j, int l, int k, int num_cells)
{
    return (j%num_cells)*num_cells*num_cells+(l%num_cells)*num_cells+(k%num_cells);
}

Lattice::Lattice(string pfile, const SSAGES::Snapshot& snapshot)
{
   int comm_size;
   MPI_Comm_size(snapshot.GetCommunicator(), &comm_size);
   int comm_rank;
   MPI_Comm_rank(snapshot.GetCommunicator(), &comm_rank);
   int world_size;
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);
   int world_rank;
   MPI_Comm_size(MPI_COMM_WORLD, &world_rank);

   // parameters

   std::map<std::string, std::string> params = readparams(pfile);
   nsep = atof(params["lattsep"].c_str());
   numcells_1d = atoi(params["lattsize"].c_str());
   ncrneigh = atoi(params["lattncrneigh"].c_str());
   map<string, bool> bmap;
   bmap["True"] = true;
   bmap["False"] = false;	 
   zperiodic = bmap[params["zperiodic"]];
   int numcells = pow(numcells_1d, 3);

   auto HMatrix = snapshot.GetHMatrix();
   lboxx = HMatrix(0,0);
   lboxy = HMatrix(1,1);
   lboxz = HMatrix(2,2);

   //cell's size
   lcellx = lboxx/numcells_1d;
   lcelly = lboxy/numcells_1d;
   lcellz = lboxz/numcells_1d; 

   double minx = 0.0;
   double miny = 0.0;
   double minz = 0.0;
   
   allnodes.resize(numcells);
   Node nod;
   for (int j=0; j< numcells_1d; ++j) 
   {
	   for (int l=0; l< numcells_1d; ++l)
	   {
		   for (int k=0; k< numcells_1d; ++k)
		   {
	         nod.pos[0] = minx + (j % numcells_1d) * lcellx;
		      nod.pos[1] = miny + (l % numcells_1d) * lcelly;
		      nod.pos[2] = minz + (k % numcells_1d) * lcellz;
            nod.number = jlk_to_number(j, l, k, numcells_1d);
	         allnodes[nod.number] = nod;
         }
      }
   }

   int leftj = comm_rank * floor(numcells_1d / comm_size);
   int rightj = (comm_rank + 1) * floor(numcells_1d / comm_size);
   if (comm_rank == comm_size - 1) 
   {
      rightj = numcells_1d;
   }
   nodes.resize((rightj - leftj) * (numcells_1d) * (numcells_1d));
   int count = 0;
   for (int j = leftj; j < rightj; ++j) 
   {
	   for (int l=0; l< numcells_1d; ++l)
	   {
		   for (int k=0; k< numcells_1d; ++k)
		   {
	         nod.pos[0] = minx + (j % numcells_1d) * lcellx;
		      nod.pos[1] = miny + (l % numcells_1d) * lcelly;
		      nod.pos[2] = minz + (k % numcells_1d) * lcellz;
            nod.number = count;
	         nodes[nod.number] = nod;
            ++count;
         }
      }
   }
}
