#include <map>
#include <string>
#include <cstdlib>
#include <iostream>
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

   // parameters

   std::map<std::string, std::string> params = readparams(pfile);
   r_sep = atoi(params["latt_sep"].c_str());
   num_cells_1d = atoi(params["num_cells_1d"].c_str());
   map<string, bool> bmap;
   int numcells = pow(num_cells_1d, 3);

   auto HMatrix = snapshot.GetHMatrix();
   lboxx = HMatrix(0,0);
   lboxy = HMatrix(1,1);
   lboxz = HMatrix(2,2);

   //cell's size
   lcellx = lboxx/num_cells_1d;
   lcelly = lboxy/num_cells_1d;
   lcellz = lboxz/num_cells_1d; 

   double minx = 0.0;
   double miny = 0.0;
   double minz = 0.0;
   
   allnodes.resize(numcells);
   Node nod;
   for (int j=0; j< num_cells_1d; ++j) 
   {
	   for (int l=0; l< num_cells_1d; ++l)
	   {
		   for (int k=0; k< num_cells_1d; ++k)
		   {
	         nod.pos[0] = minx + (j % num_cells_1d) * lcellx;
		      nod.pos[1] = miny + (l % num_cells_1d) * lcelly;
		      nod.pos[2] = minz + (k % num_cells_1d) * lcellz;
            nod.number = jlk_to_number(j, l, k, num_cells_1d);
	         allnodes[nod.number] = nod;
         }
      }
   }

   int leftj = comm_rank * floor(num_cells_1d / comm_size);
   int rightj = (comm_rank + 1) * floor(num_cells_1d / comm_size);
   if (comm_rank == comm_size - 1) 
   {
      rightj = num_cells_1d;
   }
   nodes.resize((rightj - leftj) * (num_cells_1d) * (num_cells_1d));
   int count = 0;
   for (int j = leftj; j < rightj; ++j) 
   {
	   for (int l=0; l< num_cells_1d; ++l)
	   {
		   for (int k=0; k< num_cells_1d; ++k)
		   {
	         nod.pos[0] = minx + (j % num_cells_1d) * lcellx;
		      nod.pos[1] = miny + (l % num_cells_1d) * lcelly;
		      nod.pos[2] = minz + (k % num_cells_1d) * lcellz;
            nod.number = count;
	         nodes[nod.number] = nod;
            ++count;
         }
      }
   }
}
