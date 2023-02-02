#include <map>
#include <string>
#include <cstdlib>
#include <iostream>
#include "particlesystem.h"
#include "readwrite.h"
#include "box.h"
#include "compile.h"
#include "Snapshot.h"
#include "mpi.h"

using std::map;
using std::string;
using std::cout;
using std::endl;

// Constructor for ParticleSystem object.

ParticleSystem::ParticleSystem(string pfile, const SSAGES::Snapshot& snapshot)
{
   int comm_size ;
   MPI_Comm_size(snapshot.GetCommunicator(), &comm_size);
   int comm_rank;
   MPI_Comm_rank(snapshot.GetCommunicator(), &comm_rank);
   int world_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);  

   auto n = snapshot.GetNumAtoms();
   int Ntot = 0; // the whole number of atoms in Comm
   MPI_Allreduce(&n, &Ntot, 1, MPI_INT, MPI_SUM, snapshot.GetCommunicator());

   // particle positions
   auto positions = snapshot.GetPositions(); 
   auto velocities = snapshot.GetVelocities();
   auto HMatrix = snapshot.GetHMatrix();
   auto atomID = snapshot.GetAtomIDs();
   double xpositions[n];
   double ypositions[n];
   double zpositions[n];
   double allxpositions[Ntot];
   double allypositions[Ntot];
   double allzpositions[Ntot];

   for (int i = 0; i < n; ++i) 
   {
	   xpositions[i] = positions[i][0];
	   ypositions[i] = positions[i][1];
	   zpositions[i] = positions[i][2];
	};

   int recvcounts[comm_size];
	MPI_Allgather(&n, 1, MPI_INT, &recvcounts, 1, MPI_INT, snapshot.GetCommunicator()); 
	int displs[comm_size];
	displs[0]=0;
	for (int i = 1; i < comm_size; ++i) 
   {
		displs[i] = displs[i-1] + recvcounts[i-1];
	}

   MPI_Allgatherv (&xpositions, n, MPI_DOUBLE, &allxpositions, recvcounts, displs, MPI_DOUBLE, snapshot.GetCommunicator()); 
   MPI_Allgatherv (&ypositions, n, MPI_DOUBLE, &allypositions, recvcounts, displs, MPI_DOUBLE, snapshot.GetCommunicator());
   MPI_Allgatherv (&zpositions, n, MPI_DOUBLE, &allzpositions, recvcounts, displs, MPI_DOUBLE, snapshot.GetCommunicator()); 

   pars.resize(n);
   Particle par;
   for (int i = 0; i < n; ++i) 
   {
	   par.pos[0] = positions[i][0];
		par.pos[1] = positions[i][1];
		par.pos[2] = positions[i][2];
	   pars[i]=par;
   }	
	allpars.resize(Ntot);
   for (int i = 0; i < Ntot; ++i) 
   {
		par.pos[0] = allxpositions[i];
		par.pos[1] = allypositions[i];
		par.pos[2] = allzpositions[i];
	   allpars[i]=par;
	}

   std::map<std::string, std::string> params = readparams(pfile);
   double lboxx = HMatrix(0,0);
   double lboxy = HMatrix(1,1);
   double lboxz = HMatrix(2,2);

   nsep = atof(params["stillsep"].c_str());
   map<string, bool> bmap;
   bmap["True"] = true;
   bmap["False"] = false;	 
   bool zperiodic = bmap[params["zperiodic"]];
   simbox = Box(lboxx,lboxy,lboxz,nsep,zperiodic);

   // number of surface particles
   nsurf = atoi(params["nparsurf"].c_str());

   // warning: at the moment these are used for both l=4 and l=6
   linval = atof(params["q6link"].c_str());
   nlinks = atoi(params["q6numlinks"].c_str());

   // write results and tests
   bool writeresults = bmap[params["writeresults"]];
   bool writestruct = bmap[params["writestruct"]];

   if (LOGGING) {
      cout << LOGMSG << "read " << allpars.size() << " particles" << endl
           << LOGMSG << "values for particle system: " << endl
           << LOGMSG << "lboxx "     << lboxx << endl
           << LOGMSG << "lboxy "     << lboxy << endl
           << LOGMSG << "lboxz "     << lboxz << endl
           << LOGMSG << "stillsep "  << nsep << endl
           << LOGMSG << "zperiodic " << zperiodic << endl
           << LOGMSG << "nparsurf " << nsurf << endl
           << LOGMSG << "q6link " << linval << endl
           << LOGMSG << "q6numlinks " << nlinks << endl;
   }
}
