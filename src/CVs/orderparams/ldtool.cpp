#include <iostream>
#include <string>
#include <vector>
#include "particlesystem.h"
#include "qdata.h"
#include "constants.h"

using std::vector;
using std::cout;
using std::endl;
using std::string;

// Tool for outputting classification of each particle using the LD
// method.  This reads in a file with the required parameters, and
// outputs an XYZ co-ordinate file, which can be viewed with molecular
// visualisation software, e.g. JMOL.  See README and the example for
// further information and an example parameter file.

int main(int argc, char* argv[])
{
   if (argc != 2) {
      cout << "Syntax: " << argv[0] << " paramfile" << endl;
      return 1;
   }

   // get name of input and output files
   string pfile = argv[1];
	  
   // create the particle system from name of input file
   // input file must contain the following fields:
   // filename   - name of xyz file to read positions from
   // lboxx      - x dimension of simulation box
   // lboxy      - y dimension     ""
   // lboxz      - z dimension     ""
   // stillsep   - neighbour separation in units of sigma
   // zperiodic  - box periodic or not (either "True" or "False")
   // nparsurf   - number of surface particles
   // q6link     - threshold for Sij to be considered a link
   // q6numlinks - number of links a particle needs to be xtal
   ParticleSystem psystem(pfile);

   // compute the qlm data
   // warning: at the moment the number of links, and the threshold
   // value for a link is the same for both l=4 and l=6
   // (psystem.linval and psystem.nlinks respectively)
   QData q6data(psystem, 6);
   QData q4data(psystem, 4);
	  
   // from q6data and q4 data, classify each particle as bcc, hcp
   // etc.  using Lechner Dellago approach.
   vector<LDCLASS> ldclass = classifyparticlesld(psystem, q4data, q6data);

   // strings for jmol representation
   // see constants.h form LDCLASS enum
   vector<string> jstr;
   jstr.push_back("S"); // FCC (yellow in jmol)
   jstr.push_back("P"); // HCP (orange in jmol)
   jstr.push_back("F"); // BCC (green in jmol)
   jstr.push_back("N"); // LIQUID (blue in jmol)
   jstr.push_back("B"); // ICOS (pink in jmol)
   jstr.push_back("O"); // SURFACE (red in jmol)

   // output number of particles
   cout << ldclass.size() << endl << endl;
   for (int i = 0; i != psystem.allpars.size(); ++i) {
      cout   << jstr[ldclass[i]]<< " " << psystem.allpars[i].pos[0]
           << " " << psystem.allpars[i].pos[1] << " "
           << psystem.allpars[i].pos[2] << endl;
   }

   return 0;
}
