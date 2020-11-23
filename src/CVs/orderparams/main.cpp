#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include "particlesystem.h"
#include "orderparameters.h"
#include "qdata.h"
#include "constants.h"
#include "utility.h"
#include "gtensor.h"

using std::cout;
using std::endl;
using std::string;

// Compute the values of all order parameters for a given
// configuration. Currently there are 42 different order parameters
// output, but a lot of these are 'duplicated' because there are two
// different definitions of the largest crystalline cluster in the
// system (one uses the so-called Ten-Wolde Frenkel approach, the
// other uses the Lecher-Dellago approach ). See README and the
// example for further information and an example parameter file.

int main(int argc, char* argv[])
{
   if (argc != 2) {
      cout << "Syntax: " << argv[0] << " paramfile" << endl;
      return 1;
   }

   // get name of input file
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

   // from q6 data only, classify each particle as either
   // crystalline or liquid, using TenWolde Frenkel approach
   vector<TFCLASS> tfclass = classifyparticlestf(psystem, q6data);

   // indices into particle vector (psystem.allpars) of those
   // particles in the ten-Wolde Frenkel largest cluster and those
   // in the Lechner Dellago cluster.
   vector<int> tfcnums = largestclustertf(psystem, tfclass);
   vector<int> ldcnums = largestclusterld(psystem, ldclass);

   // indices of liquid like particles that have at least one
   // neighbour in the cluster, for both ld and tf
 //  vector<int> ldliquid1nums = nparatleastone(ldclass, ldcnums, LIQUID, q6data.lneigh);
 //  vector<int> tfliquid1nums = nparatleastone(tfclass, tfcnums, LIQ, q6data.lneigh);

   // indexes of all particles (minus surface particles)
 //  vector<int> pindices = range(psystem.nsurf, psystem.allpars.size());

   // radius of gyration tensor for both clusters
  // GTensor tfgtensor(psystem, tfcnums);
 //  GTensor ldgtensor(psystem, ldcnums);

   // compute each order parameter in turn and print to stdout.
   // See orderparams.cpp for these functions.

   //////////////////////////////////////////////////////////////////
   // The following order parameters are associated in some way with
   // properties of the largest cluster.  There are two approaches
   // to determining this cluster, which I call Lecher Dellage (LD)
   // and ten-Wolde Frenkel (TF), and thus two different clusters.
   // All of the OPs are computed for both clusters.
   /////////////////////////////////////////////////////////////////
	  
   // Size of cluster by LD method
   cout << "N_ld " << csizeld(ldcnums) << endl;

   // Size of cluster by TF method
   cout << "N_tf " << csizetf(tfcnums) << endl;
	  
   // fraction of bcc pars in LD cluster
  // cout << "n_bccLD " << parfrac(ldclass, ldcnums, BCC) << endl;
	  
   // fraction of bcc pars in TF cluster	  
   //cout << "n_bccTF " << parfrac(ldclass, tfcnums, BCC) << endl;
	  
   // fraction of fcc pars in LD cluster
   //cout << "n_fccLD " << parfrac(ldclass, ldcnums, FCC) << endl;
	  
   // fraction of fcc pars in TF cluster	  
   //cout << "n_fccTF " << parfrac(ldclass, tfcnums, FCC) << endl;
	  
   // fraction of hcp pars in LD cluster
   //cout << "n_hcpLD " << parfrac(ldclass, ldcnums, HCP) << endl;
	  
   // fraction of hcp pars in TF cluster	  
   //cout << "n_hcpTF " << parfrac(ldclass, tfcnums, HCP) << endl;
	  
   // fraction of icos pars in LD cluster
  // cout << "n_icosLD " << parfrac(ldclass, ldcnums, ICOS) << endl;
	  
   // // fraction of icos pars in TF cluster	  
  // cout << "n_icosTF " << parfrac(ldclass, tfcnums, ICOS) << endl;
	 
   // average Q6 of LD cluster
  // cout << "Q6clusLD " << qavgroup(q6data, ldcnums) << endl;

   // average Q6 of TF cluster
  // cout << "Q6clusTF " << qavgroup(q6data, tfcnums) << endl;
	  	 
   // average Q4 of LD cluster
  // cout << "Q4clusLD " << qavgroup(q4data, ldcnums) << endl;

   // average Q4 of TF cluster
  // cout << "Q4clusTF " << qavgroup(q4data, tfcnums) << endl;

   // number of liquid like particles with at least one neighbour in
   // LD cluster.  Note that we could pass either q6data.lneigh or
   // q4data.lneigh, since these are identical
  // cout << "N_sLD " << ldliquid1nums.size() << endl;

   // same as above but for TF cluster
  // cout << "N_sTF " << tfliquid1nums.size() << endl;

   // total number of connections for all liquid-like particles with
   // at least one neighbour in cluster for LD cluster
  // cout << "N_lLD " << numconnections(q6data, ldliquid1nums) << endl;

   // // same as above but for TF cluster
  // cout << "N_lTF " << numconnections(q6data, tfliquid1nums) << endl;

   // average q6 of liquid-like particles with at least one
   // neighbour in cluster for LD cluster
  // cout << "Q6N_sLD " << qavgroup(q6data, ldliquid1nums) << endl;	  	  

   // same as above but for TF cluster
 //  cout << "Q6N_sTF " << qavgroup(q6data, tfliquid1nums) << endl;

   // average q4 of liquid-like particles with at least one
   // neighbour in cluster for LD cluster
 //  cout << "Q4N_sLD " << qavgroup(q4data, ldliquid1nums) << endl;

   // same as above but for LD cluster
 //  cout << "Q4N_sTF " << qavgroup(q4data, tfliquid1nums) << endl;

   // smallest eigenvalue of gyration tensor for LD cluster
 //  cout << "Rbar_g,1LD " << eigsmall(ldgtensor) << endl;	  

   // smallest eigenvalue of gyration tensor for TF cluster
 //  cout << "Rbar_g,1TF " << eigsmall(tfgtensor) << endl;

   // middle eigenvalue of gyration tensor for LD cluster
 //  cout << "Rbar_g,2LD " << eigmid(ldgtensor) << endl;

   // middle eigenvalue of gyration tensor for TF cluster
  // cout << "Rbar_g,2TF " << eigmid(tfgtensor) << endl;

   // largest eigenvalue of gyration tensor for LD cluster
 //  cout << "Rbar_g,3LD " << eiglarge(ldgtensor) << endl;
	  
   // largest eigenvalue of gyration tensor for TF cluster
 //  cout << "Rbar_g,3TF " << eiglarge(tfgtensor) << endl;

   // square of 'radius of gyration' for LD cluster
 //  cout << "Rbar_gLD " << rogsquared(ldgtensor) << endl;

   // square of 'radius of gyration' for TF cluster
 //  cout << "Rbar_gTF " << rogsquared(tfgtensor) << endl;	  
	  
   // (3,3) element of non-diagonalized gyration tensor for LD
   // cluster
 //  cout << "R_g,zLD " << element33(ldgtensor) << endl;
	  
   // (3,3) element of non-diagonalized gyration tensor for TF
   // cluster
 //  cout << "R_g,zTF " << element33(tfgtensor) << endl;

   // smallest eigenvalue of top-diagonalised gyration tensor for LD
   // cluster
 //  cout << "R_g,1LD " << eigsmalltop(ldgtensor) << endl;	  

   // smallest eigenvalue of top-diagonalised gyration tensor for TF
   // cluster
 //  cout << "R_g,1TF " << eigsmalltop(tfgtensor) << endl;	  

   // largest eigenvalue of top-diagonalised gyration tensor for LD
   // cluster
 //  cout << "R_g,2LD " << eiglargetop(ldgtensor) << endl;

   // largest eigenvalue of top-diagonalised gyration tensor for TF
   // cluster
 //  cout << "R_g,2TF " << eiglargetop(tfgtensor) << endl;	  

   //////////////////////////////////////////////////////////////////
   // These order parameter are 'global' i.e. for the entire system
   // (that is, no mention of a cluster of any kind!).  Note that we
   // exclude surface particles from the calculations.
   //////////////////////////////////////////////////////////////////
	  
   // fraction of bcc particles in entire system
 //  cout << "s_bcc " << parfrac(ldclass, pindices, BCC) << endl;
	  
   // fraction of fcc particles in entire system
 //  cout << "s_fcc " << parfrac(ldclass, pindices, FCC) << endl;
	  
   // fraction of hcp particles in entire system
 //  cout << "s_hcp " << parfrac(ldclass, pindices, HCP) << endl;

   // fraction of icosahedral particles in entire system
 //  cout << "s_icos " << parfrac(ldclass, pindices, ICOS) << endl;

   // Average q6 of all particles in system
//   cout << "Q6 " << qavgroup(q6data, pindices) << endl;

   // average q4 of all particles in system
 //  cout << "Q4 " << qavgroup(q4data, pindices) << endl;	  
}
