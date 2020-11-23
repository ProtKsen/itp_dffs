#ifndef GYRATION_H
#define GYRATION_H

#include <vector>
#include "typedefs.h"
#include "particlesystem.h"

struct GTensor
{
public:
   GTensor(const ParticleSystem& psystem, const std::vector<int>& cnums);

   // We store in this object:
   // 1) The complete gyration tensor in the normal x,y,z coordinate
   //    system.
   // 2) The eigenvalues of the complete gyration tensor
   // 3) The eigenvalues of the top 2x2 portion of the gyration
   //    tensor.
	  
   tensor gtensor;
   double fulleig[3];
   double topeig[2];
};

tensor getgytensor(const ParticleSystem&, const std::vector<int>&);

#endif
