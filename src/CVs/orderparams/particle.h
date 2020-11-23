#ifndef PARTICLE_H
#define PARTICLE_H

// Simple particle class for particle simulations e.g. MC/MD.

struct Particle
{
   double pos[3];
   double vel[3];
   double mass;
   int type; // this allows for different interaction potentials etc.
   char symbol; // for outputting e.g. jmol
};

#endif
