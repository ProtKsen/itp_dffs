#ifndef READWRITE_H
#define READWRITE_H

#include "particle.h"
#include <string>

std::map<std::string, std::string> readparams(const std::string fname);
std::vector<Particle> readxyz(const std::string fname, bool symbols = true, bool gettypes = true);
void writexyz(std::vector<Particle> pars, const std::string fname, bool writesymbols = true);

#endif
