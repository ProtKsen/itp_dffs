#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <map>
#include <vector>
#include "particle.h"

using std::vector;
using std::ifstream;
using std::ofstream;
using std::setprecision;
using std::string;
using std::cout;
using std::endl;
using std::map;
     
bool space(char c)
{
   return isspace(c);
}

bool not_space(char c)
{
   return !isspace(c);
}

// Split a string.

vector<string> split(const string& str)
{
   typedef string::const_iterator iter;
   vector<string> ret;

   iter i = str.begin();
   while (i != str.end()) {

      // ignore leading blanks
      i = find_if(i, str.end(), not_space);

      // find end of next word
      iter j = find_if(i, str.end(), space);

      // copy the characters in [i,j)
      if (i != str.end()) {
         ret.push_back(string(i,j));
      }
      i = j;
   }
   return ret;
}

// Write vector of particles to output file in the normal XYZ format.

void writexyz(vector<Particle> pars, const string fname, bool writesymbols = true)
{
   ofstream outfile(fname.c_str());
   int npar = pars.size();
   outfile << npar << endl << endl;
   outfile.precision(8); // default precision is 8 decimal places
   outfile.setf(std::iostream::fixed);
   outfile.setf(std::iostream::showpoint);

   for (vector<Particle>::iterator i = pars.begin(); i != pars.end(); i++) {
      if (writesymbols) {
         outfile << i->symbol << " ";
      }
      outfile << i->pos[0] << " " << i->pos[1] << " " << i->pos[2] << endl;
   }

   return;
}

// Read input parameters.  See 'params.out'

map<string, string> readparams(const string fname)
{
   map<string, string> params;
   ifstream infile(fname.c_str());
   string sline;
   vector<string> spline;
     
   // warning: there is no real error checking here   
   while (infile) { 

      getline(infile,sline);

      // comments # must be first char on line
      if (sline[0] != '#') {
         spline = split(sline);
         if (spline.size() == 2) {
            params[spline[0]] = spline[1];
         }
      }
   }

   return params;
}

// Read vector of particles from file in the normal XYZ format.

vector<Particle> readxyz(const string fname, bool symbols = true,
                         bool gettypes = true)
{
   vector<Particle> allpars;
   int npar = 0;
   ifstream infile(fname.c_str());
   // map to define conversion between character symbol
   // and particle type (an integer)
   map<char,int> partypes;
   partypes['O'] = 1;
   partypes['S'] = 0;
   partypes['N'] = 0;

   // check that file exists and can be read from
   if (!infile) {
      cout << "Warning: " << fname
           << " does not exist or cannot be read." << endl;
      // return empty vector of particles
      return allpars;
   }

   // read number of particles (must be top line of XYZ file)
   try {
      infile >> npar;
   } catch (...) {
      cout << "Warning: " << fname
           << " does not appear to be an XYZ file."
           << " Top line must be integer number of particles." << endl
           << " No particles found." << endl;
      // return empty vector of particles
      return allpars;
   }
   allpars.resize(npar);

   string sline;
   vector<string> spline;
   unsigned int ncols = 3 + symbols; // number of columns in XYZ file
   Particle par;
   int nread = 0;
   int lread = 1;
   int i = 0;

   // invariant : we have successfully read nread particles
   while (infile) {

      getline(infile,sline);
      ++lread;
      spline = split(sline);

      if (spline.empty()) { // we read a blank line
         continue;
      }

      // check that we read correct number of columns
      if (spline.size() != ncols) {
         cout << "Warning: corrupt XYZ file, line " << lread << " expected "
              << ncols << " columns." << endl;
         continue;
      }

      i = 0;
      try {
         if (symbols) { // note the symbol must only be a single character
            par.symbol = spline[i++][0];
            if (gettypes) 
               par.type = partypes[par.symbol];
         }
         par.pos[0] = atof(spline[i++].c_str());
         par.pos[1] = atof(spline[i++].c_str());
         par.pos[2] = atof(spline[i].c_str());

         allpars[nread] = par;
         ++nread ;
      }
      catch (...) {
         cout << "Warning: Error reading XYZ file, line "
              << lread << "." << endl;
      }
   }

   // check we read the correct number of particles
   if (nread != npar) 
      cout << "Warning: Did not read correct number of particles from "
           << fname << " (" << nread << " of " << npar << " read)"
           << endl;

   return allpars;
}
