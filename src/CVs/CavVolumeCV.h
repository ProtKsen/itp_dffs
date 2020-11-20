/**
 * This file is part of
 * SSAGES - Software Suite for Advanced General Ensemble Simulations
 *
 * Copyright 2016 Hythem Sidky <hsidky@nd.edu>
 *
 * SSAGES is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SSAGES is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SSAGES.  If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once 

#include "CollectiveVariable.h"
#include "Validator/ObjectRequirement.h"
#include "Drivers/DriverException.h"
#include "Snapshot.h"
#include "schema.h"
#include <numeric>
#include "mpi.h"
#include <ctime>
#include <vector>
#include <map>
#include "typedefs.h"
#include <math.h>
#include <complex>
#include "constants.h"
#include <boost/graph/connected_components.hpp>
#include "boost/multi_array.hpp"
#include <omp.h>
#include <queue>
#include <stdlib.h>
#include <algorithm>

namespace SSAGES
{
    struct Particle
    {
        double pos[3];
        char symbol; // for outputting e.g. jmol
    };


    class Box
    {
    public:
        Box(){ }
        Box(double lx, double ly, double lz, double ns = 1.5, bool pz = false):
        lboxx(lx), lboxy(ly), lboxz(lz), nsep(ns), nsepsq(ns * ns), periodicz(pz){}

        // functions for separation between two particles
        inline void sep(const Particle& p1, const Particle& p2, double* s) const;
        inline double sepsq(const Particle& p1, const Particle& p2) const;

        // these ones also return whether the particles are neighbours
        inline bool isneigh(double* s, double&r2) const;

        // testing whether positions are valid (in box)
        inline bool posvalid(double* pos) const;
        inline bool getvalidifnot(double* pos) const;

        inline void setdims(double lx, double ly, double lz);

    private:
        double lboxx;
        double lboxy;
        double lboxz;        
        double nsep;
        double nsepsq;
        bool periodicz;
    };


    struct ParticleSystem
    {
        // constructor: this will set the correct values for all of the
         // variables defined below.
        ParticleSystem(std::string pfile);
        Box simbox;
        unsigned int nlinks; // num links for particle to be liquid 
        double nsep;         // neighbour separation, if rij < nsep particles are neighbours
    };


    bool space(char c)    // for strings
    {
        return isspace(c);
    };

    bool not_space(char c)  // for strings
    {
        return !isspace(c);
    };


    std::vector<std::string> split(const std::string& str)  // splitting a string
    {
        typedef std::string::const_iterator iter;
        std::vector<std::string> ret;

        iter i = str.begin();
        while (i != str.end()) {
            // ignore leading blanks
            i = find_if(i, str.end(), not_space);

            // find end of next word
            iter j = find_if(i, str.end(), space);

            // copy the characters in [i,j)
            if (i != str.end()) {
                ret.push_back(std::string(i,j));
            }
        i = j;
        }
        return ret;
    };


    std::map<std::string, std::string> readparams(const std::string fname)
    {
        std::map<std::string, std::string> params;
        std::ifstream infile;
        infile.open(fname.c_str());
        std::string sline;
        std::vector<std::string> spline;
     
        // warning: there is no real error checking here   
        while (infile) 
        { 
            std::getline(infile,sline);

            // comments # must be first char on line
            if (sline[0] != '#') {
                spline = split(sline);
                if (spline.size() == 2) {
                    params[spline[0]] = spline[1];
                }
            }
        }
        infile.close();
        return params;
    };


    inline void Box::sep(const Particle& p1, const Particle& p2, double* s) const
    // distance between particles
    {
        double sepx,sepy,sepz;
        sepx = std::abs(p1.pos[0] - p2.pos[0]);  // use abs?
        sepy = std::abs(p1.pos[1] - p2.pos[1]);
        sepz = std::abs(p1.pos[2] - p2.pos[2]);
     
        if (sepx > 0.5 * lboxx) {
            sepx = lboxx - sepx;
        }
        if (sepy > 0.5 * lboxy) {
            sepy = lboxy - sepy;
        }
        if (sepz > 0.5 * lboxz) {
            sepz = lboxz - sepz;
        }    
        
        s[0] = sepx;
        s[1] = sepy;
        s[2] = sepz;
        return;
    };


    inline bool Box::isneigh(double *s, double& rsq) const   
    // are particles neighbours
    {
        if ((std::abs(s[0]) < nsep) && (std::abs(s[1]) < nsep)
            && (std::abs(s[2]) < nsep)) {
            rsq = s[0] * s[0] + s[1] * s[1] + s[2] * s[2];
            if (rsq < nsep * nsep) {     // < or <=  ???
                return true;
            }
        }
        return false;
    };


    std::vector<int> getneigh(const std::vector<Particle>& particles, const std::vector<Particle>& allparticles, const Box& simbox,
             std::vector<int>& numneigh)
    {
        std::vector<Particle>::size_type npar = particles.size();
        std::vector<Particle>::size_type npar_tot = allparticles.size();
        std::vector<int> numneigh_local(npar, 0);
        double PI=3.14159265358979323846;

        double r2;
        double sep[3];
        std::vector<Particle>::size_type i,j;
   
        for (int i = 0; i < npar; ++i) {
            for (int j = 0; j < npar_tot; ++j) {
                simbox.sep(particles[i], allparticles[j], sep);

                if (simbox.isneigh(sep, r2)) {
                    // particles i and j are neighbours
                    ++numneigh[i];
                }
            }     
            numneigh_local[i] = numneigh[i] - 1; // i and i partcile are not neighb
        } 
        return numneigh_local;
    };


    int number(int j, int l, int k, int num_cells)
    // from j,l,k to n of cell
    {
        return (j % num_cells) * num_cells * num_cells + (l % num_cells) * num_cells + (k % num_cells);
    };


    void number_to_jlk (int numb, int &j, int &l, int &k, int num_cells)
    // from n of cell to j,l,k
    {
        k = numb % num_cells;   
	    int div = numb / num_cells;  
	    l = div % num_cells;  
	    div = div / num_cells;      
	    j = div % num_cells;
    };


    void check_for_bound(int &j, int num_cells)
    {
	    if (j == num_cells) {
		    j = 0;
	    }		
	    if (j == -1) {
		    j = num_cells - 1;
	    }	
    };


    //! Collective variable on the cavity volume of a box. 
    /*!
     * Collective variable on the cavity volume of a box. 
     * 
     * \ingroup CVs
     */
    class CavVolumeCV : public CollectiveVariable
    {
    public:
        //! Constructor
        CavVolumeCV()
        {}

        //! Initialize the CV.
        void Initialize(const Snapshot& /* snapshot */) override
        {
        }

        //! Evaluate the CV. 
        /*!
         * \param snapshot Current simulation snapshot.
         */
        void Evaluate(const Snapshot& snapshot) override
        {
            int comm_size ;
            MPI_Comm_size(snapshot.GetCommunicator(), &comm_size);
            int comm_rank;
            MPI_Comm_rank(snapshot.GetCommunicator(), &comm_rank);
            int world_rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); 

            int ncell_VAP = 0; // the number of vapour cells
	        int npart_VAP=0; // the number of vapour particles

            auto timestep = snapshot.GetIteration();
            if (timestep % 1 == 0)  // 1 - check lambda on every step
            {
                // Fill empty gradient. 
			    auto n = snapshot.GetNumAtoms();
			    std::fill(grad_.begin(), grad_.end(), Vector3{0,0,0});
			    grad_.resize(n, Vector3{0,0,0});

                int Ntot = 0; // the whole number of atoms in Comm
                MPI_Allreduce(&n, &Ntot, 1, MPI_INT, MPI_SUM, snapshot.GetCommunicator());
    
                // read atoms and box 
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

                int recvcounts[size];
	            MPI_Allgather(&n, 1, MPI_INT, &recvcounts, 1, MPI_INT, snapshot.GetCommunicator()); 
	            int displs[size];
	            displs[0]=0;
	            for (int i = 1; i < size; ++i) 
                {
		            displs[i] = displs[i-1] + recvcounts[i-1];
	            }

                MPI_Allgatherv (&xpositions, n, MPI_DOUBLE, &allxpositions, recvcounts, displs, MPI_DOUBLE, snapshot.GetCommunicator()); 
                MPI_Allgatherv (&ypositions, n, MPI_DOUBLE, &allypositions, recvcounts, displs, MPI_DOUBLE, snapshot.GetCommunicator());
                MPI_Allgatherv (&zpositions, n, MPI_DOUBLE, &allzpositions, recvcounts, displs, MPI_DOUBLE, snapshot.GetCommunicator());

                std::vector<Particle> pars;
                pars.resize(n);
                Particle par;
                for (int i = 0; i < n; ++i) 
                {
		            par.pos[0] = positions[i][0];
		            par.pos[1] = positions[i][1];
		            par.pos[2] = positions[i][2];
	                pars[i]=par;
	            }
	 
	            std::vector<Particle> allpars;
	            allpars.resize(Ntot);
                for (int i = 0; i < Ntot; ++i) 
                {
		            par.pos[0] = allxpositions[i];
		            par.pos[1] = allypositions[i];
		            par.pos[2] = allzpositions[i];
	                allpars[i] = par;
	            }

                

                std::string pfile="params.out.start.txt";

                val_ = snapshot.GetVolume();
            }
            if(snapshot.GetCommunicator().rank() == 0)
                boxgrad_ = val_*Matrix3::Identity();
        }

        //! \copydoc CollectiveVariable::BuildCV()
        static CavVolumeCV* Build(const Json::Value& json, const std::string& path)
		{
			Json::ObjectRequirement validator;
			Json::Value schema;
			Json::CharReaderBuilder rbuilder;
			Json::CharReader* reader = rbuilder.newCharReader();

			reader->parse(JsonSchema::CavVolumeCV.c_str(),
			              JsonSchema::CavVolumeCV.c_str() + JsonSchema::CavVolumeCV.size(),
			              &schema, nullptr);
			validator.Parse(schema, path);

			// Validate inputs.
			validator.Validate(json, path);
			if(validator.HasErrors())
				throw BuildException(validator.GetErrors());

			return new CavVolumeCV();
		}
    };
}
