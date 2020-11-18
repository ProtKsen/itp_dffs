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

namespace SSAGES
{
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
    }

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
            // Fill empty gradient. 
			auto n = snapshot.GetNumAtoms();
			std::fill(grad_.begin(), grad_.end(), Vector3{0,0,0});
			grad_.resize(n, Vector3{0,0,0});

            val_ = snapshot.GetVolume();
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
