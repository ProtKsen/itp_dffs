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
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include "orderparams/particlesystem.h"
#include "orderparams/orderparameters.h"
#include "orderparams/qdata.h"
#include "orderparams/constants.h"
#include "orderparams/utility.h"
#include "orderparams/gtensor.h"
#include "orderparams/readwrite.h"  // readparams
#include "orderparams/box.h" // Box
#include <numeric>
#include "mpi.h"
#include <map>
#include <math.h>
#include <complex>
#include "constants.h"
#include <boost/graph/connected_components.hpp>
#include "boost/multi_array.hpp"
#include <omp.h>
#include <queue>
#include <stdlib.h>
#include <algorithm>

using std::cout;
using std::endl;
using std::string;

namespace SSAGES
{
    //! Collective variable on the volume of a box. 
    /*!
     * Collective variable on the volume of a box. 
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
            int comm_size ;
            MPI_Comm_size(snapshot.GetCommunicator(), &comm_size);
            int comm_rank;
            MPI_Comm_rank(snapshot.GetCommunicator(), &comm_rank);
            int world_rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); 
            auto timestep = snapshot.GetIteration();
            if (timestep % 1 == 0)  // 1 - check lambda on every step
            {
                auto n = snapshot.GetNumAtoms();
			    std::fill(grad_.begin(), grad_.end(), Vector3{0,0,0});
			    grad_.resize(n, Vector3{0,0,0});

                int Ntot = 0; // the whole number of atoms in Comm
                MPI_Allreduce(&n, &Ntot, 1, MPI_INT, MPI_SUM, snapshot.GetCommunicator());

                std::string pfile="params.out.start.txt";
                ParticleSystem psystem(std::string pfile, const SSAGES::Snapshot& snapshot);

                // compute the qlm data
                // warning: at the moment the number of links, and the threshold
                // value for a link is the same for both l=4 and l=6
                // (psystem.linval and psystem.nlinks respectively)
                int lval = 6;
                QData q6data(ParticleSystem psystem, const SSAGES::Snapshot& snapshot, int lval);
                lval = 4;
                QData q4data(ParticleSystem psystem, const SSAGES::Snapshot& snapshot, int lval);
    
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

                val_ = tfcnums;
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

