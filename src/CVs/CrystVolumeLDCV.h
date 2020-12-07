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


namespace SSAGES
{
    //! Collective variable on the volume of a box. 
    /*!
     * Collective variable on the volume of a box. 
     * 
     * \ingroup CVs
     */
    class CrystVolumeLDCV : public CollectiveVariable
    {
    public:
        //! Constructor
        CrystVolumeLDCV()
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
                ParticleSystem psystem(pfile, snapshot);

                // compute the qlm data
                // warning: at the moment the number of links, and the threshold
                // value for a link is the same for both l=4 and l=6
                // (psystem.linval and psystem.nlinks respectively)
                int lval = 6;
                QData q6data(psystem, snapshot, lval);
                lval = 4;
                QData q4data(psystem, snapshot, lval);
    
                // from q6data and q4 data, classify each particle as bcc, hcp
                // etc.  using Lechner Dellago approach.
                
                std::vector<LDCLASS> ldclass = classifyparticlesld(psystem, q4data, q6data);

                // indices into particle vector (psystem.allpars) of those
                // particles in the ten-Wolde Frenkel largest cluster and those
                // in the Lechner Dellago cluster.
                std::vector<int> ldcnums = largestclusterld(psystem, ldclass);

                int size_LD = csizeld(ldcnums);
                val_ = size_LD;

                MPI_Barrier(snapshot.GetCommunicator());
     
                // write results
               if (snapshot.GetCommunicator().rank() == 0)
               {
		           auto dumpfilename = snapshot.GetIteration(); 
                   std::system("mkdir -p NLDs");
                   std::string FileNLD="NLDs/NLD_"+std::to_string(snapshot.GetWalkerID())+".txt";
                   std::ofstream fout_NLD(FileNLD,std::ios_base::out | std::ios_base::app);  
                   fout_NLD << dumpfilename <<" "<< val_ << std::endl;
                   fout_NLD.close(); 
		       }           

	            MPI_Barrier(snapshot.GetCommunicator());
            }

            if(snapshot.GetCommunicator().rank() == 0)
                boxgrad_ = val_*Matrix3::Identity();
        }

        //! \copydoc CollectiveVariable::BuildCV()
        static CrystVolumeLDCV* Build(const Json::Value& json, const std::string& path)
		{
			Json::ObjectRequirement validator;
			Json::Value schema;
			Json::CharReaderBuilder rbuilder;
			Json::CharReader* reader = rbuilder.newCharReader();

			reader->parse(JsonSchema::CrystVolumeLDCV.c_str(),
			              JsonSchema::CrystVolumeLDCV.c_str() + JsonSchema::CrystVolumeLDCV.size(),
			              &schema, nullptr);
			validator.Parse(schema, path);

			// Validate inputs.
			validator.Validate(json, path);
			if(validator.HasErrors())
				throw BuildException(validator.GetErrors());

			return new CrystVolumeLDCV();
		}
    };
}

