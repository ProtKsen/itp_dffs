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
    class CrystVolumeTFCV : public CollectiveVariable
    {
    public:
        //! Constructor
        CrystVolumeTFCV()
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
                
                int lval = 6;
                QData q6data(psystem, snapshot, lval);
    
                // from q6 data only, classify each particle as either
                // crystalline or liquid, using TenWolde Frenkel approach
                std::vector<TFCLASS> tfclass = classifyparticlestf(psystem, q6data);

                // indices into particle vector (psystem.allpars) of those
                // particles in the ten-Wolde Frenkel largest cluster.
                std::vector<int> tfcnums = largestclustertf(psystem, tfclass);

                int size_TF = csizetf(tfcnums);
                val_ = size_TF;

                MPI_Barrier(snapshot.GetCommunicator());
     
                // write results
               if (snapshot.GetCommunicator().rank() == 0)
               {
		           auto dumpfilename = snapshot.GetIteration(); 
                   std::system("mkdir -p NTFs");
                   std::string FileNTF="NTFs/NTF_"+std::to_string(snapshot.GetWalkerID())+".txt";
                   std::ofstream fout_NTF(FileNTF,std::ios_base::out | std::ios_base::app);  
                   fout_NTF << dumpfilename <<" "<< val_ << std::endl;
                   fout_NTF.close(); 
		       }           

	            MPI_Barrier(snapshot.GetCommunicator());
            }

            if(snapshot.GetCommunicator().rank() == 0)
                boxgrad_ = val_*Matrix3::Identity();
        }

        //! \copydoc CollectiveVariable::BuildCV()
        static CrystVolumeTFCV* Build(const Json::Value& json, const std::string& path)
		{
			Json::ObjectRequirement validator;
			Json::Value schema;
			Json::CharReaderBuilder rbuilder;
			Json::CharReader* reader = rbuilder.newCharReader();

			reader->parse(JsonSchema::CrystVolumeTFCV.c_str(),
			              JsonSchema::CrystVolumeTFCV.c_str() + JsonSchema::CrystVolumeTFCV.size(),
			              &schema, nullptr);
			validator.Parse(schema, path);

			// Validate inputs.
			validator.Validate(json, path);
			if(validator.HasErrors())
				throw BuildException(validator.GetErrors());

			return new CrystVolumeTFCV();
		}
    };
}

