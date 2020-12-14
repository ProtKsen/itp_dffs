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

            // for testing
            /*std::string FileTest="Test_"+std::to_string(snapshot.GetWalkerID())+
                                        "_" + std::to_string(snapshot.GetCommunicator().rank()) +".txt";
            std::string FileTest_pos_local="pos_local_"+ std::to_string(timestep)+
                                        "_" + std::to_string(snapshot.GetWalkerID())+
                                        "_" + std::to_string(snapshot.GetCommunicator().rank()) + ".txt";
             std::string FileTest_pos_global="pos_global_"+ std::to_string(timestep)+
                                        "_" + std::to_string(snapshot.GetWalkerID())+
                                        "_" + std::to_string(snapshot.GetCommunicator().rank()) + ".txt";

            std::ofstream fout_test(FileTest, std::ios_base::out | std::ios_base::app);
            std::ofstream fout_pos_local(FileTest_pos_local, std::ios_base::out | std::ios_base::app);
            std::ofstream fout_pos_global(FileTest_pos_global, std::ios_base::out | std::ios_base::app);  
            fout_pos_local.precision(10);
            fout_pos_global.precision(10);
            */

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
                
                /*for  (int i = 0; i < Ntot; ++i) 
                {
                    fout_test << timestep << " numlinks " << psystem.allpars[i].pos[0] << " " <<
                    psystem.allpars[i].pos[1] << " " << psystem.allpars[i].pos[2] << " " <<
                    q6data.qlm[i][0] << " " << q6data.numlinks[i] << std::endl;
                }               
                fout_test.close();
                

                auto HMatrix = snapshot.GetHMatrix();
                fout_pos_global << Ntot << std::endl;
                fout_pos_global << HMatrix(0,0) << " " << HMatrix(1,1) << " " << HMatrix(2,2) << " " << std::endl;
                for  (int i = 0; i < Ntot; ++i) 
                {
                    fout_pos_global << "N " << psystem.allpars[i].pos[0] << " " <<
                    psystem.allpars[i].pos[1] << " " << psystem.allpars[i].pos[2] << std::endl;
                }               
                fout_pos_global.close();

                fout_pos_local << Ntot << std::endl;
                fout_pos_local << " " << std::endl;
                for  (int i = 0; i < n; ++i) 
                {
                    fout_pos_local << i << " " << psystem.allpars[i].pos[0] << " " <<
                    psystem.allpars[i].pos[1] << " " << psystem.allpars[i].pos[2] << std::endl;
                }               
                fout_pos_local.close();
                */
    
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

