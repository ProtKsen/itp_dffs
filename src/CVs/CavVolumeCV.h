/*
* calculate cavity volume by density analysis, suitable for any condensed phase
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
#include "orderparams/lattice.h" 
#include "orderparams/neighdata.h"
#include <stdio.h> 
#include <time.h>
#include <ctime>


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
            int comm_size ;
            MPI_Comm_size(snapshot.GetCommunicator(), &comm_size);

            int comm_rank;
            MPI_Comm_rank(snapshot.GetCommunicator(), &comm_rank);

            int world_rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

            auto timestep = snapshot.GetIteration();
            const auto& HMatrix = snapshot.GetHMatrix();

            val_ = 0;
            if ((timestep > 1) && (timestep % 1 == 0))  // 1 - check lambda on every step
            {
                auto n = snapshot.GetNumAtoms(); // the number of atoms for 1 proc
			    std::fill(grad_.begin(), grad_.end(), Vector3{0,0,0});
			    grad_.resize(n, Vector3{0,0,0});

                int Ntot = 0; // the whole number of atoms in Communicator
                MPI_Allreduce(&n, &Ntot, 1, MPI_INT, MPI_SUM, snapshot.GetCommunicator());

                // psystem save particle's coordiantes and some system's parameters
                std::string pfile="order_parameter.txt";
                ParticleSystem psystem(pfile, snapshot);
                
                //std::string FileTest="Test" + std::to_string(snapshot.GetWalkerID()) +"_" + std::to_string(snapshot.GetCommunicator().rank()) + ".txt";
                //std::ofstream fout_test(FileTest, std::ios_base::out | std::ios_base::app);
                //fout_test << "timestep " << snapshot.GetIteration() << std::endl;
                                

                //int start = std::clock();
                // create virtual lattice
                Lattice lattice(pfile, snapshot);
                //int stop = std::clock();
                //int delta = (stop - start);
                //fout_test << "lattice " << std::to_string(delta) << std::endl;

                // set phase of particles
                //start = std::clock();
                std::vector<CVPCLASS> cvpclass_local = classifypars(psystem);
                //stop = std::clock();
                //delta = (stop - start);
                //fout_test << "cvpclass " << std::to_string(delta) << std::endl;

                // calculate the number of liquid neighbours for every node
                //start = std::clock();
                NeighData neighdata(psystem, lattice, snapshot, cvpclass_local);
                //stop = std::clock();
                //delta = (stop - start);
                //fout_test << "neighdata " << std::to_string(delta) << std::endl;

                // set phase of nodes   
                //start = std::clock();             
                std::vector<CVNCLASS> cvnclass = classifynodes(lattice, neighdata);
                //stop = std::clock();
                //delta = (stop - start);
                //fout_test << "classifynodes " << std::to_string(delta) << std::endl;

                // find the biggest cluster
                //start = std::clock();
                std::vector<int> cnums = largestnodescluster(lattice, cvnclass);
                //stop = std::clock();
                //delta = (stop - start);
                //fout_test << "largestnodescluster " << std::to_string(delta) << std::endl;

                // write result information about lattice
                if (psystem.write_struct & snapshot.GetCommunicator().rank() >= 0)
                {
                    std::string FileOut="Structure_" + std::to_string(snapshot.GetWalkerID()) +"_" + std::to_string(snapshot.GetCommunicator().rank()) + ".txt";
                    std::ofstream fout(FileOut, std::ios_base::out | std::ios_base::app); 
                
                    fout << "timestep " << snapshot.GetIteration() << std::endl;
                    fout << "Box " << HMatrix(0,0) << " " << HMatrix(1,1) << " " << HMatrix(2,2) << std::endl;
                    fout << "id x y z num_neighbours phase is_in_clust" << std::endl;
                    for  (int i = 0; i < lattice.allnodes.size(); ++i) 
                    {
                        int in_clust = 0;
                        for (int j = 0; j < cnums.size(); ++j)
                        {                        
                            if (i == cnums[j])
                            {
                                in_clust = 1;
                            }
                        }
                        fout << i << " " << lattice.allnodes[i].pos[0] << " " << 
                                    lattice.allnodes[i].pos[1] << " " << 
                                    lattice.allnodes[i].pos[2] << " " <<
                                    neighdata.numneigh[i] << " " << 
                                    cvnclass[i] << " " << in_clust << std::endl;                    
                    } 
                }            

                // set order parameter
                val_ = cnums.size() * lattice.lcellx * lattice.lcelly * lattice.lcellz;

                // write results of order parameter
               if (psystem.write_results & snapshot.GetCommunicator().rank() >= 0)
               {
		           auto dumpfilename = snapshot.GetIteration(); 
                   std::system("mkdir -p CVs");
                   std::string FileCV="CVs/CV_"+std::to_string(snapshot.GetWalkerID())+ "_" + std::to_string(snapshot.GetCommunicator().rank()) + ".txt";
                   std::ofstream fout_CV(FileCV,std::ios_base::out | std::ios_base::app);  
                   fout_CV << snapshot.GetCommunicator().rank() << " " << dumpfilename << 
                   " "<< val_ << std::endl;
                   fout_CV.close(); 
		       }           


            }

            MPI_Barrier(snapshot.GetCommunicator()); 
            
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

