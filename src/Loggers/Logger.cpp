/**
 * This file is part of
 * SSAGES - Software Suite for Advanced General Ensemble Simulations
 *
 * Copyright 2017 Hythem Sidky <hsidky@nd.edu>
 *                Ben Sikora <bsikora906@gmail.com>
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

#include <stdexcept>
#include "Logger.h"
#include "Snapshot.h"
#include "CVs/CVManager.h"
#include "Validator/ObjectRequirement.h"
#include "Drivers/DriverException.h"
#include "schema.h"

using namespace Json;

namespace SSAGES
{
	void Logger::PreSimulation(Snapshot* /* snapshot */, const CVManager& cvmanager)
	{
		if(IsMasterRank(comm_))
		{
			if(append_)
				log_.open(filename_.c_str(), std::ofstream::out | std::ofstream::app);
			else
			{
				// Write out header.
				log_.open(filename_.c_str(), std::ofstream::out);
				log_ << "#"; 
				log_ << "Iteration "; 

				auto cvs = cvmanager.GetCVs(cvmask_);
				for(size_t i = 0; i < cvs.size() - 1; ++i)
					log_ << "cv_" + std::to_string(i) << " ";
				log_ << "cv_" + std::to_string(cvs.size() - 1) << std::endl;
			}
		}
	}

	void Logger::PostIntegration(Snapshot* snapshot, const CVManager& cvmanager)
	{
		// Get necessary info. 
		auto cvs = cvmanager.GetCVs(cvmask_);
		if(IsMasterRank(comm_))
		{
			log_.precision(8);
			log_ << snapshot->GetIteration() << " ";

			// Print out CV values.
			for(size_t i = 0; i < cvs.size() - 1; ++i)
				log_ << cvs[i]->GetValue() << " ";
			log_ << cvs.back()->GetValue() << std::endl;
		}
	}

	void Logger::PostSimulation(Snapshot* /* snapshot */, const class CVManager& /* cvmanager */)
	{
	}

	Logger* Logger::Build(const Value& json,
	                      const MPI_Comm& world,
	                      const MPI_Comm& comm,
	                      const std::string& path)
	{
		ObjectRequirement validator;
		Value schema;
		CharReaderBuilder rbuilder;
		CharReader* reader = rbuilder.newCharReader();

		reader->parse(JsonSchema::Logger.c_str(),
		              JsonSchema::Logger.c_str() + JsonSchema::Logger.size(),
		              &schema, nullptr);
		validator.Parse(schema, path);

		// Validate inputs.
		validator.Validate(json, path);
		if(validator.HasErrors())
			throw BuildException(validator.GetErrors());

		auto freq = json.get("frequency", 1).asInt();
		std::string name = "cvlog.dat";

		bool ismulti = GetNumWalkers(world, comm) > 1;
		unsigned int wid = GetWalkerID(world, comm);

		if(json["output_file"].isArray())
		{
			if(json["output_file"].size() != GetNumWalkers(world, comm))
			{
				throw BuildException({path + ": Multi-walker simulations require a separate output file for each walker."});
			}
			name = json["output_file"][wid].asString();
		}
		else if(ismulti)
			throw std::invalid_argument(path + ": Multi-walker simulations require a separate output file for each.");
		else
			name = json["output_file"].asString();

		auto* l = new Logger(freq, name, world, comm);
		l->SetAppend(json.get("append", false).asBool());

		// Load CV mask.
		std::vector<unsigned int> cvmask;
		for(auto& cv : json["cvs"])
		{
			cvmask.push_back(CVManager::LookupCV(cv, path));
		}
		l->SetCVMask(cvmask);

		return l;
	}

}
