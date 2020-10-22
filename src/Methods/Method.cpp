/**
 * This file is part of
 * SSAGES - Software Suite for Advanced General Ensemble Simulations
 *
 * Copyright 2017 Hythem Sidky <hsidky@nd.edu>
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
#include "Method.h"
#include "ABF.h"
#include "ANN.h"
#include "CFF.h"
#include "Umbrella.h"
#include "BasisFunc.h"
#include "ForwardFlux.h"
#include "Meta.h"
#include "StringMethod.h"
#include "schema.h"
#include "json/json.h"
#include "Drivers/DriverException.h"
#include "Validator/ObjectRequirement.h"
#include "CVs/CVManager.h"
#include <stdexcept>

using namespace Json;

namespace SSAGES
{
	Method* Method::BuildMethod(const Value& json,
	                            const MPI_Comm& world,
	                            const MPI_Comm& comm,
	                            const std::string& path)
	{
		ObjectRequirement validator;
		Value schema;
		CharReaderBuilder rbuilder;
		CharReader* reader = rbuilder.newCharReader();

		reader->parse(JsonSchema::Method.c_str(),
		              JsonSchema::Method.c_str() + JsonSchema::Method.size(),
		              &schema, nullptr);
		validator.Parse(schema, path);
		validator.Validate(json, path);
		if(validator.HasErrors())
			throw BuildException(validator.GetErrors());

		Method* method = nullptr;

		if(json["type"] == "ABF")
			method = ABF::Build(json, world, comm, path);
		else if(json["type"] == "ANN")
			method = ANN::Build(json, world, comm, path);
		else if(json["type"] == "BFSMethod")
			method = BFS::Build(json, world, comm, path);
		else if(json["type"] == "CFF")
			method = CFF::Build(json, world, comm, path);
		else if(json["type"] == "ForwardFlux")
			method = ForwardFlux::Build(json, world, comm, path);
		else if(json["type"] == "Metadynamics")
			method = Meta::Build(json, world, comm, path);
		else if(json["type"] == "Umbrella")
			method = Umbrella::Build(json, world, comm, path);
		else if(json["type"] == "String")
			method = StringMethod::Build(json, world, comm, path);
		else
			throw std::invalid_argument(path + ": Unknown method type specified.");

		// Load cv mask.
		std::vector<unsigned int> cvmask;
		for(auto& cv : json["cvs"])
		{
			cvmask.push_back(CVManager::LookupCV(cv, path));
		}

		method->SetCVMask(cvmask);
		return method;
	}
}
