//==============================================================================
// Pair-HMM phylogenetic tree estimator
// 
// Copyright (c) 2015 Marcin Bogusz.
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses>.
//==============================================================================

/*
 * CommandReader.hpp
 *
 *  Created on: Oct 7, 2013
 *      Author: mbogusz
 */

#ifndef COMMANDREADER_H_
#define COMMANDREADER_H_

#include "core/IParser.hpp"
#include "core/HmmException.hpp"
#include "core/Definitions.hpp"
#include <dlib/cmd_line_parser.h>

namespace EBC
{

class CommandReader
{
public:

	dlib::command_line_parser parser;

	CommandReader(int argc, char** argv);
	IParser* getParser() throw (HmmException&);
/*
	inline bool isViterbi()
	{
		return parser.option("V");
	}

	inline bool isForward()
	{
		return parser.option("F");
	}

	inline bool isMLE()
	{
		return parser.option("M");
	}

	//2 below are optional and new
	inline bool isFdist()
	{
		return parser.option("X");
	}
	inline bool isVdist()
	{
		return parser.option("Y");
	}

	inline bool isFixedAlignment()
	{
		return parser.option("fa");
	}

	vector<double> getIndelParams();

	vector<double> getSubstParams();

	bool isOutputViterbiAlignment()
	{
		return parser.option("ov");
	}

	Definitions::AlgorithmType getAlgorithmType()
		{
			if (parser.option("V"))
			{
				return Definitions::AlgorithmType::Viterbi;
			}
			if (parser.option("F"))
			{
				return Definitions::AlgorithmType::Forward;
			}
			//default;
			return Definitions::AlgorithmType::MLE;
		}

*/
	Definitions::ModelType getModelType()
	{

		if (parser.option("m0"))
		{
			return Definitions::ModelType::M0;
		}
		//default;
		return Definitions::ModelType::M0;
	}

	string getInputFileName()
	{
		return parser.option("in").argument();
	}

	Definitions::OptimizationType getOptimizationType()
	{

			return Definitions::OptimizationType::BFGS;
	}


	FileLogger::logType getLoggingLevel()
	{
		if (parser.option("lDD"))
			return FileLogger::L_DMP;
		if (parser.option("lE"))
			return FileLogger::L_ERR;
		if (parser.option("lW"))
			return FileLogger::L_WARN;
		if (parser.option("lI"))
			return FileLogger::L_INF;
		if (parser.option("lD"))
			return FileLogger::L_DBG;
		//info by default!
		return
			FileLogger::L_INF;
	}


	Definitions::SequenceType getSequenceType()
	{

			return Definitions::SequenceType::Codon;
	}


};

} /* namespace EBC */
#endif /* COMMANDREADER_H_ */
