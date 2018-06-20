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


#ifndef CORE_PAIRHMMCALCULATIONWRAPPER_HPP_
#define CORE_PAIRHMMCALCULATIONWRAPPER_HPP_

#include "core/IOptimizable.hpp"
#include "hmm/EvolutionaryPairHMM.hpp"
#include "core/OptimizedModelParameters.hpp"

namespace EBC
{

class PairHmmCalculationWrapper : public IOptimizable
{
private:
	EvolutionaryPairHMM* phmm;
	OptimizedModelParameters* modelParams;
	SubstitutionModelBase* substModel;
	IndelModel* indelModel;

public:
	PairHmmCalculationWrapper();

	virtual ~PairHmmCalculationWrapper();

	double runIteration();

	void setTargetHMM(EvolutionaryPairHMM* hmm);

	inline void setSubstModel(SubstitutionModelBase* sm){
		substModel = sm;
	}
	inline void setIndelModel(IndelModel* im){
		indelModel = im;
	}

	void setModelParameters(OptimizedModelParameters* mp);


};

}
#endif /* CORE_PAIRHMMCALCULATIONWRAPPER_HPP_ */
