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



#include "core/PairwiseEstimator.hpp"
#include "models/GTRModel.hpp"
#include "models/HKY85Model.hpp"
#include "models/AminoacidSubstitutionModel.hpp"
#include "models/NegativeBinomialGapModel.hpp"

namespace EBC
{


PairwiseEstimator::BFGS::BFGS(PairwiseEstimator* enclosing, Definitions::OptimizationType ot) : optimizationType(ot)
{
	parent = enclosing;
	paramsCount = parent->modelParams->optParamCount();
	this->initParams.set_size(paramsCount);
	this->lowerBounds.set_size(paramsCount);
	this->upperBounds.set_size(paramsCount);

	parent->modelParams->toDlibVector(initParams,lowerBounds,upperBounds);

	//cerr << "DLIB optimizer init with " << paramsCount << " parameters" << endl;
}

PairwiseEstimator::BFGS::~BFGS()
{
}

double PairwiseEstimator::BFGS::objectiveFunction(const column_vector& bfgsParameters)
{
	this->parent->modelParams->fromDlibVector(bfgsParameters);
	return parent->runIteration();
}


const column_vector PairwiseEstimator::BFGS::objectiveFunctionDerivative(const column_vector& bfgsParameters)
{
	column_vector results(this->paramsCount);
	return results;
}


void PairwiseEstimator::BFGS::optimize()
{
	using std::placeholders::_1;
	std::function<double(const column_vector&)> f_objective= std::bind( &PairwiseEstimator::BFGS::objectiveFunction, this, _1 );
	double likelihood;

	switch(optimizationType)
	{
		case Definitions::OptimizationType::BFGS:
		{
			likelihood = dlib::find_min_box_constrained(dlib::bfgs_search_strategy(),
					dlib::objective_delta_stop_strategy(1e-8),
					f_objective,
					derivative(f_objective),
					initParams,
					lowerBounds,
					upperBounds);
			break;
		}
		case Definitions::OptimizationType::BOBYQA:
		{
			likelihood = dlib::find_min_bobyqa(f_objective, initParams, parent->modelParams->optParamCount()+4,
					lowerBounds,upperBounds, 0.05, 1e-7, 20000 );
			break;
		}
	}
	this->parent->modelParams->fromDlibVector(initParams);
	//parent->modelParams->outputParameters();
	//cout  << likelihood << "\n";

}


PairwiseEstimator::PairwiseEstimator(Definitions::AlgorithmType at, Sequences* inputSeqs, Definitions::ModelType model ,std::vector<double> indel_params,
		std::vector<double> subst_params, Definitions::OptimizationType ot, bool banding, unsigned int bandPercentage,
		unsigned int rateCategories, double alpha, bool estimateAlpha, double userTime) : inputSequences(inputSeqs), gammaRateCategories(rateCategories),
		pairCount(inputSequences->getPairCount()), hmms(pairCount)
{
	maths = new Maths();
	dict = inputSequences->getDictionary();

	//Helper models
	if (model == Definitions::ModelType::GTR)
	{
		substModel = new GTRModel(dict, maths,gammaRateCategories);
	}
	else if (model == Definitions::ModelType::HKY85)
	{
		substModel = new HKY85Model(dict, maths,gammaRateCategories);
	}
	else if (model == Definitions::ModelType::LG)
	{
			substModel = new AminoacidSubstitutionModel(dict, maths,gammaRateCategories,Definitions::aaLgModel);
	}

	indelModel = new NegativeBinomialGapModel();

	estimateSubstitutionParams = subst_params.size() == 0;
	estimateIndelParams = indel_params.size() == 0;
	this->estimateAlpha = estimateAlpha;

	DEBUG("Pairwise model estimator starting");
	DEBUG("Estimate substitution parameters set to : " << estimateSubstitutionParams << " Estimate indel parameters set to : " << estimateIndelParams);
	DEBUG("Estimate alpha set to : " << estimateAlpha << " , rate categories " << gammaRateCategories << " , alpha : " << alpha);

	modelParams = new OptimizedModelParameters(substModel, indelModel,inputSequences->getSequenceCount(), pairCount, estimateSubstitutionParams,
			estimateIndelParams, estimateAlpha, userTime < 0, maths);

	if(!estimateIndelParams)
		modelParams->setUserIndelParams(indel_params);
	if(!estimateSubstitutionParams)
		modelParams->setUserSubstParams(subst_params);
	modelParams->setAlpha(alpha);



	if (userTime > 0)
	{
		//cerr << "User time " << userTime << endl;
		vector<double> times(pairCount);
		for (auto it = times.begin(); it < times.end(); it++)
		{
			*it = 2.0*(userTime/pairCount);
		}
		modelParams->setUserDivergenceParams(times);
	}


	bandFactor = bandPercentage;
	bandingEnabled = banding;

	EvolutionaryPairHMM* hmm;

	for(unsigned int i =0; i<pairCount; i++)
	{
		std::pair<unsigned int, unsigned int> idxs = inputSequences->getPairOfSequenceIndices(i);

		if (at == Definitions::AlgorithmType::Viterbi)
		{
			hmm = hmms[i] = new ViterbiPairHMM(inputSequences->getSequencesAt(idxs.first), inputSequences->getSequencesAt(idxs.second),
				substModel, indelModel);
		}
		else if (at == Definitions::AlgorithmType::Forward)
		{
			hmm = hmms[i] = new ForwardPairHMM(inputSequences->getSequencesAt(idxs.first), inputSequences->getSequencesAt(idxs.second),
					substModel, indelModel, Definitions::DpMatrixType::Full);
		}
		else
		{
			throw HmmException("Wrong algorithm type - use either Forward or viterbi\n");
		}
	}

	substModel->setObservedFrequencies(inputSequences->getElementFrequencies());
	if (estimateSubstitutionParams == false)
	{
		//set parameters and calculate the model
		substModel->setAlpha(modelParams->getAlpha());
		substModel->setParameters(modelParams->getSubstParameters());
		substModel->calculateModel();
	}

	if (estimateIndelParams == false)
	{
			//set parameters and calculate the model
		indelModel->setParameters(modelParams->getIndelParameters());
	}

	bfgs = new BFGS(this,ot);
	bfgs->optimize();
}

PairwiseEstimator::~PairwiseEstimator()
{
	delete bfgs;
	delete modelParams;
	//delete Y;
	//delete X;
	//delete M;
	//delete substModel;
	//delete indelModel;
    delete maths;
}

double PairwiseEstimator::runIteration()
{
	double result = 0;
	EvolutionaryPairHMM* hmm;

	//this->modelParams->outputParameters();
	//cerr << "iteration " << endl;

	if (estimateSubstitutionParams == true)
	{
			//set parameters and calculate the model
		substModel->setAlpha(modelParams->getAlpha());
		substModel->setParameters(modelParams->getSubstParameters());
		substModel->calculateModel();
	}

	if (estimateIndelParams == true)
	{
		//set parameters and calculate the model
		indelModel->setParameters(modelParams->getIndelParameters());
	}



	for(unsigned int i =0; i<pairCount; i++)
	{
		hmm = hmms[i];
		hmm->setDivergenceTimeAndCalculateModels(modelParams->getDivergenceTime(i));
		indelModel->setTime(modelParams->getDivergenceTime(i));
		indelModel->calculate();
		result += hmm->runAlgorithm();
		//modelParams->outputParameters();
	}
	cerr << " lnl " << result << endl;
	return result;
}

void PairwiseEstimator::outputDistanceMatrix(stringstream& ss)
{
	unsigned int count, pairCount;
	count = this->inputSequences->getSequenceCount();
	pairCount = this->inputSequences->getPairCount();

	ss << "\t" << this->inputSequences->getSequenceCount() << endl;

	for (unsigned int i = 0; i< count; i++)
	{
		ss << "S" << i << " ";
		for(unsigned int j=0; j< count; j++)
		{
			ss << this->modelParams->getDistanceBetween(i,j) << " ";
		}
		ss << endl;
	}
}

} /* namespace EBC */
