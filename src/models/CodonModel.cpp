/*
 * NucleotideSubstitutionMode.cpp
 *
 *  Created on: Jan 13, 2014
 *      Author: root
 */

#include "models/CodonModel.hpp"

namespace EBC
{

CodonModel::CodonModel(Dictionary* dict, Maths* alg, unsigned int alpha) :
	SubstitutionModelBase(dict,alg,alpha,Definitions::CodonM0ParamCount)
{
	DUMP("Codon Model: construction");
	this->parameterHiBounds[0] = Definitions::kappaHiBound;
	this->parameterLoBounds[0] = Definitions::standardLowBound;

	this->parameterHiBounds[1] = Definitions::omegaHiBound;
	this->parameterLoBounds[1] = Definitions::standardLowBound;

	this->parameters = new double[Definitions::CodonM0ParamCount];
	this->buildInitialQmatrix();
}

void CodonModel::setParameters(const vector<double>& par)
{
	for (unsigned int i = 0; i< paramsNumber; i++)
	{
		this->parameters[i] = par[i];
	}
}

void CodonModel::buildInitialQmatrix()
{
	DUMP("Codon Model: build initial qMatrix");
	bool transition;
	bool synonymous;
	unsigned int nodiff;
	unsigned int i,j;

	for(i = 0; i < dictionary->getAlphabetSize(); i++)
		for (j=0; j < dictionary->getAlphabetSize(); j++){
			if (i == j)
				continue;
			//check if stop codon
			if (cdict->getAminoacidId(i) == Definitions::stopCodonId
					|| cdict->getAminoacidId(j) == Definitions::stopCodonId)
				this->qMatrix[i * this->matrixSize + j] = 0;
			else{
				nodiff = cdict->getNumberOfDifferentPositions(i,j,synonymous, transition);
				if(nodiff > 1)
					this->qMatrix[i * this->matrixSize + j] = 0;
				else{
					if (synonymous && !transition)
						this->piCodons.emplace_back(i,j);
					else if (synonymous && transition)
						this->kCodons.emplace_back(i,j);
					else if (!synonymous && !transition)
						this->wCodons.emplace_back(i,j);
					else
						this->kwCodons.emplace_back(i,j);
				}
			}
		}



}

void CodonModel::summarize()
{
	cout << "M0 codon model summary" << endl;
	cout << "kappa " << parameters[0] << endl;
	cout << "omega " << parameters[1] << endl;

}

void CodonModel::calculateModel()
{
	DUMP("Codon Model: calculate model");

	double k = parameters[0];
	double w = parameters[1];
	double kw = k*w;
	unsigned int qpos;
	unsigned int i,j;
	double sum, meanRate = 0;

	for(auto pi : piCodons){
		qpos = pi.first * matrixSize + pi.second;
		qMatrix[qpos] = this->piFreqs[pi.second];
	}

	for(auto kappa : kCodons){
			qpos = kappa.first * matrixSize + kappa.second;
			qMatrix[qpos] = k * this->piFreqs[kappa.second];
	}

	for(auto omega : wCodons){
			qpos = omega.first * matrixSize + omega.second;
			qMatrix[qpos] = w * this->piFreqs[omega.second];
	}

	for(auto komega : kwCodons){
			qpos = komega.first * matrixSize + komega.second;
			qMatrix[qpos] = k* w * this->piFreqs[komega.second];
	}

	for (i=0; i< this->matrixSize; i++)
		{
			qMatrix[i*matrixSize+i] = 0;
			sum = 0.0;
			for (j=0; j < this->matrixSize; j++)
			{
				sum -= qMatrix[i*matrixSize + j];
			}
			qMatrix[i*matrixSize + i] = sum;
			meanRate -= sum*piFreqs[i];
		}
		for (i=0; i< this->matrixSize; i++)
		{
				for (j=0; j < this->matrixSize; j++)
				{
					qMatrix[i*matrixSize + j] /= meanRate;
				}
		}

		//TODO Should I multiply it by Pis ?

		this->doEigenDecomposition();

	DUMP("Codon Model: calculated");

}

CodonModel::~CodonModel()
{
	if (this->parameters)
		delete[] this->parameters;
}

} /* namespace EBC */

