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


#include <core/PMatrixDouble.hpp>

namespace EBC
{

PMatrixDouble::PMatrixDouble(SubstitutionModelBase* m) : PMatrix(m)
{
	this->fastPairGammaPt = new double[matrixFullSize];
	this->fastLogPairGammaPt = new double[matrixFullSize];
	this->sitePatterns = new double*[matrixSize+1];
	for (int i =0; i<= matrixSize; i++ )
	{
		sitePatterns[i] = new double[matrixSize+1];
	}

}

PMatrixDouble::~PMatrixDouble()
{
	delete [] fastPairGammaPt;
	delete [] fastLogPairGammaPt;

	for (int i =0; i<= matrixSize; i++ )
	{
			delete[] sitePatterns[i];
	}
	delete[] sitePatterns;
}

void PMatrixDouble::calculatePairSitePatterns()
{
	//includes gaps - does not discard missing data!
	for (int i =0; i<= matrixSize; i++ )
		for (int j =0; j<= matrixSize; j++ )
		{
			if (i == j && i == matrixSize)
				continue;
			if (i == matrixSize)
			{
				sitePatterns[i][j]  = getLogEquilibriumFreq(j);
				if (std::isnan(sitePatterns[i][j]))
					getLogEquilibriumFreq(j);
			}
			else if (j == matrixSize)
			{
				sitePatterns[i][j]  = getLogEquilibriumFreq(i);
				if (std::isnan(sitePatterns[i][j]))
					getLogEquilibriumFreq(i);
			}
			else
			{
				sitePatterns[i][j] = getLogPairTransition(i,j);
				if (std::isnan(sitePatterns[i][j]))
					getLogPairTransition(i,j);
			}
		}
	sitePatterns[matrixSize][matrixSize] = 0;
}

void PMatrixDouble::calculate()
{
	if (time != 0)
	{
		std::fill(fastPairGammaPt, fastPairGammaPt+matrixFullSize, 0);

		for(unsigned int i = 0; i< rateCategories; i++)
		{
			if (ptMatrices[i] != NULL)
				delete [] ptMatrices[i];
			ptMatrices[i] = this->model->calculatePt(time, i);
			for (int j=0; j< matrixFullSize; j++)
			{
				fastPairGammaPt[j] += ptMatrices[i][j] * model->gammaFrequencies[i];
				fastLogPairGammaPt[j] = fastPairGammaPt[j] < Definitions::almostZero ? Definitions::minMatrixLikelihood : log(fastPairGammaPt[j]);
			}

		}

		calculatePairSitePatterns();
	}
	else
		throw HmmException("PMatrixDouble : attempting to calculate p(t) with t set to 0");
}


double PMatrixDouble::getPairTransition(array<unsigned int, 2>& nodes)
{
	return getPairTransition(nodes[0],nodes[1]);
}

double PMatrixDouble::getPairTransition(unsigned int xi, unsigned int yi)
{
	return model->getEquilibriumFrequencies(xi) * fastPairGammaPt[xi*matrixSize+yi];
}

double PMatrixDouble::getLogPairTransition(unsigned int xi, unsigned int yi)
{
	return model->getLogEquilibriumFrequencies(xi) + fastLogPairGammaPt[xi*matrixSize+yi];
}

void PMatrixDouble::summarize()
{
	cout << "P(t) matrix summary :" << endl;
	cout << "Divergence time : " << time << endl;

	cout << "Pairwise Site patterns " << endl;
	for (int i =0; i<= matrixSize; i++ )
	{
		for (int j =0; j<= matrixSize; j++ )
		{
			cout << sitePatterns[i][j] << "\t\t";
		}
		cout << endl;
	}

	cout << endl << "Averaged gamma cat P(t) matrix" << endl;
	for (int i =0; i< matrixSize; i++ )
	{
		for (int j =0; j<  matrixSize; j++ )
		{
			cout << fastPairGammaPt[(i*matrixSize)+j] << "\t\t";
		}
		cout << endl;
	}
	cout << endl;

}

} /* namespace EBC */


