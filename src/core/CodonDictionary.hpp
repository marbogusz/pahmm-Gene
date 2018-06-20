//==============================================================================
// Pair-HMM dN/dS and indel rate estimator
// 
// Copyright (c) 2015-2018 Marcin Bogusz.
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


#ifndef CODONDICTIONARY_H_
#define CODONDICTIONARY_H_

#include <vector>
#include <string>
#include <map>
#include "core/SequenceElement.hpp"
#include "core/Definitions.hpp"
#include "core/Dictionary.hpp"

using namespace std;

namespace EBC
{

	class CodonDictionary : public Dictionary
	{
	protected:
		int reducedGenCode[64];

		bool isPurine(char base);
		bool isPyramidine(char base);

		void calculateFreqeuencies();

	public:
		CodonDictionary(Definitions::FrequencyScheme fs);

		vector<SequenceElement*>* translate(string &sequence, bool disregardIndels = false);

		//void setGeneticCode(const int[]);

		int getAminoacidId(unsigned int codonId);

		unsigned int getNumberOfDifferentPositions(unsigned int codon1, unsigned int codon2,
				bool& synonymous, bool& transition);
	};
}


#endif /* DICTIONARY_H_ */
