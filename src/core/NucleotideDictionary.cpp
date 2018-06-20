//==============================================================================
// Pair-HMM dN/dS and indel rate estimator
// 
// Copyright (c) 2015-2018 Marcin Bogusz.
// 
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

#include "core/NucleotideDictionary.hpp"
#include "core/Definitions.hpp"
#include "core/HmmException.hpp"
#include <algorithm>

namespace EBC
{

NucleotideDictionary::NucleotideDictionary(Definitions::FrequencyScheme fs) : Dictionary(fs)
{
	gapId = 4;
	this->setAlphabet(Dictionary::nucleotides,4);
	this->translator.insert(std::make_pair("U",new SequenceElement(false, 0, NULL, "U")));

	elementFrequencies = new double[alphabetSize];

	for (unsigned int i=0; i<alphabetSize; i++)
				this->elementFrequencies[i] = 0;

	if (fs == Definitions::FrequencyScheme::Equal){
		setEqualFrequencies();
	}



}

}//Namespace definition
