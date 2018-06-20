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


#ifndef AADICTIONARY_H_
#define AADICTIONARY_H_

#include <vector>
#include <string>
#include <map>
#include "core/SequenceElement.hpp"
#include "core/Definitions.hpp"
#include "core/Dictionary.hpp"

using namespace std;

namespace EBC
{
	class AminoacidDictionary : public Dictionary
	{
	public:
		AminoacidDictionary(Definitions::FrequencyScheme fs);
	};
}


#endif /* DICTIONARY_H_ */
