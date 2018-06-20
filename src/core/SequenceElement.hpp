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


#ifndef SEQUENCEELEMENT_H_
#define SEQUENCEELEMENT_H_

#include <vector>
#include <string>
using namespace std;


namespace EBC
{

class SequenceElement
{
//FIXME - this is potentially slow - rework

protected:
	bool isGap;
	unsigned char matrixIndex;
	string symbol;
	//vector<short> alternativeIndexes;
public:
	SequenceElement(bool, unsigned char, short*, string smbl);

	inline bool isIsGap() const
	{
		return isGap;
	}

	inline unsigned char getMatrixIndex()
	{
		return matrixIndex;
	}

	inline string& getSymbol()
	{
		return symbol;
	}
};

} /* namespace EBC */
#endif /* SEQUENCEELEMENT_H_ */
