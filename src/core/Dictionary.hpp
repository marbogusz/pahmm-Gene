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



#ifndef DICTIONARY_H_
#define DICTIONARY_H_

#include <vector>
#include <string>
#include <map>
#include "core/SequenceElement.hpp"
#include "core/Definitions.hpp"

using namespace std;

namespace EBC
{
	class Dictionary
	{
	protected:
		unsigned short alphabetSize;
		
		unsigned char gapId;

		Definitions::FrequencyScheme fScheme;

		//equilibruim frequencies based on various strategies
		double* elementFrequencies;	

		vector<string> alphabet;
		
		map<string,SequenceElement*> translator;

		void setEqualFrequencies();

		bool freqsCalculated;

		unsigned int simpleCount;

		virtual void calculateFreqeuencies();

	public:

		static const string nucleotides[5];
		static const string aminoacids[21];
		static const string codons[65];
		static const int geneticCode[64];
		static const string gapChar;

		virtual vector<SequenceElement*>* translate(string &sequence, bool disregardIndels = false);

		virtual unsigned short getAlphabetSize();

		Dictionary(Definitions::FrequencyScheme fs);

		virtual ~Dictionary();

		SequenceElement* getSequenceElement(string& symbol);

		virtual string& getSymbolAt(unsigned char i);

		virtual void outputAlphabet();

		virtual double* getElementFrequencies();

		inline unsigned char getGapID()
		{
			return gapId;
		}

	protected:
		virtual void setAlphabet(const string alphabet[], unsigned short size);

	};
}


#endif /* DICTIONARY_H_ */
