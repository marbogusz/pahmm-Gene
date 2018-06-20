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

#include "core/Dictionary.hpp"
#include "core/Definitions.hpp"
#include "core/HmmException.hpp"
#include <algorithm>

namespace EBC
{


//FIXME  - U and T equivalence!!!

Dictionary::Dictionary(Definitions::FrequencyScheme fs) : fScheme(fs), freqsCalculated(false)
{
	simpleCount = 0;
}


Dictionary::~Dictionary()
{
	delete[] elementFrequencies;
}

void Dictionary::setAlphabet(const string dict[], unsigned short size)
{
	unsigned short i;

	for (i=0; i<=size;i++){
		alphabet.push_back(dict[i]);
	}

	for(unsigned short i=0; i<=size; i++)
	{
		//this->alphabet.push_back(string(1,dict[i]));
		this->translator.insert(std::make_pair(alphabet[i],new SequenceElement(i==gapId, i, NULL, alphabet[i])));
	}

	//alphabet size does not include gap e.g. size is 4 for nucleotides
	this->alphabetSize = size;


}

void Dictionary::setEqualFrequencies()
{
	for(unsigned int i = 0; i < alphabetSize; i++)
		elementFrequencies[i] = 1.0/alphabetSize;
}

void Dictionary::calculateFreqeuencies()
{
	freqsCalculated = true;
	for (unsigned int i=0; i < alphabetSize; i++){
		elementFrequencies[i] /= simpleCount;
	}
}

double* Dictionary::getElementFrequencies()
{
	if (!freqsCalculated)
		calculateFreqeuencies();
	return elementFrequencies;
}

void Dictionary::outputAlphabet()
{
	cout << "Model dictionary: " << endl;
	for(auto sym : alphabet){
		cout << sym << ",";
	}
}

string& Dictionary::getSymbolAt(unsigned char i)
{
	return (alphabet[i]);
}

SequenceElement* Dictionary::getSequenceElement(string& symbol)
{
	return translator[symbol];
}

vector<SequenceElement*>* Dictionary::translate(string& sequence, bool disregardIndels)
{

	vector<SequenceElement*> *translatedVector = new vector<SequenceElement*>(sequence.size());
	unsigned short currentEl;
	SequenceElement* se = NULL;


	unsigned int pos = 0;
	for(string::iterator it = sequence.begin(); it < sequence.end(); it++)
	{
		string cstrg = string(1, *it);
		se = getSequenceElement(cstrg);
		elementFrequencies[se->getMatrixIndex()] += 1;
		(*translatedVector)[pos] = se;
		pos++;
		simpleCount++;

	}
	return translatedVector;

}

const string Dictionary::nucleotides[] = {"T", "C", "A", "G", "-"};

//FIXME - T and U equivalence needs to be addressed
//const string Dictionary::UracilSymbol = "U";
//const unsigned int Dictionary::UracilId = 0;

//const char Dictionary::nucleotides[] = {'A', 'C', 'G', 'T'};
const string Dictionary::aminoacids[] = {"A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","-"};
									//	  0   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18  19  20

const string Dictionary::codons[] = {"TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG","TAT", "TAC", "TAA", "TAG", "TGT", "TGC", "TGA", "TGG",
									 "CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG","CAT", "CAC", "CAA", "CAG", "CGT", "CGC", "CGA", "CGG",
									 "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA", "ACG","AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG",
									 "GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA", "GCG","GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG", "---"};


//aminoacid Id to codon conversion
const int Dictionary::geneticCode[] = {13,13,10,10,15,15,15,15,18,18,-1,-1, 4, 4,-1,17,
									   10,10,10,10,14,14,14,14, 8, 8, 5, 5, 1, 1, 1, 1,
                                        9, 9, 9,12,16,16,16,16, 2, 2,11,11,15,15, 1, 1,
                                       19,19,19,19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7 };


const string Dictionary::gapChar = "-";

unsigned short Dictionary::getAlphabetSize()
{
	return this->alphabetSize;
}

}//Namespace definition
