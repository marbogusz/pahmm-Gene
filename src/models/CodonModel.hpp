/*
 * NucleotideSubstitutionModel.hpp
 *
 *  Created on: Jan 13, 2014
 *      Author: root
 */

#ifndef CODON_MODEL_H_
#define CODON_MODEL_H_

#include "core/Dictionary.hpp"
#include "core/CodonDictionary.hpp"
#include "core/Definitions.hpp"
#include "models/SubstitutionModelBase.hpp"
#include "core/Maths.hpp"
#include <cmath>

namespace EBC
{

class CodonModel : public EBC::SubstitutionModelBase
{

protected:


	vector<pair<unsigned int, unsigned int> > piCodons;
	vector<pair<unsigned int, unsigned int> > kCodons;
	vector<pair<unsigned int, unsigned int> > wCodons;
	vector<pair<unsigned int, unsigned int> > kwCodons;

	//now we can check for syn/nonsyn/ti/tv
	CodonDictionary* cdict = dynamic_cast<CodonDictionary*>(dictionary);

	//params[0] - kappa
	//params[1] - omega

	void buildInitialQmatrix();

public:

	CodonModel(Dictionary*, Maths*i, unsigned int);

	virtual ~CodonModel();

	void calculateModel();

	void setParameters(const vector<double>&);

	void summarize();

};

} /* namespace EBC */
#endif /* MODEL_H_ */
