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


#include "core/CommandReader.hpp"
#include "core/FileParser.hpp"
#include <sstream>
#include <cstring>
#include <cstdlib>

using namespace std;

namespace EBC
{

CommandReader::CommandReader(int argc, char** argv)
{
	try
	{
		/*
		parser.add_option("V", "Run Viterbi algorithm using user parameters");
		parser.add_option("F", "Run Forward algorithm");
		parser.add_option("M", "Run MLE");
		parser.add_option("X", "Run Forward pair dist est with specified parameters");
		parser.add_option("Y", "Run Viterbi pair dist est with specified parameters");
		parser.add_option("fa", "Fixed Alignment");
		*/
		parser.add_option("in","This option takes one argument which specifies the name of the file we want to analyze",1);
		/*
		parser.add_option("rev", "REV Substitution Model");
		parser.add_option("hky", "HKY85 Substitution Model");
		parser.add_option("lg", "Le & Gasquel AA Substitution Model");
		*/
		parser.add_option("m0", "M0 codon model");
		/*
		parser.add_option("i","indel parameters (NB probability and rate)",2);
		//FIXME - remove
		parser.add_option("d","evolutionary distance",1);
		parser.set_group_name("Miscellaneous Options");

		parser.add_option("b","Toggle banding, default is no banding");
		parser.add_option("o","Set optimizer, 0- BFGS, 1- BOBYQA default is 0",1);
		parser.add_option("param_rev","GTR model parameters",5);
		parser.add_option("param_hky","HKY85 model parameters",1);
		parser.add_option("param_m0","M0 model parameters",2);
		parser.add_option("ov","Output viterbi alignment for estimated parameters");
		parser.add_option("rateCat", "Specify gamma rate categories, default is 5",1);
		parser.add_option("initAlpha", "Specify initial alpha parameter, default is 0.5",1 );
		parser.add_option("estimateAlpha", "Specify to estimate alpha 0|1, default is 1",1 );
		parser.add_option("pi", "Specify element frequency scheme, default is empirical count",1 );
*/

		parser.add_option("h","Display this help message.");
		parser.add_option("lE", "log error");
		parser.add_option("lW", "log warning");
		parser.add_option("lI", "log info");
		parser.add_option("lD", "log debug");
		parser.add_option("lDD", "log dump");

		parser.parse(argc,argv);

		//const char* one_time_opts[] = {"V", "F", "M","X", "Y", "in", "i","d" ,"h","b","o", "ov"};
		const char* one_time_opts[] = {"in","m0"};
		parser.check_one_time_options(one_time_opts);
/*
		parser.check_incompatible_options("V", "F");
		parser.check_incompatible_options("X", "Y");
		parser.check_incompatible_options("rev", "hky");
		//parser.check_incompatible_options("d", "F");



		const char* f_sub_opts[] = {"b","o","ov"};
		const char* rev_sub_opts[] = {"param_rev"};
		const char* hky_sub_opts[] = {"param_hky"};
		const char* m0_sub_opts[] = {"param_m0"};
		parser.check_sub_options("F", f_sub_opts);
		parser.check_sub_options("rev", rev_sub_opts);
		parser.check_sub_options("hky", hky_sub_opts);
		parser.check_sub_options("m0", m0_sub_opts);

		parser.check_option_arg_range("param_hky", 0.0000001, 20.0);
		parser.check_option_arg_range("param_rev", 0.0, 10.0);
		parser.check_option_arg_range("i", 0.0, 1.0);
		parser.check_option_arg_range("d", 0.0000001, 3.5);
		parser.check_option_arg_range("initAlpha", 0.0000001, 1000.0);
*/
		if (!parser.option("m0") && !parser.option("in"))
		{
		    cout << "Usage: paHMM-Gene --m0 --in input_file\n";
		    parser.print_options();
			throw HmmException("Specify the model and file name\n");
		}

		/*
		parser.check_option_arg_range("o", 0, 1);
		parser.check_option_arg_range("estimateAlpha", 0, 1);
		parser.check_option_arg_range("rateCat", 0, 1000);
		 */
		if (parser.option("h"))
		{
			// display all the command line options
		    cout << "Usage: paHMM-Gene --m0 --in input_file\n";
		    parser.print_options();
		}
	}
	catch (exception& e)
	{
	        throw HmmException(e.what());
	}
}

/*
vector<double> CommandReader::getSubstParams()
{
	int i;
	vector<double> vec;
if (parser.option("m0"))
	{
			if (parser.option("param_m0"))
			{
				for (i=0; i< 2; i++)
				{
					DEBUG("m0 parameter " << i <<  ": " << parser.option("param_m0").argument(i));
					vec.push_back(atof(parser.option("param_m0").argument(i).c_str()));
				}
			}
	}
	else throw HmmException("Model not specified");

	return vec;
}

*/
IParser* CommandReader::getParser() throw (HmmException&)
{
	if (parser.option("in"))
	{
		return new FileParser((string(parser.option("in").argument())).c_str());
	}
	else
		throw HmmException("input file not specified");
}

} /* namespace EBC */
