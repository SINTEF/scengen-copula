// main routine that generates scenarios for a copula given in a file

#include <iostream>
#include <fstream>
#include <ctime> // needed on gcc-win
//#include <tclap/CmdLine.h> // processing of command-line arguments

// testing only!
#include <boost/numeric/ublas/symmetric.hpp>

#include "copula-info.hpp"
#include "cop2Dsample.hpp"
#include "copula-sample.hpp"
#include "margins.hpp"

using namespace std;
using namespace Copula2D;
using namespace CopulaScen;
using namespace CopulaDef;
using namespace MarginDistrib;

#include <boost/program_options.hpp>
#include <sstream>
namespace prOpt = boost::program_options; // short-cut name


int main(int argc, char *argv[]) {

	/*
		The code should work as follows: by default, it outputs the copula.
		If we provide info about the margins, it transforms them and outputs
		the final distribution. If we, in addition, use "--output-cop [filename]",
		it will also output copula to the provided file (might have a default
		value). We should also have a switch to say that we want the copula as
		ranks, instead of U(0,1) values. (And if this is used with margins, it
		should imply "--output-cop".)
	**/

	// variables whose values is read from the command line
	int nSc = 0;                // number of scenarios to generate
	string copType;             // type of the copula
	string copParamsF;          // file with copula parameters
	string margType;            // type of the margins
	string margParamsF;         // file with parameters for the margins
	string outputFName;         // output file name
	bool transfMargins = false; // transform margins to target distrib.?
	bool outputCopula = false;  // output copula, in addition to the distrib.?
	string outCopFName;         // output file name for the copula
	bool copAsRanks = false;    // output copula in terms of ranks?
	int numCandPts = 1;         // minimal number of candidate scenarios
	int randSeed;               // random seed
	bool writeProbAllocData = false; // write data for the prob-alloc model?

	stringstream helpHeader;
	helpHeader << "Copula generation code by Michal Kaut" << endl << endl
	           << "Usage: " << argv[0] << " sample|normal nmb-scens [options]"
	           << endl;

	// NEW: re-implementing using boost::program_options
	try{
		// generic options; allowed only from the command line
		prOpt::options_description genOpt("Generic options");
		std::string configFName; // config file name
		genOpt.add_options()
			("version,v", "print version string")
			("help,h",    "produce help message")
			("config", prOpt::value<string>(&configFName), "config file name")
			("all-help",  "help including hidden options")
			;

		// main options - allowed from both command line and config. file
		prOpt::options_description confOpt("Main options");
		confOpt.add_options()
			("output,o", prOpt::value<string>(&outputFName)
			             ->default_value("cop-gen.out"), "output file")
			("rand-seed,s", prOpt::value<int>(&randSeed)
			                ->default_value(-1), "random seed (-1 means random)")
			("save-cop", prOpt::bool_switch(&outputCopula),
			             "save copula, in addition to the result")
			("cop-output", prOpt::value<string>(&outCopFName),
			               "output file for the copula\n(implies --save-cop)")
			("cop-as-ranks,r", prOpt::bool_switch(&copAsRanks),
			                   "output copula as ranks, instead of U(0,1)")
			;

		// method-specific options - allowed from both sources
		prOpt::options_description methodOpt("Method-specific options");
		methodOpt.add_options()
			("cop-type,c", prOpt::value<string>(&copType), "copula type/family")
			("input,i", prOpt::value<string>(&copParamsF)
			            ->default_value("cop-params.dat"),
			            "input file: copula parameters")
			("marg-type,m", prOpt::value<string>(&margType),
			               "type of the marginal distributions")
			("marg-par,d", prOpt::value<string>(&margParamsF),
			               "input file: params of marg. distrib.")
			;

		// hidden options - allowed everywhere, but not shown to the user
		prOpt::options_description hiddenOpt("Hidden options");
		hiddenOpt.add_options()
			("nmb-scen", prOpt::value<int>(&nSc), "number of scenarios")
			("num-cand-pts", prOpt::value<int>(&numCandPts)->default_value(1),
			                 "min number of cand. scenarios")
			("out-prob-alloc", prOpt::bool_switch(&writeProbAllocData),
			                   "write data for prob-alloc model")
			;

		prOpt::options_description visibleOpt;
		visibleOpt.add(genOpt).add(confOpt).add(methodOpt);

		prOpt::options_description cmdLineOpt; // processed from command line
		cmdLineOpt.add(genOpt).add(confOpt).add(methodOpt).add(hiddenOpt);

		prOpt::options_description confFileOpt; // processed from config. file
		confFileOpt.add(confOpt).add(methodOpt).add(hiddenOpt);

		prOpt::positional_options_description posArg;
		posArg.add("nmb-scen", 1);

		// process the command line parameters
		prOpt::variables_map optV;
		prOpt::store(prOpt::command_line_parser(argc, argv).
		             options(cmdLineOpt).positional(posArg).run(), optV);
		// a simpler syntax without positional options:
		//prOpt::store(prOpt::parse_command_line(argc, argv, allOpt), optV);
		prOpt::notify(optV);

		// process the config file
		ifstream confF(configFName.c_str());
		if (!confF && optV.count("config")) {
			throw ios_base::failure("the specified config file does not exist");
		}
		store(prOpt::parse_config_file(confF, confOpt), optV);
		prOpt::notify(optV);

		// deal with options that need extra action
		if (optV.count("version")) {
			cout << "This code was built at " << __TIME__ << ", " << __DATE__
			     << "." << endl;
			return 0;
		}
		if (optV.count("help")) {
			cout << helpHeader.str() << visibleOpt;
			return 0;
		}
		if (optV.count("all-help")) {
			visibleOpt.add(hiddenOpt); // add the hidden options
			cout << helpHeader.str() << visibleOpt;
			return 0;
		}
		if (!optV.count("nmb-scen")) {
			throw ios_base::failure(helpHeader.str());
		}
		if (randSeed < 0) {
			randSeed = time(NULL);
		}
		if (optV.count("cop-output")) {
			outputCopula = true;
		}
		if (optV.count("marg-type") || optV.count("marg-par")) {
			// some info about marginal distributions -> process margins
			transfMargins = true;
		}
		else {
			// no info about marginal distributions -> output copula
			outputCopula = true;
		}
	}
	catch(exception& e) {
		cerr << e.what() << endl;
		exit(1);
	}

	srand(randSeed); // set the random seed


	// -------------------------------------------------------------------
	// setup copula objects + generate the copula scenarios

	CopInfoBy2D::Ptr p2tgCop;
	MarginsInfo::Ptr p2tgMargins;

	if (copType == "sample") {
		cout << "copula of type 'sample'" << endl;

		// create a new object of the specific class
		CopInfoData * p2tgCopData = new CopInfoData(copParamsF);
		p2tgCopData->setup_2d_targets();

		p2tgCop.reset(p2tgCopData); // p2tgCop takes over the pointer

		if (margType == "") {
			margType = "sample"; // default for sample cop. is sample margins
		}
		if (margType == "sample") {
			// set the margins here, where we have the CopInfoData * pointer
			// !using default for the second param -> no post-processing!
			p2tgMargins.reset(new SampleMargins(p2tgCopData->data_vals()));
		}
	}
	if (copType == "normal") {
		cout << "normal copula" << endl;

		// create a new object of the specific class
		CopInfoNormal * p2tgCopNormal = new CopInfoNormal(copParamsF);
		p2tgCopNormal->setup_2d_targets();

		p2tgCop.reset(p2tgCopNormal); // p2tgCop takes over the pointer

		if (margType == "") {
			margType = "normal"; // default for sample copula
		}
	}

	CopulaSample copSc(p2tgCop, nSc, numCandPts);
	copSc.gen_sample();

	if (outputCopula) {
		string fName;
		if (outCopFName.size() > 0) {
			fName = outCopFName;
		} else {
			if (!transfMargins) {
				fName = outputFName;
			} else {
				// file name not given -> append "_cop" to outputFName
				size_t lastDot = outputFName.find_last_of('.');
				fName = outputFName.substr(0, lastDot) + "_cop"
				        + outputFName.substr(lastDot);
			}
		}
		try {
			copSc.print_as_txt(fName.c_str(), !copAsRanks);
		}
		catch(exception& e) {
			cerr << "Error: There was some problem when saving the copula!"
			     << endl << e.what() << endl;
		}
	}


	// -------------------------------------------------------------------
	// transform the margins to the target distributions
	if (transfMargins) {

		// the sample margins have been set up together with the copula
		if (margType == "sample" && copType != "sample") {
			throw std::logic_error("this combination is not supported");
		}

		if (margType == "normal") {
			try {
				// !using default for the second param -> no post-processing!
				p2tgMargins.reset(new NormalMargins(margParamsF));
			}
			catch(exception & e) {
				cerr << "Error while reading the target means and std. devs "
				     << "from file " << margParamsF << endl << e.what() << endl;
				throw; // re-throw the exception
			}
		}

		MatrixI copRanks;
		copSc.get_result_ranks(copRanks); // get the results in terms of ranks
		MatrixD resValues;
		p2tgMargins->get_margin_distr(copRanks, resValues);

		try {
			std::ofstream oFile;
			oFile.open(outputFName.c_str(), std::ios::out);
			if (!oFile) {
				throw ios_base::failure("could not open output file "
				                         + outputFName);
			}
			// output transposed, to get margins in columns instead of rows
			oFile << trans(resValues) << endl;
			oFile.close();
		}
		catch(exception & e) {
			cerr << "Error while writing the output to file " << outputFName;
			throw; // re-throw the exception
		}
	}


	// -------------------------------------------------------------------
	if (writeProbAllocData) {
		// write the AMPL/GMP data file
		copSc.write_gmp_data();
	}

	return 0;
}
