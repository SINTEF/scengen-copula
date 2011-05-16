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

using namespace std;
using namespace Copula2D;
using namespace CopulaScen;
using namespace CopulaDef;

#include <boost/program_options.hpp>
#include <sstream>
namespace prOpt = boost::program_options; // short-cut name


int main(int argc, char *argv[]) {

	// variables whose values is read from the command line
	int nSc = 0;               // number of scenarios to generate
	std::string copType;       // type of the copula
	std::string copParamsF;    // file with copula parameters
	std::string margParamsF;   // file with parameters for the margins
	std::string outputFName;   // output file name
	int numCandPts = 1;        // minimal number of candidate scenarios
	bool writeProbAllocData = false; // write data for the prob-alloc model?
	bool outputRanks = false;  // output margins as ranks?
	bool outputCopula = false; // output margins as values from U(0,1)?
	int randSeed;              // random seed

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
			("config,c", prOpt::value<string>(&configFName), "config file name")
			("all-help",  "help including hidden options")
			;

		// main options - allowed from both command line and config. file
		prOpt::options_description confOpt("Main options");
		confOpt.add_options()
			("output,o", prOpt::value<string>(&outputFName)
			             ->default_value("cop-gen.out"), "output file")
			("rand-seed,r", prOpt::value<int>(&randSeed)
			                ->default_value(-1), "random seed (-1 means random)")
			("output-ranks", prOpt::bool_switch(&outputRanks),
			                 "output margins as ranks")
			("output-copula", prOpt::bool_switch(&outputCopula),
			                  "output margins as values from U(0,1)")
			;

		// method-specific options - allowed from both sources
		prOpt::options_description methodOpt("Method-specific options");
		methodOpt.add_options()
			("tg-dist,t", prOpt::value<string>(&copParamsF)
			             ->default_value("target-dist.dat"),
			              "file: target distrib. (for cop=sample)")
			("cop-par,p", prOpt::value<string>(&copParamsF)
			              ->default_value("cop-params.dat"),
			              "file: copula params (for cop=normal)")
			("marg-par,m", prOpt::value<string>(&margParamsF)
			               ->default_value("marg-params.dat"),
			               "file: margin params (for cop=normal)")
			;

		// hidden options - allowed everywhere, but not shown to the user
		prOpt::options_description hiddenOpt("Hidden options");
		hiddenOpt.add_options()
			("nmb-scen", prOpt::value<int>(&nSc), "number of scenarios")
			("cop-type", prOpt::value<string>(&copType), "copula type/family")
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
		posArg.add("cop-type", 1);
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
	}
	catch(exception& e) {
		cerr << e.what() << endl;
		exit(1);
	}

	srand(randSeed); // set the random seed


	CopInfoBy2D::Ptr p2tgCop;

	if (copType == "sample") {
		cout << "copula of type 'sample'" << endl;

		// create a new object of the specific class
		CopInfoData * p2tgCopData = new CopInfoData(copParamsF);
		p2tgCopData->setup_2d_targets();

		p2tgCop.reset(p2tgCopData); // p2tgCop takes over the pointer
	}
	if (copType == "normal") {
		cout << "normal copula" << endl;

		// create a new object of the specific class
		CopInfoNormal * p2tgCopNormal = new CopInfoNormal(copParamsF);
		p2tgCopNormal->setup_2d_targets();

		p2tgCop.reset(p2tgCopNormal); // p2tgCop takes over the pointer
	}

	CopulaSample copSc(p2tgCop, nSc, numCandPts);

	copSc.gen_sample();
	copSc.print_as_txt(outputFName.c_str(), true);

	if (writeProbAllocData) {
		// write the AMPL/GMP data file
		copSc.write_gmp_data();
	}

	// cleaning
//	while (p2copData.size() > 0) {
//		delete p2copData.back(); // deletes the object
//		p2copData.pop_back(); // removes it from the list
//	}

	return 0;
}
