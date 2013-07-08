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


// output level
#ifdef NDEBUG
	OutputLevel defOutLvl = TrInfo;    // default value for release code
#else
	OutputLevel defOutLvl = TrDetail2; // default value for debug code
#endif
OutputLevel outLvl = defOutLvl; // can be changed later


/**
	TO-DO list
	- \todo Add shuffling of the outputs, either by shuffling the first margin
	        (this does not work now), or by shuffling the results at the end,
	        either in the class or in the output methods (probably in the class,
	        to ensure consistency in case of multiple outputs).
	- \todo Add a generic class with mixed margins
	- \todo Add margins controlled by moments
	- \todo Add a possibility to input incomplete correlation matrix
	        - this could be used to generate values for multiple periods at once,
	          including auto-correlations!
**/


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
	string copType = "";        // type of the copula
	string copParamsF = "";     // file with copula parameters
	string margType = "";       // type of the margins
	string margParamsF = "";    // file with parameters for the margins
	string outputFName = "";    // output file name
	bool transfMargins = false; // transform margins to target distrib.?
	bool outputCopula = false;  // output copula, in addition to the distrib.?
	string outCopFName = "";    // output file name for the copula
	bool copAsRanks = false;    // output copula in terms of ranks?
	int numCandPts = 1;         // minimal number of candidate scenarios
	bool sorted1stMarg = false; // should the first margin be sorted?
	int randSeed;               // random seed
	bool writeProbAllocData = false; // write data for the prob-alloc model?
	int outLvlInt;   // output level as an integer (for easier processing)

	stringstream helpHeader;
	helpHeader << "Copula generation code by Michal Kaut" << endl << endl
	           << "Usage: " << argv[0] << " --cop-type <value> [--input <file>]"
	           << " [options] nmb-scens" << endl
	           << " Help: " << argv[0] << " --help" << endl;

	// NEW: re-implementing using boost::program_options
	try{
		// generic options; allowed only from the command line
		prOpt::options_description genOpt("Generic options");
		std::string configFName; // config file name
		genOpt.add_options()
			("version,v", "print version string")
			("help,h",    "produce help message")
			("config", prOpt::value<string>(&configFName), "config file name")
			("out-lvl,l", prOpt::value<int>(&outLvlInt)
			              ->default_value(static_cast<int>(defOutLvl)),
			              "level of output")
			("all-help",  "help including hidden options")
			;

		// main options - allowed from both command line and config. file
		prOpt::options_description confOpt("Main options");
		confOpt.add_options()
			("output,o", prOpt::value<string>(&outputFName)
			             ->default_value("cop-gen.out"), "output file")
			("rand-seed,s", prOpt::value<int>(&randSeed)
			                ->default_value(-1), "random seed (-1 means time-based)")
			("save-cop", prOpt::bool_switch(&outputCopula),
			             "save copula, in addition to the result")
			("cop-output", prOpt::value<string>(&outCopFName),
			               "copula output file (implies --save-cop)")
			("cop-as-ranks,r", prOpt::bool_switch(&copAsRanks),
			                   "output copula as ranks, instead of U(0,1)")
			("sort-marg-1", prOpt::bool_switch(&sorted1stMarg),
			                "sort the output after the 1st marg")
			;

		// method-specific options - allowed from both sources
		prOpt::options_description methodOpt("Method-specific options");
		methodOpt.add_options()
			("cop-type,c", prOpt::value<string>(&copType)->required(),
			               "copula type/family")
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
			("nmb-scen", prOpt::value<int>(&nSc)->required(),
			             "number of scenarios")
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

		// process the config file
		if (optV.count("config")) {
			ifstream confF(configFName.c_str());
			if (!confF) {
				throw ios_base::failure("the specified config file does not exist");
			}
			store(prOpt::parse_config_file(confF, confOpt), optV);
		}

		/*
			At this point, we have processed both the command line and config file.
			If we now called 'prOpt::notify(optV)', it would throw an exception
			if any of the required parameters were missing .. even if we called
			the program with "--help"; see stackoverflow.com/questions/5395503/.
			This means that we have to check all the special options first,
			before checking the correctness of the input.
		*/
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

		// now check the input - with throw an exception on errors
		prOpt::notify(optV);

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

		// post-processing
		outLvl = static_cast<OutputLevel>(outLvlInt);
	}
	catch(exception& e) {
		cerr << e.what() << endl;
		cout << endl << "For help, call: " << argv[0] << " --help" << endl;
		exit(1);
	}
	catch(...) {
		cerr << "Unknown error!" << endl;
		exit(2);
	}

	srand(randSeed); // set the random seed


	// -------------------------------------------------------------------
	// setup copula objects + generate the copula scenarios

	CopInfoBy2D::Ptr p2tgCop;
	MarginsInfo::Ptr p2tgMargins;

	CopNameMapT copNameMap;    ///< convert copula name to type
	MargNameMapT margNameMap;  ///< convert margins name to type

	// first check if the given copula type exists
	make_cop_name_map(copNameMap);
	if (copNameMap.count(copType) == 0) {
		cerr << "Unknown copula type `" << copType << "' .. aborting!" << endl;
		cout << "Known copula types are:";
		for (auto cIt = copNameMap.begin(); cIt != copNameMap.end(); ++ cIt)
			cout << " " << cIt->first;
		cout << endl;
		exit(1);
	}
	// copula type is OK -> it is safe to use copNameMap[copType]
	try {
	switch (copNameMap[copType]) {
	case CopTypeID::sample: // "sample"
		MSG (TrInfo, "copula of type 'sample'");
		p2tgCop = boost::make_shared<CopInfoData>(copParamsF);
		if (margType == "") {
			margType = "sample"; // default for sample cop. is sample margins
		}
		if (margType == "sample") {
			// for this, we need CopInfoData::data_vals(), which is specific to
			// this class, i.e. does not exist in CopInfoBy2D -> use casting
			//! using default for the second param -> no post-processing!
			CopInfoData * p2tgCopData = static_cast<CopInfoData *> (p2tgCop.get());
			p2tgMargins = boost::make_shared<SampleMargins>(p2tgCopData->data_vals());
		}
		break;
	case CopTypeID::normal: // "normal"
		MSG (TrInfo, "copula of type 'normal'");
		if (margType == "") {
			margType = "normal"; // default for sample copula
		}
		// create a new object of the specific class
		p2tgCop = boost::make_shared<CopInfoNormal>(copParamsF);
		break;
	case CopTypeID::indep: // independent margins
		MSG (TrInfo, "copula of type 'independent'");
		// create a new object of the specific class
		p2tgCop = boost::make_shared<CopInfoIndep>(copParamsF);
		break;
	case CopTypeID::mixed: // generic mixed of 2D copulas
		cout << "copula of type 'mixed 2D copulas'" << endl;
		// create a new object of the specific class
		p2tgCop = boost::make_shared<CopInfoGen2D>(copParamsF);
		break;
	default:
		cerr << "ERROR: file " << __FILE__ << ", line " << __LINE__
		     << " .. should never be here!" << endl;
		exit(1);
	}
	}
	catch(exception& e) {
		cerr << "Error: There was some problem initializing the copula!" << endl
		     << "       The error message was: " << e.what() << endl;
		exit(1);
	}

	// We need the margins first after we have generated the copula, but
	// it is better to check here, so we can find an input error early..
	if (transfMargins) {
		// check that we have a known class of margins
		make_marg_name_map(margNameMap);
		if (margNameMap.count(margType) == 0) {
			cerr << "Unknown margins type `" << margType << "' .. aborting!" << endl;
			cout << "Known margins types are:";
			for (auto mIt = margNameMap.begin(); mIt != margNameMap.end(); ++ mIt)
				cout << " " << mIt->first;
			cout << endl;
			exit(1);
		}

		// margins type is OK -> it is safe to use margNameMap[margType]
		try {
		switch (margNameMap[margType]) {
		case MargTypeID::sample: // "sample"
			MSG (TrInfo, "margins of type 'sample'");
			if (copNameMap[copType] != CopTypeID::sample) {
				throw std::logic_error("Sample margins are supported only for "
				                       "sample copula!");
			}
			{
				// for this, we need CopInfoData::data_vals(), which is specific to
				// this class, i.e. does not exist in CopInfoBy2D -> use casting
				//! using default for the second param -> no post-processing!
				CopInfoData * p2tgCopData = static_cast<CopInfoData *> (p2tgCop.get());
				p2tgMargins = boost::make_shared<SampleMargins>(p2tgCopData->data_vals());
			}
			break;
		case MargTypeID::normal: // "normal"
			MSG (TrInfo, "margins of type 'normal'");
			try {
				//! !using default for the second param -> no post-processing!
				p2tgMargins = boost::make_shared<NormalMargins>(margParamsF);
			}
			catch(exception & e) {
				cerr << "Error while reading the target means and std. devs "
				     << "from file " << margParamsF << endl;// << e.what() << endl;
				throw; // re-throw the exception
			}
			break;
		case MargTypeID::fixed: // generic mixed of 2D copulas
			MSG (TrInfo, "margins of type 'fixed'");
			// create a new object of the specific class
			throw std::logic_error("not yet implemented");
			break;
		case MargTypeID::mixed: // generic mixed of 2D copulas
			MSG (TrInfo, "margins of type 'mixed'");
			// create a new object of the specific class
			throw std::logic_error("not yet implemented");

/* splitting ifstream by lines into a stream:
			std::ifstream file("plop");
			std::string   line;

			while(std::getline(file, line))
			{
				 std::stringstream   linestream(line);
				 int                 val1;
				 int                 val2;

				 // Read the integers using the operator >>
				 linestream >> val1 >> val2;
			}
*/

			break;
		default:
			cerr << "ERROR: file " << __FILE__ << ", line " << __LINE__
				  << " .. should never be here!" << endl;
			exit(1);
		}
		}
		catch(exception& e) {
			cerr << "Error: There was some problem initializing the margins!" << endl
				  << "       The error message was: " << e.what() << endl;
			exit(1);
		}
	}


	CopulaSample copSc(p2tgCop, nSc, numCandPts);
	copSc.gen_sample();
	// the first margin is sorted (see CopulaSample::gen_sample())
	//  - if we do not want this, we have to shuffle it manually
	if (!sorted1stMarg)
		copSc.shuffle_results();

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
		TRACE (TrDetail, ""); // empty line, for easier reading

/*
		//! TEMP - testing the new mixed class with triangular data
		{
			DimT nMargs;
			DimT i;
			UnivarMargin::Ptr p2marg;

			MatrixD margSpec;
			std::ifstream margSpecF("test_triang.txt");
			if (!margSpecF) {
				throw std::ios_base::failure("Could not open input file `"
				                             "test_triang.txt" "'!");
			}
			margSpecF >> nMargs;
			if (nMargs != p2tgCop->dim())
				throw std::runtime_error ("dimension inconsistency between "
				                          "margins and copula specifications");
			margSpec.resize(nMargs, 3);
			margSpecF >> margSpec;

			p2tgMargins = boost::make_shared<MixedMargins>(nMargs);
			MixedMargins * p2MixedMargs
				= dynamic_cast<MixedMargins *>(p2tgMargins.get());

			for (i = 0; i < nMargs; ++i) {
				p2marg = boost::make_shared<MarginTriang>
				         (margSpec(i,0), margSpec(i,1), margSpec(i,2), true);
				p2MixedMargs->attach_margin(p2marg, i);
			}
		}
*/

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
			oFile << trans(resValues); // no end-line needed here
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
