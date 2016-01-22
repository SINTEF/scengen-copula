/// main function for the forecast-error-based generator
/**
	\note The main function cannot be in \c cop-gen_fc-err.cpp, since
	      this file is used also with other build targets.
**/

#include "cop-gen_fc-err.hpp"
#include "copula-sample.hpp"
#include "margins.hpp"

//#include <boost/numeric/ublas/assignment.hpp>
#include <iostream>
#include <fstream>
//#include <sstream>
//#include <ctime> // needed on gcc-win
#include <boost/program_options.hpp>

using namespace Copula2D;
using namespace FcErr_Gen;
namespace prOpt = boost::program_options; // short-cut name
using std::string;
using std::cout;
using std::endl;

// output level
#ifdef NDEBUG
	OutputLevel defOutLvl = TrInfo;    // default value for release code
#else
	OutputLevel defOutLvl = TrDetail2; // default value for debug code
#endif
OutputLevel outLvl = defOutLvl; // can be changed later


int main(int argc, char *argv[]) {

	string histDataFName = "";  // file with all historical data
	string forecastFName = "";  // file with the current forecast
	string outputFName = "";    // output file name
	DimT nSc = 0;               // number of scenarios to generate
	DimT nPer = 0;              // number of periods to generate
	DimT nVar = 0;              // number of variables/margins
	string branchingStr;        // list of branching factors as a string
	VectorI branching;          // vector of branching factors
	string curValStr;           // current values as a comma-separated list
	int perVarDt = 1;           // gen. 2D copulas up to t ± dt, for given var
	int intVarDt = 0;           // gen. 2D copulas up to t ± dt, between vars.
	HistDataSort hDataSort = HistDataSort::fCastTimeAsc;  // for the input file
	bool fcastInclCur = false;  // forecast starts with current values

	int randSeed;               // random seed
	int outLvlInt;              // output level as an integer

	MatrixD histData;
	MatrixD forecast;
	VectorD curVal;

	DimT i;

	std::stringstream helpHeader;
	helpHeader << "Scenario-generation code using historical forecast errors"
	           << endl << "author: Michal Kaut" << endl << endl
	           << "Usage: " << argv[0] << " --hist-data <file> --forecast <file>"
	           << " [options]" << endl
	           << " Help: " << argv[0] << " --help" << endl;

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
			("hist-data,d", prOpt::value<string>(&histDataFName)->required(),
			                "file with historical data")
			("forecast,f", prOpt::value<string>(&forecastFName)->required(),
			               "file with current forecast")
			("scens,s", prOpt::value<DimT>(&nSc), "number of scenarios")
			("periods,t", prOpt::value<DimT>(&nPer), "number of periods")
			("branching,b", prOpt::value<string>(&branchingStr),
			                "branches per per. (comma-separated list)")
			("cur-val,c", prOpt::value<string>(&curValStr),
			              "current values (comma-separated list)")
			("per-var-dt", prOpt::value<int>(&perVarDt)->default_value(1),
			               "gen. 2D copulas up to t±dt, for given var")
			("int-var-dt", prOpt::value<int>(&intVarDt)->default_value(0),
			               "gen. 2D copulas up to t±dt, between vars")
			("fcast-incl-cur", prOpt::bool_switch(&fcastInclCur),
			                   "forecast starts with current values")
			;

		// hidden options - allowed everywhere, but not shown to the user
		prOpt::options_description hiddenOpt("Hidden options");
		hiddenOpt.add_options()
			("out-lvl,l", prOpt::value<int>(&outLvlInt)
			              ->default_value(static_cast<int>(defOutLvl)),
			              "level of output")
			("rand-seed,s", prOpt::value<int>(&randSeed)->default_value(-1),
			                "random seed (-1 means time-based)")
			;

		prOpt::options_description visibleOpt;
		visibleOpt.add(genOpt).add(confOpt);

		prOpt::options_description cmdLineOpt; // processed from command line
		cmdLineOpt.add(genOpt).add(confOpt).add(hiddenOpt);

		prOpt::options_description confFileOpt; // processed from config. file
		confFileOpt.add(confOpt).add(hiddenOpt);

		//prOpt::positional_options_description posArg;
		//posArg.add("nmb-scen", 1);

		// process the command line parameters
		prOpt::variables_map optV;
		//prOpt::store(prOpt::command_line_parser(argc, argv).
		//             options(cmdLineOpt).positional(posArg).run(), optV);
		// a simpler syntax without positional options:
		prOpt::store(prOpt::parse_command_line(argc, argv, cmdLineOpt), optV);

		// process the config file
		if (optV.count("config")) {
			std::ifstream confF(configFName.c_str());
			if (!confF) {
				throw std::ios_base::failure
					("the specified config file does not exist");
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

		// NB: optV.count("fcast-incl-cur") = 1 always!
		if (fcastInclCur && optV.count("cur-val") > 0)
			throw std::runtime_error("options --fcast-incl-cur and --curval "
			                         "are mutually exclusive.");

		// read the input files
		std::ifstream inFileStr(histDataFName);
		if (!inFileStr) {
			throw std::ios_base::failure("Could not open input file `"
			                             + histDataFName + "'!");
		}
		inFileStr >> histData;
		inFileStr.close();
		TRACE(TrDetail, "histData: [" << histData.size1() << ", "
		                << histData.size2() << "]");
		//
		inFileStr.open(forecastFName);
		if (!inFileStr) {
			throw std::ios_base::failure("Could not open input file `"
			                             + forecastFName + "'!");
		}
		inFileStr >> forecast;
		inFileStr.close();
		if (fcastInclCur) {
			TRACE(TrDetail, "forecast incl. cur.: [" << forecast.size1() << ", "
			                << forecast.size2() << "]");
			DISPLAY_NL(forecast);
			curVal = ublas::row(forecast, 0);
			for (i = 1; i < forecast.size1(); ++i)
				ublas::row(forecast, i - 1) = ublas::row(forecast, i);
			forecast.resize(forecast.size1() - 1, forecast.size2(), true);
			DISPLAY_NL(forecast);
		}
		TRACE(TrDetail, "forecast: [" << forecast.size1() << ", "
		                << forecast.size2() << "]");
		// dimension checks
		nVar = forecast.size2();
		DimT T = histData.size2() / nVar - 1;
		if (histData.size2() != nVar * (T + 1))
			throw std::length_error
				("Inconsistent dimension of input matrices.");
		if (optV.count("periods")) {
			// number of periods given -> check that we have enough data
			if (T < nPer)
				throw std::length_error
					("Not enough columns in the historical-data matrix.");
			if (forecast.size1() < nPer)
				throw std::length_error
					("Not enough rows in the forecast matrix.");
		} else {
			// number of periods not given -> imply from other data
			// (nPer = min(dt in forecast, dt in histData)
			nPer = forecast.size1();
			if (T < nPer)
				nPer = T;
		}

		// post-processing - keep it here, so we can throw input-errors
		outLvl = static_cast<OutputLevel>(outLvlInt);
		srand(randSeed >= 0 ? randSeed : time(nullptr)); // set the random seed
		if (optV.count("branching")) {
			if (optV.count("scens"))
				throw std::runtime_error
					("'branching' and 'scens' are mutually exclusive");
			// parse the string into VectorI
			// problem: ublas::vector does not have a push_back funtion!
			std::stringstream brStrStr(branchingStr);
			std::vector<DimT> brVector;
			// src: http://stackoverflow.com/a/1894955/842693
			while (brStrStr >> i) {
				brVector.push_back(i);
				if (brStrStr.peek() == ',')
					brStrStr.ignore();
			}
			// compare the length of the vector to nPer (if provided)
			if (optV.count("periods")) {
				// specified number of periods - has to have enough data
				if (brVector.size() < nPer) {
					throw std::length_error
						("not enough periods in the branching list");
				}
			} else  {
				assert (nPer > 0
				        && "no --periods -> nPer is set from the forecast");
				if (brVector.size() > nPer) {
					// not enough periods in either forecast or hist. data
					throw std::length_error("branching list contains more "
					                        "periods than we have data for");
				}
				if (brVector.size() < nPer) {
					// branching contains fewer periods than the forecast
					nPer = brVector.size();
				}
			}
			// now we have nPer > 0 && nPer <= length of the forecast
			branching.resize(nPer);
			for (i = 0; i < nPer; ++i)
				branching(i) = brVector[i];
			DBGSHOW(TrDetail, branching);
		} else {
			// no 'branching' -> need scenarios
			if (!optV.count("scens"))
				throw std::runtime_error("needs either --scens or --branching");
		}
		assert(nPer > 0 && "should have the number of periods at this point");
		// dimension check:
		if (histData.size2() < nVar * (nPer+1))
			throw std::length_error
				("Inconsistent dimension - not enough columns in hist. data");

	}
	catch(std::exception & e) {
		std:: cerr << "Input error: " << e.what() << endl;
		cout << endl << "For help, call: " << argv[0] << " --help" << endl;
		exit(1);
	}
	catch(...) {
		std::cerr << "Unknown error!" << endl;
		exit(2);
	}

	DBGSHOW(TrDetail, nVar);
	DBGSHOW(TrDetail, nPer);
	DBGSHOW(TrDetail, nSc);

	FcErrTreeGen scenGen(nVar, histData, hDataSort, perVarDt, intVarDt);
	ScenTree scTree;

	if (nSc > 0) {
		// two-stage tree
		INFO("Generating a 2-stage tree with " << nSc << " scenarios and "
		)
		scenGen.gen_2stage_tree(forecast, nSc, scTree);
	} else {
		// multi-stage tree
		assert (branching.size() > 0);
		scenGen.gen_reg_tree(forecast, branching, scTree);
	}
	if (curValStr.size() > 0) {
		// parse the string into std::vector
		// problem: ublas::vector does not have a push_back funtion!
		std::stringstream curValStrStr(curValStr);
		std::vector<double> curValV;
		double x;
		// src: http://stackoverflow.com/a/1894955/842693
		while (curValStrStr >> x) {
			curValV.push_back(x);
			if (curValStrStr.peek() == ',')
				curValStrStr.ignore();
		}
		DBGSHOW(TrDetail, curValStr);
		DBGSHOW(TrDetail, curValV);
		// now copy it to the scenario-tree object
		if (curValV.size() == nVar) {
			scTree.set_root_values(curValV);
		} else {
			std::cerr << "Warning: wrong size of the current-value vector "
			             "- ignoring it!" << std::endl;
		}
	}

	scTree.display_per_scen();
	scTree.display_per_var();

	scTree.make_gnuplot_charts(&forecast);

	return 0;
}
