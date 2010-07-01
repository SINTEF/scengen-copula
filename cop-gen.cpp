// main routine that generates scenarios for a copula given in a file

#include <iostream>
#include <fstream>
//#include <ctime>
#include <tclap/CmdLine.h> // processing of command-line arguments

#include "cop2Dsample.hpp"
#include "copula-sample.hpp"

using namespace std;
using namespace Copula2D;
using namespace CopulaScen;


int main(int argc, char *argv[]) {

	int i, j;

	// variables whose values is read from the command line
	int nSc = 0;             // number of scenarios to generate
	std::string tgDistFName; // input file name
	std::string outputFName; // output file name

	// parameters for processing of command-line arguments
	try {
		TCLAP::CmdLine cmd("Copula-generation alg. by Michal Kaut", ' ', "0.1");
		TCLAP::UnlabeledValueArg<int> argNScen("nscen", "number of scenarios",
																					 true, 0, "int", cmd);
		TCLAP::ValueArg<int> argRSeed ("r", "rseed", "random seed",
																	 false, time(NULL), "int", cmd);
		TCLAP::ValueArg<std::string> argTgDistFName ("i", "tgdist",
																		             "file with target distrib.",
																	               false, "target-dist.dat",
																	               "file name", cmd);
		/*
		TCLAP::ValueArg<std::string> argTgCopFName ("d", "tgcop",
																		            "file with target copula",
																	              false, "target-cop.dat",
																	              "file name", cmd); */
		TCLAP::ValueArg<std::string> argOutputFName ("o", "outfile",
																		             "output file name",
																	               false, "out_cop-gen.txt",
																	               "file name", cmd);

		// parse the arguments
		cmd.parse( argc, argv );
		nSc = argNScen.getValue();
		srand(argRSeed.getValue()); cerr << "srand = " << argRSeed.getValue() << endl;
		tgDistFName = argTgDistFName.getValue();
		outputFName = argOutputFName.getValue();

	} catch (TCLAP::ArgException &e) {
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
	}

	// read the input file
	ifstream tgDistF(tgDistFName.c_str());
	if (!tgDistF) {
		cerr << "Problem with the input file!" << endl;
		exit(1);
	}
	int nSamples, nVars;
	tgDistF >> nSamples >> nVars;
	// need std::vector for each margin
	std::vector< std::vector<double> > histDist(nVars);
	for (i = 0; i < nVars; i++) {
		histDist[i].resize(nSamples);
	}
	for (i = 0; i < nSamples; i++) {
		for (j = 0; j < nVars; j++) {
			tgDistF >> histDist[j][i];
		}
	}

	// NB: at the moment, we only have classes for 2D copulas,
	//     so the multi-variate case has to be handled manually :-(
	int nCopulas = nVars * (nVars - 1) / 2;
	CopulaSample copSc(nVars, nSc);
	//
	std::vector< Cop2DData<Vector> * > p2copData(nCopulas);
	Cop2DData<Vector> * p2cop2Ddata;
	int c = 0;
	for (i = 0; i < nVars; i++) {
		for (j = i+1; j < nVars; j++) {
			p2cop2Ddata = new Cop2DData<Vector>(&histDist[i], &histDist[j],
																					nSamples, nSc);
			assert (c < nCopulas && "sanity check");
			p2copData[c] = p2cop2Ddata;
			copSc.attach_tg_2Dcop(p2copData[c], i, j);
			c++;
		}
	}
	assert (c == nCopulas && "sanity check");

	copSc.gen_sample();
	copSc.print_as_txt(outputFName.c_str());

	return 0;
}
