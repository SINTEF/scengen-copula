// main routine that generates scenarios for a copula given in a file

#include <iostream>
#include <fstream>
#include <ctime> // needed on gcc-win
#include <tclap/CmdLine.h> // processing of command-line arguments

#include "copula-info.hpp"
#include "cop2Dsample.hpp"
#include "copula-sample.hpp"

using namespace std;
using namespace Copula2D;
using namespace CopulaScen;


int main(int argc, char *argv[]) {

	int i, j;

	/*{
		// testing the array view/slicing
		typedef boost::multi_array_types::index_range IRange;

		boost::multi_array<double, 2> A(boost::extents[2][3]);

		boost::multi_array<double, 1>::array_view<1>::type Acol1
			= A[ boost::indices[IRange()][1] ];
		boost::multi_array<double, 1> Acol2
			= A[ boost::indices[IRange()][1] ];
		boost::multi_array<double, 1>::array_view<1>::type Acol3
			= Acol2[ boost::indices[IRange()] ];

		cout << "A[0,1] at addr " << &A[0][1]
				 << "; Acol1[0] at addr. " << &Acol1[0] // same as A[0,1]
				 << "; Acol2[0] at addr. " << &Acol2[0] // different!
				 << "; Acol3[0] at addr. " << &Acol3[0] << endl; // same as Acol2!
	}*/

	// variables whose values is read from the command line
	int nSc = 0;             // number of scenarios to generate
	std::string tgDistFName; // input file name
	std::string outputFName; // output file name
	int numCandPts = 1;      // minimal number of candidate scenarios
	bool writeProbAllocData; // should we write data for the prob-alloc model?

	// parameters for processing of command-line arguments
	// value-arguments are shown in the reversed order!
	try {
		TCLAP::CmdLine cmd("Copula-generation alg. by Michal Kaut", ' ', "0.1");
		TCLAP::UnlabeledValueArg<int> argNScen("nscen", "number of scenarios",
		                                       true, 0, "int", cmd);
		TCLAP::ValueArg<int> argRSeed ("r", "rseed", "random seed",
		                               false, time(NULL), "int", cmd);
		TCLAP::ValueArg<std::string> argOutputFName ("o", "outfile",
		                                             "output file name",
		                                             false, "out_cop-gen.txt",
		                                             "file name", cmd);
		TCLAP::ValueArg<std::string> argTgDistFName ("i", "tgdist",
		                                             "file with target distrib.",
		                                             false, "target-dist.dat",
		                                             "file name", cmd);
		TCLAP::ValueArg<int> argNumCandPts ("c", "numcand",
		                                    "min number of candidate scenarios",
		                                    false, 1, "integer >= 1", cmd);
		TCLAP::SwitchArg argWriteProbAllocData ("", "proballoc",
		                                        "write data for prob-alloc model?",
		                                        cmd, false);

		// parse the arguments
		cmd.parse( argc, argv );
		nSc = argNScen.getValue();
		srand(argRSeed.getValue());
		tgDistFName = argTgDistFName.getValue();
		outputFName = argOutputFName.getValue();
		numCandPts = argNumCandPts.getValue();
		writeProbAllocData = argWriteProbAllocData.getValue();

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

	TMatrixI histRanks(nVars, nSamples);
	get_ranks_or_rows(histDist, histRanks);


	CopulaDef::CopInfoData tgCopInfo(histRanks);

	// NB: at the moment, we only have classes for 2D copulas,
	//     so the multi-variate case has to be handled manually :-(
	int nCopulas = nVars * (nVars - 1) / 2;
	CopulaSample copSc(nVars, nSc, numCandPts);
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
	copSc.print_as_txt(outputFName.c_str(), true);

	if (writeProbAllocData) {
		// write the AMPL/GMP data file
		copSc.attach_tg_cop_info(&tgCopInfo);
		copSc.write_gmp_data();
	}

	/*
	TVectorD uV(nVars);
	cout << endl << "the target cdf at the scenario points:" << endl;
	for (int s = 0; s < nSc; ++s) {
		for (i = 0; i < nVars; ++i) {
			uV[i] = rank2U01(copSc.tmp_get_res(i, s), nSc);
		}
		cout << "scen " << s << ": " << tgCopInfo.cdf(uV) << endl;
	}
	*/

	return 0;
}
