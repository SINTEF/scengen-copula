#include <iostream>
#include <fstream>
#include <ctime> // needed by gcc-win

#include <ql/math/distributions/normaldistribution.hpp>
#include <boost/lexical_cast.hpp>

#include "cop2Dsample.hpp"
#include "copula-sample.hpp"

using namespace std;
using namespace Copula2D;
using namespace CopulaScen;

int main(int argc, char *argv[]) {
	// set the random seed
	srand ( time(NULL) );
	//srand(0); // for debugging

	int i, j, k;

	const int defaultNmbSamples = 10;
	int nS = (argc > 1 ? atoi(argv[1]) : defaultNmbSamples);

	ifstream tgCorrF("tg_corrs.txt");
	if (!tgCorrF) {
		cerr << "Problem opening the input file 'tg_corrs.txt'!" << endl;
		exit(1);
	}
	int nVars;
	tgCorrF >> nVars >> i;

	TMatrixD tgCorrs(nVars, nVars);
	for (i = 0; i < nVars; ++i) {
		for (j = 0; j < nVars; ++j) {
			tgCorrF >> tgCorrs[i][j];
			if (j == i) {
				if (!isEq(tgCorrs[i][j], 1.0)) {
					cerr << "Error in the target correlations!" << endl;
					exit(1);
				}
			}
		}
	}

	int nmb2Dcop = (nVars * (nVars - 1)) / 2;
	std::vector< Cop2DNormal<Vector> * > tg2Dcopulas(nmb2Dcop);

	// minimal length of the candidate list - passed as a second arg.
	int minNmbCand = (argc <= 2 ? 1 : boost::lexical_cast<int>(argv[2]));
	//
	CopulaSample copScens(nVars, nS, minNmbCand);

	k = 0;
	for (i = 0; i < nVars; ++i) {
		for (j = i+1; j < nVars; ++j) {
			tg2Dcopulas[k] = new Cop2DNormal<Vector> (tgCorrs[i][j], nS);
			copScens.attach_tg_2Dcop(tg2Dcopulas[k], i, j);
		}
	}

	copScens.gen_sample();
	copScens.print_as_txt("out_normal-cop.txt", true);


	// -----------------------------------------------------------------------
	// finished the copula .. now get the real margins
	ifstream tgMargF("tg_margins.txt");
	if (tgMargF) {
		tgMargF >> i;
		if (i != nVars) {
			cerr << "Error: 'tg_margins.txt' disagrees with 'tg_corrs.txt'!"
			     << endl;
		} else {
			TMatrixD normScens(nS, nVars);
			TMatrixD normScFixMV(nS, nVars); // scenarios with fixed mean-var
			string margT;
			double tgMean, tgStD, tgSkew, tgKurt;
			ofstream normScenF("out_normal-scens.txt");
			if (normScenF) {
				normScenF << nS << endl << nVars << endl;
				for (i = 0; i < nVars; ++i) {
					tgMargF >> margT >> tgMean >> tgStD;
					if (margT[0] == 'm') {
						tgMargF >> tgSkew >> tgKurt;
						if (!isEq(tgSkew, 0) || !isEq(tgKurt, 3)) {
							cerr << "ERROR: the distribution in 'tg_margins.txt' is not normal!"
								  << endl;
							exit(1);
						}
					}
					// transform the copula to the correct margin
					// also compute the actual mean and std. dev.
					double mean = 0.0;
					double stD = 0.0;
					QuantLib::InverseCumulativeNormal normCdf(tgMean, tgStD);
					for (j = 0; j < nS; ++j) {
						normScens[j][i] = normCdf((copScens.tmp_get_res(i, j) + 0.5)
						                          / static_cast<double>(nS));
						mean += normScens[j][i];
						stD += pow(normScens[j][i], 2);
					}
					mean /= nS;
					stD = sqrt(stD / nS - pow(mean, 2));
					// now we can compute the mean-var-fixed values
					double fixMult = tgStD / stD;
					for (j = 0; j < nS; ++j) {
						normScFixMV[j][i] = fixMult * (normScens[j][i] - mean)
						                    + tgMean;
					}
				}
				normScenF << normScens;
				normScenF.close();
			}
			ofstream normScenFixedF("out_normal-scens_mv-fixed.txt");
			if (normScenF) {
				normScenFixedF << nS << endl << nVars << endl << normScFixMV;
				normScenFixedF.close();
			}
		}
	}


	return 0;
}
