#include <iostream>
#include <fstream>

#include "copula-sample.hpp"

using namespace std;
using namespace CopulaScen;


// --------------------------------------------------------------------------
// CONSTRUCTORS AND DESTRUCTORS
CopulaSample::CopulaSample(int const dim, int const S,
                           double const numCandPtsRel)
: haveSc4Marg(dim, false),
  p2copInfo(NULL),
  p2tgInfo(boost::extents[dim][dim]), p2sample2D(boost::extents[dim][dim]),
  p2prob(NULL), sample(dim),
  minNumCandPtsRel(minNumCandPtsRel),
  nVar(dim), nSc(S)
{
	int i, j;
	for (i = 0; i < nVar; i++) {
		for (j = 0; j < nVar; j++) {
			p2tgInfo[i][j] = NULL;
			p2sample2D[i][j] = NULL;
		}
		sample[i].resize(nSc);
	}
}

double CopulaSample::gen_new_margin(int const marg)
{
	int i, j, s;
	int tg, iR;

	assert (haveSc4Marg[marg] == false && "the margin has not been generated");

	/// vector of target specifications for all the involved 2D copulas
	std::vector<Cop2DInfo const*> tg2Dcopulas;
	tg2Dcopulas.reserve(nVar);

	/// list of already existing margins that are being matched
	std::vector<int> oldMargins;
	oldMargins.reserve(nVar);

	double CdfDistEps = DblEps; // DblEps should give 'optimal' discretization

	// find the copulas to be matched
	// for the moment, assume that the new variable is the second index!!!
	for (i = 0; i < nVar; i++) {
		if (i != marg && haveSc4Marg[i] && p2tgInfo[i][marg]) {
			// store the margin in the list of copulas (i, marg) to match
			oldMargins.push_back(i);
			tg2Dcopulas.push_back(p2tgInfo[i][marg]);
			// create a new copula-2D-sample object
			stringstream cop2DId; // using stringstream to get simple conversions
			cop2DId << "sample_" << i << "_" << marg;
			p2sample2D[i][marg] = new Cop2DSample(nSc, p2tgInfo[i][marg],
			                                      cop2DId.str());
			#ifndef NDEBUG
				cout << "Created new Cop2DSample with id = " << cop2DId.str() << endl;
			#endif
			// initialize the 2D-sample generator with scenarios of the known margin
			p2sample2D[i][marg]->set_scen_of_marg_1(sample[i], p2prob);
		}
	}
	int nTg2Dcops = tg2Dcopulas.size();

	double dist = 0.0; // init not really needed - just to avoid a gcc warning..
	double minDist;
	double totDist = 0.0;
	double scProb;

	/// \todo rename bestRows to bestScen
	/// \todo fix the numbering: rows OR scens !!!!

	/*
		The logic in the following is as follows:
		- loop through the ranks of the new margin (sequentially)
		- for each rank, loop through all scenarios and compute the distance
		  (error) of putting the rank to the given scenario,
		  based on the values of already assigned margins
		- store the best scenarios
		- after the loop, chose one of the candidates
	*/

	unsigned minNumCandScens
		= static_cast<unsigned>(ceil(minNumCandPtsRel * nSc)) + 1;
	Copula2D::Cop2DSample::CandList candScens(minNumCandScens, CdfDistEps);
	//IVector bestScens;     ///< list of scenarios that minimize the distance
	//bestScens.reserve(10); // this should be enough to prevent reallocations(?)
	std::vector<bool> scenUsed(nSc, false);

	std::vector<Vector> prevRowCdf(nTg2Dcops); // cdf of the previous row
	std::vector<Vector> colCdfDist(nTg2Dcops); // dist. of using the columns
	for (tg = 0; tg < nTg2Dcops; tg++) {
		prevRowCdf[tg] = Vector(nSc, 0);
		colCdfDist[tg] = Vector(nSc);
	}

	// compute a safe upper bound for the distance:
	double maxMinDist = 0.0;
	for (tg = 0; tg < nTg2Dcops; tg++) {
		maxMinDist += p2sample2D[oldMargins[tg]][marg]->cdfDist(0.0, 1.0);
	}
	maxMinDist *= nSc;

	// loop over ranks of the newly added margin
	for (iR = 0; iR < nSc; iR++) {
		// init minDist & clean the bestRows array
		minDist = maxMinDist;
		//bestScens.clear();
		candScens.clear();

		for (tg = 0; tg < nTg2Dcops; tg++) {
			// get the dist for the rows
			/// \todo check this!
			p2sample2D[oldMargins[tg]][marg]->cdf_dist_of_row(iR, prevRowCdf[tg],
			                                                  colCdfDist[tg]);
		}

		// brute-force approach - this will be SLOW
		for (s = 0; s < nSc; s++) {
			// compute the distance for putting rank iR into scenario s
			if (!scenUsed[s]) {
				dist = 0.0;
				for (tg = 0; tg < nTg2Dcops; tg++) {
					int rankInTgMargin = sample[oldMargins[tg]][s]; // CHECK!
					dist += colCdfDist[tg][rankInTgMargin];
				}

				candScens.insert_cand(s, dist);
/*
				if (dist < minDist - CdfDistEps) {
					// j is a new best row -> delete current candidates and store j
					bestScens.clear();
					bestScens.push_back(s);
					minDist = dist;
				} else {
					if (dist < minDist + CdfDistEps)  {
						// j is just as good as the best known value -> add to list
						bestScens.push_back(s);
						// We could update the minDist here, but that could cause the
						// list to 'slide' downwards, so we end-up with bad solutions
						// On the other hand, if we do not update, then we on the next
						// change of minDist remove also points that are in the range..
						// The best would probably be to keep all the distances as well,
						// keep updating minDist and then clean the list at the end.
						//minDist = min(dist, minDist); // causes problems!
					}
				}
*/
			}
		}
		// now we should have a couple..
		assert (candScens.get_num_cand() > 0
						&& "we should have found something...");
		candScens.get_rand_cand(s, minDist);
/*		assert (bestScens.size() > 0 && "we should have found something...");

		/// \todo Check!
		if (bestScens.size() == 1) {
			// only have one candidate row -> easy
			s = bestScens[0];
		} else {
			// more candidates -> choose randomly, using my macro for random int
			s = bestScens[irand(bestScens.size())];
		}
*/
		#ifndef NDEBUG
			cout << "CopulaSample::gen_new_margin(" << marg
			     << "): new link: iR=" << iR << ", s=" << s
			     << "; dist = " << minDist << " (min dist = "
			     << candScens.get_best_cand_value() << ")" << endl << endl;
			cout.flush();
		#endif
		scProb = (p2prob == NULL ? 1.0 / nSc : p2prob[s]);
		totDist += scProb * minDist;
		sample[marg][s] = iR;
		scenUsed[s] = true;

		/// \todo Check!
		for (tg = 0; tg < nTg2Dcops; tg++) {
			i = oldMargins[tg];
			j = sample[i][s]; // the row-rank for the given target copula

			// update prevRowCdf:
			for (int jj = j; jj < nSc; jj++) {
				prevRowCdf[tg][jj] += scProb;
			}

			// update p2sample2D
			bool linkAdded = p2sample2D[i][marg]->add_link(j, iR);
			assert (linkAdded && "This should always succeed...");

			// update also the [marg][i] copula sample? !!!
		}
	}
	haveSc4Marg[marg] = true;

	return dist;
}


/// \todo Write something here!!!
void CopulaSample::attach_tg_2Dcop(Cop2DInfo const* p2cop,
                                   int const i, int const j,
                                   bool makeTranspTg)
{
	p2tgInfo[i][j] = p2cop;
	if (makeTranspTg) {
		p2tgInfo[j][i] = new Copula2D::Cop2DInfTr(p2cop);
	}
}



/// the main routine; returns the KS-distance
/// \todo Do something here!
double CopulaSample::gen_sample()
{
	int marg;
	int s;
	double totDist = 0.0;

	// initialize the first margin
	marg = 0;
	for (s = 0; s < nSc; s++) {
		/// \todo Once this works, test it with $nSc - s$ or random numbers
		sample[marg][s] = s; // alternative is to use random numbers
	}
	haveSc4Marg[marg] = true;

	for (marg = 1; marg < nVar; marg++) {
		totDist += gen_new_margin(marg);
	}

	return totDist;
}


void CopulaSample::print_as_txt(string const fName, bool const scaleTo01,
                                int const sortByMarg)
{
	std::ofstream oFile;
	oFile.open(fName.c_str(), std::ios::out);
	if (!oFile) {
		cerr << "WARNING: could not open output file " << fName << "!" << endl;
	} else {
		int marg, s;
		if (sortByMarg >= 0) {
			cerr << "Sorting of output by margins si not implemented yet!" << endl;
		}
		if (scaleTo01) {
			for (s = 0; s < nSc; s++) {
				for (marg = 0; marg < nVar - 1; marg++) {
					oFile << (static_cast<double>(sample[marg][s]) + 0.5) / nSc << "\t";
				}
				oFile << (static_cast<double>(sample[nVar - 1][s]) + 0.5) / nSc << endl;
			}
		} else {
			for (s = 0; s < nSc; s++) {
				for (marg = 0; marg < nVar - 1; marg++) {
					oFile << sample[marg][s] << "\t";
				}
				oFile << sample[nVar - 1][s] << endl;
			}
		}
		oFile.close();
	}

	/*
	// TMP !!!
	#ifndef NDEBUG
	Vector iScOfR(nSc);
	Vector jScOfR(nSc);
	Vector tgScPr(nSc);
	for (int j = 1; j < nVar; ++j) {
		for (int i = 0; i < j; ++i) {
			if (j != i && haveSc4Marg[i] && p2sample2D[i][j]) {
				assert (p2tgInfo[i][j]);
				cout << "summary for 2D-copula (" << i << "," << j << ")" << endl;

				for (int s = 0; s < nSc; ++s) {
					iScOfR[sample[i][s]] = s;
					jScOfR[sample[j][s]] = s;
					tgScPr[s] = 0.0; // init
				}
				double sumTgPr = 0.0;
				for (int iR = 0; iR < nSc; ++iR) {
					int s = iScOfR[iR];
					int jR = sample[j][s];
					assert (jScOfR[jR] == s && "sanity check");
					double u = (iR + 0.5) / (double) nSc;
					double v = (jR + 0.5) / (double) nSc;
					double tgCdf = p2tgInfo[i][j]->cdf(u, v);
					double tgProb = tgCdf;
					// now have to subtract all probs included in the cdf:
					for (int iR2 = 0; iR2 < iR; ++iR2) {
						int s2 = iScOfR[iR2];
						int jR2 = sample[j][s2];
						if (jR2 < jR) {
							//assert (tgScPr[s2] > 0.0 && "should have been assigned");
							tgProb -= tgScPr[s2];
							//assert (tgProb > 0.0 && "sanity check");
						}
					}
					sumTgPr += tgProb;
					tgScPr[s] = tgProb;
				}
				cout << "sum of prob. = " << sumTgPr << endl;

				for (int s = 0; s < nSc; ++s) {
					double scCdf = p2sample2D[i][j]->cdfOfR(sample[i][s], sample[j][s]);
					//double u = (sample[i][s] + 0.5) / nSc;
					//double v = (sample[j][s] + 0.5) / nSc;
					double u = rank2U01(sample[i][s], nSc);
					double v = rank2U01(sample[j][s], nSc);
					double u0 = (sample[i][s] == 0 ? 0.0 : (sample[i][s] - 0.5) / nSc);
					double v0 = (sample[j][s] == 0 ? 0.0 : (sample[j][s] - 0.5) / nSc);
					double tgCdf = p2tgInfo[i][j]->cdf(u, v);
					double scProb = (p2prob ? p2prob[s] : 1.0 / (double) nSc);
					cout << "\tscen " << s << ": (" << sample[i][s] << "," << sample[j][s]
					     << ") - sc-cdf = " << scCdf << ", tg-cdf = " << tgCdf
					     << "; sc-prob = " << scProb << ", tg-prob = " << tgScPr[s]
					     << " ?=? " << sqrt(p2tgInfo[i][j]->cdf(u, v) - p2tgInfo[i][j]->cdf(u0, v)
					        - p2tgInfo[i][j]->cdf(u, v0) + p2tgInfo[i][j]->cdf(u0, v0)) << endl;
				}
			}
		}
	}
	#endif
*/

}


void CopulaSample::write_gmp_data(string const fName)
{
	int i, s;
	std::ofstream oFile;
	oFile.open(fName.c_str(), std::ios::out);
	if (!oFile) {
		cerr << "WARNING: could not open output file " << fName << "!" << endl;
	} else {
		oFile << "# this file was written by the copula-generation code" << endl
		      << endl
		      << "param nVars := " << nVar << ";" << endl
		      << "param nScens := " << nSc << ";" << endl
		      << endl
		      << "param wAvgErr := " << 1.0 << ";" << endl
		      << "param wMaxErr := " << 1.0 << ";" << endl
		      << "param wAvgDev := " << 1.0 << ";" << endl
		      << endl
		      << "param rank (tr)" << endl << "\t:";
		for (i = 0; i < nVar; ++i) {
			oFile << "\t" << i;
		}
		oFile << " :=";
		for (s = 0; s < nSc; ++s) {
			oFile << endl << "\t" << s;
			for (i = 0; i < nVar; ++i) {
				oFile << "\t" << sample[i][s];
			}
		}
		oFile << endl << ";" << endl
		      << endl;

		double tgScProb;
		TVectorD uScen(nVar);
		oFile << "param tgCdf :=";
		for (s = 0; s < nSc; ++s) {
			for (i = 0; i < nVar; ++i) {
				uScen[i] = rank2U01(sample[i][s], nSc);
			}
			if (p2copInfo) {
				// have the multivar target info
				tgScProb = p2copInfo->cdf(uScen);
			} else {
				// no multivar target -> take the average of all the 2D copulas!
				tgScProb = 0.0;
				int nTgCops = 0;
				for (i = 0; i < nVar; ++i) {
					for (int j = 0; j < nVar; ++j) {
						if (j != i && p2tgInfo[i][j]) {
							tgScProb += p2tgInfo[i][j]->cdf(uScen[0], uScen[1]);
							nTgCops ++;
						}
					}
				}
				tgScProb /= nTgCops;
			}
			oFile << endl << "\t" << s << "\t" << tgScProb;
		}
		oFile << endl << ";" << endl;

		oFile << endl << "end;" << endl;
		oFile.close();
	}
}
