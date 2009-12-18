#include <iostream>
#include <fstream>

#include "copula-sample.hpp"

using namespace std;
using namespace Copula2D;


// --------------------------------------------------------------------------
// CONSTRUCTORS AND DESTRUCTORS
CopulaSample::CopulaSample(int const dim)
: haveSc4Marg(dim, false), nVar(dim), p2tgInfo(boost::extents[dim][dim])
{
	int i, j;
	for (i = 0; i < nVar; i++) {
		for (j = 0; j < nVar; j++) {
			p2tgInfo[i][j] = NULL;
		}
	}
}

double CopulaSample::gen_new_margin(int const marg)
{
	int i, j, s;

	std::vector<Cop2DInfo const*> tg2Dcopulas;
	tg2Dcopulas.reserve(nVar);

	double CdfDistEps = DblEps; // DblEps should give 'optimal' discretization

	// find the copulas to be matched
	// for the moment, assume that the new variable is the second index!!!
	for (i = 0; i < nVar; i++) {
		if (i != marg && haveSc4Marg[i] && p2tgInfo[i][marg]) {
			// store the margin in the list of copulas (i, marg) to match
			tg2Dcopulas.push_back(i);
			// initialize the 2D-sample generator with scenarios of the known margin
			p2sample2D->set_scen_of_i(sample[i], p2prob);
		}
	}
	int nTg2Dcops = tg2Dcopulas.size();

	double dist;
	double minDist;

	/// \todo rename bestRows to bestScen
	/// \todo fix the numbering: rows OR scens !!!!

	IVector bestRows;     ///< all rows ('j') that minimize the distance
	bestRows.reserve(10); // this should be enough to prevent reallocations(?)
	std::vector<bool> scenUsed(nSc, false);

	std::vector<Vector> prevColCdf(nTg2Dcops);
	std::vector<Vector> rowCdfDist(nTg2Dcops);
	for (tg = 0; tg < nTg2Dcops; tg++) {
		prevColCdf[tg] = Vector(nSc, 0);
		rowCdfDist[tg] = Vector(nSc);
	}

	for (iR = 0; iR < nSc; iR++) {
		minDist = N * cdfDist(0.0, 1.0); // we should always find something better
		bestRows.clear();

		for (tg = 0; tg < nTg2Dcops; tg++) {
			// get the dist for the rows
			tg2Dcopulas->cdf_dist_of_col(i, prevColCdf[tg], rowCdfDist[tg]);
		}

		// brute-force approach - this will be SLOW
		for (s = 0; s < nSc; s++) {
			// compute the distance for putting rank iR into scenario s
			if (!scenUsed) {
				dist = 0;
				for (tg = 0; tg < nTg2Dcops; tg++) {
					dist += rowCdfDist[tg][s];
				}

				if (dist < minDist - CdfDistEps) {
					// j is a new best row -> delete current candidates and store j
					bestRows.clear();
					bestRows.push_back(s);
					minDist = dist;
				} else {
					if (dist < minDist + CdfDistEps)  {
						// j is just as good as the best known value -> add to list
						bestRows.push_back(j);
						// We could update the minDist here, but that could cause the
						// list to 'slide' downwards, so we end-up with bad solutions
						// On the other hand, if we do not update, then we on the next
						// change of minDist remove also points that are in the range..
						// The best would probably be to keep all the distances as well,
						// keep updating minDist and then clean the list at the end.
						//minDist = min(dist, minDist); // causes problems!
					}
				}
			}
		}
		// now we should have a couple..
		assert (bestRows.size() > 0 && "we should have found something...");

		/// \todo adapt this for the CopulaSample class
		/*
		if (bestRows.size() == 1) {
			// only have one candidate row -> easy
			j = bestRows[0];
		} else {
			// more candidates -> choose randomly, using my macro for random int
			j = bestRows[irand(bestRows.size())];
		}
		totDist += probCol * minDist;
		i2jC[i] = j;
		j2iC[j] = i;

		// update prevColCdf:
		for (jj = j; jj < N; jj++) {
			prevColCdf[jj] += probCol;
		}
		*/
	}
}
