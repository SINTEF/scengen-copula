#include <iostream>
#include <fstream>

#include "cop2Dsample.hpp"

using namespace std;
using namespace Copula2D;


// --------------------------------------------------------------------------
// CONSTRUCTORS AND DESTRUCTORS

Cop2DSample::Cop2DSample(int const nSamples, Cop2DInfo const *const p2TgCop)
: N(nSamples), i2jC(nSamples, -1), j2iC(nSamples, -1),
  p2prob(NULL), scenOfRow(0), //cumProb(N),
  evalPtPos(1.0), copEvalPts(N),
  p2tgInfo(p2TgCop)
{
	for (int i = 0; i < N; i++) {
		copEvalPts[i] = (i + evalPtPos) / N; // discretization points
		//cumProb[i] = (i + 1.0) / N;          // cummulative probabilities
	}
}


// --------------------------------------------------------------------------
// PUBLIC METHODS

void Cop2DSample::set_scen_of_i(IVector const &iScenVect,
																double const *p2scProb)
{
	if ((int) scenOfRow.size() != N) {
		assert ((int) scenOfRow.size() == 0
						&& "number of samples should never change");
		scenOfRow.resize(N);
	}
	int i, s;

	for (s = 0; s < N; s++) {
		scenOfRow[iScenVect[s]] = s;
	}

	p2prob = p2scProb;
	double cumPr = 0.0;
	for (i = 0; i < N; i++) {
		s = scenOfRow[i]; // scenario of row 'i'
		// remember that copEvalPts is sorted by rows, not by scenarios
		copEvalPts[i] = cumPr + evalPtPos * p2prob[s]; // point used in the cdf
		cumPr += p2prob[s];
		//cumProb[i] = cumPr;
	}
	assert (fabs(cumPr - 1.0) < DblEps && "probabilities must sum up to 1!");
}


double Cop2DSample::cdf(int const i, int const j) const
{
	return prob_in_box(0, 0, i, j);
}


void Cop2DSample::print_as_txt(string const fName, bool const sortByScen)
{
	std::ofstream oFile;
	oFile.open(fName.c_str(), std::ios::out);
	if (!oFile) {
		cerr << "WARNING: could not open output file " << fName << "!" << endl;
	} else {
		int i, s;
		for (s = 0; s < N; s++) {
			i = (sortByScen ? scenOfRow[s] : s);
			oFile << i << "\t" << i2jC[i] << endl;
		}
	}
}


// --------------------------------------------------------------------------
// PRIVATE METHODS


void Cop2DSample::update_j2iC_from_i2jC()
{
	int i, j;

	for (i = 0; i < N; i++) {
		j = i2jC[i];
		assert (j < N && "bound check: unassigned link are designated by -1");
		if(j >= 0) { // j < 0 means unassigned
			j2iC[j] = i;
		}
	}
}

void Cop2DSample::update_i2jC_from_j2iC()
{
	int i, j;

	for (j = 0; j < N; j++) {
		i = j2iC[j];
		assert (i < N && "bound check: unassigned link are designated by -1");
		if(i >= 0) { // i < 0 means unassigned
			i2jC[i] = j;
		}
	}
}


bool Cop2DSample::have_valid_cop()
{
	int i;
	for (i = 0; i < N; i++) {
		if (i2jC[i] < 0 || j2iC[i2jC[i]] != i) {
			return false;
		}
	}
	return true;
}


double Cop2DSample::prob_in_box(int const i0, int const j0,
										            int const i1, int const j1) const
{
	int i, j;
	double prob = 0.0;

	if (p2prob) {
		// we have custom probabilities
		cerr << "ERROR: non-equiprobable scenarios not implemented." << endl;
	} else {
		// no custom probabilities -> assume equiprobable
		int n = 0;
		if ((i1 - i0) <= (j1 - j0)) {
			// fewer i's than j's -> do it column-wise
			for (i = i0; i < i1; i++) {
				j = i2jC[i];
				assert (j >= 0 && j < N && "bound check; if it fails, check i2jC");
				if ((j >= j0) && (j < j1)) {
					n++;
				}
			}
		} else {
			// fewer j's than i's -> do it row-wise
			for (j = j0; j < j1; j++) {
				i = j2iC[j];
				assert (i >= 0 && i < N && "bound check; if it fails, check j2iC");
				if ((i >= i0) && (i < i1)) {
					n++;
				}
			}
		}
		prob = n / static_cast<double>(N); // alt: n * 1.0 / N
	}

	return prob;
}


/*
	If one wants a direct comparison to the MIP method, one should probably
	round the tgRCdf values - but it gives worse results (bigger KS-dist.)!!

	\todo Deal with equal values -> choose randomly between them!
*/
double Cop2DSample::gen_heur()
{
	int i, j, jj;
	double tgVal;
	double dist;
	double minDist;

	assert ((!p2prob || (int) scenOfRow.size() == N)
					&& "if we have probabilities, we must have column-scens as well");

	/// minimum cdf-improvement for a new best solution
	/// This controls how much better (in terms of cdf-distance) must a row be
	/// to be considered a new best solution. Small values give best samples,
	/// bigger values give worse solutions, but different at each run.
	double CdfDistEps = DblEps; // DblEps should give 'optimal' discretization

	Vector prevColCdf(N, 0.0);
	IVector bestRows;     ///< all rows ('j') that minimize the distance
	bestRows.reserve(10); // this should be enough to prevent reallocations(?)

	double probCol;       ///< probability of column 'i'
	double totDist = 0.0;
	for (i = 0; i < N; i++) {
		// init
		probCol = (p2prob ? p2prob[scenOfRow[i]] : 1.0 / N);
		minDist = N * cdfDist(0.0, 1.0); // we should always find something better
		bestRows.clear();

		/*
			The error/distance of putting the point at (i,j) is
			dist(i,j) = sum_k^{j-1} d(F(i-1,k), T(i,k))
			            + sum_j^{n-1} d(F(i-1,k) + 1, T(i,k))
			where: d() is the distance measure  = cdfDist()
			       F() is the sample (rank) Cdf -> prevColRCdf[] is F(i-1,*)
			       T() is the target (rank) Cdf = tgRankCdf[][]
			We use the fact that for j>0 we have the following:
			dist(i,j) = dist(i,j-1) + d(F(i-1,j-1), T(i,j-1))
			                        - d(F(i-1,j-1) + 1, T(i,j-1))
			It means that we have to treat j=0 separately to get dist(i,0)
		*/
		j = 0;
		// compute dist(i, 0):
		dist = 0;
		for (jj = 0; jj < N; jj++) {
			// RCdf(i,jj) = prevColRCdf(jj) + 1 for all jj - we have point (i, 0)
			dist += cdfDist(prevColCdf[jj] + evalPtPos * probCol, tgCdf(i, jj));
		}
		if (j2iC[j] < 0) {
			// j=0 is available, i.e. it is a valid choice
			bestRows.push_back(j);
			minDist = dist;
		}

		for (j = 1; j < N; j++) {
			// using the recursive formula shown above
			tgVal = tgCdf(i, j-1);
			dist += cdfDist(prevColCdf[j-1], tgVal)
			        - cdfDist(prevColCdf[j-1] + evalPtPos * probCol, tgVal);
			if (j2iC[j] < 0) {
				// j is an available row
				if (dist < minDist - CdfDistEps) {
					// j is a new best row -> delete current candidates and store j
					bestRows.clear();
					bestRows.push_back(j);
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
	}
	assert (have_valid_cop() && "should have a valid copula sample now");

	#ifndef NDEBUG
		cout << endl << "Cop2DSample::gen_heur() finished with distance/error = "
				 << totDist << endl;
	#endif

	return totDist;
}


/// \todo make this into a private method and call it also from the
///       \a gen_heur() method !!!
void Cop2DSample::cdf_dist_of_col(int const i, Vector const prevColCdf,
											            Vector rowCdfDist, bool rowFree[])
{
	assert (rowFree == NULL && "handling of rowFree is not yet implented...");

	double probCol = (p2prob ? p2prob[scenOfRow[i]] : 1.0 / N);

	int j = 0;
	// compute dist(i, 0):
	double dist = 0;
	for (jj = 0; jj < N; jj++) {
		// RCdf(i,jj) = prevColRCdf(jj) + 1 for all jj - we have point (i, 0)
		dist += cdfDist(prevColCdf[jj] + evalPtPos * probCol, tgCdf(i, jj));
	}
	rowCdfDist[0] = dist;

	for (j = 1; j < N; j++) {
		// using the recursive formula shown above
		tgVal = tgCdf(i, j-1);
		dist += cdfDist(prevColCdf[j-1], tgVal)
						- cdfDist(prevColCdf[j-1] + evalPtPos * probCol, tgVal);
		rowCdfDist[j] = dist;
	}
}
