#include <iostream>
#include <fstream>

#include "cop2Dsample.hpp"

using namespace std;
using namespace Copula2D;


// --------------------------------------------------------------------------
// CONSTRUCTORS AND DESTRUCTORS

Cop2DSample::Cop2DSample(int const nSamples, Cop2DInfo * const p2TgCop,
                         string const id)
: N(nSamples), i2jC(nSamples, -1), j2iC(nSamples, -1),
  tgCdfOfR(p2TgCop->get_cdf_grid()),
  p2prob(NULL), scenOfMarg1R(0), //cumProb(N),
  evalPtPos(0.5), copEvalPts(N),
  sampleId(id),
  p2tgInfo(p2TgCop)
{
	/// \todo do we need this here? it is repeated at set_scen_of_marg_1()!
	for (int i = 0; i < N; i++) {
		copEvalPts[i] = (i + evalPtPos) / N; // discretization points
		//cumProb[i] = (i + 1.0) / N;          // cummulative probabilities
	}
	// avoiding numerical issues:
	if (evalPtPos >= 1.0 - DblEps) {
		copEvalPts[N-1] = 1.0;
	}

	// initialize the cdf grid of the target info object
	p2TgCop->init_cdf_grid(nSamples, evalPtPos);

	// fill the matrix of target rank-cdf values
	//fill_tgCdfOfR();
}


// --------------------------------------------------------------------------
// PUBLIC METHODS

/*
void Cop2DSample::set_scen_of_i(IVector const &scenOfColR,
																double const *p2scProb)
{
	if ((int) scenOfMarg1R.size() != N) {
		assert ((int) scenOfMarg1R.size() == 0
						&& "number of samples should never change");
		scenOfMarg1R.resize(N);
	}
	int i, s;
	double scPr;

	for (i = 0; i < N; i++) {
		scenOfMarg1R[i] = scenOfColR[i];
	}

	p2prob = p2scProb;
	double cumPr = 0.0;
	for (i = 0; i < N; i++) {
		s = scenOfMarg1R[i]; // scenario of col 'i'
		// remember that copEvalPts is sorted by ranks, not by scenarios
		scPr = (p2prob ? p2prob[s] : 1.0 / N);
		copEvalPts[i] = cumPr + evalPtPos * scPr; // point used in the cdf
		cumPr += scPr;
		//cumProb[i] = cumPr;
	}
	assert (fabs(cumPr - 1.0) < DblEps && "probabilities must sum up to 1!");
}
*/

void Cop2DSample::set_scen_of_marg_1(IVector const &margScen,
                                     double const *p2scProb)
{
	if ((int) scenOfMarg1R.size() != N) {
		assert ((int) scenOfMarg1R.size() == 0
		        && "number of samples should never change");
		scenOfMarg1R.resize(N);
	}
	int i, s;
	double scPr;

	#ifndef NDEBUG
		// extra checks - see the assert below
		for (i = 0; i < N; i++) {
			scenOfMarg1R[i] = -1;
		}
	#endif

	for (s = 0; s < N; s++) {
		i = margScen[s]; // scenario s includes rank i
		assert (scenOfMarg1R[i] < 0 && "each rank should be set only once!");
		scenOfMarg1R[i] = s;
	}

	p2prob = p2scProb;
	double cumPr = 0.0;
	for (i = 0; i < N; i++) {
		s = scenOfMarg1R[i]; // scenario of col 'i'
		scPr = (p2prob ? p2prob[s] : 1.0 / N);
		// remember that copEvalPts is sorted by ranks, not by scenarios
		copEvalPts[i] = cumPr + evalPtPos * scPr; // point used in the cdf
		cumPr += scPr;
		//cumProb[i] = cumPr;
	}
	assert (fabs(cumPr - 1.0) < DblEps && "probabilities must sum up to 1!");
}



double Cop2DSample::cdfOfR(int const i, int const j) const
{
	return grid_prob_below(i, j);
}


void Cop2DSample::print_as_txt(string const fName, bool const scaleTo01,
                               bool const sortByScen)
{
	std::ofstream oFile;
	oFile.open(fName.c_str(), std::ios::out);
	if (!oFile) {
		cerr << "WARNING: could not open output file " << fName << "!" << endl;
	} else {
		int i, r1;
		for (r1 = 0; r1 < N; r1++) {
			i = (sortByScen ? scenOfMarg1R[r1] : r1);
			if (scaleTo01) {
				oFile << i << "\t" << rank2U01(i2jC[i]) << endl;
			} else {
				oFile << i << "\t" << i2jC[i] << endl;
			}
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


/*
void Cop2DSample::fill_tgCdfOfR()
{
	int i, j;
	for (i = 0; i < N; i++) {
		double u = copEvalPts[i];
		for (j = 0; j < N; j++) {
			tgCdfOfR[i][j] = p2tgInfo->cdf(u, copEvalPts[j]);
		}
	}
}
*/


int Cop2DSample::grid_pts_below (int const iR, int const jR) const
{
	assert (!p2prob && "does not make sense for non-equiprobable case!");
	int i, j;

	int n = 0;
	if (iR <= jR) {
		// fewer i's than j's -> do it column-wise
		for (i = 0; i <= iR; ++i) {
			j = i2jC[i];
			if (j < 0 || j >= N) {
				cout << " (i2jC[" << i << "] = " << j << ") ";
			}
			assert (j >= 0 && j < N && "bound check; if it fails, check i2jC");
			if (j <= jR) {
				//cout << " +pt_c(" << i << ", " << j << ") ";
				n++;
			}
		}
	} else {
		// fewer j's than i's -> do it row-wise
		for (j = 0; j <= jR; ++j) {
			i = j2iC[j];
			assert (i >= 0 && i < N && "bound check; if it fails, check j2iC");
			if (i <= iR) {
				//cout << " +pt_r(" << i << ", " << j << ") ";
				n++;
			}
		}
	}
	return n;
}


double Cop2DSample::grid_prob_below (int const iR, int const jR) const
{
	// in the equiprobable case, it should be faster to count the points
	// below and divide by the total number of points
	if (!p2prob) {
		return static_cast<double>(grid_pts_below(iR, jR)) / (double) N;
	}

	int i, j;

	double p = 0.0;
	if (iR <= jR) {
		// fewer i's than j's -> do it column-wise
		for (i = 0; i <= iR; ++i) {
			j = i2jC[i];
			assert (j >= 0 && j < N && "bound check; if it fails, check i2jC");
			if (j <= jR) {
				p += p2prob[scenOfMarg1R[i]];
			}
		}
	} else {
		// fewer j's than i's -> do it row-wise
		for (j = 0; j <= jR; ++j) {
			i = j2iC[j];
			assert (i >= 0 && i < N && "bound check; if it fails, check j2iC");
			if (i <= iR) {
				p += p2prob[scenOfMarg1R[i]];
			}
		}
	}
	return p;
}


double Cop2DSample::grid_prob_box(int const i0, int const j0,
                                  int const i1, int const j1) const
{
	int i, j;
	double prob = 0.0;
	const double eqPrProb = 1.0 / (double) N;

	if ((i1 - i0) <= (j1 - j0)) {
		// fewer i's than j's -> do it column-wise
		for (i = i0; i <= i1; ++i) {
			j = i2jC[i];
			assert (j >= 0 && j < N && "bound check; if it fails, check i2jC");
			if ((j >= j0) && (j <= j1)) {
				prob += (p2prob ? p2prob[scenOfMarg1R[i]] : eqPrProb);
			}
		}
	} else {
		// fewer j's than i's -> do it row-wise
		for (j = j0; j <= j1; j++) {
			i = j2iC[j];
			assert (i >= 0 && i < N && "bound check; if it fails, check j2iC");
			if ((i >= i0) && (i <= i1)) {
				prob += (p2prob ? p2prob[scenOfMarg1R[i]] : eqPrProb);
			}
		}
	}

	return prob;
}


/*
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
*/

/*
	If one wants a direct comparison to the MIP method, one should probably
	round the tgRCdfOfR values - but it gives worse results (bigger KS-dist.)!!

	\todo Deal with equal values -> choose randomly between them!
*/
double Cop2DSample::gen_heur()
{
	int i, j, jj;
	double tgVal;
	double dist;
	double minDist;

	assert ((!p2prob || (int) scenOfMarg1R.size() == N)
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
		probCol = (p2prob ? p2prob[scenOfMarg1R[i]] : 1.0 / N);
		minDist = N * cdfDist(0.0, 1.0); // we should always find something better
		bestRows.clear();

		/*
			The error/distance of putting the point at (i,j) is
			dist(i,j) = sum_{k=0}^{j-1} d(F(i-1,k), T(i,k))
			            + sum_{k=j}^{n-1} d(F(i-1,k) + 1/n, T(i,k))
			where: d() is the distance measure  = cdfDist()
			     : F() is the sample (rank) Cdf -> prevColRCdf[] is F(i-1,*)
			     : T() is the target (rank) Cdf = tgRankCdf[][]
			     : all scenarios have prob. 1/n
			We use the fact that for j>0 we have the following:
			dist(i,j) = dist(i,j-1) + d(F(i-1,j-1), T(i,j-1))
			                        - d(F(i-1,j-1) + 1/n, T(i,j-1))
			It means that we have to treat j=0 separately to get dist(i,0)
		*/
		j = 0;
		// compute dist(i, 0):
		dist = 0;
		for (jj = 0; jj < N; jj++) {
			// RCdf(i,jj) = prevColRCdf(jj) + 1 for all jj - we have point (i, 0)
			/* comment and edit 2010-01-26:
			   evalPtPos should not matter for the actual cdf - we just take the
			   previous one and increase it by the column probability; the value
			   is taken to be valid in the whole box. evalPtPos comes into play
			   through tgCdfOfR(), where it is used in computing the evaluation points.
			*/
			ECHO ("Cop2DSample::gen_heur(): i=" << i << ",j=0 ; jj=" << jj
			      << "; prevColCdf[" << jj << "]=" << prevColCdf[jj]
			      << ", tgCdfOfR(" << i << "," << jj << ")=" << tgCdfOfR(i, jj));
			dist += cdfDist(prevColCdf[jj] + probCol, tgCdfOfR(i, jj));
		}
		ECHO ("Cop2DSample::gen_heur(): dist(" << i << "," << j << ") = " << dist);
		if (j2iC[j] < 0) {
			// j=0 is available, i.e. it is a valid choice
			bestRows.push_back(j);
			minDist = dist;
		}

		for (j = 1; j < N; j++) {
			// using the recursive formula shown above
			tgVal = tgCdfOfR(i, j-1);
			dist += cdfDist(prevColCdf[j-1], tgVal)
			        - cdfDist(prevColCdf[j-1] + probCol, tgVal);
			ECHO ("Cop2DSample::gen_heur(): dist(" << i << "," << j << ") = "
			      << dist);
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
		ECHO ("Cop2DSample::gen_heur(): new link: (" << i << ", " << j
					 << "); dist = " << minDist << endl);

		// update prevColCdf:
		for (jj = j; jj < N; jj++) {
			prevColCdf[jj] += probCol;
		}
	}
	assert (have_valid_cop() && "should have a valid copula sample now");

	ECHO (endl << "Cop2DSample::gen_heur() finished with distance/error = "
				<< totDist);

	return totDist;
}


// ------------------------------------------------------------------------

void Cop2DSample::CandList::add_item_to_list(int const index,
                                             double const error,
                                             bool const sortList)
{
	if (list.size() + 1 > list.capacity()) {
		// adding a number will lead to re-allocation -> double the capacity
		list.reserve(2 * list.capacity());
	}
	IdxValPair tmpPair(index, error);
	list.push_back(tmpPair);
	if (sortList) {
		sort_list();
	}
}

void Cop2DSample::CandList::insert_cand(int const index, double const error)
{
	if (list.size() < minNumCand) {
		// we do not have enough candidates -> accept without conditions
		add_item_to_list(index, error, true);
	} else {
		// we alsready have enough candidates
		// -> add the new one only if it is good enough
		if (error > list.back().value + maxDiff) {
			return;
		}
		// add the new item; this also sorts the list
		add_item_to_list(index, error, true);

		// remove values from the end of the list that are no longer good enough
		// compare to the last item on the minimal-length list
		while (list.size() >= minNumCand
		       && list.back().value > list[minNumCand-1].value + maxDiff) {
			list.pop_back(); // removes the last (i.e. worst) element
		}
	}
	#ifndef NDEBUG
		cout << "DEBUG: CandList::insert_cand - added item (" << index << ", "
		     << error << ") to the canditate list" << endl;
	#endif
}

// ------------------------------------------------------------------------


/// \todo make this into a private method and call it also from the
///       \a gen_heur() method !!!
void Cop2DSample::cdf_dist_of_col(int const i, Vector const &prevColCdf,
                                  Vector &rowCdfDist, bool rowFree[])
{
	assert (rowFree == NULL && "handling of rowFree is not yet implented...");

	double probCol = (p2prob ? p2prob[scenOfMarg1R[i]] : 1.0 / N);
	int j, jj;
	double tgVal;

	// Using the recursive formula
	// -> have to treat (i,0) separately to start the recursion
	double dist = 0;
	for (j = 0; j < N; j++) {
		if (j == 0) {
			for (jj = 0; jj < N; jj++) {
				// point/link (i, 0) -> Cdf(i,jj) = prevColRCdf(jj) + prob[i]
				dist += cdfDist(prevColCdf[jj] + probCol, tgCdfOfR(i, jj));
				ECHO ("Cop2DSample::cdf_dist_of_col(" << i << "): j=0; jj=" << jj
				      << "; prevColCdf[" << jj << "]=" << prevColCdf[jj]
				      << ", tgCdfOfR(" << i << "," << jj << ")="
				      << tgCdfOfR(i, jj));
			}
		} else {
			// using the recursive formula shown above
			// -> have to evaluate all j's, even if they won't be stored!
			tgVal = tgCdfOfR(i, j-1);
			dist += cdfDist(prevColCdf[j-1], tgVal)
			        - cdfDist(prevColCdf[j-1] + probCol, tgVal);
		}
		ECHO ("Cop2DSample::cdf_dist_of_col(" << i
		      << "): dist(" << i << "," << j << ") = " << dist);
		if (rowFree == NULL || rowFree[j])
			rowCdfDist[j] = dist;
		else
			rowCdfDist[j] = DblInf;
	}
}


/// \todo make this into a private method and call it also from the
///       \a gen_heur() method ???
void Cop2DSample::cdf_dist_of_row(int const j, Vector const &prevRowCdf,
                                  Vector &colCdfDist, bool colFree[])
{
	assert (colFree == NULL && "handling of colFree is not yet implented...");

	double probRow = (p2prob ? p2prob[scenOfMarg1R[j]] : 1.0 / N); /// \todo CHECK !
	int i, ii;
	double tgVal;

	// Using the recursive formula
	// -> have to treat (i,0) separately to start the recursion
	double dist = 0;
	for (i = 0; i < N; i++) {
		if (i == 0) {
			for (ii = 0; ii < N; ii++) {
				// point/link (0, j) -> Cdf(ii,j) = prevRowRCdf(ii) + prob[j]
				dist += cdfDist(prevRowCdf[ii] + probRow, tgCdfOfR(ii, j));
				ECHO ("Cop2DSample::cdf_dist_of_row(" << j << "): i=0; ii=" << ii
				      << "; prevRowCdf[" << ii << "]=" << prevRowCdf[ii]
				      << ", tgCdfOfR(" << ii << "," << j << ")="
				      << tgCdfOfR(ii, j));
			}
		} else {
			// using the recursive formula shown above
			// -> have to evaluate all i's, even if they won't be stored!
			tgVal = tgCdfOfR(i-1, j);
			dist += cdfDist(prevRowCdf[i-1], tgVal)
			        - cdfDist(prevRowCdf[i-1] + probRow, tgVal);
		}
		ECHO ("Cop2DSample::cdf_dist_of_row(" << j
		      << "): dist(" << i << "," << j << ") = " << dist);
		if (colFree == NULL || colFree[i])
			colCdfDist[i] = dist;
		else
			colCdfDist[i] = DblInf;
	}
}


/// \todo do this; return false if this breaks another link
bool Cop2DSample::add_link(int const colR, int const rowR)
{
	if (i2jC[colR] >= 0 || j2iC[rowR] >= 0) {
		cout.flush();
		cerr << "Error in add_link(" << colR << "," << rowR
		     << ") - link already exists!" << endl;
		cerr.flush();
		return false;
	}

	ECHO (sampleId << " : new link (" << colR << "," << rowR << ")");

	i2jC[colR] = rowR;
	j2iC[rowR] = colR;

	return true;
}
