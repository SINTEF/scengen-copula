#include <iostream>
#include <fstream>
//#include <cmath>
//#include <algorithm>
//#include <deque>
//#include <vector>
//#include <cassert>

#include "copula-info.hpp"

using namespace std;
using namespace CopulaDef;


CopInfoData::CopInfoData(std::string & tgFName)
: CopInfoBy2D()
{
	try {
		read_tg_file(tgFName); // fills hData
	}
	catch(exception& e) {
		cerr << e.what() << endl;
		exit(1);
	}

	fill_ranks_etc(); // fills hRanks and hU01
}


CopInfoData::CopInfoData(UIMatrix const & hDataMat)
: CopInfoBy2D(hDataMat.size1()),
  nPts(hDataMat.size2()),
  hData(nVars), hRanks(nVars, nPts), hU01(nVars, nPts)
{
	for (int i = 0; i < nVars; ++i) {
		hData[i] = ublas::row(hDataMat, i);
	}

	fill_ranks_etc(); // fills hRanks and hU01
}

/*
CopInfoData::CopInfoData(TMatrixI & hRanksMat, bool const fillU01Data)
: CopulaInfo(),
  hRanks(hRanksMat), hU01(hRanksMat.num_rows(), hRanksMat.num_cols()),
  nVars(hRanksMat.num_rows()), nPts(hRanksMat.num_cols())
{
	if (fillU01Data) {
		int i, j;
		for (i = 0; i < nVars; ++i) {
			for (j = 0; j < nPts; ++j) {
				hU01[i][j] = rank2U01(hRanks[i][j], nPts);
			}
		}
	}
}
*/


/**
	Note that we have two ways of deciding whether a given hRanks[i][j] is below
	the specified value u[i]:
	# check if <code>hRanks[i][j] <= u012Rank(u[i], nPts)</code>
	# check if <code>rank2U01(hRanks[i][j], nPts) <= u[i]</code>
	This is not equivalent, since the transformations are not inverses: we can
	have <code>u012Rank(z) = N</code> and <code>rank2U01(N) = z + eps</code>!

	This function uses the comparison of U(0,1) values, i.e. the second option
	from above. This makes it equivalent to the 2D version.
**/
double CopInfoData::cdf(TVectorD const u) const
{
	assert (u.size() == (unsigned) nVars);

	int i, j;

	// if any of the margins is zero, the cdf is zero
	for (i = 0; i < nVars; ++i) {
		if (isEq(u[i], 0.0)) {
			return 0.0;
		}
	}

	// If nVars-1 margins are equal to 1.0, return the remaining margin
	// This can be turned off, if it slows the code down!
	{
		int numNon1s = 0;
		double lastNon1 = 0.0; // init just to avoid gcc's warning
		i = 0;
		while (i < nVars && numNon1s <= 1) {
			if (!isEq(u[i], 1.0)) {
				numNon1s ++;
				lastNon1 = u[i];
			}
			i++;
		}
		if (numNon1s == 0) { return 1.0; }
		if (numNon1s == 1) { return lastNon1; }
	}

	assert (hU01.size1() == hRanks.size1()
	        && hU01.size2() == hRanks.size2()
	        && isEq(hU01(0,0), rank2U01(hRanks(0,0), nPts))
	        && "sanity check");

	// Count the number of data points below the given point
	int nPtsBelow = 0;
	for (j = 0; j < nPts; ++j) {
		i = 0;
		while (i < nVars && rank2U01(hRanks(i,j), nPts) <= u[i]) {
			i++;
		}
		if (i == nVars) {
			nPtsBelow += 1; // sample point 'j' is below uR in all dimensions;
		}
	}

	return (double) nPtsBelow / (double) nPts;
}


void CopInfoData::read_tg_file(string tgFName)
{
	// read the input file
	std::ifstream tgDistF(tgFName.c_str());
	if (!tgDistF) {
		throw ios_base::failure("Could not open input file " + tgFName + "!");
	}
	tgDistF >> nPts >> nVars;

	hData.resize(nVars);
	for (int i = 0; i < nVars; ++i) {
		hData[i].resize(nPts);
	}
	// read from the file - margins are in columns -> must transpose!
	for (int j = 0; j < nPts; ++j) {
		for (int i = 0; i < nVars; ++i) {
			tgDistF >> hData[i](j);
		}
	}
}


void CopInfoData::fill_ranks_etc()
{
	assert (nVars > 0 && nPts > 0 && hData.size1() * hData.size2() > 0
	        && "data must be ready at this point");

	get_ranks_or_rows(hData, hRanks);

	hU01.resize(nVars, nPts);
	for (int i = 0; i < nVars; ++i) {
		for (int j = 0; j < nPts; ++j) {
			hU01(i,j) = rank2U01(hRanks(i,j), nPts);
		}
	}
}


void CopInfoData::setup_2d_targets(bool const makeTranspTgs)
{
	int i, j;
	p2Info2D.resize(boost::extents[nVars][nVars]);
	for (i = 0; i < nVars; i++) {
		for (j = i+1; j < nVars; j++) {
			p2Info2D[i][j].reset(new Copula2D::Cop2DData(hData[i], hData[j], nPts));
			if (makeTranspTgs) {
				assert (p2Info2D[i][j] == NULL);
				// the .get() method returns a raw pointer from a smart one
				p2Info2D[j][i].reset(new Copula2D::Cop2DInfTr(p2Info2D[i][j].get()));
			}
		}
	}
}
