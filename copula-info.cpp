#include <iostream>
//#include <fstream>
//#include <cmath>
//#include <algorithm>
//#include <deque>
//#include <vector>
//#include <cassert>

#include "copula-info.hpp"

using namespace std;
using namespace CopulaDef;


CopInfoData::CopInfoData(TMatrixI & ranksMat, bool const fillU01Data)
: CopulaInfo(),
  ranks(ranksMat), u01Data(ranksMat.num_rows(), ranksMat.num_cols()),
  nVars(ranksMat.num_rows()), nPts(ranksMat.num_cols())
{
	if (fillU01Data) {
		int i, j;
		for (i = 0; i < nVars; ++i) {
			for (j = 0; j < nPts; ++j) {
				u01Data[i][j] = rank2U01(ranks[i][j], nPts);
			}
		}
	}
}


/**
	Note that we have two ways of deciding whether a given ranks[i][j] is below
	the specified value u[i]:
	# check if <code>ranks[i][j] <= u012Rank(u[i], nPts)</code>
	# check if <code>rank2U01(ranks[i][j], nPts) <= u[i]</code>
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
		double lastNon1;
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

	assert (u01Data.num_cols() == ranks.num_cols()
	        && u01Data.num_rows() == ranks.num_rows()
	        && isEq(u01Data[0][0], rank2U01(ranks[0][0], nPts))
	        && "sanity check");

	// Count the number of data points below the given point
	int nPtsBelow = 0;
	for (j = 0; j < nPts; ++j) {
		i = 0;
		while (i < nVars && rank2U01(ranks[i][j], nPts) <= u[i]) {
			i++;
		}
		if (i == nVars) {
			nPtsBelow += 1; // sample point 'j' is below uR in all dimensions;
		}
	}

	return (double) nPtsBelow / (double) nPts;
}
