//#include <cmath>
#include "common.hpp"
#include "ranker.h"

//bool isEq(double const x, double const y) { return abs(x - y) < DblEps; }
bool isEq(double const x, double const y)
{
	// warning: abs() is integer, this needs std::abs() or fabs()
	return std::abs(x - y) < DblEps;
}

void get_ranks(const std::vector<double> &inputVect, std::vector<int> &ranks) {
	// using the rank() from "ranker.h"
	// this gives numbers from 1 to N, so we have to subtract 1!
	rank(inputVect, ranks);
	unsigned len = ranks.size();
	for (unsigned i = 0; i < len; i++) {
		ranks[i]--;
	}
}

/*
void get_ranks(double const inputVect[], std::vector<int> ranks); {
	size_t len = sizeof(inputVect) / sizeof(double);
	rank(inputVect, len, ranks);
	for (unsigned i = 0; i < len; i++) {
		ranks[i]--;
	}
}
*/

void get_ranks_or_rows(TMatrixD const & valMat, TMatrixI & rankMat)
{
	assert (valMat.num_rows() == rankMat.num_rows()
	     && valMat.num_cols() == rankMat.num_cols()
	     && "dimension check");

	int nR = valMat.num_rows();
	int nC = valMat.num_cols();
	int r, c;

	std::vector<double> valV(nC);
	std::vector<int> rankV(nC);
	for (r = 0; r < nR; ++r) {
		// make a temp. copy of the row
		for (c = 0; c < nC; ++c) {
			valV[c] = valMat[r][c];
		}
		rank (valV, rankV); // this computes the ranks
		// copy to the target matrix; subtract 1 to get ranks from zero
		for (c = 0; c < nC; ++c) {
			rankMat[r][c] = rankV[c] - 1;
		}
	}
}

void get_ranks_or_rows(std::vector< std::vector<double> > const & valMat,
                       TMatrixI & rankMat)
{
	assert (valMat.size() == rankMat.num_rows()
	     && valMat[0].size() == rankMat.num_cols()
	     && "dimension check");

	int nR = rankMat.num_rows();
	int nC = rankMat.num_cols();
	int r, c;

	std::vector<int> rankV(nC);
	for (r = 0; r < nR; ++r) {
		rank (valMat[r], rankV); // this computes the ranks
		// copy to the target matrix; subtract 1 to get ranks from zero
		for (c = 0; c < nC; ++c) {
			rankMat[r][c] = rankV[c] - 1;
		}
	}
}
