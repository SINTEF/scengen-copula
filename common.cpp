//#include <cmath>
#include <ostream>

#include "common.hpp"
#include "ranker.h"

using std::cout;
using std::cerr;
using std::endl;


// ----------------------------------------------------------------------------
template<class T>
std::ostream & operator<< (std::ostream & os, ublas::vector<T> & v)
{
	typename ublas::vector<T>::iterator it;
	for (it = v.begin(); it != v.end(); ++it) {
		os << *it << " ";
	}
	return os;
}

/// stream output for ublas matrices
template<class T>
std::ostream & operator<< (std::ostream & os, ublas::matrix<T> & M)
{
	typename ublas::matrix<T>::iterator1 it1;
	typename ublas::matrix<T>::iterator2 it2;
	os << endl;
	for (it1 = M.begin1(); it1 != M.end1(); ++it1) {
		for (it2 = it1.begin(); it2 != it1.end(); ++it2) {
			os << *it2 << "\t";
		}
		os << endl;
	}
	return os;
}

/// stream input for ublas vectors
template<class T>
std::istream & operator>> (std::istream & is, ublas::vector<T> & v)
{
	if (v.size() == 0) {
		// vector not allocated -> assume the size is in the file!
		DimT vLen;
		is >> vLen;
		v.resize(vLen);
	}
	typename ublas::vector<T>::iterator it;
	for (it = v.begin(); it != v.end(); ++it) {
		is >> *it;
	}
	return is;
}

/// stream input for ublas matrices
template<class T>
std::istream & operator>> (std::istream & is, ublas::matrix<T> & M)
{
	if (M.size1() * M.size2() == 0) {
		// matrix not allocated -> assume sizes are the first numbers in 'is'!
		DimT nRows, nCols;
		is >> nRows >> nCols;
		M.resize(nRows, nCols);
	}
	typename ublas::matrix<T>::iterator1 it1;
	typename ublas::matrix<T>::iterator2 it2;
	for (it1 = M.begin1(); it1 != M.end1(); ++it1) {
		for (it2 = it1.begin(); it2 != it1.end(); ++it2) {
			is >> *it2;
		}
	}
	return is;
}

// list of instances to compile ("explicit instantiation")
template std::ostream & operator<< (std::ostream &, UVector &);
template std::ostream & operator<< (std::ostream &, UIVector &);
template std::ostream & operator<< (std::ostream &, UMatrix &);
template std::ostream & operator<< (std::ostream &, UIMatrix &);
//
template std::istream & operator>> (std::istream &, UVector &);
template std::istream & operator>> (std::istream &, UIVector &);
template std::istream & operator>> (std::istream &, UMatrix &);
template std::istream & operator>> (std::istream &, UIMatrix &);


// ----------------------------------------------------------------------------
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


void get_ranks_or_rows(std::vector<UVector> const & valMat, UIMatrix & rankMat)
{
	unsigned nR = valMat.size();
	unsigned nC = valMat[0].size(); // valMat is a vector of vectors, not a matrix!

	// resize rankMat if needed; complain if resizing non-empty matrix
	if (rankMat.size1() != nR || rankMat.size2() != nC) {
		if (rankMat.size1() + rankMat.size2() > 0) {
			cerr << "Warning from get_ranks_or_rows(): wrong size of output matrix"
			     << " - resizing!" << endl;
			rankMat.resize(valMat.size(), valMat[0].size());
		}
	}

	unsigned r, c;

	std::vector<double> valV(nC);
	std::vector<int> rankV(nC);
	for (r = 0; r < nR; ++r) {
		// make a temp. copy of the row
		for (c = 0; c < nC; ++c) {
			valV[c] = valMat[r](c);
		}
		rank (valV, rankV); // this computes the ranks - needs std::vector<>!
		// copy to the target matrix; subtract 1 to get ranks from zero
		for (c = 0; c < nC; ++c) {
			rankMat(r,c) = rankV[c] - 1;
		}
	}
}


void get_ranks_or_rows(UMatrix const & valMat, UIMatrix & rankMat)
{
	// resize rankMat if needed; complain if resizing non-empty matrix
	if (valMat.size1() != rankMat.size1() || valMat.size2() != rankMat.size2()) {
		if (rankMat.size1() + rankMat.size2() > 0) {
			cerr << "Warning from get_ranks_or_rows(): wrong size of output matrix"
			     << " - resizing!" << endl;
		}
		rankMat.resize(valMat.size1(), valMat.size2());
	}

	int nR = valMat.size1();
	int nC = valMat.size2();
	int r, c;

	std::vector<double> valV(nC);
	std::vector<int> rankV(nC);
	for (r = 0; r < nR; ++r) {
		// make a temp. copy of the row
		for (c = 0; c < nC; ++c) {
			valV[c] = valMat(r,c);
		}
		rank (valV, rankV); // this computes the ranks - needs std::vector<>!
		// copy to the target matrix; subtract 1 to get ranks from zero
		for (c = 0; c < nC; ++c) {
			rankMat(r,c) = rankV[c] - 1;
		}
	}
}


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

// '<<' operator for an matrix
/// \todo .. make it work for the template
/// \todo .. add '>>' as well
std::ostream & operator<<(std::ostream & output, const TMatrixD & M) {
	for (unsigned i = 0; i < M.num_rows(); ++i) {
		for (unsigned j = 0; j < M.num_cols(); ++j) {
			if (j > 0)
				output << "\t";
			output << M[i][j];
		}
		output << std::endl;
	}
	return output;  // for multiple << operators.
}
