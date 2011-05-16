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

#include <boost/numeric/ublas/symmetric.hpp>

/// stream input for ublas matrices
template<class T>
std::istream & operator>> (std::istream & is, ublas::symmetric_matrix<T> & M)
{
	// using >> for a standard matrix
	typename ublas::matrix<T> tmpM(M.size1(), M.size2());
	is >> tmpM;
	cout << endl << "tmpM = " << tmpM << endl;
	DimT N = tmpM.size1();
	assert (tmpM.size2() == N && "must be a square matrix");

	M.resize(N);
	for (DimT i = 0; i < N; ++i) {
		for (DimT j = 0; j <= i; ++j) {
			M(i, j) = tmpM(i, j);
			assert (isEq(tmpM(j, i), M(i, j)) && "must be symmetric");
		}
	}
	return is;
}

// list of instances to compile ("explicit instantiation")
template std::ostream & operator<< (std::ostream &, VectorD &);
template std::ostream & operator<< (std::ostream &, VectorI &);
template std::ostream & operator<< (std::ostream &, MatrixD &);
template std::ostream & operator<< (std::ostream &, MatrixI &);
//
template std::istream & operator>> (std::istream &, VectorD &);
template std::istream & operator>> (std::istream &, VectorI &);
template std::istream & operator>> (std::istream &, MatrixD &);
template std::istream & operator>> (std::istream &, MatrixI &);
//
template std::istream & operator>> (std::istream &,
                                    ublas::symmetric_matrix<double> &);

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


void get_ranks_or_rows(std::vector<VectorD> const & valMat, MatrixI & rankMat)
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


void get_ranks_or_rows(MatrixD const & valMat, MatrixI & rankMat)
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
