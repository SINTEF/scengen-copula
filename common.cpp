#include "common.hpp"
#include "ranker.h"

#include <boost/numeric/ublas/symmetric.hpp>
#include <ostream>
//#include <cmath>

using std::cout;
using std::cerr;
using std::endl;


// ----------------------------------------------------------------------------
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
		if (is.fail()) {
			// stream is in a failed state _before_ the read!
			is.clear();       // reset the state
			string nonNumber;
			is >> nonNumber;  // read the offending string
			*it = std::numeric_limits<T>::signaling_NaN();
		}
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
			if (is.fail()) {
				// stream is in a failed state _before_ the read!
				is.clear();       // reset the state
				string nonNumber;
				is >> nonNumber;  // read the offending string
				*it2 = std::numeric_limits<T>::signaling_NaN();
			}
		}
	}
	return is;
}

/// stream input for ublas symmetric matrices
template<class T>
std::istream & operator>> (std::istream & is, ublas::symmetric_matrix<T> & M)
{
	// using >> for a standard matrix
	typename ublas::matrix<T> tmpM(M.size1(), M.size2());
	is >> tmpM;
	DimT N = tmpM.size1();
	assert (tmpM.size2() == N && "must be a square matrix");

	M.resize(N);
	for (DimT i = 0; i < N; ++i) {
		for (DimT j = 0; j <= i; ++j) {
			M(i, j) = tmpM(i, j);
			assert (isEq(tmpM(j, i), M(i, j)) && "must be symmetric");
		}
	}
	TRACE (TrDetail2, "symmetric matrix read from a stream: " << endl << M);
	return is;
}

// list of instances to compile ("explicit instantiation")
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
	// the rank function support four different methods of resolving ties:
	// - average - rank of equal values is the average of their ranks [default]
	//           - can be non-integral; will be casted to the target type
	// - min,max - rank of equal values is the min/max rank
	// - random  - rank of equal values is not equal (depends on order?)
	//! \todo do we always want the tie-resolving variant ("random")?
	rank(inputVect, ranks, "random");
	auto len = ranks.size();
	for (unsigned i = 0; i < len; i++) {
		ranks[i]--;
	}
}


void get_ranks_or_rows(std::vector<VectorD> const & valMat, MatrixI & rankMat)
{
	DimT nR = valMat.size();
	DimT nC = valMat[0].size(); // valMat is a vector of vectors, not a matrix!
	DimT r, c;

	// resize rankMat if needed; complain if resizing non-empty matrix
	if (rankMat.size1() != nR || rankMat.size2() != nC) {
		if (rankMat.size1() + rankMat.size2() > 0) {
			cerr << "Warning from get_ranks_or_rows(): wrong size of output matrix"
			     << " - resizing!" << endl;
			rankMat.resize(valMat.size(), valMat[0].size());
		}
	}

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

	DimT nR = valMat.size1();
	DimT nC = valMat.size2();
	DimT r, c;

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


template<typename T>
double vec_mean(ublas::vector<T> const & v)
{
	double mu = 0.0;
	for (auto it = v.begin(); it != v.end(); ++it)
		mu += *it;
	return mu / v.size();
}
template double vec_mean(VectorD const &); // compile this for VectorD


template<typename T>
double vec_std_dev(ublas::vector<T> const & v, double const mean,
                   bool const unbiased)
{
	DimT N = v.size();
	double stD = 0.0;
	for (auto it = v.begin(); it != v.end(); ++it)
		stD += pow(*it, 2);
	stD /= N;            // std = E[X^2]
	stD -= pow(mean, 2); // std = E[X^2] - (EX)^2 = Var(X)
	if (stD < 0)
		throw std::domain_error("negative variance in std_dev()!");
	if (unbiased)
		stD *= static_cast<double>(N) / (N - 1.0); // unbiased (instead of MLE)
	return sqrt(stD);
}
template double vec_std_dev(VectorD const &, double const, bool const);


void fix_mean_std(VectorD & sample, double const mean, double const stD,
                  bool const unbiasedStD)
{
	double curMean = vec_mean(sample);
	if (stD < 0.0) { // this means std. dev. not given -> fix only the mean
		// ublas does NOT have "vector + scalar" overloaded operators!
		// -> we add by creating a temp. vector with the same value at all el.
		sample += ublas::scalar_vector<double>(sample.size(), mean - curMean);
	} else {
		double curStD = vec_std_dev(sample, curMean, unbiasedStD);
		if (curStD < DblEps) { // zero variance -> problems
			if (stD < DblEps) { // the target is also zero -> just fix the mean
				// again creating a temp. scalar vector to do "vector += constant"
				sample += ublas::scalar_vector<double>(sample.size(),
				                                       mean - curMean);
				return;
			}
		}
		double scaling = stD / curStD;
		sample *= scaling;
		// again creating a temp. scalar vector to do "vector += constant"
		sample += ublas::scalar_vector<double>(sample.size(),
		                                       mean - scaling * curMean);
	}
}
