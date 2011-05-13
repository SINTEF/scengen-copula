#ifndef COMMON_HPP
#define COMMON_HPP

#include <iostream>
#include <cassert>
#include <string>

#ifdef NDEBUG
	// disable range-checking for boost - must come before other boost headers
	#define BOOST_DISABLE_ASSERTS
#endif
#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp> // matrix rows and columns
//#include <boost/numeric/ublas/io.hpp> // input/output methods
// new matrix types definitions - everything should converts to these!
namespace ublas = boost::numeric::ublas; // shortcut name
typedef ublas::vector<double> UVector;
typedef ublas::matrix<double> UMatrix;
typedef ublas::vector<int> UIVector;
typedef ublas::matrix<int> UIMatrix;
typedef boost::shared_ptr<UIMatrix> UIMatrixPtr;
typedef UVector::size_type DimT;
DimT const MaxVecLen = std::numeric_limits<UVector::size_type>::max();

/// \name input-output routines for ublas objects
///@{
	/// stream output for ublas vectors
	template<class T>
	std::ostream & operator<< (std::ostream & os, ublas::vector<T> & v);

	/// stream output for ublas matrices
	template<class T>
	std::ostream & operator<< (std::ostream & os, ublas::matrix<T> & M);

	/// stream input for ublas vectors
	template<class T>
	std::istream & operator>> (std::istream & is, ublas::vector<T> & v);

	/// output ublas vectors
	template<class T>
	std::istream & operator>> (std::istream & is, ublas::matrix<T> & M);
///@}


// at a couple of places, we use the boost multi_array for more flexibility
#include <boost/multi_array.hpp>


/// Random number generator - integer from 0 to Max-1
/**
	 rand() is the only random number generator in C++ <br>
	 Note that rand() is from 0 to RAND_MAX inclusive!
**/
#define irand(Max) (rand() % (Max))

/// Random number generator - double from [0, 1)
// note that RAND_MAX is casted to double before we add 1,
// in case that RAND_MAX = MAXINT
#define frand() ((double) rand() / ((double) RAND_MAX + 1))

const double DblEps = 1e-9; // for comparing doubles
bool isEq(double const x, double const y);

#include <limits>
const double DblInf = sqrt(std::numeric_limits<double>::max());

// Debug messages
// Note that the definition does not end with ";", so we need one when used!
#ifndef NDEBUG
	#define ECHO(message) cout << "DEBUG: " << message << endl; cout.flush()
#else
	#define ECHO(message)
#endif

/// compute ranks of a given vector
/**
	ranks go from 0 to N-1
**/
void get_ranks(const std::vector<double> &inputVect, std::vector<int> &ranks);

/// compute ranks of a given vector
/**
	ranks go from 0 to N-1
**/
//void get_ranks(double const inputVect[], std::vector<int> ranks);

/// compute ranks of rows of a matrix (i.e. taking one row at a time)
void get_ranks_or_rows(std::vector<UVector> const & valMat, UIMatrix & rankMat);

/// compute ranks of rows of a matrix (i.e. taking one row at a time)
void get_ranks_or_rows(UMatrix const & valMat, UIMatrix & rankMat);

/// compute ranks of rows of a matrix (i.e. taking one row at a time)
//void get_ranks_or_rows(TMatrixD const & valMat, TMatrixI & rankMat);

/// compute ranks of rows of a matrix (i.e. taking one row at a time)
/// delete this version when it is no longer needed!
//void get_ranks_or_rows(std::vector< std::vector<double> > const & valMat,
//                       TMatrixI & rankMat);


/// position of the mass inside the discretization intervals
//double const discrPt = 0.5; // 0.5 means in the middle

/// convert rank {-1,0,..,N-1} into number from [0,1]
/**
	This is used for computing probabilities of scenarios, so it places the
	points in the top-right corners of the squares (to get the whole area).

	\param[in] r input rank value
	\param[in] N size of the grid (max. rank)
	\return value between zero and one (including the two points)
	For compatibility with \c u012Rank(), we allow r = -1, which returns 0.0 !
**/
inline double rank2U01(double const r, int const N) {
	return (r + 1.0) / (double) N;
}

/// convert unif[0,1] values to ranks {-1,0,..,N-1}
/**
	With N scenarios, we have N possible ranks, which means N intervals - but
	those intervals have N+1 border points. To deal with this, the formula
	implies <code>u012Rank(0.0, N) = -1</code> for all N. This is also consistent
	with the fact that z=0.0 should mean zero probability, while each rank has
	some probability.
	\warning Note that this means that all codes using this function must be
	         prepared to get -1 in return!
**/
inline int u012Rank(double const z, int const N) {
	// subtract epsilon, to avoid numerical issues with ceil()
	return static_cast<int>(ceil(N * (z - DblEps))) - 1;
}


// print vectors using cout
template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
	copy(v.begin(), v.end(), std::ostream_iterator<T>(std::cout, " "));
	return os;
}

#endif
