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


#include <boost/multi_array.hpp>
typedef std::vector<double> Vector;
typedef std::vector<int> IVector;
typedef boost::multi_array<double, 2> Matrix;
typedef boost::multi_array<int, 2> IMatrix;

/// template for a vector class based on boost
template <typename T>
class TVectorGen : public boost::multi_array<T, 1> {
	public:
		TVectorGen (int const len)
		: boost::multi_array<T, 1>(boost::extents[len])
		{}

		/// copy-constructor from the base type, using the base-type constructor
		TVectorGen (boost::multi_array<T, 1> const &X)
		: boost::multi_array<T, 1>(X)
		{}

		/// copy-constructor from the array view, using the base-type constructor
		TVectorGen (boost::detail::multi_array::multi_array_view<T, 1ul> const &X)
		: boost::multi_array<T, 1>(X)
		{}

/*
		/// assignment using the parent boost class
		TVectorGen & operator= (boost::multi_array<T, 1> &X) {
			if (this != &X) {
				this->boost::multi_array<T, 1>::operator= (X);
			}
			return *this;
		}
*/
};
typedef TVectorGen<double> TVectorD;
typedef TVectorGen<int> TVectorI;


/// matrix class using the boost libraries
/**
	Note that this allows to have methods for column and rows that return
	boost vectors !
**/
template <typename T>
class TMatrixGen : public boost::multi_array<T, 2> {
	private:
		typedef boost::multi_array_types::index_range IRange;
		typename boost::multi_array<T, 2>::index_gen indices;

	public:
		TMatrixGen (int const nR, int const nC)
			: boost::multi_array<T, 2>(boost::extents[nR][nC])
			{}

		/// get a vector consisting of a column of the matrix
		/**
			!!! This version creates a new objects, allocating new data !!!
			An alternative is to return \c boost::multi_array<T,1>
		*/
		TVectorGen<T> col(int const c) const {
			return operator[](indices[IRange()][c]);
		}

		/// get a view consisting of a column of the matrix
		/**
			This version does not copy data, i.e. the resulting object points
			to the data structures of the matrix.
			\todo Is it possible to do with a vector?
		*/
		boost::detail::multi_array::multi_array_view<T, 1> col(int const c) {
			return operator[](indices[IRange()][c]);
		}

		/// get a vector consisting of a row of the matrix
		/**
			!!! This version creates a new objects, allocating new data !!!
			An alternative is to return \c boost::multi_array<T,1>
		*/
		TVectorGen<T> row(int const r) const {
			return operator[](indices[r][IRange()]);
		}

		/// get a view consisting of a row of the matrix
		/**
			This version does not copy data, i.e. the resulting object points
			to the data structures of the matrix.
		*/
		boost::detail::multi_array::multi_array_view<T, 1> row(int const r) {
			return operator[](indices[r][IRange()]);
		}

		unsigned num_rows() const { return this->shape()[0]; }
		unsigned num_cols() const { return this->shape()[1]; }

};
typedef TMatrixGen<double> TMatrixD;
typedef TMatrixGen<int> TMatrixI;
typedef boost::multi_array<double, 1>::array_view<1>::type TMatSliceD;
typedef boost::multi_array<int, 1>::array_view<1>::type TMatSliceI;
typedef boost::shared_ptr<TMatrixD> TMatrixDPtr; // smart pointer to TMatrixD
typedef boost::shared_ptr<TMatrixI> TMatrixIPtr; // smart pointer to TMatrixI

/// '<<' operator for an matrix.
/// \todo .. make it work for the template
/// \todo .. add '>>' as well
std::ostream & operator<<(std::ostream & output, const TMatrixD & M);

// use this to get arrays starting at -1: boost::extents[ExtRange(-1,N)]
//typedef boost::multi_array_types::extent_range ExtRange;


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
void get_ranks_or_rows(TMatrixD const & valMat, TMatrixI & rankMat);

/// compute ranks of rows of a matrix (i.e. taking one row at a time)
/// delete this version when it is no longer needed!
void get_ranks_or_rows(std::vector< std::vector<double> > const & valMat,
                       TMatrixI & rankMat);


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
