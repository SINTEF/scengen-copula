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
#include <boost/make_shared.hpp>

// ---------------------------------------------------------------------------
// vector and matrix definitions, using boost::numeric::ublas
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp> // matrix rows and columns

namespace ublas = boost::numeric::ublas; // shortcut name
typedef ublas::vector<double> VectorD;
typedef ublas::matrix<double> MatrixD;
typedef ublas::vector<int> VectorI;
typedef ublas::matrix<int> MatrixI;
typedef boost::shared_ptr<MatrixI> MatrixIPtr;
typedef VectorD::size_type DimT;
DimT const MaxVecLen = std::numeric_limits<VectorD::size_type>::max();

/// \name input-output routines for ublas objects
/**
	These methods exist in <boost/numeric/ublas/io.hpp>, but they use
	a different format (put a whole matrix on one line, ..). We use
	ublas::vector_expression<> and ublas::matrix_expression<> instead of
	ublas::vector<> and ublas::matrix<>, so we can print results of operations
	like trans() etc.

	The implementations of the output operations are here in the header file,
	to avoid explicit instantiation in the .cpp file. For example, we would need
	to add a new type if we wanted to print a transposed matrix, since trans()
	returns a special (and rather cryptic) type..
**/
///@{
	/// stream output for ublas vectors (\a VT is the vector class)
	template<class VT>
	std::ostream & operator<< (std::ostream & os,
	                           const ublas::vector_expression<VT> & v)
	{
		// Note: the vector_expression<VT> has only one method "()", which
		// returns "& VT" or "const & VT" - a ref. to the included vector object.
		for (auto it = v().begin(); it != v().end(); ++it) {
			os << *it << " ";
		}
		return os;
	}

	/// stream output for ublas matrices (\a MT is the matrix class)
	template<class MT>
	std::ostream & operator<< (std::ostream & os,
	                           const ublas::matrix_expression<MT> & M)
	{
		// Note: the matrix_expression<MT> has only one method "()", which
		// returns "& MT" or "const & MT" - a ref. to the included matrix object.
		for (auto it1 = M().begin1(); it1 != M().end1(); ++it1) {
			for (auto it2 = it1.begin(); it2 != it1.end(); ++it2) {
				os << *it2 << "\t";
			}
			os << std::endl;
		}
		return os;
	}

	/// stream input for ublas vectors
	template<class T>
	std::istream & operator>> (std::istream & is, ublas::vector<T> & v);

	/// stream input for ublas matrices
	template<class T>
	std::istream & operator>> (std::istream & is, ublas::matrix<T> & M);

	/// stream input for ublas symmetric matrices
	template<class T>
	std::istream & operator>> (std::istream & is,
	                           ublas::symmetric_matrix<T> & M);
///@}

/// \name other functions related to the vector and matrix objects
///@{
	template<typename T>
	double vec_mean(ublas::vector<T> const & v);

	template<typename T>
	double vec_std_dev(ublas::vector<T> const & v, double const mean,
	                   bool const unbiased = false);
///@}


// ---------------------------------------------------------------------------
// at a couple of places, we use the boost multi_array for more flexibility
#include <boost/multi_array.hpp>
template <typename T>
class Array2D : public boost::multi_array<T, 2> {
	private:
		typedef boost::multi_array_types::index_range IRange;
		typedef boost::multi_array_types::extent_range ERange;
		typename boost::multi_array<T, 2>::index_gen indices;

	public:
		Array2D (DimT const nR, DimT const nC)
			: boost::multi_array<T, 2>(boost::extents[nR][nC]) {}

		Array2D (DimT const nR, DimT const nC, int const i0)
			: boost::multi_array<T, 2>(boost::extents[ERange(i0,nR)]
			                                         [ERange(i0,nC)]) {}
};


// ---------------------------------------------------------------------------
// in addition, we also use std::vector at a couple of places
/// print std::vector to a stream
template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
	copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " "));
	return os;
}


// ---------------------------------------------------------------------------
// generic macros and difinitions
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


// ---------------------------------------------------------------------------
// misc. functions for the copula-generation code
/// compute ranks of a given vector
/**
	ranks go from 0 to N-1
**/
void get_ranks(const std::vector<double> &inputVect, std::vector<int> &ranks);

/// compute ranks of rows of a matrix (i.e. taking one row at a time)
void get_ranks_or_rows(std::vector<VectorD> const & valMat, MatrixI & rankMat);

/// compute ranks of rows of a matrix (i.e. taking one row at a time)
void get_ranks_or_rows(MatrixD const & valMat, MatrixI & rankMat);

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

/// fix mean and standard deviation (if given) of a sample (in-place)
/**
	\param          sample  the sample to fix
	\param[in]        mean  the target mean
	\param[in]         stD  the target std. deviation; -1 means do not fix
	\param[in] unbiasedStD  should we use the unbiased std. dev. estimator?

	\note At the moment, this assumes equiprobable sample values
**/
void fix_mean_std(VectorD & sample, double const mean, double const stD = -1.0,
                  bool const unbiasedStD = false);


// ---------------------------------------------------------------------------
// Output from the code (for info and debug)
enum OutputLevel {
	TrNone    =  0, // no output whatever happens
	TrFatal   =  1, // report fatal errors (errors that stop execution)
	TrError   =  2, // report all errors
	TrWarn    =  3, // as above + warnings
	TrInfo    =  4, // as above + info
	TrInfo2   =  5, // as above + more info
	TrInfo3   =  6, // as above + more info
	TrDetail  =  7, // as above + more detailed info
	TrDetail2 =  8, // as above + more detailed info
	TrDetail3 =  9, // as above + more detailed info
	TrAll     = 10  // absolute all messages
};
extern OutputLevel outLvl;
// Note that the definitions do not end with ";", so we need one when used!
// 'internal' macro, used by the others
#define ECHO(message) std::cout << message << std::endl; std::cout.flush()
// this will produce message in all profiles (debug and release)
#define MSG(lvl, txt) if (lvl <= outLvl) { ECHO(txt); }
// this will produce message only in debug profiles -> can be used in parts
// where we do not waste time for checking outLvl in the release versions
#ifndef NDEBUG
	#define TRACE(lvl, txt) if (lvl <= outLvl) { ECHO(txt); }
#else
	#define TRACE(lvl, txt)
#endif
// this is for temporary tracing
#define DBGOUT(txt) ECHO("DEBUG " << __FILE__ << "/l." << __LINE__ << ": " << txt)
// macros for displaying variable values (stand-alone and one for use in cout)
// see http://hamersun.blogspot.com/2010/09/
#define DISPLAY(var) std::cout << #var " = " << var << std::endl
#define VALUE(var) #var " = " << var


#endif
