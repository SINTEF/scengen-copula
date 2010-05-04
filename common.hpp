#ifndef COMMON_HPP
#define COMMON_HPP


// matrix types definitions
#ifdef NDEBUG
	// disable range-checking for boost arrays
	// must come before #include <boost/multi_array.hpp>
	#define BOOST_DISABLE_ASSERTS
#endif
#include <boost/multi_array.hpp>
typedef std::vector<double> Vector;
typedef std::vector<int> IVector;
typedef boost::multi_array<double, 2> Matrix;
typedef boost::multi_array<int, 2> IMatrix;

// use this to get arrays starting at -1: boost::extents[Range(-1,N)]
typedef boost::multi_array_types::extent_range Range;


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

const double DblEps = 1e-8; // for comparing doubles
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

#endif
